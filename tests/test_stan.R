# tests/test_stan.R
#
# Runs the Stan (NUTS) sampler on the real simulated dataset, printing and
# plotting the same diagnostic battery as test_cpp.R -- same data, same
# hyperparameters/initial values, so Stan is directly comparable by eye and
# by metric to the other four samplers.
#
# Runs N_chains chains with dispersed initializations (see
# make_chain_inits() below, identical to test_cpp.R) and passes all of them
# to print_and_plot_diagnostics(), which computes R_hat/ESS across chains
# and plots only the chain selected via `plot_chain`.
#
# Unlike test_cpp.R, chains are NOT parallelized via foreach/makeCluster
# here: rstan::sampling()'s own chains=/cores= arguments already run the
# N_chains chains in parallel internally (one process per chain), which is
# more efficient than spawning N_chains separate rstan::sampling(chains=1)
# calls -- so sample_stan() takes the full list of chain_inits and returns
# all N_chains results in one call.

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

# See test_cpp.R for why debug=FALSE matters here (elapsed_time feeds
# directly into ESS/second diagnostics further down).
pkgload::load_all("../PoissonLTDM", debug = FALSE)

source("../PoissonLTDM/R/sampler_stan.R")
source("plot_diagnostics.R")

printf <- function(...) cat(paste(sprintf(...), "\n"))

method_grid <- c("amh_montoril", "pg_as", "sir_laplace", "sir_collapsed", "stan")
Tt_grid <- c(200, 400, 800, 1600)
function_grid <- c("constant", "linear", "quadratic", "sinusoidal")

# ---- Choose method and data data here ----
method <- method_grid[5]
f <- function_grid[3]
Tt <- Tt_grid[4]
replica <- 1     # one of: 1,...,200
# ------------------------------------------

# Load the data
source_name <- sprintf("%s_%s_%s", f, Tt, replica)
data <- readRDS(paste("../data/simulated/", source_name, ".rds", sep = ""))
y <- data$y
theta1_true <- data$theta
changepoints <- c(0.25, 0.5, 0.75)*Tt

# Compute the seed (based on pattern)
method_idx <- match(method, method_grid)
Tt_idx <- match(Tt, Tt_grid)
f_idx <- match(f, function_grid)
seed_base <- method_idx*1e5 + f_idx*1e4 + Tt_idx*1e3 + replica*10

# ---- Parameters and initialization ----
# General simulation parameters
N_chains <- 3        # number of chains
plot_chain <- 1      # which chain (1..N_chains) to plot / trace in detail
N <- 11000           # number of iterations (Stan: total, i.e. warmup + post-warmup)
burnin <- 1000
plots <- TRUE
compute_rhat <- TRUE
compute_ess <- TRUE
verbose <- FALSE     # TRUE: show all Stan output; FALSE: suppressed/invisible

# Parallel execution of the N_chains chains.
#
# Stan chains are parallelized natively via rstan::sampling(chains=,
# cores=), not via foreach/makeCluster (see file header). n_cores = NULL
# (default) picks the number of PHYSICAL cores automatically, same
# resolve_n_cores() logic as test_cpp.R. Set n_cores explicitly (e.g.
# n_cores = N_chains) to override.
n_cores <- NULL

# Prior hyperparameters
# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- 0
sigma2_01 <- 100

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 0
sigma2_02 <- 100

# phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
nu_01  <- 2
eta_01 <- 0.01

# phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
nu_02  <- 2
eta_02 <- 0.0001

printf("Running %s for %s, seed_base=%d, N=%d, burnin=%d",
       method, source_name, seed_base, N, burnin)


# ---- Chain initialization ----
#
# Reference values derived from the data (method-of-moments style), used
# only as the CENTER around which each chain's starting point is dispersed.
# Not used directly as an initial value in any chain. Identical scheme to
# test_cpp.R, so Stan's initialization is directly comparable to the other
# samplers.
theta1_ref <- log(y + 0.5)
theta2_ref <- c(diff(theta1_ref), 0)
theta_01_ref <- theta1_ref[1]
theta_02_ref <- theta2_ref[1]
W1_ref <- var(diff(theta1_ref))
W2_ref <- var(diff(theta2_ref))

make_chain_inits <- function(chain_id, seed) {
    set.seed(seed)

    theta_01_init <- theta_01_ref + rnorm(1, 0, sd = 2 * sqrt(abs(theta_01_ref) + 1))
    theta_02_init <- theta_02_ref + rnorm(1, 0, sd = 2 * sqrt(abs(theta_02_ref) + 1))
    W1_init <- W1_ref * exp(rnorm(1, 0, sd = 1))
    W2_init <- W2_ref * exp(rnorm(1, 0, sd = 1))

    theta1_init <- numeric(Tt) # Stan explores theta1 via HMC starting from
                                # this initial value -- kept at zero, as in
                                # sir_laplace/sir_collapsed, since Stan's
                                # own warmup adaptation (step size, mass
                                # matrix) handles finding the typical set
                                # regardless of this starting point
    theta2_init <- numeric(Tt)

    list(chain_id = chain_id, seed = seed,
         theta_01 = theta_01_init, theta_02 = theta_02_init,
         W1 = W1_init, W2 = W2_init,
         theta1 = theta1_init, theta2 = theta2_init)
}

chain_inits <- lapply(1:N_chains, function(k) make_chain_inits(k, seed_base + (k - 1)))


# ---- Physical core count (used when n_cores is not set explicitly) ----
# Same resolve_n_cores() as test_cpp.R.
resolve_n_cores <- function(n_cores, N_chains) {
    if (!is.null(n_cores)) return(n_cores)

    phys <- NA_integer_
    if (Sys.info()[["sysname"]] == "Linux" && nzchar(Sys.which("lscpu"))) {
        out <- tryCatch(
            system("lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#'", intern = TRUE),
            error = function(e) character(0)
        )
        if (length(out) > 0) phys <- length(unique(out)) # distinct (core,socket) pairs
    }

    if (is.na(phys)) phys <- parallel::detectCores(logical = FALSE)

    if (is.na(phys)) {
        logi <- parallel::detectCores(logical = TRUE)
        phys <- max(1, logi %/% 2)
        printf("Physical core count unavailable from OS; falling back to floor(logical/2) = %d", phys)
    }
    min(phys, N_chains) # no point requesting more workers than chains
}
n_cores_used <- resolve_n_cores(n_cores, N_chains)
printf("Running %d chains in parallel on %d cores (Stan native chains=/cores=)", N_chains, n_cores_used)


# ---- Prepare Stan model ----
options(mc.cores = n_cores_used)
rstan::rstan_options(auto_write = FALSE)

if (file.exists("../cache/poisson_ltdm.rds")) {
    model <- readRDS("../cache/poisson_ltdm.rds")
} else {
    printf("Building the model")
    stan_file <- "../PoissonLTDM/inst/stan/poisson_ltdm.stan"
    model <- rstan::stan_model(file = stan_file, model_name = "PoissonLTDM")
    saveRDS(model, file = "../cache/poisson_ltdm.rds")
}


# NOTE: elapsed_time below uses wall-clock time (Sys.time()), not
# proc.time()'s user.self field -- consistent with test_cpp.R, and
# necessary here too since rstan::sampling(cores=) runs chains in separate
# processes.
start_time <- Sys.time()

stan_out <- sample_stan(
    model      = model,
    y          = y,
    N          = N,
    burnin     = burnin,
    seed       = seed_base,
    N_chains   = N_chains,
    n_cores    = n_cores_used,
    chain_inits = chain_inits,
    mu_01      = mu_01,
    sigma2_01  = sigma2_01,
    mu_02      = mu_02,
    sigma2_02  = sigma2_02,
    nu_01      = nu_01,
    eta_01     = eta_01,
    nu_02      = nu_02,
    eta_02     = eta_02,
    verbose    = verbose)

elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
printf("Total wall-clock time (%d chains, parallel, %d cores): %.2f s",
       N_chains, n_cores_used, elapsed_time)

results <- stan_out$results

if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 400) t_obs <- c(75, 100, 200, 300)
if (Tt == 800) t_obs <- c(100, 300, 500, 700)
if (Tt == 1600) t_obs <- c(400, 800, 1200, 1600)

gc(full = TRUE) # Release unused memory

print_and_plot_diagnostics(results, y, changepoints,
                            theta1_true = theta1_true, theta2_true = NULL,
                            t_obs = t_obs, burnin = burnin, elapsed_time = elapsed_time,
                            plot_chain = plot_chain,
                            nu_01 = nu_01, eta_01 = eta_01,
                            nu_02 = nu_02, eta_02 = eta_02, plots = plots,
                            compute_rhat = compute_rhat, compute_ess = compute_ess)
