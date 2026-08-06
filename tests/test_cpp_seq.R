# tests/test_cpp.R
#
# Runs the chosen method using the Rcpp-compiled implementation
# (PoissonLTDM/src/), on the real simulated dataset, printing and plotting
# the same diagnostic battery as test_prototype_R.R -- same data, same
# hyperparameters/initial values, so the two are directly comparable by
# eye. For a numeric R-vs-C++ comparison, see test_validation.R instead.
#
# Runs N_chains independent chains with dispersed initializations (see
# make_chain_inits() below) and passes all of them to
# print_and_plot_diagnostics(), which computes R_hat/ESS across chains and
# plots only the chain selected via `plot_chain`.

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

# IMPORTANT: load_all() compiles C++ in DEBUG mode by default (-g -O0),
# silently overriding any -O3 in src/Makevars, causing a ~4-5x slowdown of
# every sampler and of rhat_ess_fast(). Always use debug=FALSE here, since
# this script's elapsed_time is used directly in ESS/second diagnostics
# further down. Using pkgload::load_all() directly (not devtools::load_all())
# because devtools::load_all()'s `...` does not reliably forward debug= to
# pkgload in this environment (devtools 2.5.2 / pkgload 1.5.3), producing a
# "must be used" warning even though the argument is valid. recompile is
# left at its default (FALSE): load_all() already detects changed .cpp/.h
# files via timestamp and recompiles only what's needed, so repeated runs
# of this script without source changes stay fast. If a .cpp edit doesn't
# seem to take effect, rerun once with recompile = TRUE to force a full rebuild.
pkgload::load_all("../PoissonLTDM", debug = FALSE)

source("plot_diagnostics.R")

printf <- function(...) cat(paste(sprintf(...), "\n"))

method_grid <- c("amh_montoril", "pg_as", "sir_laplace", "sir_collapsed")
Tt_grid <- c(200, 400, 800, 1600)
function_grid <- c("constant", "linear", "quadratic", "sinusoidal")

# ---- Choose method and data data here ----
method <- method_grid[1]
f <- function_grid[1]
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
N <- 80000          # number of iterations
burnin <- 5000
verbose <- FALSE
print_every <- 1000
plots <- TRUE
compute_rhat <- TRUE
compute_ess <- TRUE


# Adaptive Metropolis hyperparameters (for amh_montoril)
ac_ref <- 0.44             # acceptance ratio target

# Particle Gibbs (pg_as) hyperparameter
K <- 30  # Number of particles

# SIR Laplace and SIR Collapsed  hyperparameters
M_is <- 3             # Number of particles - IS for W1 integrated likelihood
M_irls_max <- 20
tol <- 1e-4

# Only for SIR Collapsed
R_prerun <- 3000      # pre-run iterations to calibrate CE proposals (phi1 and phi2)
M_sir <- 3            # Number of particles - SIR of theta1

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

printf("Running %s (Rcpp) for %s, seed_base=%d, N=%d, burnin=%d", 
       method, source_name, seed_base, N, burnin)


# ---- Chain initialization ----
#
# Reference values derived from the data (method-of-moments style), used
# only as the CENTER around which each chain's starting point is dispersed.
# Not used directly as an initial value in any chain.
theta1_ref <- log(y + 0.5)
theta2_ref <- c(diff(theta1_ref), 0)
theta_01_ref <- theta1_ref[1]
theta_02_ref <- theta2_ref[1]
W1_ref <- var(diff(theta1_ref))
W2_ref <- var(diff(theta2_ref))

# For amh_montoril and pg_as, theta1 is part of the Markov chain state
# itself (MH / Particle Gibbs reference trajectory) -- its initial value
# genuinely matters and must be dispersed across chains for R_hat to be
# informative. For sir_laplace and sir_collapsed, theta1 is redrawn via
# SIR/IS every iteration, so zero-initialization is sufficient; only the
# hyperparameters (theta_01, theta_02, W1, W2) need dispersion there.
theta1_dispersed <- method %in% c("amh_montoril", "pg_as")

make_chain_inits <- function(chain_id, seed) {
    set.seed(seed)

    theta_01_init <- theta_01_ref + rnorm(1, 0, sd = 2 * sqrt(abs(theta_01_ref) + 1))
    theta_02_init <- theta_02_ref + rnorm(1, 0, sd = 2 * sqrt(abs(theta_02_ref) + 1))
    W1_init <- W1_ref * exp(rnorm(1, 0, sd = 1))
    W2_init <- W2_ref * exp(rnorm(1, 0, sd = 1))

    if (theta1_dispersed) {
        theta1_init <- theta1_ref + rnorm(Tt, 0, sd = 2 * sqrt(W1_ref))
    } else {
        theta1_init <- numeric(Tt)
    }
    theta2_init <- numeric(Tt) # Chan (exact) in all four samplers -- dispersion
                                # is inherited from theta_02_init/W2_init

    list(chain_id = chain_id, seed = seed,
         theta_01 = theta_01_init, theta_02 = theta_02_init,
         W1 = W1_init, W2 = W2_init,
         theta1 = theta1_init, theta2 = theta2_init,
         theta1_tilde = rep(0, Tt),    # sir_laplace / sir_collapsed: IRLS
                                       # warm-start only, no dispersion needed
         varsigma2 = rep(0.02, Tt))    # amh_montoril: adaptation hyperparameter,
                                       # not a model state -- identical across chains
}

chain_inits <- lapply(1:N_chains, function(k) make_chain_inits(k, seed_base + (k - 1)))


# amh_montoril_cpp / pg_as_cpp / sir_laplace_cpp / sir_collapsed_cpp are
# already available here -- compiled and exported by the load_all() call
# above (PoissonLTDM/src/*.cpp), no separate sourceCpp() needed.

run_one_chain <- function(init) {
    printf("Chain %d/%d (seed=%d)", init$chain_id, N_chains, init$seed)
    set.seed(init$seed)

    if (method == "pg_as") {
        pg_as_cpp(y = y, K = K, N = N,
                  mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                  nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                  theta1 = init$theta1, theta2 = init$theta2,
                  theta_01 = init$theta_01, theta_02 = init$theta_02,
                  W1 = init$W1, W2 = init$W2, verbose = verbose, print_every = print_every)
    } else if (method == "amh_montoril") {
        amh_montoril_cpp(y = y, N = N,
                          mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                          nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                          theta1 = init$theta1, theta2 = init$theta2,
                          theta_01 = init$theta_01, theta_02 = init$theta_02,
                          W1 = init$W1, W2 = init$W2, ac_ref = ac_ref, varsigma2 = init$varsigma2,
                          verbose = verbose, print_every = print_every)
    } else if (method == "sir_laplace") {
        sir_laplace_cpp(y = y, N = N,
                         mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                         nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                         theta1 = init$theta1, theta2 = init$theta2,
                         theta_01 = init$theta_01, theta_02 = init$theta_02,
                         W1 = init$W1, W2 = init$W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol,
                         verbose = verbose, print_every = print_every)
    } else if (method == "sir_collapsed") {
        sir_collapsed_cpp(y = y, R_prerun = R_prerun, N = N,
                           mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                           nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                           theta1 = init$theta1, theta2 = init$theta2, theta1_tilde = init$theta1_tilde,
                           theta_01 = init$theta_01, theta_02 = init$theta_02, W1 = init$W1, W2 = init$W2,
                           M_is_lik = M_is, M_sir_theta1 = M_sir, M_irls_max = M_irls_max, tol = tol,
                           verbose = verbose, print_every = print_every)
    } else {
        stop(sprintf("Unknown method: '%s'.", method))
    }
}

start_time <- proc.time()
results <- lapply(chain_inits, run_one_chain)
tf <- proc.time()
elapsed_time <- (tf - start_time)[[1]]
printf("Total CPU time (%d chains): %.2f s", N_chains, elapsed_time)

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
                            nu_02 = nu_02, eta_02 = eta_02, ac_ref=ac_ref, plots = plots,
                            compute_rhat = compute_rhat, compute_ess = compute_ess)
