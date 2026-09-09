# tests/test_cpp.R
#
# Runs the chosen method using the Rcpp-compiled implementation
# (PoissonLTDM/src/), on the real simulated dataset, printing and plotting
# the same diagnostic battery as test_prototype_R.R -- same data, same
# hyperparameters/initial values, so the two are directly comparable by
# eye. For a numeric R-vs-C++ comparison, see test_R_vs_cpp.R instead.
#
# Runs N_chains independent chains with dispersed initializations (see
# make_chain_inits() below) and passes all of them to
# print_and_plot_diagnostics(), which computes R_hat/ESS across chains and
# plots only the chain selected via `plot_chain`.
#
# Chains run sequentially or in parallel depending on `parallel_chains`
# (see "Parallel execution" below); the parallel path uses foreach +
# doParallel + doRNG (same packages already used by simulation.R).

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
method <- method_grid[4]
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

# Full method grid used ONLY for the seed formula below -- index must
# match simulation_run.R's methods_grid_ref / calibration_run.R's
# method_grid exactly, including "stan" (which this script doesn't run),
# so amh_montoril/pg_as/sir_laplace/sir_collapsed keep the SAME index
# here as in production regardless of which subset this script exercises.
method_grid_seed_ref <- c("amh_montoril", "pg_as", "sir_laplace", "sir_collapsed", "stan")

# Compute the seed -- SAME pattern as simulation_run.R's production
# task_grid$seed formula:
#   seed_base(m, g, tau, r) = m*1e7 + g*1e6 + tau*1e5 + 10000 + r*10
# where m/g/tau are 1-based indices into method_grid_seed_ref/
# function_grid/Tt_grid. seed_base here IS the seed of that (method, f,
# Tt, replica)'s production chain 1; chain k below reuses the existing
# seed_base + (k - 1) convention (same as calibration_run.R) to reproduce
# every one of that task's production chains exactly (see simulation_run.R's
# header comment for the full derivation and the collision-avoidance
# rationale against calibration's seeds).
#
# calibration_run.R uses a related but distinct formula, replacing the
# fixed "10000" term with config_idx*1e3 (config_idx = the 1-based row
# index of that method's candidate N/burnin/K configuration in
# calibration_run.R's method_configs table -- not applicable outside
# calibration):
#   seed_base_calib(m, g, tau, c, r) = m*1e7 + g*1e6 + tau*1e5 + c*1e3 + r*10
# To instead reproduce a specific calibration chain's seed here, replace
# "+ 10000" below with "+ config_idx * 1e3", using that config's
# config_idx from calibration_run.R's method_configs.
method_idx <- match(method, method_grid_seed_ref)
Tt_idx <- match(Tt, Tt_grid)
f_idx <- match(f, function_grid)
seed_base <- method_idx*1e7 + f_idx*1e6 + Tt_idx*1e5 + 10000 + replica*10

# ---- Parameters and initialization ----
# General simulation parameters
N_chains <- 3        # number of chains
plot_chain <- 1      # which chain (1..N_chains) to plot / trace in detail
N <- 11000          # number of iterations
burnin <- 1000
verbose <- FALSE
print_every <- 1000
plots <- TRUE
compute_rhat <- TRUE
compute_ess <- TRUE

# Parallel execution of the N_chains chains.
#
# parallel_chains = TRUE runs the chains concurrently via foreach %dorng%
# (doParallel backend, doRNG for per-chain reproducible RNG streams --
# same pattern used in simulation.R). FALSE keeps the original sequential
# lapply(), useful when verbose = TRUE (interleaved worker output is hard
# to read) or when debugging a single chain.
#
# n_cores = NULL (default) picks the number of PHYSICAL cores automatically
# (see resolve_n_cores() below, which uses `lscpu` on Linux -- more
# reliable than parallel::detectCores(logical = FALSE) on this hardware).
# Hyperthreaded logical cores rarely help much for this kind of dense,
# branch-heavy MCMC/SMC code, so physical-core count is a better default
# than all logical cores. Set n_cores explicitly (e.g. n_cores = N_chains)
# to override.
parallel_chains <- TRUE
n_cores <- NULL


# Adaptive Metropolis hyperparameters (for amh_montoril)
ac_ref <- 0.44             # acceptance ratio target

# Particle Gibbs (pg_as) hyperparameter
K <- 100  # Number of particles

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

# ---- Physical core count (used when n_cores is not set explicitly) ----
#
# parallel::detectCores(logical = FALSE) is NOT reliable on this hardware:
# on the user's i7-4510U (2 physical cores, hyperthreaded to 4 logical), it
# returned 4 -- same as detectCores(logical = TRUE) -- silently counting
# hyperthreads as physical cores instead of falling back or erroring.
# `lscpu -p=CORE,SOCKET` reads the same /proc or /sys topology info but
# parses it correctly here (verified against `lscpu`'s own "Thread(s) per
# core" / "CPU(s)" summary), so it's used as the primary source on Linux;
# detectCores(logical = FALSE), then floor(logical/2), are fallbacks only
# for non-Linux systems or if lscpu isn't installed (e.g. some cluster
# images).
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

# NOTE: elapsed_time below uses wall-clock time (Sys.time()), not
# proc.time()'s user.self field. proc.time() only measures CPU time of the
# CURRENT process; under parallel_chains = TRUE the chains run in worker
# processes, so proc.time() on the main process would not capture their
# work and would badly understate elapsed_time -- which feeds directly
# into the ESS/second diagnostics in print_and_plot_diagnostics(). Wall
# clock is correct in both the sequential and parallel cases.
start_time <- Sys.time()

if (parallel_chains) {
    n_cores_used <- resolve_n_cores(n_cores, N_chains)
    printf("Running %d chains in parallel on %d cores", N_chains, n_cores_used)

    cl <- parallel::makeCluster(n_cores_used)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    `%dorng%` <- doRNG::`%dorng%`

    # Each worker needs the compiled samplers and the globals run_one_chain()
    # closes over (method, y, N, hyperparameters, etc.). .export handles the
    # globals; load_all() (debug=FALSE, see note at top of file) puts
    # amh_montoril_cpp/pg_as_cpp/sir_laplace_cpp/sir_collapsed_cpp in scope
    # on each worker.
    current_wd <- getwd()
    parallel::clusterExport(cl, "current_wd")
    parallel::clusterEvalQ(cl, {
        setwd(current_wd)
        pkgload::load_all("../PoissonLTDM", debug = FALSE)
    })

    # foreach auto-detects and exports the globals referenced inside the
    # loop body (run_one_chain and everything it closes over: method, y, N,
    # hyperparameters, printf, etc.) -- no need to list them manually.
    results <- foreach::foreach(init = chain_inits) %dorng% {
        run_one_chain(init)
    }

    #parallel::stopCluster(cl)
} else {
    results <- lapply(chain_inits, run_one_chain)
}

elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
printf("Total wall-clock time (%d chains, %s): %.2f s",
       N_chains, if (parallel_chains) sprintf("parallel, %d cores", n_cores_used) else "sequential",
       elapsed_time)

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
