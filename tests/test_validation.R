# tests/test_validation.R
#
# Runs the chosen algorithm's R prototype AND its Rcpp port with IDENTICAL
# inputs (same seed, same data, same hyperparameters/initial values) and
# checks that they agree numerically. This is the formal regression test:
# every algorithm here was hand-validated this way during development
# (R_prototypes/*.R match their PoissonLTDM/src/*.cpp counterparts to
# ~1e-14, once both use the same natural-order resampling algorithm --
# see R_prototypes/utils.R's sample_one_from_logw()/systematic_resample()
# for why that matters).
#
# N/R_prerun below are the same production-scale values used in
# test_prototype_R.R/test_cpp.R -- correctness doesn't need that many
# iterations to show up, so reduce N (and R_prerun, for sir_collapsed)
# here if a faster check is preferred.

library(Rcpp)
suppressMessages(library(testthat))

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

# ---- Choose the algorithm here ----
algorithm <- "pg_as" # one of: "pg_as", "amh_montoril", "sir_laplace", "sir_collapsed"
tolerance <- 1e-10    # both implementations solve the same tridiagonal
                       # system via different code paths (sparse Cholesky
                       # in R vs a direct LDL^T solve in C++), so expect
                       # floating-point-noise-level differences, not exact
                       # bit-identity

printf <- function(...) cat(paste(sprintf(...), "\n"))

# Load the data
functions_grid <- c("constant", "linear", "quadratic", "sinusoidal")
f <- 3
Tt <- 1600
replica <- 1
source_name <- sprintf("%s_%s_%s", functions_grid[f], Tt, replica)
data <- readRDS(paste("../data/simulated/", source_name, ".rds", sep = ""))
y <- data$y
Tt_grid <- c(200, 400, 800, 1600)
method <- 2
tau <- match(Tt, Tt_grid)
seed <- method * 1e5 + f * 1e4 + tau * 1e3 + replica

printf("Validating %s: R prototype vs Rcpp, seed=%d", algorithm, seed)

# ---- Parameters and initialization ----
# General simulation parameters
N <- 10000          # number of iterations
burnin <- 1000
# Adaptive Metropolis hyperparameters
# (for montoril)
varsigma2_scal <- 0.02     # initial varsigma2
ac_ref <- 0.44             # acceptance ratio target
# Particle Gibbs (sir_pg) hyperparameter
K <- 50  # Number of particles
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
# Initialization
W2 <- 0.01
W1 <- 0.01
theta_01 <- 0
theta_02 <- 0
theta1_tilde_scal <- 0 # sir_laplace and sir_collapsed

theta1 <- numeric(Tt)
theta2 <- numeric(Tt)
theta1_tilde <- rep(theta1_tilde_scal, Tt)
varsigma2 <- rep(varsigma2_scal, Tt)

# Source both implementations
source(file.path("..", "R_prototypes", paste0(algorithm, ".R")))
sourceCpp(file.path("..", "PoissonLTDM", "src", paste0(algorithm, ".cpp")))

# Build the argument list once -- _r() and _cpp() have identical signatures,
# so the same list can be do.call()'d against both.
if (algorithm == "pg_as") {
    args <- list(y = y, K = K, N = N,
                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                 theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                 W1 = W1, W2 = W2, verbose = FALSE)
} else if (algorithm == "amh_montoril") {
    args <- list(y = y, N = N,
                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                 theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                 W1 = W1, W2 = W2, ac_ref = ac_ref, varsigma2 = varsigma2, verbose = FALSE)
} else if (algorithm == "sir_laplace") {
    args <- list(y = y, N = N,
                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                 theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                 W1 = W1, W2 = W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol, verbose = FALSE)
} else if (algorithm == "sir_collapsed") {
    args <- list(y = y, R_prerun = R_prerun, N = N,
                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                 theta1 = theta1, theta2 = theta2, theta1_tilde = theta1_tilde,
                 theta_01 = theta_01, theta_02 = theta_02, W1 = W1, W2 = W2,
                 M_is_lik = M_is, M_sir_theta1 = M_sir, M_irls_max = M_irls_max, tol = tol, verbose = FALSE)
} else {
    stop(sprintf("Unknown algorithm: '%s'.", algorithm))
}

fn_r   <- get(paste0(algorithm, "_r"))
fn_cpp <- get(paste0(algorithm, "_cpp"))

set.seed(seed)
result_r <- do.call(fn_r, args)

set.seed(seed)
result_cpp <- do.call(fn_cpp, args)

# ---- Compare every history/summary field common to both results ----
common_fields <- intersect(names(result_r), names(result_cpp))
printf("Comparing %d common field(s): %s", length(common_fields), paste(common_fields, collapse = ", "))

all_passed <- TRUE
for (field in common_fields) {
    a <- result_r[[field]]
    b <- result_cpp[[field]]
    if (is.logical(a)) a <- as.numeric(a)
    if (is.logical(b)) b <- as.numeric(b)

    max_diff <- max(abs(a - b))
    test_result <- tryCatch({
        expect_equal(a, b, tolerance = tolerance)
        "PASS"
    }, error = function(e) {
        all_passed <<- FALSE
        "FAIL"
    })
    printf("  [%s] %-16s max abs diff = %.3e", test_result, field, max_diff)
}

printf("")
if (all_passed) {
    printf("ALL CHECKS PASSED: %s R prototype matches the Rcpp port (tolerance = %.0e).", algorithm, tolerance)
} else {
    printf("SOME CHECKS FAILED: %s R prototype and Rcpp port disagree beyond tolerance -- investigate.", algorithm)
}
