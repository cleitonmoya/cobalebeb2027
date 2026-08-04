# tests/test_validation.R
#
# Runs the chosen method's R prototype AND its Rcpp port with IDENTICAL
# inputs (same seed, same data, same hyperparameters/initial values) and
# checks that they agree numerically. This is the formal regression test:
# every method here was hand-validated this way during development
# (R_prototypes/*.R match their PoissonLTDM/src/*.cpp counterparts to
# ~1e-14, once both use the same natural-order resampling method --
# see R_prototypes/utils.R's sample_one_from_logw()/systematic_resample()
# for why that matters).
#
# N/R_prerun below are the same production-scale values used in
# test_prototype_R.R/test_cpp.R -- correctness doesn't need that many
# iterations to show up, so reduce N (and R_prerun, for sir_collapsed)
# here if a faster check is preferred.
#
# NOTE on large-N/large-T validation failures (R prototype vs Rcpp port)
# 
# For large (T, N) -- e.g. T = 1600, N = 10000 -- this test can FAIL even
# though R and C++ are algorithmically equivalent. Confirmed root cause
# (pg_as, constant_1600_1, seed = 214001): both chains match to fp noise
# (~1e-8) for thousands of iterations, then diverge abruptly at a single
# iteration (n = 8699, t = 401, ancestor k = 46), where the resampling
# draw landed only 9.35e-10 from the cumulative-weight boundary. R's
# vectorized rnorm()/exp()/cumsum() vs C++'s std::exp() + sequential sum
# are mathematically equivalent but not bit-identical (IEEE 754 is
# non-associative); the ~1-ULP-per-call error compounds over T = 1600
# sequential steps within one Gibbs iteration until it crosses a
# resampling threshold. Once one ancestor index flips, R and C++ track
# different but EQUALLY VALID realizations of the same correct PG-AS
# transition, and the divergence propagates forward permanently.
#
# With T*K*N ~ 8e8 threshold checks, such a crossing is near-certain
# eventually -- a matter of when, not if. Smaller T/N (T <= 800,
# N <= 5000) pass to machine precision and remain a valid correctness
# check. For large-scale runs, compare posterior summaries (means, ESS)
# of two independent chains instead of requiring sample-path identity.
# ----------------------------------------------------------------------

library(Rcpp)
library(testthat) # provide expect_equal

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

printf <- function(...) cat(paste(sprintf(...), "\n"))

method_grid <- c("amh_montoril", "pg_as", "sir_laplace", "sir_collapsed")
Tt_grid <- c(200, 400, 800, 1600)
function_grid <- c("constant", "linear", "quadratic", "sinusoidal")

# ---- Choose method and data data here ----
method <- method_grid[2]
f <- function_grid[1]
Tt <- Tt_grid[4] 
replica <- 1     # one of: 1,...,200

tolerance <- 1e-10      # both implementations solve the same tridiagonal
                        # system via different code paths (sparse Cholesky
                        # in R vs a direct LDL^T solve in C++), so expect
                        # floating-point-noise-level differences, not exact
                        # bit-identity
# ------------------------------------------

# Load the data
source_name <- sprintf("%s_%s_%s", f, Tt, replica)
data <- readRDS(paste("../data/simulated/", source_name, ".rds", sep = ""))
y <- data$y
theta1_true <- data$theta

# Compute the seed (based on pattern)
method_idx <- match(method, method_grid)
Tt_idx <- match(Tt, Tt_grid)
f_idx <- match(f, function_grid)
seed <- method_idx*1e5 + f_idx*1e4 + Tt_idx*1e3 + replica

printf("Validating %s: R prototype vs Rcpp, source=%s, seed=%d", method, source_name, seed)

# ---- Parameters and initialization ----
# General simulation parameters
N <- 10000            # number of iterations
burnin <- 1000
verbose <- TRUE
print_every <- 1000

# Adaptive Metropolis hyperparameters (for montoril)
ac_ref <- 0.44        # acceptance ratio target

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
theta1 <- numeric(Tt)
theta2 <- numeric(Tt)
theta1_tilde <- rep(0, Tt)    # sir_laplace and sir_collapsed
varsigma2 <- rep(0.02, Tt)     # amh_montoril - initial varsigma2

# Source both implementations
source(file.path("..", "R_prototypes", paste0(method, ".R")))
sourceCpp(file.path("..", "PoissonLTDM", "src", paste0(method, ".cpp")))

# Build the argument list once -- _r() and _cpp() have identical signatures,
# so the same list can be do.call()'d against both.
if (method == "pg_as") {
    args <- list(y = y, K = K, N = N,
                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                 theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                 W1 = W1, W2 = W2, verbose = FALSE)
} else if (method == "amh_montoril") {
    args <- list(y = y, N = N,
                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                 theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                 W1 = W1, W2 = W2, ac_ref = ac_ref, varsigma2 = varsigma2, verbose = FALSE)
} else if (method == "sir_laplace") {
    args <- list(y = y, N = N,
                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                 theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                 W1 = W1, W2 = W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol, verbose = FALSE)
} else if (method == "sir_collapsed") {
    args <- list(y = y, R_prerun = R_prerun, N = N,
                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                 theta1 = theta1, theta2 = theta2, theta1_tilde = theta1_tilde,
                 theta_01 = theta_01, theta_02 = theta_02, W1 = W1, W2 = W2,
                 M_is_lik = M_is, M_sir_theta1 = M_sir, M_irls_max = M_irls_max, tol = tol, verbose = FALSE)
} else {
    stop(sprintf("Unknown method: '%s'.", method))
}

fn_r   <- get(paste0(method, "_r"))
fn_cpp <- get(paste0(method, "_cpp"))

set.seed(seed)
time_r <- system.time(result_r <- do.call(fn_r, args))

set.seed(seed)
time_cpp <- system.time(result_cpp <- do.call(fn_cpp, args))

printf("R prototype time:  %.3f s", time_r[["user.self"]])
printf("Rcpp time:         %.3f s", time_cpp[["user.self"]])
printf("Speedup (R / Cpp): %.1fx", time_r[["user.self"]] / time_cpp[["user.self"]])

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
    printf("ALL CHECKS PASSED: %s R prototype matches the Rcpp port (tolerance = %.0e).", method, tolerance)
} else {
    printf("SOME CHECKS FAILED: %s R prototype and Rcpp port disagree beyond tolerance -- investigate.", method)
}
