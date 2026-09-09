# tests/test_prototype_R.R
#
# Runs the chosen method using the pure-R prototype implementation
# (R_prototypes/), on the simulated dataset, printing and plotting the
# diagnostic battery (see plot_diagnostics.R). This is the R-only
# counterpart to test_cpp.R -- same data, same diagnostics, so the two are
# directly comparable by eye. For a numeric R-vs-C++ comparison, see
# test_R_vs_cpp.R instead.

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

source("plot_diagnostics.R")

printf <- function(...) cat(paste(sprintf(...), "\n"))

method_grid <- c("amh_montoril", "pg_as", "sir_laplace", "sir_collapsed")
Tt_grid <- c(200, 400, 800, 1600)
function_grid <- c("constant", "linear", "quadratic", "sinusoidal")

# ---- Choose method and data data here ----
method <- method_grid[3]
f <- function_grid[1]
Tt <- Tt_grid[4] 
replica <- 1     # one of: 1,...,200     
# ------------------------------------------

# Load the data
source_name <- sprintf("%s_%s_%s", f, Tt, replica)
data <- readRDS(paste("../data/simulated/", source_name, ".rds", sep = ""))
y <- data$y
theta1_true <- data$theta

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
# function_grid/Tt_grid. This reproduces the EXACT seed of that
# (method, f, Tt, replica)'s production chain 1 (see simulation_run.R's
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

seed <- method_idx*1e7 + f_idx*1e6 + Tt_idx*1e5 + 10000 + replica*10
set.seed(seed)
printf("Running %s (R prototype) for %s, seed=%d", method, source_name, seed)


# ---- Parameters and initialization ----
# General simulation parameters
N <- 10000          # number of iterations
burnin <- 1000
verbose <- TRUE
print_every <- 1000
plots <- FALSE

# Adaptive Metropolis hyperparameters (for montoril)
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
theta1 <- numeric(Tt)
theta2 <- numeric(Tt)
theta1_tilde <- rep(0, Tt)    # sir_laplace and sir_collapsed
varsigma2 <- rep(0.02, Tt)     # amh_montoril - initial varsigma2
  
# Source the R prototype for the chosen method
source(file.path("..", "R_prototypes", paste0(method, ".R")))

start_time <- proc.time()

if (method == "pg_as") {
    result <- pg_as_r(y = y, K = K, N = N,
                       mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                       nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                       theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                       W1 = W1, W2 = W2,
                       verbose = verbose, print_every = print_every)
} else if (method == "amh_montoril") {
    result <- amh_montoril_r(y = y, N = N,
                              mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                              nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                              theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                              W1 = W1, W2 = W2, ac_ref = ac_ref, varsigma2 = varsigma2,
                              verbose = verbose, print_every = print_every)
} else if (method == "sir_laplace") {
    result <- sir_laplace_r(y = y, N = N,
                             mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                             nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                             theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                             W1 = W1, W2 = W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol,
                             verbose = verbose, print_every = print_every)
} else if (method == "sir_collapsed") {
    result <- sir_collapsed_r(y = y, R_prerun = R_prerun, N = N,
                               mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                               nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                               theta1 = theta1, theta2 = theta2, theta1_tilde = theta1_tilde,
                               theta_01 = theta_01, theta_02 = theta_02, W1 = W1, W2 = W2,
                               M_is_lik = M_is, M_sir_theta1 = M_sir, M_irls_max = M_irls_max, tol = tol,
                               verbose = verbose, print_every = print_every)
} else {
    stop(sprintf("Unknown method: '%s'.", method))
}

elapsed_time <- (proc.time() - start_time)[[1]]
printf("Total CPU time: %.2f s", elapsed_time)


if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 400) t_obs <- c(75, 100, 200, 300)
if (Tt == 800) t_obs <- c(100, 300, 500, 700)
if (Tt == 1600) t_obs <- c(400, 800, 1200, 1600)

print_and_plot_diagnostics(result, y,
                            theta1_true = theta1_true, theta2_true = NULL,
                            t_obs = t_obs, burnin = burnin, elapsed_time = elapsed_time,
                            nu_01 = nu_01, eta_01 = eta_01,
                            nu_02 = nu_02, eta_02 = eta_02,
                            ac_ref = ac_ref, plots=plots)
