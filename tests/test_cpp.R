# tests/test_cpp.R
#
# Runs the chosen algorithm using the Rcpp-compiled implementation
# (PoissonLTDM/src/), on the real simulated dataset, printing and plotting
# the same diagnostic battery as test_prototype_R.R -- same data, same
# hyperparameters/initial values, so the two are directly comparable by
# eye. For a numeric R-vs-C++ comparison, see test_validation.R instead.

library(Rcpp)
library(coda)

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

source("plot_diagnostics.R")

# ---- Choose the algorithm here ----
algorithm <- "pg_as" # one of: "pg_as", "amh_montoril", "sir_laplace", "sir_collapsed"

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
set.seed(seed)
printf("Running %s (Rcpp) for %s, seed=%d", algorithm, source_name, seed)
if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 400) t_obs <- c(75, 100, 200, 300)
if (Tt == 800) t_obs <- c(100, 300, 500, 700)
if (Tt == 1600) t_obs <- c(400, 800, 1200, 1600)
theta1_true <- data$theta

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

# Compile the chosen algorithm's C++ source (utils.h is picked up locally
# via #include, so it must sit alongside the .cpp files in PoissonLTDM/src/)
sourceCpp(file.path("..", "PoissonLTDM", "src", paste0(algorithm, ".cpp")))

start_time <- proc.time()
if (algorithm == "pg_as") {
    result <- pg_as_cpp(y = y, K = K, N = N,
                         mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                         nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                         theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                         W1 = W1, W2 = W2,
                         verbose = TRUE, print_every = 1000)
} else if (algorithm == "amh_montoril") {
    result <- amh_montoril_cpp(y = y, N = N,
                                mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                                nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                                theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                                W1 = W1, W2 = W2, ac_ref = ac_ref, varsigma2 = varsigma2,
                                verbose = TRUE, print_every = 1000)
} else if (algorithm == "sir_laplace") {
    result <- sir_laplace_cpp(y = y, N = N,
                               mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                               nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                               theta1 = theta1, theta2 = theta2, theta_01 = theta_01, theta_02 = theta_02,
                               W1 = W1, W2 = W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol,
                               verbose = TRUE, print_every = 1000)
} else if (algorithm == "sir_collapsed") {
    result <- sir_collapsed_cpp(y = y, R_prerun = R_prerun, N = N,
                                 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                                 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                                 theta1 = theta1, theta2 = theta2, theta1_tilde = theta1_tilde,
                                 theta_01 = theta_01, theta_02 = theta_02, W1 = W1, W2 = W2,
                                 M_is_lik = M_is, M_sir_theta1 = M_sir, M_irls_max = M_irls_max, tol = tol,
                                 verbose = TRUE, print_every = 1000)
} else {
    stop(sprintf("Unknown algorithm: '%s'.", algorithm))
}
elapsed_time <- (proc.time() - start_time)[[1]]
printf("Total CPU time: %.2f s", elapsed_time)

print_and_plot_diagnostics(result, y,
                            theta1_true = theta1_true, theta2_true = NULL,
                            t_obs = t_obs, burnin = burnin, elapsed_time = elapsed_time,
                            nu_01 = nu_01, eta_01 = eta_01,
                            nu_02 = nu_02, eta_02 = eta_02)
