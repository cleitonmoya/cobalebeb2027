# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))


rm(list = ls())     # clear the environment
set.seed(42)

source("../PoissonLTDM/R/utils.R")
source("../PoissonLTDM/R/sampler_sir_collapsed.R")

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

# Load the data
data <- readRDS("../data/poisson_pol2_200.rds")
y <- data$y
Tt <- length(y)
theta1_true <- data$theta

# Simulation parameters
N <- 10000            # Gibbs iterations
burnin <- 1000
R_prerun   <- 3000    # pre-run iterations to calibrate CE proposals (phi1 and phi2)
M_is <- 3             # Number of particles - IS for W1 integrated likelihood
M_sir <- 3            # Number of particles - SIR of theta1
M_irls_max <- 20
tol <- 1e-4

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
theta1_tilde <- numeric(Tt)

res <- sample_sir_collapsed(
	        y            = y,
			N            = N,
			burnin       = burnin,
			R_prerun     = R_prerun,
			M_is         = M_is,
			M_sir        = M_sir,
			M_irls_max   = M_irls_max,
			tol          = tol, 
			mu_01        = mu_01,
			sigma2_01    = sigma2_01,
			mu_02        = mu_02,
			sigma2_02    = sigma2_02,
			nu_01        = nu_01,
			eta_01       = eta_01,
			nu_02        = nu_02,
			eta_02       = eta_02,
			W1           = W1,
			W2           = W2,
			theta_01     = theta_01,
			theta_02     = theta_02,
			theta1       = theta1,
			theta2       = theta2,
			theta1_tilde = theta1_tilde)


theta_01_hist  <- res$theta_01_hist
theta_02_hist  <- res$theta_02_hist
W1_hist        <- res$W1_hist
W2_hist        <- res$W2_hist
theta1_hist    <- res$theta1_hist
theta2_hist    <- res$theta2_hist
itr_irls_hist  <- res$itr_irls_hist
ess_is_hist    <- res$ess_is_hist
ess_sir_hist   <- res$ess_sir_hist
accepted1_hist <- res$accepted1_hist
accepted2_hist <- res$accepted2_hist

#####
theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)

printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))

loglik <- sum(dpois(y, lambda_mean, log=TRUE))
printf("Log-likelihood: %.2f", loglik)

# y, theta1_true, theta1_mean ####
x <- 1:Tt
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
ylim_range <- range(theta1_mean, theta1_true)
plot(x, theta1_mean, type="l", col="red", lwd=2, ylim=ylim_range,
	 xlab="t", ylab="", main="theta_t1")
lines(x, theta1_true, col="blue", lwd=2)
legend("topright", legend=expression(hat(theta)[t1], theta[t1]),
		   col=c("red","blue"), lwd=2, bty="n")

# theta2_true, theta2_mean
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(x, theta2_mean, type="l", col="red", lwd=2,
	 xlab="t", ylab="", main="theta_t2")
legend("topright", legend=expression(hat(theta)[t2]), col="red", lwd=2, bty="n")

# Traceplots for theta_t1 ####
t_obs <- c(50, 100, 150, 175)
par(mfrow = c(2, 2))
for (t in t_obs) {
	plot(theta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
	abline(v=burnin, col="red")
}

# Effective sample size for IS ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_is_hist, type="l", main="Effective Sample Size - IS W1")
abline(v=burnin, col="red")

# Effective sample size of theta1 SIR ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_sir_hist, type="l", main="Effective Sample Size - SIR theta1")
abline(v=burnin, col="red")