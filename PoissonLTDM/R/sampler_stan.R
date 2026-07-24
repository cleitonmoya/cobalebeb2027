# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: NUTS (Stan)
# Author: Cleiton Moya de Almeida

graphics.off()    # close the plots
rm(list = ls())   # clear the environment
cat("\014")       # clear the console
options(error = function() traceback(2))
set.seed(42)

library(rstan)

# RStan configuration
options(mc.cores = 1)   
rstan_options(auto_write = FALSE)

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))


# Load or compile the model
if (file.exists("poisson_ltdm.rds")) {
  model <- readRDS("poisson_ltdm.rds")
  printf("Model loaded")
} else {
  printf("Building the model")
  model <- stan_model(file = "poisson_pol2.stan", model_name = "PoissonLTDM")
  saveRDS(model, file = "poisson_ltdm.rds")
}


# Mapping stan data
stan_data <- list(
    Tt        = Tt,
    y         = y,
    mu_01     = mu_01,
    sigma2_01 = sigma2_01,
    mu_02     = mu_02,
    sigma2_02 = sigma2_02,
    nu_01     = nu_01,
    eta_01    = eta_01,
    nu_02     = nu_02,
    eta_02    = eta_02
)

init_list <- list(
  list(
    theta_01 = theta_01,
    theta_02 = theta_02,
    theta1   = theta1,
    theta2   = theta2,
    phi1     = phi1,
    phi2     = phi2
  )
)

fit <- sampling(
    object = model,
    data   = stan_data,
    chains = 1,
    iter   = N,
    warmup = burnin,
    thin   = 1,
    seed   = 42,
    init   = init_list
)

warmup_time <- get_elapsed_time(fit)[1, 1]
sample_time <- get_elapsed_time(fit)[1, 2]
elapsed_time <- warmup_time + sample_time

# Extract the samples
samples <- extract(fit, permuted = FALSE, inc_warmup = TRUE)
theta_01_hist <- samples[, 1, "theta_01"]
theta_02_hist <- samples[, 1, "theta_02"]
W1_hist <- samples[, 1, "W1"]
W2_hist <- samples[, 1, "W2"]
theta1_hist <- samples[, 1, paste0("theta1[", 1:Tt, "]")]
theta2_hist <- samples[, 1, paste0("theta2[", 1:Tt, "]")]

sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
ac_hist <- sampler_params[[1]][, "accept_stat__"]

# Simulation summary ####
printf("Total elapsed CPU time: %.0f s", elapsed_time)
printf("Mean acception ratio of theta1: %.2f", mean(ac_hist))

# Posterior mean
theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)

printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))


# Log-likelihood
loglik <- sum(dpois(y, lambda_mean, log=TRUE))
printf("Log-likelihood: %.2f", loglik)


# Effective sample size
printf("Effective Sample Size:")
ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)

ess_theta1 <- effectiveSize(mcmc(theta1_hist[-(1:burnin),]))
ess_theta2 <- effectiveSize(mcmc(theta2_hist[-(1:burnin),]))
printf("\ttheta1 (mean): %.2f", mean(ess_theta1))
printf("\ttheta2 (mean): %.2f", mean(ess_theta2))

# Effective sample size per second
printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)

ess_sec_theta1 <- ess_theta1/elapsed_time
ess_sec_theta2 <- ess_theta2/elapsed_time
printf("\ttheta1 (mean): %.2f", mean(ess_sec_theta1))
printf("\ttheta2 (mean): %.2f", mean(ess_sec_theta2))


# Geweke diagnostic: Z test for two mean difference
#   H0: segments same means -> chain has converged
printf("Geweke convergence diagnostic")
z_w1 <- unname(geweke.diag(W1_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
z_w2 <- unname(geweke.diag(W2_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
printf("\tz_w1: %.2f", z_w1)
printf("\tz_w2: %.2f", z_w2)

# Percent of instants in the H_0 rejection region:
z_theta1 <- unname(geweke.diag(theta1_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z_theta2 <- unname(geweke.diag(theta2_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z1_out <- sum((z_theta1 < -1.96) | (z_theta1 > 1.96))/Tt
z2_out <- sum((z_theta2 < -1.96) | (z_theta2 > 1.96))/Tt
printf("\tPercent of theta1 out: %.3f", z1_out)
printf("\tPercent of theta2 out: %.3f", z2_out)


#####
# Plots
# y, lambda_true, lambda_estimated
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson Local Trend Polynomial Model")
points(x, y, pch = 20)
lines(x, lambda_mean, col="red", lwd=2)
if (theta1_present) {
  lines(x, lambda_true, col="blue", lwd=2)
  legend("topright",
         legend = expression(y[t], lambda[t], hat(lambda)[t]),
         col = c("black", "blue", "red"),
         lty = c(NA, 1, 1),
         lwd = c(NA, 2, 2),
         pch = c(20, NA, NA),
         bty = "n")
} else {
  legend("topright",
         legend = expression(y[t], hat(lambda)[t]),
         col = c("black", "red"),
         lty = c(NA, 1),
         lwd = c(NA, 2),
         pch = c(20, NA),
         bty = "n")
}


# y, theta1_true, theta1_mean ####
x <- 1:Tt
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
if (theta1_present) {
  ylim_range <- range(theta1_mean, theta1_true)
} else {
  ylim_range <- range(theta1_mean)
}
plot(x, theta1_mean, type="l", col="red", lwd=2, ylim=ylim_range,
     xlab="t", ylab="", main="theta_t1")
if (theta1_present) {
  lines(x, theta1_true, col="blue", lwd=2)
  legend("topright", legend=expression(hat(theta)[t1], theta[t1]),
         col=c("red","blue"), lwd=2, bty="n")
} else {
  legend("topright", legend=expression(hat(theta)[t1]),
         col="red", lwd=2, bty="n")
}


# y, theta2_true, theta2_mean ####
#y_range <- range(theta2_true, theta2_mean)
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
if (theta2_present) {
  ylim_range <- range(theta2_mean, theta2_true)
} else {
  ylim_range <- range(theta2_mean)
}
plot(x, theta2_mean, type="l", col="red", lwd=2, ylim=ylim_range,
     xlab="t", ylab="", main="theta_t2")
if (theta2_present) {
  lines(x, theta2_true, col="blue", lwd=2)
  legend("topright", legend=expression(hat(theta)[t2], theta[t2]),
         col=c("red","blue"), lwd=2, bty="n")
} else {
  legend("topright", legend=expression(hat(theta)[t2]),
         col="red", lwd=2, bty="n")
}


# Posterior distribution of theta_t1 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
  hist(theta1_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
       xlab = bquote(theta[.(t) * "," * 1]),
       main = bquote("Posterior of " * theta[.(t) * "," * 1]))
  lines(density(theta1_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of theta_t2 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
  hist(theta2_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
       xlab = bquote(theta[.(t) * "," * 2]),
       main = bquote("Posterior of " * theta[.(t) * "," * 2]))
  lines(density(theta2_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of W1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W1")
lines(density(W1_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W2")
lines(density(W2_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W1 and W2 ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W1_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W1")
plot(W2_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W2")


# Traceplots for theta_t1 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
  plot(theta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
}


# Traceplots for theta_t2 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
  plot(theta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
}


# Traceplot of acceptance ratio of theta_t1 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta[t*1]))


# Effective sample size ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ess_theta1, type="l", main=expression("Effective sample of " * theta[t1]), xlab="t")
plot(ess_theta2, type="l", main=expression("Effective sample of " * theta[t2]), xlab="t")

par(mfrow = c(1, 2), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
hist(ess_theta1)
hist(ess_theta2)


# Geweke diagnostic ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(z_theta1, type="l", main=expression("Geweke diagnostic for " * theta[t1]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")

plot(z_theta2, type="l", main=expression("Geweke diagnostic for " * theta[t2]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")


# ACF for theta1 e theta2 ####
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
for (t in t_obs) {
  acf(theta1_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*1]))
  acf(theta2_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*2]))
}


# Prior vs posterior for phi2 ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
curve(dgamma(x, shape=nu_02, rate=eta_02), from=0, to=max(1/W2_hist[-(1:burnin)]),
      main="phi2 prior vs. posterior", col="red", lwd=2)
lines(density(1/W2_hist[-(1:burnin)]), col="blue", lwd=2)
legend("topright", legend=c("Prior","Posterior"), col=c("red","blue"), lwd=2)

