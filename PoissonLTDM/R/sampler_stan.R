# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: NUTS (Stan)
# Author: Cleiton Moya de Almeida

graphics.off()    # close the plots
rm(list = ls())   # clear the environment
cat("\014")       # clear the console
set.seed(42)

library(rstan)

# Print auxiliary function
printf <- function(...) {
    x = paste(sprintf(...),"\n")
    return(cat(x))
}

# RStan configuration
options(mc.cores = parallel::detectCores())  # paralelizar chains
rstan_options(auto_write = TRUE)             # recompila só se o .stan mudar

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
theta2_present <- FALSE
#t_obs <- c(50, 100, 150, 200)
source <- "poisson_sin_200"
t_obs <- c(500, 1000, 1500, 2000)
source <- "poisson_sin_2000"

data <- readRDS(paste("../data/", source, ".rds", sep=""))
y <- data$y
Tt <- length(y)

if (theta2_present) {
    theta1_true <- data$theta1
    theta2_true <- data$theta2

} else {
    theta1_true <- data$theta
}
lambda_true <- exp(theta1_true)


#####
# Complile once and store the object
#model <- stan_model(file = "poisson_pol2.stan")
#saveRDS(model, file = "poisson_pol2.rds")

#####
model <- readRDS("poisson_pol2.rds")

stan_data <- list(
    T         = Tt,
    y         = y,
    mu_01     = y[1],       # theta_01 ~ N(mu_01, sigma2_01)
    sigma2_01 = 10,
    nu_01     = 0.01,       # phi_01   ~ Gamma(nu_01, eta_01)
    eta_01    = 0.01,
    mu_02     = 0.1,        # theta_02 ~ N(mu_02, sigma2_02)
    sigma2_02 = 1,
    nu_02     = 0.01,       # phi_01   ~ Gamma(nu_02, eta_02)
    eta_02    = 0.01
)

burnin <- 200

fit <- sampling(
    object = model,
    data   = stan_data,
    chains = 1,
    iter   = 5000,
    warmup = burnin,
    thin   = 1,
    seed   = 42
)
warmup_time <- get_elapsed_time(fit)[1, 1]
sample_time <- get_elapsed_time(fit)[1, 2]
elapsed_time <- warmup_time + sample_time


# Results ####
samples <- extract(fit, inc_warmup = TRUE)
sampler_params <- get_sampler_params(fit, inc_warmup = TRUE)
ac_hist <- sampler_params[[1]][, "accept_stat__"]
ac_mean <- mean(ac_hist)

W1_samples   <- samples$W1
W2_samples   <- samples$W2

theta1_samples <- samples$theta1   # matrix: n_iter x T
theta1_mean  <- apply(theta1_samples[-(1:burnin), ], 2, mean)

theta2_samples <- samples$theta2   # matrix: n_iter x T
theta2_mean  <- apply(theta2_samples[-(1:burnin), ], 2, mean)

lambda_samples <- samples$lambda_hat
lambda_mean  <- apply(lambda_samples[-(1:burnin), ], 2, mean)

W1_mean  <- mean(W1_samples)
W2_mean  <- mean(W2_samples)

printf("Execution time: %.0f s", elapsed_time)
printf("Mean acception ratio: %.2f", ac_mean)
print(fit, pars = c("W1", "W2",
                    "theta1[10]", "theta1[50]", "theta1[100]", "theta1[150]",
                    "theta2[10]", "theta2[50]", "theta2[100]", "theta2[150]"),
      probs = c(0.25, 0.5, 0.75))
printf("W1 mean: %.5f", W1_mean)
printf("W2 mean: %.5f", W2_mean)

# Effective sample size / elapsed time ####
s <- summary(fit)$summary
ess_bulk <- s[, "n_eff"]
ess_per_sec <- ess_bulk / elapsed_time

ess_per_sec_w1 <- ess_per_sec[["W1"]]

printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_per_sec[["W1"]])
printf("\tW2: %.2f", ess_per_sec[["W2"]])

for (t in t_obs) {
    printf("\ttheta1[%d]: %.2f", t, ess_per_sec[[sprintf("theta1[%d]", t)]])
}

for (t in t_obs) {
    printf("\ttheta2[%d]: %.2f", t, ess_per_sec[[sprintf("theta2[%d]", t)]])
}


# Plots #####
# y, lambda_true, lambda_mean
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson local level model")
points(x, y, pch = 20)
lines(x, lambda_mean, col="red", lwd=2)
lines(x, lambda_true, col="blue", lwd=2)
legend("topright",
       legend = expression(y[t], lambda[t], hat(lambda)[t]),
       col = c("black", "blue", "red"),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2),
       pch = c(20, NA, NA),
       bty = "n")


# theta1_true, theta1_mean ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(x, theta1_mean, type="l", col="blue", lwd="2")
lines(x, theta1_true)
legend("topright",
       legend = expression(theta[t1], hat(theta)[t1]),
       col = c("black", "blue"),
       lty = c(1, 1),
       lwd = c(1, 2),
       bty = "n")

# theta2_true, theta2_mean ####
    par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
    plot(x, theta2_mean, type="l", col="blue", lwd="2")
    #lines(x, theta2_true)
    legend("topright",
           legend = expression(theta[t2], hat(theta)[t2]),
           col = c("black", "blue"),
           lty = c(1, 1),
           lwd = c(1, 2),
           bty = "n")



# Posterior distribution of theta_t1 #####
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    hist(theta1_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[1][.(t)]), main = bquote("Posterior of " * theta[1][.(t)]))
    lines(density(theta1_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}

# Posterior distribution of theta_t2 #####
par(mfrow = c(2, 2))
for (t in c(10, 50, 100, 150)) {
    hist(theta2_samples[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[2][.(t)]), main = bquote("Posterior of " * theta[2][.(t)]))
    lines(density(theta2_samples[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of W1 #####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_samples[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = bquote(theta[1][.(t)]), main ="Posterior of W1")
lines(density(W1_samples[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2 #####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_samples[-(1:burnin)], breaks = 50, freq = FALSE,
     xlab = bquote(theta[2][.(t)]), main ="Posterior of W2")
lines(density(W2_samples[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W1 and W2 #####
plot(W1_samples, type="l", xlab="n", ylab="W1", main="Traceplot of W1")
plot(W2_samples, type="l", xlab="n", ylab="W2", main="Traceplot of W2")


# Traceplots for theta_t1 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta1_samples[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
}

# Traceplots for theta_t2 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta2_samples[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
}


# Traceplot of acceptance ratio ####
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
plot(ac_hist, type="l", xlab="n", ylab="ratio",
     main=expression("Acceptance ratio of " * theta))
