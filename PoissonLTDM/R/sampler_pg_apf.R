# Poisson - 2nd Order Polynomial Dynamic Model
# Particle Gibbs (PG) with Backward Sampling + Component-wise Gibbs for theta_t2
# Strategy: - Auxiliary Particle Filter (APF) for theta_t1 (scalar state);
#           - theta2 sampled via Preicison Matrix (Chan Method)
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269-342.
# Author: Cleiton Moya de Almeida



# Auxiliary functions
printf <- function(...) cat(paste(sprintf(...), "\n"))





# Prior Hyperparameters

# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- 0
sigma2_01 <- 100

# theta_02 ~ N(mu_02, sigma2_02)
mu_02     <- 0
sigma2_02 <- 100

# phi1 = 1/W1 ~ Gamma(nu_01, eta_01)
nu_01 <- 2
eta_01  <- 0.01

# phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
nu_02 <- 2
eta_02  <- 0.0001


# Simulation parameters
N <- 10000      # number of Gibbs iterations
K <- 50          # number of particles
burnin <- 2500

# Initial Values
theta1 <- numeric(Tt)
theta2 <- numeric(Tt)
theta_01 <- 0
theta_02 <- 0
W1 <- 0.01
W2 <- 0.01


# Auxiliary variables
W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)
ess_smc <- numeric(N)


#####
# FIXED SPARSE STRUCTURES FOR CHAN METHOD ####
start_time = proc.time() # execution time

# Base for the prior Precision Matrix K
sub_diag_base <- rep(-1, Tt-1)
main_diag_base <- c(rep(2, Tt-1), 1)
K0 <- bandSparse(n=Tt, k=c(0, -1),
                 diagonals=list(main_diag_base, sub_diag_base),
                 symmetric = TRUE)

# diagonal mask
# @x: slot of the Sparce matrix (S4 object) that contains the non-zero values
diag_pattern <- bandSparse(n=Tt, k=c(0, -1),
                           diagonals=list(rep(TRUE, Tt), rep(FALSE, Tt-1)),
                           symmetric=TRUE)
idx_diag <- which(diag_pattern@x) # index of subpattern@x which is non-zero

# subdiagonal mask
sub_pattern <- bandSparse(n=Tt, k=c(0, -1),
                          diagonals=list(rep(FALSE, Tt), rep(TRUE, Tt-1)),
                          symmetric=TRUE)
idx_sub <- which(sub_pattern@x)

# Initial symbolic Cholesky factor
Ch02_factor <- Cholesky(K0, perm = FALSE, LDL = TRUE)

# Work precision matrix (static)
P2_matrix <- K0

# Sparse matrix building time
time1 <- proc.time()
building_time <- (time1 - start_time)[[1]]
printf("Sparse structures building: %.4f s", building_time)



chan_sample_theta2 <- function(theta1, phi1, phi2, theta_02) {

    z <- diff(theta1)   # z_t = theta1[t+1] - theta1[t], t=1,...,T-1

    diag_obs <- c(rep(phi1, Tt-1), 0)
    P2_matrix@x[idx_diag] <- (main_diag_base*phi2) + diag_obs
    P2_matrix@x[idx_sub]  <- -phi2

    Ch2_factor <- update(Ch02_factor, P2_matrix)

    b <- numeric(Tt)
    b[1:(Tt-1)] <- z * phi1
    b[1] <- b[1] + theta_02 * phi2

    theta2_hat <- as.numeric(Matrix::solve(Ch2_factor, b, system="A"))
    d <- Matrix::diag(Ch2_factor)
    u <- rnorm(Tt)
    w <- u/sqrt(d)
    x <- as.vector(Matrix::solve(Ch2_factor, w, system="Lt"))

    return(theta2_hat + x)
}


# Gibbs sampling
start_time = proc.time() # execution time
for (n in 1:N) {

    if (n %% 1000 == 0) {
        time <- proc.time()
        elapsed_time <- (time - start_time)[[1]]
        printf("Iteration %d / %d | Elapsed CPU time: %.0f s", n, N, elapsed_time)
    }

    # Sample theta_02
    sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
    mu_02_bar <- sigma2_02_bar * (
        (theta1[1] - theta_01)/W1 + theta2[1]/W2 +
        mu_02/sigma2_02
    )
    theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))

    # Sample theta_01
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar *
        (mu_01 / sigma2_01 + (theta1[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

    # Sample W1
    dif1   <- theta1 - c(theta_01, theta1[-Tt])
    diffs1 <- dif1 - c(theta_02, theta2[-Tt])
    nu_01_bar <- nu_01 + Tt / 2
    eta_01_bar  <- eta_01 + 0.5 * sum(diffs1^2)
    phi1 <- rgamma(1, shape = nu_01_bar, rate = eta_01_bar)
    W1 <- 1 / phi1
    sd_W1 <- sqrt(W1)

    # Sample W2
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    nu_02_bar <- nu_02 + Tt / 2
    eta_02_bar  <- eta_02 + 0.5 * sum(diffs2^2)
    phi2 <- rgamma(1, shape = nu_02_bar, rate = eta_02_bar)
    W2 <- 1 / phi2
    sd_W2 <- sqrt(W2)


    #
    # Conditional SMC for theta_t1
    #
    theta_1_k   <- matrix(0, Tt, K)
    log_w_tilde <- matrix(0, Tt, K)

    # t = 1
    theta_1_k[1, ] <- rnorm(K, mean = theta_01 + theta_02, sd = sd_W1)
    theta_1_k[1, K] <- theta1[1]

    log_w_1 <- log_p_yt(y[1], theta_1_k[1, ])
    log_w_tilde[1, ] <- log_w_1 - logsumexp(log_w_1) # normalizing

    # t = 2, ..., T
    for (t in 2:Tt) {
        # Predictor (auxiliary variable)
        theta_hat_t1_k = theta_1_k[t - 1, ] + theta2[t - 1]

        # Auxiliary weights
        log_lambda_k = y[t] * theta_hat_t1_k - exp(theta_hat_t1_k)

        # First stage resampling
        log_aux <- log_w_tilde[t - 1, ] + log_lambda_k
        A <- sample(1:K, K, replace = TRUE, prob = exp(log_aux - max(log_aux)))
        A[K] <- K   # reference path

        # Propagation
        theta_t1_k <- rnorm(K, mean=theta_1_k[t-1, A] + theta2[t-1], sd=sd_W1)
        theta_t1_k[K] <- theta1[t]
        theta_1_k[t, ] <- theta_t1_k

        # Updated weights
        log_w_t <- log_p_yt(y[t], theta_1_k[t, ]) - log_lambda_k[A]
        log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)
    }

    ess_smc[n] <- 1 / sum(exp(2 * (log_w_tilde[Tt, ])))


    # backward sampling for theta1
    # Backward weights: w_t^k * N(theta1[t+1] | theta_t1^k + theta2[t], W1)
    k_final <- sample(1:K, 1, prob = exp(log_w_tilde[Tt, ] - max(log_w_tilde[Tt, ])))
    theta1[Tt] <- theta_1_k[Tt, k_final]

    for (t in (Tt - 1):1) {
        log_bw <- log_w_tilde[t, ] +
            dnorm(theta1[t+1],
                mean = theta_1_k[t, ] + theta2[t],
                sd = sd_W1,
                log = TRUE
            )
        log_bw <- log_bw - max(log_bw)
        bw <- exp(log_bw)
        bw <- bw / sum(bw)
        b  <- sample(1:K, 1, prob = bw)
        theta1[t] <- theta_1_k[t, b]
    }


    # Sample theta2 - Chan method
    theta2 <- chan_sample_theta2(theta1, phi1, phi2, theta_02)

    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1
    theta2_hist[n, ] <- theta2
}

#### Simulation summary

# Execution time
end_time <- proc.time()
sampling_time <- (end_time - time1)[[1]]
elapsed_time <- (end_time - start_time)[[1]]
printf("Sampling: %.2f s", sampling_time)
printf("Total CPU time: %.0f s", elapsed_time)

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
plot(W1_hist, type="l", xlab="n", ylab="W", main="Traceplot of W1")
abline(v=burnin, col="red")
plot(W2_hist, type="l", xlab="n", ylab="W", main="Traceplot of W2")
abline(v=burnin, col="red")


# Traceplots for theta_t1 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
    abline(v=burnin, col="red")
}


# Traceplots for theta_t2 ####
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
    abline(v=burnin, col="red")
}


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


# ACF for theta1 e theta2
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
for (t in t_obs) {
    acf(theta1_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*1]))
    acf(theta2_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*2]))
}


# Prior vs posterior for phi2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8) # bottom left, top, right
curve(dgamma(x, shape=nu_02, rate=eta_02), from=0, to=max(1/W2_hist[-(1:burnin)]),
      main="phi2 prior vs. posterior", col="red", lwd=2)
lines(density(1/W2_hist[-(1:burnin)]), col="blue", lwd=2)
legend("topright", legend=c("Prior","Posterior"), col=c("red","blue"), lwd=2)


# Effective sample size ####
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_smc, type="l", main="Effective Sample Size - SMC")
abline(v=burnin, col="red")
