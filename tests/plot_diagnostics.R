# tests/plot_diagnostics.R
#
# Shared diagnostics/plotting routine for the Poisson 2nd-order polynomial
# DLM samplers. Used by both test_prototype_R.R (pure R) and test_cpp.R
# (Rcpp), so the same battery of plots/summary stats is produced regardless
# of which implementation was run -- only test_validation.R (which compares
# the two numerically) does not use this file.
#
# `result` must contain at least: theta1_hist, theta2_hist, theta_01_hist,
# theta_02_hist, W1_hist, W2_hist. Algorithm-specific extra diagnostics are
# added when the corresponding fields are present in `result` (ess_smc for
# pg_as; ac_hist for amh_montoril; ess_is/itr_irls for sir_laplace;
# accepted_hist/accepted2_hist/itr_irls/ess_is_hist/ess_sir_hist/ce1_*/ce2_*
# for sir_collapsed).

library(coda)

print_and_plot_diagnostics <- function(result, y, theta1_true = NULL, theta2_true = NULL,
                                        t_obs, burnin, elapsed_time,
                                        nu_01 = NULL, eta_01 = NULL,
                                        nu_02 = NULL, eta_02 = NULL,
                                        ac_ref = NULL) {

    printf <- function(...) cat(paste(sprintf(...), "\n"))
    
    
    theta1_hist <- result$theta1_hist
    theta2_hist <- result$theta2_hist
    theta_01_hist <- result$theta_01_hist
    theta_02_hist <- result$theta_02_hist
    W1_hist <- result$W1_hist
    W2_hist <- result$W2_hist

    N <- nrow(theta1_hist)
    Tt <- ncol(theta1_hist)
    theta1_present <- !is.null(theta1_true)
    theta2_present <- !is.null(theta2_true)

    
    # Summary stats
    theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
    theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
    lambda_mean <- exp(theta1_mean)

    printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
    printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))
    printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
    printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))

    loglik <- sum(dpois(y, lambda_mean, log = TRUE))
    printf("Log-likelihood: %.2f", loglik)

    ess_theta01 <- effectiveSize(mcmc(theta_01_hist[-(1:burnin)]))
    ess_theta02 <- effectiveSize(mcmc(theta_02_hist[-(1:burnin)]))
    ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
    ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
    ess_theta1 <- effectiveSize(mcmc(theta1_hist[-(1:burnin), ]))
    ess_theta2 <- effectiveSize(mcmc(theta2_hist[-(1:burnin), ]))
    printf("Effective Sample Size:")
    printf("\ttheta_01: %.2f", ess_theta01)
    printf("\ttheta_02: %.2f", ess_theta02)
    printf("\tW1: %.0f", ess_w1)
    printf("\tW2: %.0f", ess_w2)
    printf("\ttheta1 (mean): %.2f", mean(ess_theta1))
    printf("\ttheta2 (mean): %.2f", mean(ess_theta2))

    if (!missing(elapsed_time) && !is.null(elapsed_time)) {
        printf("Effective Sample Size / second:")
        printf("\tW1: %.2f", ess_w1 / elapsed_time)
        printf("\tW2: %.2f", ess_w2 / elapsed_time)
        printf("\ttheta1 (mean): %.2f", mean(ess_theta1 / elapsed_time))
        printf("\ttheta2 (mean): %.2f", mean(ess_theta2 / elapsed_time))
    }

    printf("Geweke convergence diagnostic")
    z_w1 <- unname(geweke.diag(W1_hist[-(1:burnin)], frac1 = 0.1, frac2 = 0.5)[[1]])
    z_w2 <- unname(geweke.diag(W2_hist[-(1:burnin)], frac1 = 0.1, frac2 = 0.5)[[1]])
    printf("\tz_w1: %.2f", z_w1)
    printf("\tz_w2: %.2f", z_w2)

    z_theta1 <- unname(geweke.diag(theta1_hist[-(1:burnin), ], frac1 = 0.1, frac2 = 0.5)[[1]])
    z_theta2 <- unname(geweke.diag(theta2_hist[-(1:burnin), ], frac1 = 0.1, frac2 = 0.5)[[1]])
    printf("\tPercent of theta1 out: %.3f", sum(abs(z_theta1) > 1.96) / Tt)
    printf("\tPercent of theta2 out: %.3f", sum(abs(z_theta2) > 1.96) / Tt)

    # Algorithm-specific summary lines
    if (!is.null(result$ac_hist)) {
        printf("Mean acceptance ratio of theta1: %.2f", mean(result$ac_hist))
    }
    if (!is.null(result$accepted_hist)) {
        printf("W1 MH acceptance rate: %.3f", mean(result$accepted_hist))
    }
    if (!is.null(result$accepted2_hist)) {
        printf("W2 MH acceptance rate: %.3f", mean(result$accepted2_hist))
    }
    if (!is.null(result$ce1_shape)) {
        printf("CE Gamma proposal for phi1: shape=%.4f, rate=%.4f", result$ce1_shape, result$ce1_rate)
    }
    if (!is.null(result$ce2_shape)) {
        printf("CE Gamma proposal for phi2: shape=%.4f, rate=%.4f", result$ce2_shape, result$ce2_rate)
    }

    # Plots

    x <- 1:Tt
    par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    plot(x, y, type = "l", xlab = "t", ylab = "", col = "gray",
         main = "Poisson Local Trend Polynomial Model")
    points(x, y, pch = 20)
    lines(x, lambda_mean, col = "red", lwd = 2)
    if (theta1_present) lines(x, exp(theta1_true), col = "blue", lwd = 2)
    legend("topright",
           legend = if (theta1_present) expression(y[t], lambda[t], hat(lambda)[t]) else expression(y[t], hat(lambda)[t]),
           col = if (theta1_present) c("black", "blue", "red") else c("black", "red"),
           lty = if (theta1_present) c(NA, 1, 1) else c(NA, 1),
           lwd = if (theta1_present) c(NA, 2, 2) else c(NA, 2),
           pch = if (theta1_present) c(20, NA, NA) else c(20, NA),
           bty = "n")

    par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    ylim_range <- if (theta1_present) range(theta1_mean, theta1_true) else range(theta1_mean)
    plot(x, theta1_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
         xlab = "t", ylab = "", main = "theta_t1")
    if (theta1_present) {
        lines(x, theta1_true, col = "blue", lwd = 2)
        legend("topright", legend = expression(hat(theta)[t1], theta[t1]), col = c("red", "blue"), lwd = 2, bty = "n")
    }

    par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    ylim_range <- if (theta2_present) range(theta2_mean, theta2_true) else range(theta2_mean)
    plot(x, theta2_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
         xlab = "t", ylab = "", main = "theta_t2")
    if (theta2_present) {
        lines(x, theta2_true, col = "blue", lwd = 2)
        legend("topright", legend = expression(hat(theta)[t2], theta[t2]), col = c("red", "blue"), lwd = 2, bty = "n")
    }

    par(mfrow = c(2, 2))
    for (t in t_obs) {
        hist(theta1_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
             xlab = bquote(theta[.(t) * "," * 1]), main = bquote("Posterior of " * theta[.(t) * "," * 1]))
        lines(density(theta1_hist[-(1:burnin), t]), col = "blue", lwd = 2)
    }

    par(mfrow = c(2, 2))
    for (t in t_obs) {
        hist(theta2_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
             xlab = bquote(theta[.(t) * "," * 2]), main = bquote("Posterior of " * theta[.(t) * "," * 2]))
        lines(density(theta2_hist[-(1:burnin), t]), col = "blue", lwd = 2)
    }

    par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    hist(W1_hist[-(1:burnin)], breaks = 50, freq = FALSE, main = "Posterior of W1")
    lines(density(W1_hist[-(1:burnin)]), col = "blue", lwd = 2)

    par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    hist(W2_hist[-(1:burnin)], breaks = 50, freq = FALSE, main = "Posterior of W2")
    lines(density(W2_hist[-(1:burnin)]), col = "blue", lwd = 2)

    par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    plot(W1_hist[-(1:burnin)], type = "l", xlab = "n", ylab = "W", main = "Traceplot of W1")
    plot(W2_hist[-(1:burnin)], type = "l", xlab = "n", ylab = "W", main = "Traceplot of W2")

    par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    plot(theta_01_hist, type = "l", main = "Traceplot of theta_01", xlab = "", ylab = "")
    plot(theta_02_hist, type = "l", main = "Traceplot of theta_02", xlab = "", ylab = "")

    par(mfrow = c(2, 2))
    for (t in t_obs) {
        plot(theta1_hist[, t], type = "l", main = bquote(theta[.(t) * "," * 1]), xlab = "", ylab = "")
        abline(v = burnin, col = "red")
    }

    par(mfrow = c(2, 2))
    for (t in t_obs) {
        plot(theta2_hist[, t], type = "l", main = bquote(theta[.(t) * "," * 2]), xlab = "", ylab = "")
        abline(v = burnin, col = "red")
    }

    par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    plot(ess_theta1, type = "l", main = expression("Effective sample of " * theta[t1]), xlab = "t")
    plot(ess_theta2, type = "l", main = expression("Effective sample of " * theta[t2]), xlab = "t")

    par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
    plot(z_theta1, type = "l", main = expression("Geweke diagnostic for " * theta[t1]), xlab = "t", ylab = "Z score")
    abline(h = c(-1.96, 1.96), col = "red")
    plot(z_theta2, type = "l", main = expression("Geweke diagnostic for " * theta[t2]), xlab = "t", ylab = "Z score")
    abline(h = c(-1.96, 1.96), col = "red")

    if (!is.null(nu_02) && !is.null(eta_02)) {
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        curve(dgamma(x, shape = nu_02, rate = eta_02), from = 1e-6, to = max(1 / W2_hist[-(1:burnin)]),
              main = "phi2 prior vs. posterior", col = "red", lwd = 2)
        lines(density(1 / W2_hist[-(1:burnin)]), col = "blue", lwd = 2)
        legend("topright", legend = c("Prior", "Posterior"), col = c("red", "blue"), lwd = 2)
    }
    if (!is.null(nu_01) && !is.null(eta_01) && !is.null(result$ce1_shape)) {
        # only sir_collapsed collapses W1 too -- otherwise phi1's posterior
        # isn't directly comparable to a fixed conjugate prior curve
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        curve(dgamma(x, shape = nu_01, rate = eta_01), from = 1e-6, to = max(1 / W1_hist[-(1:burnin)]),
              main = "phi1 prior vs. posterior", col = "red", lwd = 2)
        lines(density(1 / W1_hist[-(1:burnin)]), col = "blue", lwd = 2)
        legend("topright", legend = c("Prior", "Posterior"), col = c("red", "blue"), lwd = 2)
    }

    # ---- Algorithm-specific plots ----

    if (!is.null(result$ess_smc)) { # pg_as
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        plot(result$ess_smc, type = "l", main = "Effective Sample Size - SMC")
        abline(v = burnin, col = "red")
    }

    if (!is.null(result$ac_hist) && !is.null(ac_ref)) { # amh_montoril
        roll_mean <- function(v, k) {
            n <- length(v); out <- rep(NA, n)
            for (i in k:n) out[i] <- mean(v[(i - k + 1):i])
            out
        }
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        plot(rowMeans(result$ac_hist), type = "l", xlab = "n", ylab = "ratio",
             main = expression("Acceptance ratio of " * theta[t * 1] * " (mean over t)"))
        abline(h = ac_ref, col = "blue", lty = 2)
        abline(v = burnin, col = "red")
    }

    if (!is.null(result$ess_is) && !is.null(result$itr_irls)) { # sir_laplace / sir_collapsed
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        plot(result$ess_is, type = "l", main = "Effective Sample Size - IS", xlab = "n")
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        plot(result$itr_irls, type = "l", main = "IRLS iterations per Gibbs step", xlab = "n")
    }

    if (!is.null(result$ess_is_hist)) { # sir_collapsed
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        plot(result$ess_is_hist, type = "l", main = "ESS - IS Integrated Likelihood (W1)", xlab = "n")
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        plot(result$ess_sir_hist, type = "l", main = "ESS - SIR theta1", xlab = "n")
    }

    return(NULL)
}
