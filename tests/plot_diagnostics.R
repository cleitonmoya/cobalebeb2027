# tests/plot_diagnostics.R
#
# Shared diagnostics/plotting routine for the Poisson 2nd-order polynomial
# DLM samplers. Used by both test_prototype_R.R (pure R) and test_cpp.R
# (Rcpp), so the same battery of plots/summary stats is produced regardless
# of which implementation was run -- only test_validation.R (which compares
# the two numerically) does not use this file.
#
# `result_list` must be a list of N_chains elements, each containing at
# least: theta1_hist, theta2_hist, theta_01_hist, theta_02_hist, W1_hist,
# W2_hist. Algorithm-specific extra diagnostics are added when the
# corresponding fields are present in result_list[[plot_chain]] (ess_smc
# for pg_as; ac_hist for amh_montoril; ess_is/itr_irls for sir_laplace;
# accepted_hist/accepted2_hist/itr_irls/ess_is_hist/ess_sir_hist/ce1_*/ce2_*
# for sir_collapsed).
#
# R_hat and ESS bulk/tail are computed via PoissonLTDM::metrics_convergence()
# / metrics_convergence_by_chain() (R/metrics.R + src/convergence.cpp),
# which wrap a C++ (Rcpp + RcppArmadillo + FFTW3 + OpenMP) reimplementation
# of the R package posterior's rhat()/ess_bulk()/ess_tail() (Vehtari et al.
# 2021). Output was validated to floating-point precision against
# posterior::summarise_draws() (this file no longer depends on the
# posterior package itself), avoiding the R-level per-variable loop
# overhead that made the posterior-based version take 90-100+ seconds per
# metric at Tt=1600.
#
# R_hat (metrics_convergence()) is computed across all chains in
# result_list. ESS bulk/tail (metrics_convergence_by_chain()) are computed
# PER CHAIN (not pooled): R_hat answers whether the chains agree with each
# other; ESS answers how much independent information a given chain
# carries. A chain can have high individual ESS while still sitting in a
# different mode from the others -- pooling ESS across chains would mask
# that, since R_hat (computed separately) is what is supposed to catch it.
# Only the chain selected via `plot_chain` is used for the plots.
#
# pkgload::load_all("../PoissonLTDM", debug = FALSE) must be called before
# this file is sourced (done in test_cpp.R / test_prototype_R.R), so that
# metrics_convergence(), metrics_convergence_by_chain(), and
# rhat_ess_fast() are available.

print_and_plot_diagnostics <- function(result_list, y, changepoints, 
                                        theta1_true = NULL, theta2_true = NULL,
                                        t_obs, burnin, elapsed_time,
                                        plot_chain = 1,
                                        nu_01 = NULL, eta_01 = NULL,
                                        nu_02 = NULL, eta_02 = NULL,
                                        ac_ref = NULL,
                                        plots = TRUE,
                                        compute_rhat = TRUE, 
                                        compute_ess = TRUE) {

    printf <- function(...) cat(paste(sprintf(...), "\n"))

    N_chains <- length(result_list)
    if (plot_chain < 1 || plot_chain > N_chains) {
        stop(sprintf("plot_chain = %d out of range (1..%d).", plot_chain, N_chains))
    }

    # ---- Chain selected for plotting / detailed trace ----
    result <- result_list[[plot_chain]]

    theta_01_hist <- result$theta_01_hist
    theta_02_hist <- result$theta_02_hist
    W1_hist <- result$W1_hist
    W2_hist <- result$W2_hist

    N <- nrow(result$theta1_hist)
    Tt <- ncol(result$theta1_hist)
    theta1_present <- !is.null(theta1_true)
    theta2_present <- !is.null(theta2_true)

    # ---- Memory: slice post-burnin once per chain, up front ----
    #
    # theta1_hist/theta2_hist are N x Tt matrices, potentially large enough
    # (depending on N, Tt) to dominate RAM when kept in full for every
    # chain through the whole function. Downstream code needs the
    # post-burnin slice repeatedly (summary stats, metrics_convergence()
    # calls); slicing once here and dropping the pre-burnin portion from
    # result_list (when plots=FALSE) lets it be garbage-collected
    # immediately instead of staying retained for the whole function.
    #
    # Only safe to drop the pre-burnin portion when plots=FALSE: the trace
    # plots further down (theta1_hist[, t] with abline(v = burnin)) need the
    # full history to show the burnin cutoff. When plots=TRUE, behavior is
    # unchanged (full history kept, as before).
    post_idx0 <- (burnin + 1):N
    for (c in 1:N_chains) {
        result_list[[c]]$theta1_hist_post <- result_list[[c]]$theta1_hist[post_idx0, ]
        result_list[[c]]$theta2_hist_post <- result_list[[c]]$theta2_hist[post_idx0, ]
        if (!plots) {
            result_list[[c]]$theta1_hist <- NULL
            result_list[[c]]$theta2_hist <- NULL
        }
    }
    gc()

    result <- result_list[[plot_chain]]
    theta1_hist_post <- result$theta1_hist_post
    theta2_hist_post <- result$theta2_hist_post
    if (plots) {
        # Full pre-burnin history, only kept/needed for the trace plots below
        theta1_hist <- result$theta1_hist
        theta2_hist <- result$theta2_hist
    }

    # Summary stats (plotted/selected chain only)
    theta1_mean <- colMeans(theta1_hist_post)
    theta2_mean <- colMeans(theta2_hist_post)
    lambda_mean <- exp(theta1_mean)

    printf("---- Chain %d/%d (selected for plots/summary) ----", plot_chain, N_chains)
    printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
    printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))
    printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
    printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))

    loglik <- sum(dpois(y, lambda_mean, log = TRUE))
    printf("Log-likelihood: %.2f", loglik)


    # ---- R_hat and ESS bulk/tail via PoissonLTDM::metrics_convergence() ----
    if (compute_rhat) {
        conv_theta01 <- metrics_convergence(result_list, "theta_01_hist", burnin)
        conv_theta02 <- metrics_convergence(result_list, "theta_02_hist", burnin)
        conv_w1 <- metrics_convergence(result_list, "W1_hist", burnin)
        conv_w2 <- metrics_convergence(result_list, "W2_hist", burnin)
        
        rhat_theta01 <- conv_theta01$rhat
        rhat_theta02 <- conv_theta02$rhat
        rhat_w1 <- conv_w1$rhat
        rhat_w2 <- conv_w2$rhat
        
        printf("R_hat (across %d chains):", N_chains)
        printf("\ttheta_01: %.4f", rhat_theta01)
        printf("\ttheta_02: %.4f", rhat_theta02)
        printf("\tW1: %.4f", rhat_w1)
        printf("\tW2: %.4f", rhat_w2)
        
        conv_theta1 <- metrics_convergence(result_list, "theta1_hist_post", 0)
        conv_theta2 <- metrics_convergence(result_list, "theta2_hist_post", 0)
        
        rhat_theta1 <- conv_theta1$rhat   # length Tt
        rhat_theta2 <- conv_theta2$rhat   # length Tt
        
        printf("\ttheta1 (mean): %.4f", mean(rhat_theta1))
        printf("\ttheta1 (max): %.4f", max(rhat_theta1))
        printf("\ttheta2 (mean): %.4f", mean(rhat_theta2))
        printf("\ttheta2 (max): %.4f", max(rhat_theta2))
    }
    
    

    # ---- ESS bulk/tail (PER CHAIN, all N_chains) ----
    #
    # metrics_convergence_by_chain() computes ESS separately per chain (not
    # pooled) -- see file header for why. For theta1_hist/theta2_hist, we
    # reuse the already-sliced *_hist_post fields (burnin=0 here, since the
    # slicing already happened above) instead of re-slicing the full
    # history again.
    if (compute_ess) {
        ess01 <- metrics_convergence_by_chain(result_list, "theta_01_hist", burnin)
        ess02 <- metrics_convergence_by_chain(result_list, "theta_02_hist", burnin)
        essw1 <- metrics_convergence_by_chain(result_list, "W1_hist", burnin)
        essw2 <- metrics_convergence_by_chain(result_list, "W2_hist", burnin)
        
        ess_theta01 <- ess01$ess_bulk
        ess_theta01_tail <- ess01$ess_tail
        ess_theta02 <- ess02$ess_bulk
        ess_theta02_tail <- ess02$ess_tail
        ess_w1 <- essw1$ess_bulk
        ess_w1_tail <- essw1$ess_tail
        ess_w2 <- essw2$ess_bulk
        ess_w2_tail <- essw2$ess_tail
        
        ess_theta1_by_t <- metrics_convergence_by_chain(result_list, "theta1_hist_post", 0)
        ess_theta2_by_t <- metrics_convergence_by_chain(result_list, "theta2_hist_post", 0)
        ess_theta1_bychain <- ess_theta1_by_t$ess_bulk             # Tt x N_chains
        ess_theta1_tail_bychain <- ess_theta1_by_t$ess_tail
        ess_theta2_bychain <- ess_theta2_by_t$ess_bulk
        ess_theta2_tail_bychain <- ess_theta2_by_t$ess_tail
        
        printf("Effective Sample Size (bulk / tail) per chain:")
        printf("\ttheta_01: %s", paste(sprintf("%.2f/%.2f", ess_theta01, ess_theta01_tail), collapse = "  "))
        printf("\t  CV(bulk): %.1f%%", sd(ess_theta01) / mean(ess_theta01) * 100)
        printf("\ttheta_02: %s", paste(sprintf("%.2f/%.2f", ess_theta02, ess_theta02_tail), collapse = "  "))
        printf("\t  CV(bulk): %.1f%%", sd(ess_theta02) / mean(ess_theta02) * 100)
        printf("\tW1: %s", paste(sprintf("%.0f/%.0f", ess_w1, ess_w1_tail), collapse = "  "))
        printf("\t  CV(bulk): %.1f%%", sd(ess_w1) / mean(ess_w1) * 100)
        printf("\tW2: %s", paste(sprintf("%.0f/%.0f", ess_w2, ess_w2_tail), collapse = "  "))
        printf("\t  CV(bulk): %.1f%%", sd(ess_w2) / mean(ess_w2) * 100)
        printf("\ttheta1 (mean over t): %s",
               paste(sprintf("%.2f/%.2f", colMeans(ess_theta1_bychain), colMeans(ess_theta1_tail_bychain)), collapse = "  "))
        printf("\t  CV(bulk): %.1f%%", sd(colMeans(ess_theta1_bychain)) / mean(colMeans(ess_theta1_bychain)) * 100)
        printf("\ttheta1 (min over t): %s",
               paste(sprintf("%.2f/%.2f", apply(ess_theta1_bychain, 2, min), apply(ess_theta1_tail_bychain, 2, min)), collapse = "  "))
        printf("\t  CV(bulk): %.1f%%", sd(apply(ess_theta1_bychain, 2, min)) / mean(apply(ess_theta1_bychain, 2, min)) * 100)
        printf("\ttheta2 (mean over t): %s",
               paste(sprintf("%.2f/%.2f", colMeans(ess_theta2_bychain), colMeans(ess_theta2_tail_bychain)), collapse = "  "))
        printf("\t  CV(bulk): %.1f%%", sd(colMeans(ess_theta2_bychain)) / mean(colMeans(ess_theta2_bychain)) * 100)
        printf("\ttheta2 (min over t): %s",
               paste(sprintf("%.2f/%.2f", apply(ess_theta2_bychain, 2, min), apply(ess_theta2_tail_bychain, 2, min)), collapse = "  "))
        printf("\t  CV(bulk): %.1f%%", sd(apply(ess_theta2_bychain, 2, min)) / mean(apply(ess_theta2_bychain, 2, min)) * 100)
        
        
        # Per-t ESS for the selected (plot_chain) chain -- used below in the ESS plot
        ess_theta1 <- ess_theta1_bychain[, plot_chain]
        ess_theta1_tail <- ess_theta1_tail_bychain[, plot_chain]
        ess_theta2 <- ess_theta2_bychain[, plot_chain]
        ess_theta2_tail <- ess_theta2_tail_bychain[, plot_chain]
        
        if (!missing(elapsed_time) && !is.null(elapsed_time)) {
            printf("Effective Sample Size / second (bulk / tail, chain %d, elapsed_time = total for all chains):", plot_chain)
            printf("\tW1: %.2f / %.2f", ess_w1[plot_chain]/elapsed_time, ess_w1_tail[plot_chain]/elapsed_time)
            printf("\tW2: %.2f / %.2f", ess_w2[plot_chain]/elapsed_time, ess_w2_tail[plot_chain]/elapsed_time)
            printf("\ttheta1 (mean): %.2f / %.2f", mean(ess_theta1 / elapsed_time), mean(ess_theta1_tail / elapsed_time))
            printf("\ttheta1 (min.): %.2f / %.2f", min(ess_theta1 / elapsed_time), min(ess_theta1_tail / elapsed_time))
            printf("\ttheta2 (mean): %.2f / %.2f", mean(ess_theta2 / elapsed_time), mean(ess_theta2_tail / elapsed_time))
            printf("\ttheta2 (min.): %.2f / %.2f", min(ess_theta2 / elapsed_time), min(ess_theta2_tail / elapsed_time))
        }
    }
    

    # Algorithm-specific summary lines (selected chain only)
    if (!is.null(result$ac_hist)) {
        printf("Mean acceptance ratio of theta1: %.2f", mean(result$ac_hist))
        if (is.null(dim(result$ac_hist))) { 
            ac_mean <- mean(result$ac_hist[changepoints])
        } else {
            ac_mean <- mean(result$ac_hist[, changepoints])
        }
        printf("Mean acceptance ratio of theta1 at changepoints: %.2f", ac_mean)
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

    # Plots (selected chain only, `result`/`plot_chain`)
    if (plots) {
        
        # Observations and theta1
        x <- 1:Tt
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        plot(x, y, type = "l", xlab = "t", ylab = "", col = "gray",
             main = sprintf("Poisson Local Trend Polynomial Model (chain %d)", plot_chain))
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
        
        # theta1, theta2 
        par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
        ylim_range <- if (theta1_present) range(theta1_mean, theta1_true) else range(theta1_mean)
        plot(x, theta1_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
             xlab = "t", ylab = "", main = sprintf("theta_t1 (chain %d)", plot_chain))
        if (theta1_present) {
            lines(x, theta1_true, col = "blue", lwd = 2)
            legend("topright", legend = expression(hat(theta)[t1], theta[t1]), col = c("red", "blue"), lwd = 2, bty = "n")
        }
        
        ylim_range <- if (theta2_present) range(theta2_mean, theta2_true) else range(theta2_mean)
        plot(x, theta2_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
             xlab = "t", ylab = "", main = sprintf("theta_t2 (chain %d)", plot_chain))
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
        
        
        if (compute_rhat) {
            par(mfrow = c(3, 1), mar = c(4, 4, 2, 2), cex = 0.8)
            plot(rhat_theta1, type = "l", main = expression("R-hat of " * theta[t1]), xlab = "t")
            abline(h = 1.01, col = "red", lty = 2)
            abline(v=changepoints, col="green")
            
            plot(rhat_theta2, type = "l", main = expression("R-hat of " * theta[t2]), xlab = "t")
            abline(h = 1.01, col = "red", lty = 2)
            abline(v=changepoints, col="green")
            
            ylim_range <- if (theta1_present) range(theta1_mean, theta1_true) else range(theta1_mean)
            plot(x, theta1_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
                 xlab = "t", ylab = "", main = sprintf("theta_t1 (chain %d)", plot_chain))
            if (theta1_present) {
                lines(x, theta1_true, col = "blue", lwd = 2)
                legend("topright", legend = expression(hat(theta)[t1], theta[t1]), col = c("red", "blue"), lwd = 2, bty = "n")
            }
        }
        
        if (compute_ess) {
            par(mfrow = c(3, 1), mar = c(4, 4, 2, 2), cex = 0.8)
            plot(ess_theta1, type = "l", main = expression("Effective sample of " * theta[t1]), xlab = "t")
            abline(h = 400, col = "red", lty = 2)
            abline(v=changepoints, col="green")
            
            plot(ess_theta2, type = "l", main = expression("Effective sample of " * theta[t2]), xlab = "t")
            abline(h = 400, col = "red", lty = 2)
            abline(v=changepoints, col="green")
            
            ylim_range <- if (theta1_present) range(theta1_mean, theta1_true) else range(theta1_mean)
            plot(x, theta1_mean, type = "l", col = "red", lwd = 2, ylim = ylim_range,
                 xlab = "t", ylab = "", main = sprintf("theta_t1 (chain %d)", plot_chain))
            if (theta1_present) {
                lines(x, theta1_true, col = "blue", lwd = 2)
                legend("topright", legend = expression(hat(theta)[t1], theta[t1]), col = c("red", "blue"), lwd = 2, bty = "n")
            }
            
        }
        
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
        
        if (!is.null(result$ac_hist) && !is.null(ac_ref)) { # amh_montoril (matrix) / stan (vector)
            roll_mean <- function(v, k) {
                n <- length(v); out <- rep(NA, n)
                for (i in k:n) out[i] <- mean(v[(i - k + 1):i])
                out
            }
            # amh_montoril's ac_hist is an N x Tt matrix (one MH acceptance
            # rate per (iteration, t) -- rowMeans() gives the mean over t
            # per iteration. stan's ac_hist (accept_stat__) is already a
            # length-N vector, one value per iteration -- nothing to
            # average over, plot it directly. Same is.null(dim(...))
            # distinction already used for the console printout above.
            ac_series <- if (is.null(dim(result$ac_hist))) result$ac_hist else rowMeans(result$ac_hist)
            par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
            plot(ac_series, type = "l", xlab = "n", ylab = "ratio",
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
    }
    
    # ---- Return computed diagnostics (no recomputation needed by callers) ----
    #
    # Every value printed to the console above is captured here, unchanged --
    # nothing new is computed. NULL where the corresponding compute_*/optional
    # block did not run (e.g. ess_w1_tail is NULL when compute_ess = FALSE;
    # ac_ratio_mean is NULL for methods other than amh_montoril). rhat_max and
    # ess_bulk_min are convenience aggregates (worst-case R_hat across all
    # scalar and time-indexed parameters; worst-case ESS bulk across all
    # chains/parameters) for callers that only need a one-number pass/fail
    # summary, e.g. calibration_phase.R.
    out <- list(
        # ---- Fit summary (selected/plot_chain only) ----
        W1_mean = mean(W1_hist[-(1:burnin)]),
        W1_median = median(W1_hist[-(1:burnin)]),
        W2_mean = mean(W2_hist[-(1:burnin)]),
        W2_median = median(W2_hist[-(1:burnin)]),
        loglik = loglik,

        # ---- R_hat (pooled across chains) ----
        rhat_theta01 = if (compute_rhat) rhat_theta01 else NULL,
        rhat_theta02 = if (compute_rhat) rhat_theta02 else NULL,
        rhat_w1      = if (compute_rhat) rhat_w1 else NULL,
        rhat_w2      = if (compute_rhat) rhat_w2 else NULL,
        rhat_theta1       = if (compute_rhat) rhat_theta1 else NULL,       # length Tt
        rhat_theta1_mean  = if (compute_rhat) mean(rhat_theta1) else NULL,
        rhat_theta1_max   = if (compute_rhat) max(rhat_theta1) else NULL,
        rhat_theta2       = if (compute_rhat) rhat_theta2 else NULL,       # length Tt
        rhat_theta2_mean  = if (compute_rhat) mean(rhat_theta2) else NULL,
        rhat_theta2_max   = if (compute_rhat) max(rhat_theta2) else NULL,
        rhat_max     = if (compute_rhat) {
            max(rhat_theta01, rhat_theta02, rhat_w1, rhat_w2, rhat_theta1, rhat_theta2)
        } else NULL,

        # ---- ESS bulk/tail, PER CHAIN (all N_chains) ----
        ess_theta01 = if (compute_ess) ess_theta01 else NULL,               # length N_chains
        ess_theta01_tail = if (compute_ess) ess_theta01_tail else NULL,
        ess_theta01_cv   = if (compute_ess) sd(ess_theta01) / mean(ess_theta01) * 100 else NULL,
        ess_theta02 = if (compute_ess) ess_theta02 else NULL,
        ess_theta02_tail = if (compute_ess) ess_theta02_tail else NULL,
        ess_theta02_cv   = if (compute_ess) sd(ess_theta02) / mean(ess_theta02) * 100 else NULL,
        ess_w1 = if (compute_ess) ess_w1 else NULL,
        ess_w1_tail = if (compute_ess) ess_w1_tail else NULL,
        ess_w1_cv   = if (compute_ess) sd(ess_w1) / mean(ess_w1) * 100 else NULL,
        ess_w2 = if (compute_ess) ess_w2 else NULL,
        ess_w2_tail = if (compute_ess) ess_w2_tail else NULL,
        ess_w2_cv   = if (compute_ess) sd(ess_w2) / mean(ess_w2) * 100 else NULL,

        ess_theta1_bychain = if (compute_ess) ess_theta1_bychain else NULL,           # Tt x N_chains
        ess_theta1_tail_bychain = if (compute_ess) ess_theta1_tail_bychain else NULL,
        ess_theta1_mean_over_t = if (compute_ess) colMeans(ess_theta1_bychain) else NULL,           # length N_chains
        ess_theta1_mean_over_t_tail = if (compute_ess) colMeans(ess_theta1_tail_bychain) else NULL,
        ess_theta1_mean_over_t_cv = if (compute_ess) {
            sd(colMeans(ess_theta1_bychain)) / mean(colMeans(ess_theta1_bychain)) * 100
        } else NULL,
        ess_theta1_min_over_t = if (compute_ess) apply(ess_theta1_bychain, 2, min) else NULL,       # length N_chains
        ess_theta1_min_over_t_tail = if (compute_ess) apply(ess_theta1_tail_bychain, 2, min) else NULL,
        ess_theta1_min_over_t_cv = if (compute_ess) {
            sd(apply(ess_theta1_bychain, 2, min)) / mean(apply(ess_theta1_bychain, 2, min)) * 100
        } else NULL,

        ess_theta2_bychain = if (compute_ess) ess_theta2_bychain else NULL,
        ess_theta2_tail_bychain = if (compute_ess) ess_theta2_tail_bychain else NULL,
        ess_theta2_mean_over_t = if (compute_ess) colMeans(ess_theta2_bychain) else NULL,
        ess_theta2_mean_over_t_tail = if (compute_ess) colMeans(ess_theta2_tail_bychain) else NULL,
        ess_theta2_mean_over_t_cv = if (compute_ess) {
            sd(colMeans(ess_theta2_bychain)) / mean(colMeans(ess_theta2_bychain)) * 100
        } else NULL,
        ess_theta2_min_over_t = if (compute_ess) apply(ess_theta2_bychain, 2, min) else NULL,
        ess_theta2_min_over_t_tail = if (compute_ess) apply(ess_theta2_tail_bychain, 2, min) else NULL,
        ess_theta2_min_over_t_cv = if (compute_ess) {
            sd(apply(ess_theta2_bychain, 2, min)) / mean(apply(ess_theta2_bychain, 2, min)) * 100
        } else NULL,

        ess_bulk_min = if (compute_ess) {
            min(ess_theta01, ess_theta02, ess_w1, ess_w2,
                ess_theta1_bychain, ess_theta2_bychain)
        } else NULL,

        # ---- ESS / second (chain = plot_chain, elapsed_time = total for all chains) ----
        # Guarded by the same condition as the console printout above:
        # elapsed_time may be missing/NULL if the caller omits it.
        ess_sec_w1_bulk = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            ess_w1[plot_chain] / elapsed_time
        } else NULL,
        ess_sec_w1_tail = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            ess_w1_tail[plot_chain] / elapsed_time
        } else NULL,
        ess_sec_w2_bulk = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            ess_w2[plot_chain] / elapsed_time
        } else NULL,
        ess_sec_w2_tail = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            ess_w2_tail[plot_chain] / elapsed_time
        } else NULL,
        ess_sec_theta1_mean_bulk = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            mean(ess_theta1 / elapsed_time)
        } else NULL,
        ess_sec_theta1_mean_tail = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            mean(ess_theta1_tail / elapsed_time)
        } else NULL,
        ess_sec_theta1_min_bulk  = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            min(ess_theta1 / elapsed_time)
        } else NULL,
        ess_sec_theta1_min_tail  = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            min(ess_theta1_tail / elapsed_time)
        } else NULL,
        ess_sec_theta2_mean_bulk = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            mean(ess_theta2 / elapsed_time)
        } else NULL,
        ess_sec_theta2_mean_tail = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            mean(ess_theta2_tail / elapsed_time)
        } else NULL,
        ess_sec_theta2_min_bulk  = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            min(ess_theta2 / elapsed_time)
        } else NULL,
        ess_sec_theta2_min_tail  = if (compute_ess && !missing(elapsed_time) && !is.null(elapsed_time)) {
            min(ess_theta2_tail / elapsed_time)
        } else NULL,

        # ---- Algorithm-specific summary values (selected chain only) ----
        ac_ratio_mean = if (!is.null(result$ac_hist)) mean(result$ac_hist) else NULL,          # amh_montoril
        ac_ratio_at_changepoints = if (!is.null(result$ac_hist)) {
            if (is.null(dim(result$ac_hist))) mean(result$ac_hist[changepoints])
            else mean(result$ac_hist[, changepoints])
        } else NULL,
        w1_mh_acceptance_rate = if (!is.null(result$accepted_hist)) mean(result$accepted_hist) else NULL,   # sir_collapsed
        w2_mh_acceptance_rate = if (!is.null(result$accepted2_hist)) mean(result$accepted2_hist) else NULL, # sir_collapsed
        ce1_shape = if (!is.null(result$ce1_shape)) result$ce1_shape else NULL,  # sir_collapsed
        ce1_rate  = if (!is.null(result$ce1_shape)) result$ce1_rate else NULL,
        ce2_shape = if (!is.null(result$ce2_shape)) result$ce2_shape else NULL,
        ce2_rate  = if (!is.null(result$ce2_shape)) result$ce2_rate else NULL
    )

    return(out)
}
