# Poisson - 2nd Order Polynomial Dynamic Model
#
# Strategy:
#  - theta1: Importance Sampling
#      Importance density: Laplace (normal) approximation
#  - theta2: Precision sampling (Chan)
#
# This file only defines sir_laplace_r(). It has no side effects (no data
# loading, no plotting) -- see tests/ for scripts that call it, either
# standalone (test_prototype_R.R) or against the Rcpp port (test_cpp.R,
# test_R_vs_cpp.R). Signature mirrors sir_laplace_cpp() 1:1.
#
# Author: Cleiton Moya de Almeida

library(Matrix)

# Resolve utils.R relative to THIS file's own location, not the caller's
# working directory -- so this file works whether it's source()'d from its
# own folder (standalone use) or from elsewhere (e.g. tests/), and even
# when nested inside another source() call.
.proto_dir <- local({
    ofile <- NULL
    for (i in rev(seq_along(sys.frames()))) {
        f <- sys.frame(i)$ofile
        if (!is.null(f)) { ofile <- f; break }
    }
    if (is.null(ofile)) "." else dirname(ofile)
})
source(file.path(.proto_dir, "utils.R"))

sir_laplace_r <- function(y, N,
                           mu_01, sigma2_01, mu_02, sigma2_02,
                           nu_01, eta_01, nu_02, eta_02,
                           theta1, theta2, theta_01, theta_02, W1, W2,
                           M_irls_max, M_is, tol,
                           verbose = TRUE, print_every = 1000) {

    printf <- function(...) cat(paste(sprintf(...), "\n"))
    Tt <- length(y)
    Ttp1 <- Tt + 1

    chan_smoothing_theta2 <- make_chan_theta2_smoother_ext(Ttp1)
    chan_smoothing_theta1 <- make_chan_theta1_smoother_ext(Ttp1)

    W1_hist <- numeric(N)
    W2_hist <- numeric(N)
    theta_01_hist <- numeric(N)
    theta_02_hist <- numeric(N)
    theta1_hist <- matrix(0, N, Tt)
    theta2_hist <- matrix(0, N, Tt)
    ess_is <- numeric(N)
    itr_irls <- numeric(N)
    theta1_tilde <- numeric(Tt)

    start_time <- proc.time()
    for (n in 1:N) {

        if (verbose && (n %% print_every == 0)) {
            elapsed_time <- (proc.time() - start_time)[[1]]
            printf("Iteration %d / %d | Elapsed CPU time: %.0f s", n, N, elapsed_time)
        }

        # Sample phi1 (conjugated gamma)
        phi1 <- gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1, theta_02, theta2, Tt)
        W1 <- 1 / phi1

        # Sample phi2 (conjugated gamma)
        phi2 <- gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt)
        W2 <- 1 / phi2

        #
        # Importance Sampling for theta_t1
        #

        # Local approximation (IRLS)
        theta1_tilde_old <- theta1_tilde
        for (j in 1:M_irls_max) {
            f_t <- exp(-theta1_tilde)         # observational variance
            z_t <- theta1_tilde + f_t * y - 1 # pseudo-observation

            res <- chan_smoothing_theta1(z_t, 1 / f_t, phi1, mu_01, sigma2_01, theta_02, theta2)
            theta1_tilde <- res$theta1_hat[-1] # drop node 0 (theta_01)

            if (max(abs(theta1_tilde - theta1_tilde_old)) < tol) break
            theta1_tilde_old <- theta1_tilde
        }
        itr_irls[n] <- j

        # Importance Sampling step
        if (M_is > 1) {
            log_w <- numeric(M_is)
            trajectories <- matrix(0, M_is, Tt)
            theta_01_trajectories <- numeric(M_is)
            for (i in 1:M_is) {
                draw_prop <- chan_sample_from_build(res, Ttp1)
                theta_01_trajectories[i] <- draw_prop[1]
                theta1_prop <- draw_prop[-1]
                trajectories[i, ] <- theta1_prop

                log_p <- sum(y * theta1_prop - exp(theta1_prop))
                log_g <- sum(-0.5 * log(2 * pi * f_t) - 0.5 * (theta1_prop - z_t)^2 / f_t)
                log_w[i] <- log_p - log_g
            }

            log_w <- log_w - logsumexp(log_w)
            ess_is[n] <- 1 / sum(exp(2 * log_w))

            idx <- sample_one_from_logw(log_w)
            theta_01 <- theta_01_trajectories[idx]
            theta1   <- trajectories[idx, ]
        } else {
            draw_prop <- chan_sample_from_build(res, Ttp1)
            theta_01  <- draw_prop[1]
            theta1    <- draw_prop[-1]
        }

        # (theta_02, theta2) jointly via extended block
        build2 <- chan_smoothing_theta2(theta1, phi1, phi2, mu_02, sigma2_02, theta_01)
        draw2  <- chan_sample_from_build(build2, Ttp1)
        theta_02 <- draw2[1]
        theta2   <- draw2[-1]

        theta_01_hist[n] <- theta_01
        theta_02_hist[n] <- theta_02
        W1_hist[n] <- W1
        W2_hist[n] <- W2
        theta1_hist[n, ] <- theta1
        theta2_hist[n, ] <- theta2
    }

    list(theta1_hist = theta1_hist,
         theta2_hist = theta2_hist,
         theta_01_hist = theta_01_hist,
         theta_02_hist = theta_02_hist,
         W1_hist = W1_hist,
         W2_hist = W2_hist,
         ess_is = ess_is,
         itr_irls = itr_irls)
}
