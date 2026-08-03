# Poisson - 2nd Order Polynomial Dynamic Model
#
# Model:
#  y_t ~ Poisson(exp{theta_t1})
#  theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1, omega_t1 ~ N(0, W1)
#  theta_t2 =                 theta_{t-1,1} + omega_t2, omega_t2 ~ N(0, W2)
#
# MCMC:
#  theta_t1: Adaptive Metropolis within Gibbs (Roberts & Rosenthal, 2009)
#            with continuous Robbins-Monro update of the \log\sigma_t
#            (Andrieu & Thoms 2008, Eq. 20/22) instead of batch update.
#  theta2: Chan
#
# This file only defines amh_montoril_r(). It has no side effects (no data
# loading, no plotting) -- see tests/ for scripts that call it, either
# standalone (test_prototype_R.R) or against the Rcpp port (test_cpp.R,
# test_validation.R). Signature mirrors amh_montoril_cpp() 1:1.
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

amh_montoril_r <- function(y, N,
                            mu_01, sigma2_01, mu_02, sigma2_02,
                            nu_01, eta_01, nu_02, eta_02,
                            theta1, theta2, theta_01, theta_02, W1, W2,
                            ac_ref, varsigma2,
                            verbose = TRUE, print_every = 1000) {

    printf <- function(...) cat(paste(sprintf(...), "\n"))
    Tt <- length(y)
    Ttp1 <- Tt + 1

    chan_smoothing_theta2 <- make_chan_theta2_smoother_ext(Ttp1)

    theta1_hist <- matrix(nrow = N, ncol = Tt)
    theta2_hist <- matrix(nrow = N, ncol = Tt)
    W1_hist <- numeric(N)
    W2_hist <- numeric(N)
    theta_01_hist <- numeric(N)
    theta_02_hist <- numeric(N)
    ac_hist <- matrix(0, nrow = N, ncol = Tt)

    start_time <- proc.time()
    for (n in 1:N) {

        if (verbose && (n %% print_every == 0)) {
            elapsed_time <- (proc.time() - start_time)[[1]]
            printf("Iteration %d / %d | mean acception rate of theta_t: %.2f | Elapsed CPU time: %.0f s",
                   n, N, mean(ac_hist[1:(n - 1), ]), elapsed_time)
        }

        # Sample theta_01
        theta_01 <- gibbs_sample_theta01(mu_01, sigma2_01, theta1[1], theta_02, W1)

        # Sample phi1
        phi1 <- gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1, theta_02, theta2, Tt)
        W1 <- 1 / phi1

        # Sample phi2
        phi2 <- gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt)
        W2 <- 1 / phi2

        # Sample theta_t1 (random walking Metropolis step)
        for (t in 1:Tt) {

            if (t < Tt) {
                if (t == 1) {
                    res <- sample_theta_t1_mh(theta1[t], theta_01, theta1[t + 1],
                                           theta2[t], theta_02,
                                           y[t], W1, varsigma2[t], final_t = FALSE)
                    theta1[t] <- res$theta_t1
                } else {
                    res <- sample_theta_t1_mh(theta1[t], theta1[t - 1], theta1[t + 1],
                                           theta2[t], theta2[t - 1],
                                           y[t], W1, varsigma2[t], final_t = FALSE)
                    theta1[t] <- res$theta_t1
                }
            } else {
                res <- sample_theta_t1_mh(theta1[t], theta1[t - 1], NULL,
                                       theta2[t], theta2[t - 1],
                                       y[t], W1, varsigma2[t], final_t = TRUE)
                theta1[t] <- res$theta_t1
            }

            ac_hist[n, t] <- res$ac
        }

        # Adaptive stage of varsigma2
        delta <- min(0.01, 1 / sqrt(n))
        ls <- log(varsigma2) / 2 + delta * (ac_hist[n, ] - ac_ref)
        varsigma2 <- exp(2 * ls)

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
         ac_hist = ac_hist)
}
