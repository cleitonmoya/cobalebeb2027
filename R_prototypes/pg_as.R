# Poisson - 2nd Order Polynomial Dynamic Model
# Particle Gibbs with Ancestral Sampling (PG-AS) + Component-wise Gibbs for theta_t2
# Strategy: - Bootstrap Particle Filter with Ancestral Sampling for theta_t1
#             (scalar state), reconstructed by simple Backward Tracking;
#           - theta2 sampled via Precision Matrix (Chan Method)
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269-342.
#
# PERFORMANCE NOTE: the resampling of the K-1 non-reference particles uses
# SYSTEMATIC resampling (Kitagawa, 1996; Doucet, de Freitas & Gordon, 2001)
# instead of multinomial resampling: a single runif() draw generates all K-1
# indices via one monotone sweep through the cumulative weights
# (findInterval), instead of K-1 independent draws inside sample(). This is
# a standard, unbiased SMC resampling scheme with LOWER variance than
# multinomial resampling (not merely a speed shortcut). The single-index
# Ancestral Sampling draw (a_ref) is unchanged -- it was already a single
# draw, so there is nothing to gain there.
#
# This file only defines pg_as_r(). It has no side effects (no data
# loading, no plotting) -- see tests/ for scripts that call it, either
# standalone (test_prototype_R.R) or against the Rcpp port (test_cpp.R,
# test_validation.R). Signature mirrors pg_as_cpp() 1:1.
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

pg_as_r <- function(y, K, N,
                     mu_01, sigma2_01, mu_02, sigma2_02,
                     nu_01, eta_01, nu_02, eta_02,
                     theta1, theta2, theta_01, theta_02, W1, W2,
                     verbose = TRUE, print_every = 100) {

    printf <- function(...) cat(paste(sprintf(...), "\n"))
    Tt <- length(y)
    Ttp1 <- Tt + 1

    # Auxiliary variables
    W1_hist <- numeric(N)
    W2_hist <- numeric(N)
    theta_01_hist <- numeric(N)
    theta_02_hist <- numeric(N)
    theta1_hist <- matrix(0, N, Tt)
    theta2_hist <- matrix(0, N, Tt)
    ess_smc <- numeric(N)

    # Chan method (theta_02, theta2 EXTENDED, (T+1)-dimensional block)
    chan_smoothing_theta2 <- make_chan_theta2_smoother_ext(Ttp1)

    start_time <- proc.time()
    for (n in 1:N) {

        if (verbose && (n %% print_every == 0)) {
            elapsed_time <- (proc.time() - start_time)[[1]]
            printf("Iteration %d / %d | Elapsed CPU time: %.0f s", n, N, elapsed_time)
        }

        # Sample W1
        phi1 <- gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1, theta_02, theta2, Tt)
        W1 <- 1 / phi1
        sd_W1 <- sqrt(W1)

        # Sample W2
        phi2 <- gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt)
        W2 <- 1 / phi2

        # Sample theta_01 (conjugated Normal)
        # theta_01 is NOT part of the particle system: it has no likelihood
        # term (no y_0), so its exact full conditional is available in
        # closed form.
        theta_01 <- gibbs_sample_theta01(mu_01, sigma2_01, theta1[1], theta_02, W1)

        #
        # Conditional SMC for theta_t1 (Bootstrap PF + Ancestral Sampling)
        #
        theta_1_k   <- matrix(0, Tt, K)
        log_w_tilde <- matrix(0, Tt, K)
        A_hist      <- matrix(0, Tt, K)   # ancestor indices, for Backward Tracking

        # t = 1
        theta_1_k[1, ] <- rnorm(K, mean = theta_01 + theta_02, sd = sd_W1)
        theta_1_k[1, K] <- theta1[1]

        log_w_1 <- log_p_yt(y[1], theta_1_k[1, ])
        log_w_tilde[1, ] <- log_w_1 - logsumexp(log_w_1)

        # t = 2, ..., T
        for (t in 2:Tt) {
            A <- systematic_resample(log_w_tilde[t - 1, ], K - 1)
            theta_t1_k <- numeric(K)
            theta_t1_k[1:(K - 1)] <- rnorm(K - 1,
                                            mean = theta_1_k[t - 1, A] + theta2[t - 1],
                                            sd = sd_W1)

            theta_t1_k[K] <- theta1[t]
            log_as <- log_w_tilde[t - 1, ] +
                dnorm(theta1[t],
                    mean = theta_1_k[t - 1, ] + theta2[t - 1],
                    sd = sd_W1,
                    log = TRUE
                )
            log_as <- log_as - max(log_as)
            a_ref  <- sample_one_from_logw(log_as)

            A_hist[t, ] <- c(A, a_ref)
            theta_1_k[t, ] <- theta_t1_k

            log_w_t <- log_p_yt(y[t], theta_1_k[t, ])
            log_w_tilde[t, ] <- log_w_t - logsumexp(log_w_t)
        }

        ess_smc[n] <- 1 / sum(exp(2 * (log_w_tilde[Tt, ])))

        # Backward tracking for theta1
        k_final <- sample_one_from_logw(log_w_tilde[Tt, ])
        theta1[Tt] <- theta_1_k[Tt, k_final]

        b <- k_final
        for (t in (Tt - 1):1) {
            b <- A_hist[t + 1, b]
            theta1[t] <- theta_1_k[t, b]
        }

        # (theta_02, theta2) jointly via extended block
        build2 <- chan_smoothing_theta2(theta1, phi1, phi2, mu_02, sigma2_02, theta_01)
        draw2 <- chan_sample_from_build(build2, Ttp1)
        theta_02 <- draw2[1]
        theta2 <- draw2[-1]

        # Store the results
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
         ess_smc = ess_smc)
}
