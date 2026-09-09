# Poisson - 2nd Order Polynomial Dynamic Model
# Gibbs with Collapsed Samplers for W1 AND W2
#
# theta_01/theta1: EXTENDED-state joint block ((T+1)-dimensional), sampled
#                  together via chan_smoothing_theta1 - still needs SIR (the
#                  Poisson likelihood on theta1[1..T] makes it only a
#                  Laplace/IRLS approximation; theta_01 itself needs no
#                  approximation, no likelihood on it).
# theta_02/theta2: EXTENDED-state joint block ((T+1)-dimensional), sampled
#                  together via chan_smoothing_theta2 - exact, no SIR
#                  needed (no Poisson likelihood involved).
#
# This file only defines sir_collapsed_r(). It has no side effects (no data
# loading, no plotting) -- see tests/ for scripts that call it, either
# standalone (test_prototype_R.R) or against the Rcpp port (test_cpp.R,
# test_R_vs_cpp.R). Signature mirrors sir_collapsed_cpp() 1:1.
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

sir_collapsed_r <- function(y, R_prerun, N,
                             mu_01, sigma2_01, mu_02, sigma2_02,
                             nu_01, eta_01, nu_02, eta_02,
                             theta1, theta2, theta1_tilde,
                             theta_01, theta_02, W1, W2,
                             M_is_lik, M_sir_theta1, M_irls_max, tol,
                             verbose = TRUE, print_every = 1000) {

    printf <- function(...) cat(paste(sprintf(...), "\n"))
    Tt <- length(y)
    Ttp1 <- Tt + 1
    theta1_star <- theta1 # internal name (matches the algorithm's own notation)
    phi1 <- 1 / W1
    phi2 <- 1 / W2

    #####
    # Static Chan-method structures

    # T-dimensional, ANCHORED: only needed for log_det_K0 and for the
    # dedicated W2-marginal-likelihood block below.
    log_det_K0 <- chan_log_det_K0(Tt)
    w2lik_obj <- chan_build_static_objects(Tt)
    P_w2lik_matrix  <- w2lik_obj$K0
    Ch_w2lik_factor <- w2lik_obj$Ch0_factor
    idx_diag        <- w2lik_obj$idx_diag
    idx_sub         <- w2lik_obj$idx_sub
    main_diag_base  <- w2lik_obj$main_diag_base

    # (T+1)-dimensional EXTENDED blocks
    chan_smoothing_theta2 <- make_chan_theta2_smoother_ext(Ttp1)
    chan_smoothing_theta1 <- make_chan_theta1_smoother_ext(Ttp1)

    #####
    # Nested helper functions (closures over Tt, Ttp1, mu_01, sigma2_01,
    # chan_smoothing_theta1, log_det_K0, P_w2lik_matrix, etc. above)

    run_irls <- function(theta1_tilde, theta2, theta_02, phi1, y, tol, M_irls_max) {
        for (j in 1:M_irls_max) {
            f_t   <- exp(-theta1_tilde)
            phi_V <- 1 / f_t
            z_t   <- theta1_tilde + f_t * y - 1

            res <- chan_smoothing_theta1(z_t, phi_V, phi1, mu_01, sigma2_01, theta_02, theta2)
            theta1_tilde_new <- res$theta1_hat[-1]
            if (any(!is.finite(theta1_tilde_new))) break

            if (max(abs(theta1_tilde_new - theta1_tilde)) < tol) {
                theta1_tilde <- theta1_tilde_new
                break
            }
            theta1_tilde <- theta1_tilde_new
        }
        return(list(res = res, theta1_tilde = theta1_tilde, itr = j,
             f_t = f_t, phi_V = phi_V))
    }

    is_log_lik <- function(irls_res, y, phi1, theta2, theta_02, M_is_lik) {
        res     <- irls_res$res
        eta_hat <- res$theta1_hat
        ch      <- res$ch
        W1 <- 1 / phi1

        log_det_H <- 2 * as.numeric(determinant(ch, logarithm = TRUE, sqrt = TRUE)$modulus)
        th_lag2_fixed <- c(theta_02, theta2[-Tt])
        log_norm_H <- -Ttp1 / 2 * log(2 * pi) + 0.5 * log_det_H

        log_w <- numeric(M_is_lik)
        u_mat <- matrix(rnorm(M_is_lik * Ttp1), nrow = M_is_lik)
        d_ch  <- Matrix::diag(ch)

        for (i in 1:M_is_lik) {
            u <- u_mat[i, ]
            w <- u / sqrt(d_ch)
            x <- as.vector(solve(ch, w, system = "Lt"))
            draw_i <- eta_hat + x
            theta_01_i <- draw_i[1]
            th         <- draw_i[-1]

            log_py <- sum(y * th - exp(th))
            log_prior_01 <- -0.5 * log(2 * pi * sigma2_01) - (theta_01_i - mu_01)^2 / (2 * sigma2_01)
            th_lag1 <- c(theta_01_i, th[-Tt])
            eps     <- th - th_lag1 - th_lag2_fixed
            log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
            log_q <- log_norm_H - 0.5 * sum(u^2)

            log_w[i] <- log_py + log_prior_01 + log_prior_th - log_q
        }

        log_lik <- logsumexp(log_w) - log(M_is_lik)
        w_norm  <- exp(log_w - logsumexp(log_w))
        ess_is  <- 1 / sum(w_norm^2)

        return(list(log_lik = log_lik, ess_is = ess_is))
    }

    sir_theta1 <- function(irls_res, y, phi1, theta2, theta_02, M_sir_theta1) {
        eta_hat <- irls_res$res$theta1_hat
        ch      <- irls_res$res$ch
        W1 <- 1 / phi1

        log_det_H     <- 2 * as.numeric(determinant(ch, logarithm = TRUE, sqrt = TRUE)$modulus)
        th_lag2_fixed <- c(theta_02, theta2[-Tt])
        log_norm_H    <- -Ttp1 / 2 * log(2 * pi) + 0.5 * log_det_H

        log_w  <- numeric(M_sir_theta1)
        draws  <- matrix(0, M_sir_theta1, Ttp1)
        d_ch   <- Matrix::diag(ch)

        for (i in 1:M_sir_theta1) {
            u  <- rnorm(Ttp1)
            w <- u / sqrt(d_ch)
            x  <- as.vector(Matrix::solve(ch, w, system = "Lt"))
            draw_i <- eta_hat + x
            draws[i, ] <- draw_i
            theta_01_i <- draw_i[1]
            th         <- draw_i[-1]

            log_py       <- sum(y * th - exp(th))
            log_prior_01 <- -0.5 * log(2 * pi * sigma2_01) - (theta_01_i - mu_01)^2 / (2 * sigma2_01)
            th_lag1      <- c(theta_01_i, th[-Tt])
            eps          <- th - th_lag1 - th_lag2_fixed
            log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
            log_q        <- log_norm_H - 0.5 * sum(u^2)

            log_w[i] <- log_py + log_prior_01 + log_prior_th - log_q
        }

        w_norm  <- exp(log_w - logsumexp(log_w))
        ess_sir <- 1 / sum(w_norm^2)
        idx     <- sample_one_from_logw(log_w)

        list(theta_01 = draws[idx, 1], theta1 = draws[idx, -1], ess = ess_sir)
    }

    log_marginal_lik_w2 <- function(theta1, phi1, phi2, theta_02) {
        z <- diff(theta1)
        diag_obs <- c(rep(phi1, Tt - 1), 0)
        P_w2lik_matrix@x[idx_diag] <- (main_diag_base * phi2) + diag_obs
        P_w2lik_matrix@x[idx_sub]  <- -phi2
        ch <- update(Ch_w2lik_factor, P_w2lik_matrix)

        b <- numeric(Tt)
        b[1:(Tt - 1)] <- z * phi1
        b[1] <- b[1] + theta_02 * phi2
        theta2_hat <- as.numeric(Matrix::solve(ch, b, system = "A"))

        log_pz <- -0.5 * (Tt - 1) * log(2 * pi / phi1) - 0.5 * phi1 * sum((z - theta2_hat[1:(Tt - 1)])^2)

        diffs2 <- theta2_hat - c(theta_02, theta2_hat[-Tt])
        log_p_theta2 <- -0.5 * Tt * log(2 * pi / phi2) + 0.5 * log_det_K0 - 0.5 * phi2 * sum(diffs2^2)

        log_det_H <- 2 * as.numeric(determinant(ch, logarithm = TRUE, sqrt = TRUE)$modulus)
        log_q <- -0.5 * Tt * log(2 * pi) + 0.5 * log_det_H

        list(log_lik = log_pz + log_p_theta2 - log_q)
    }

    calibrate_ce_gamma <- function(samples) {
        m1    <- mean(samples)
        m_log <- mean(log(samples))
        rhs   <- log(m1) - m_log

        c_hat <- if (rhs <= 0.5772) {
            (3 - rhs + sqrt((rhs - 3)^2 + 24 * rhs)) / (12 * rhs)
        } else {
            1 / rhs
        }

        for (k in 1:100) {
            f_val <- log(c_hat) - digamma(c_hat) - rhs
            fp    <- 1 / c_hat - trigamma(c_hat)
            c_new <- c_hat - f_val / fp
            if (!is.finite(c_new) || c_new <= 0) break
            if (abs(c_new - c_hat) < 1e-10) { c_hat <- c_new; break }
            c_hat <- c_new
        }
        d_hat <- c_hat / m1
        return(list(shape = c_hat, rate = d_hat))
    }

    mh_w1_collapsed <- function(phi1_cur, log_lik_cur, irls_cur,
                                ce_params, nu_01, eta_01,
                                theta2, theta_02,
                                theta1_tilde, y, tol, M_irls_max, M_is_lik) {

        phi1_prop <- rgamma(1, shape = ce_params$shape, rate = ce_params$rate)

        irls_prop    <- run_irls(theta1_tilde, theta2, theta_02,
                                 phi1_prop, y, tol, M_irls_max)
        res_lik_prop <- is_log_lik(irls_prop, y, phi1_prop, theta2, theta_02, M_is_lik)
        log_lik_prop <- res_lik_prop$log_lik

        log_prior <- function(phi) dgamma(phi, shape = nu_01, rate = eta_01, log = TRUE)
        log_prop  <- function(phi) dgamma(phi, shape = ce_params$shape, rate = ce_params$rate, log = TRUE)

        log_alpha <- (log_lik_prop + log_prior(phi1_prop) + log_prop(phi1_cur)) -
            (log_lik_cur  + log_prior(phi1_cur)  + log_prop(phi1_prop))

        if (log(runif(1)) < log_alpha) {
            return(list(phi1 = phi1_prop, log_lik = log_lik_prop,
                        irls_res = irls_prop, accepted = TRUE))
        } else {
            return(list(phi1 = phi1_cur, log_lik = log_lik_cur,
                        irls_res = irls_cur, accepted = FALSE))
        }
    }

    mh_w2_collapsed <- function(phi2_cur, log_lik2_cur,
                                ce2_params, nu_02, eta_02,
                                theta1_star, phi1, theta_02) {

        phi2_prop <- rgamma(1, shape = ce2_params$shape, rate = ce2_params$rate)

        res_prop <- log_marginal_lik_w2(theta1_star, phi1, phi2_prop, theta_02)
        log_lik2_prop <- res_prop$log_lik

        log_prior <- function(phi) dgamma(phi, shape = nu_02, rate = eta_02, log = TRUE)
        log_prop  <- function(phi) dgamma(phi, shape = ce2_params$shape, rate = ce2_params$rate, log = TRUE)

        log_alpha <- (log_lik2_prop + log_prior(phi2_prop) + log_prop(phi2_cur)) -
            (log_lik2_cur  + log_prior(phi2_cur)  + log_prop(phi2_prop))

        if (!is.finite(log_alpha)) {
            return(list(phi2 = phi2_cur, log_lik2 = log_lik2_cur, accepted = FALSE))
        }

        if (log(runif(1)) < log_alpha) {
            return(list(phi2 = phi2_prop, log_lik2 = log_lik2_prop, accepted = TRUE))
        } else {
            return(list(phi2 = phi2_cur, log_lik2 = log_lik2_cur, accepted = FALSE))
        }
    }

    #####
    # Histories
    W1_hist <- numeric(N)
    W2_hist <- numeric(N)
    theta_01_hist <- numeric(N)
    theta_02_hist <- numeric(N)
    theta1_hist <- matrix(0, N, Tt)
    theta2_hist <- matrix(0, N, Tt)
    accepted_hist  <- logical(N)
    accepted2_hist <- logical(N)
    itr_irls <- numeric(N)
    ess_is_hist <- numeric(N)
    ess_sir_hist <- numeric(N)

    #####
    # PRE-RUN (simple, non-collapsed Gibbs, purely to calibrate the CE
    # Gamma proposals for phi1 and phi2)
    if (verbose) printf("Starting pre-run (%d iterations)...", R_prerun)
    phi1_prerun <- numeric(R_prerun)
    phi2_prerun <- numeric(R_prerun)
    time_prerun <- proc.time()

    for (r in 1:R_prerun) {

        phi1 <- gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1_star, theta_02, theta2, Tt)
        W1 <- 1 / phi1
        phi1_prerun[r] <- phi1

        phi2 <- gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt)
        W2 <- 1 / phi2
        phi2_prerun[r] <- phi2

        irls_out <- run_irls(theta1_tilde, theta2, theta_02, phi1, y, tol, M_irls_max)
        theta1_tilde <- irls_out$theta1_tilde
        res_sir <- sir_theta1(irls_out, y, phi1, theta2, theta_02, M_sir_theta1)
        theta_01    <- res_sir$theta_01
        theta1_star <- res_sir$theta1

        build2 <- chan_smoothing_theta2(theta1_star, phi1, phi2, mu_02, sigma2_02, theta_01)
        draw2  <- chan_sample_from_build(build2, Ttp1)
        theta_02 <- draw2[1]
        theta2   <- draw2[-1]
    }

    if (verbose) {
        elapsed_prerun <- (proc.time() - time_prerun)[[1]]
        printf("Pre-run done in %.0f s", elapsed_prerun)
    }

    #####
    # CE Calibration
    ce_params  <- calibrate_ce_gamma(phi1_prerun)
    ce2_params <- calibrate_ce_gamma(phi2_prerun)

    if (verbose) {
        printf("CE Gamma proposal for phi1: shape = %.4f, rate = %.4f (mean phi1 = %.6f, mean W1 = %.6f)",
               ce_params$shape, ce_params$rate,
               ce_params$shape / ce_params$rate,
               ce_params$rate / ce_params$shape)
        printf("CE Gamma proposal for phi2: shape = %.4f, rate = %.4f (mean phi2 = %.6f, mean W2 = %.6f)",
               ce2_params$shape, ce2_params$rate,
               ce2_params$shape / ce2_params$rate,
               ce2_params$rate / ce2_params$shape)
    }

    #####
    # Main Gibbs loop
    irls_cur <- run_irls(theta1_tilde, theta2, theta_02, phi1, y, tol, M_irls_max)
    theta1_tilde <- irls_cur$theta1_tilde

    res_lik <- is_log_lik(irls_cur, y, phi1, theta2, theta_02, M_is_lik)
    log_lik_cur <- res_lik$log_lik

    if (verbose) printf("Starting main MCMC (%d iterations)...", N)
    start_time <- proc.time()

    for (n in 1:N) {

        if (verbose && (n %% print_every == 0)) {
            elapsed <- (proc.time() - start_time)[[1]]
            acc_rate  <- mean(accepted_hist[1:(n - 1)])
            acc_rate2 <- mean(accepted2_hist[1:(n - 1)])
            printf("Iter %d / %d | Elapsed: %.0f s | W1 accept: %.2f | W2 accept: %.2f",
                   n, N, elapsed, acc_rate, acc_rate2)
        }

        # Collapsed MH for W1
        mh_res <- mh_w1_collapsed(
            phi1_cur = phi1, log_lik_cur = log_lik_cur, irls_cur = irls_cur,
            ce_params = ce_params, nu_01 = nu_01, eta_01 = eta_01,
            theta2 = theta2, theta_02 = theta_02, theta1_tilde = theta1_tilde,
            y = y, tol = tol, M_irls_max = M_irls_max, M_is_lik = M_is_lik
        )
        phi1 <- mh_res$phi1
        W1 <- 1 / phi1
        log_lik_cur <- mh_res$log_lik
        irls_cur <- mh_res$irls_res
        itr_irls[n] <- irls_cur$itr
        accepted_hist[n] <- mh_res$accepted

        # Collapsed MH for W2
        res_lik2_cur <- log_marginal_lik_w2(theta1_star, phi1, phi2, theta_02)
        log_lik2_cur <- res_lik2_cur$log_lik

        mh2_res <- mh_w2_collapsed(
            phi2_cur = phi2, log_lik2_cur = log_lik2_cur, ce2_params = ce2_params,
            nu_02 = nu_02, eta_02 = eta_02,
            theta1_star = theta1_star, phi1 = phi1, theta_02 = theta_02
        )
        phi2 <- mh2_res$phi2
        W2 <- 1 / phi2
        accepted2_hist[n] <- mh2_res$accepted

        # (theta_02, theta2) jointly via extended block, BEFORE (theta_01, theta1)
        build2 <- chan_smoothing_theta2(theta1_star, phi1, phi2, mu_02, sigma2_02, theta_01)
        draw2  <- chan_sample_from_build(build2, Ttp1)
        theta_02 <- draw2[1]
        theta2   <- draw2[-1]

        # (theta_01, theta1) jointly (SIR, extended, given accepted phi1 and
        # the fresh theta2)
        res_sir <- sir_theta1(irls_cur, y, phi1, theta2, theta_02, M_sir_theta1)
        theta_01     <- res_sir$theta_01
        theta1_star  <- res_sir$theta1
        theta1_tilde <- irls_cur$theta1_tilde
        ess_sir_hist[n] <- res_sir$ess

        # Update irls_cur and log_lik_cur for next iteration
        irls_cur <- run_irls(theta1_tilde, theta2, theta_02, phi1, y, tol, M_irls_max)
        theta1_tilde <- irls_cur$theta1_tilde

        res_lik <- is_log_lik(irls_cur, y, phi1, theta2, theta_02, M_is_lik)
        log_lik_cur <- res_lik$log_lik
        ess_is_hist[n] <- res_lik$ess_is

        theta_01_hist[n] <- theta_01
        theta_02_hist[n] <- theta_02
        W1_hist[n] <- W1
        W2_hist[n] <- W2
        theta1_hist[n, ] <- theta1_star
        theta2_hist[n, ] <- theta2
    }

    list(theta1_hist = theta1_hist,
         theta2_hist = theta2_hist,
         theta_01_hist = theta_01_hist,
         theta_02_hist = theta_02_hist,
         W1_hist = W1_hist,
         W2_hist = W2_hist,
         accepted_hist = accepted_hist,
         accepted2_hist = accepted2_hist,
         itr_irls = itr_irls,
         ess_is_hist = ess_is_hist,
         ess_sir_hist = ess_sir_hist,
         ce1_shape = ce_params$shape,
         ce1_rate = ce_params$rate,
         ce2_shape = ce2_params$shape,
         ce2_rate = ce2_params$rate)
}
