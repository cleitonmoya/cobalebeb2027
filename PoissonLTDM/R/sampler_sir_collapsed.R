# Poisson - 2nd Order Polynomial Dynamic Model
#
# MCMC: Gibbs with SIR for theta1 and Collapsed Samplers for W1 and W2
#
# Strategy:
#   - Pre-run: standard Gibbs (R_prerun iterations) to collect phi1 and phi2 samples
#   - CE calibration: fit Gamma(c_hat, d_hat) proposals for both phi1 and phi2
#   - Gibbs sampling:
#       * theta_02, theta_01: Conjugated Normal
#       * W1: Collapsed MH with CE proposal marginalizing theta1;
#             Integrated likelihod: IS with Laplace approx.
#       * W2: Collapsed MH with CE proposal marginalizing theta2;
#             Integrate likelihood: Gaussian    
#       * theta1: SIR + Chan sampler
#       * theta2: Chan sampler
#
# Author: Cleiton Moya de Almeida


# IRLS: build Laplace approximation around theta1_tilde
run_irls <- function(theta1_tilde, theta2, theta_01, theta_02, phi1,
                     y, tol, M_irls_max, chan_smoothing_theta1) {
    for (j in 1:M_irls_max) {
        f_t   <- exp(-theta1_tilde)
        phi_V <- 1 / f_t
        z_t   <- theta1_tilde + f_t * y - 1   # Poisson pseudo-observation

        theta1_build <- chan_smoothing_theta1(z_t, phi_V, phi1, theta_01, theta_02, theta2)
        theta1_tilde_new <- theta1_build$theta1_hat
        if (any(!is.finite(theta1_tilde_new))) break

        if (max(abs(theta1_tilde_new - theta1_tilde)) < tol) {
            theta1_tilde <- theta1_tilde_new
            break
        }
        theta1_tilde <- theta1_tilde_new
    }
  
    return(list(theta1_build = theta1_build, 
                theta1_tilde = theta1_tilde, 
                itr = j,
                f_t = f_t, 
                phi_V = phi_V))
}


# IS estimator of log p(y | phi1, theta2)  [integrated over theta1]
# Approximate: the Poisson likelihood makes p(theta1|y,phi1,theta2) non-Gaussian,
# so this uses IS around the Laplace mode (computed via IRLS).
is_log_lik <- function(irls_res, y, phi1, theta2, theta_01, theta_02, M_is) {
    
    theta1_build <- irls_res$theta1_build
    eta_hat <- theta1_build$theta1_hat
    ch <- theta1_build$ch
    
    W1 <- 1/phi1

    log_det_H <- 2 * as.numeric(Matrix::determinant(ch, logarithm = TRUE)$modulus)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])
    log_norm_H <- -Tt / 2 * log(2 * pi) + 0.5 * log_det_H

    log_w <- numeric(M_is)
    u_mat <- matrix(rnorm(M_is * Tt), nrow = M_is)
    d_ch  <- Matrix::diag(ch)

    for (i in 1:M_is) {
        u <- u_mat[i, ]
        w <- u / sqrt(d_ch)
        x <- as.vector(Matrix::solve(ch, w, system = "Lt"))
        th <- eta_hat + x

        log_py <- sum(y * th - exp(th))
        th_lag1 <- c(theta_01, th[-Tt])
        eps     <- th - th_lag1 - th_lag2_fixed
        log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
        log_q <- log_norm_H - 0.5 * sum(u^2)

        log_w[i] <- log_py + log_prior_th - log_q
    }

    log_lik <- logsumexp(log_w) - log(M_is)
    w_norm  <- exp(log_w - logsumexp(log_w))
    ess_is  <- 1 / sum(w_norm^2)

    return(list(log_lik = log_lik, ess_is = ess_is))
}



# EXACT integrated likelihood log p(z | phi2, phi1, theta_02) [integrated over theta2]
# Exact: z|theta2 is Gaussian valuated via the
# Chan sparse Cholesky factor at the conditional mode theta2_hat:
#
#   log p(z|phi2) = log p(z|theta2_hat) + log p(theta2_hat|phi2) - log q(theta2_hat)
#
# where q(.) ~ N(theta2_hat, H^{-1})
log_marginal_lik_w2 <- function(theta1, phi1, phi2, theta_02, 
                                chan_smoothing_theta2, log_det_K0) {
  
    build <- chan_smoothing_theta2(theta1, phi1, phi2, theta_02)
    theta2_hat <- build$theta2_hat
    ch <- build$ch
    z <- build$z

    log_pz <- -0.5*(Tt-1)*log(2*pi/phi1) - 0.5*phi1*sum((z - theta2_hat[1:(Tt-1)])^2)

    diffs2 <- theta2_hat - c(theta_02, theta2_hat[-Tt])
    log_p_theta2 <- -0.5*Tt*log(2*pi/phi2) + 0.5*log_det_K0 - 0.5*phi2*sum(diffs2^2)

    log_det_H <- 2 * as.numeric(Matrix::determinant(ch, logarithm = TRUE)$modulus)
    log_q <- -0.5*Tt*log(2*pi) + 0.5*log_det_H

    list(log_lik = log_pz + log_p_theta2 - log_q, build = build)
}


# Cross-entropy calibration: fit Gamma(c, d) to samples via moment matching
calibrate_ce_gamma <- function(samples) {
    m1 <- mean(samples)
    m_log <- mean(log(samples))
    rhs <- log(m1) - m_log

    c_hat <- if (rhs <= 0.5772) {
        (3 - rhs + sqrt((rhs - 3)^2 + 24 * rhs)) / (12 * rhs)
    } else {
        1 / rhs
    }

    for (k in 1:100) {
        f_val <- log(c_hat) - digamma(c_hat) - rhs
        fp <- 1 / c_hat - trigamma(c_hat)
        c_new <- c_hat - f_val / fp
        if (!is.finite(c_new) || c_new <= 0) break
        if (abs(c_new - c_hat) < 1e-10) { c_hat <- c_new; break }
        c_hat <- c_new
    }
    d_hat <- c_hat / m1
    return(list(shape = c_hat, rate = d_hat))
}


# Collapsed MH step for W1 (independence chain with CE proposal, approximate IS likelihood)
mh_phi1_collapsed <- function(phi1_cur, log_lik1_cur, irls_cur,
                            ce1_params, nu_01, eta_01,
                            theta2, theta_01, theta_02,
                            theta1_tilde, y, tol, M_irls_max, M_is,
                            chan_smoothing_theta1) {

    phi1_prop <- rgamma(1, shape = ce1_params$shape, rate = ce1_params$rate)

    irls_prop    <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                             phi1_prop, y, tol, M_irls_max, chan_smoothing_theta1)
    res_lik_prop <- is_log_lik(irls_prop, y, phi1_prop, theta2, theta_01, theta_02, M_is)
    log_lik_prop <- res_lik_prop$log_lik

    log_prior <- function(phi) dgamma(phi, shape = nu_01, rate = eta_01, log = TRUE)
    log_prop  <- function(phi) dgamma(phi, shape = ce1_params$shape, rate = ce1_params$rate, log = TRUE)

    log_alpha <- (log_lik_prop + log_prior(phi1_prop) + log_prop(phi1_cur)) -
        (log_lik1_cur  + log_prior(phi1_cur)  + log_prop(phi1_prop))

    if (log(runif(1)) < log_alpha) {
        return(list(phi1 = phi1_prop, log_lik = log_lik_prop,
                    irls_res = irls_prop, accepted = TRUE))
    } else {
        return(list(phi1 = phi1_cur, log_lik = log_lik1_cur,
                    irls_res = irls_cur, accepted = FALSE))
    }
}


# Collapsed MH step for W2 (independence chain with CE proposal, EXACT likelihood)
mh_phi2_collapsed <- function(phi2_cur, log_lik2_cur, theta2_build_cur,
                            ce2_params, nu_02, eta_02,
                            theta1, phi1, theta_02,
                            chan_smoothing_theta2,
                            log_det_K0) {

    phi2_prop <- rgamma(1, shape = ce2_params$shape, rate = ce2_params$rate)

    res_prop <- log_marginal_lik_w2(theta1, phi1, phi2_prop, theta_02, 
                                    chan_smoothing_theta2, log_det_K0)
    log_lik2_prop <- res_prop$log_lik

    log_prior <- function(phi) dgamma(phi, shape = nu_02, rate = eta_02, log = TRUE)
    log_prop  <- function(phi) dgamma(phi, shape = ce2_params$shape, rate = ce2_params$rate, log = TRUE)

    log_alpha <- (log_lik2_prop + log_prior(phi2_prop) + log_prop(phi2_cur)) -
        (log_lik2_cur  + log_prior(phi2_cur)  + log_prop(phi2_prop))

    if (!is.finite(log_alpha)) {
        return(list(phi2 = phi2_cur, log_lik2 = log_lik2_cur,
                    build = theta2_build_cur, accepted = FALSE))
    }

    if (log(runif(1)) < log_alpha) {
        return(list(phi2 = phi2_prop, log_lik2 = log_lik2_prop,
                    build = res_prop$build, accepted = TRUE))
    } else {
        return(list(phi2 = phi2_cur, log_lik2 = log_lik2_cur,
                    build = theta2_build_cur, accepted = FALSE))
    }
}


# Sampling Importance Resampling for theta1
sir_theta1 <- function(irls_res, y, phi1, theta2, theta_01, theta_02, M_sir) {
    
    eta_hat <- irls_res$theta1_build$theta1_hat
    ch <- irls_res$theta1_build$ch
    W1 <- 1/phi1

    log_det_H <- 2*as.numeric(Matrix::determinant(ch, logarithm = TRUE)$modulus)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])
    log_norm_H <- -(Tt/2)*log(2*pi) + 0.5*log_det_H

    log_w <- numeric(M_sir)
    draws <- matrix(0, M_sir, Tt)
    d_ch <- Matrix::diag(ch)

    for (i in 1:M_sir) {
        u <- rnorm(Tt)
        w <- u / sqrt(d_ch)
        x <- as.vector(Matrix::solve(ch, w, system = "Lt"))
        th <- eta_hat + x
        draws[i, ] <- th

        log_py <- sum(y * th - exp(th))
        th_lag1 <- c(theta_01, th[-Tt])
        eps <- th - th_lag1 - th_lag2_fixed
        log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
        log_q <- log_norm_H - 0.5 * sum(u^2)

        log_w[i] <- log_py + log_prior_th - log_q
    }

    w_norm <- exp(log_w - logsumexp(log_w))
    ess_sir <- 1 / sum(w_norm^2)
    idx <- sample.int(M_sir, size = 1, prob = w_norm)

    list(theta1 = draws[idx, ], ess = ess_sir)
}

# PRE-RUN: standard Gibbs to collect phi1 and phi2 samples for CE calibration
collapsed_prerun <- function(R_prerun, y, Tt,
                             M_irls_max, tol,
                             mu_01, sigma2_01, mu_02, sigma2_02,
                             nu_01, eta_01, nu_02, eta_02,
                             W1, W2, theta_01, theta_02, theta1, theta2, 
                             theta1_tilde,
                             chan_smoothing_theta1, chan_smoothing_theta2) {
  
  phi1 <- 1/W1
  phi2 <- 1/W2
  
  phi1_prerun <- numeric(R_prerun)
  phi2_prerun <- numeric(R_prerun)
  
  for (r in 1:R_prerun) {
    
    # Sample theta_02 (conjugated Normal)
    theta_02 <- gibbs_sample_theta02(mu_02, sigma2_02, theta_01, theta1[1],
                                     theta2[1], W1, W2)
    
    # Sample theta_01 (conjugated Normal)
    theta_01 <- gibbs_sample_theta01(mu_01, sigma2_01, theta1[1],
                                     theta_02, W1)
    
    # Sample phi1 (conjugated Gamma)
    phi1 <- gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1,
                              theta_02, theta2, Tt)
    phi1_prerun[r] <- phi1
    W1 <- 1/phi1
    
    
    # Sample phi2 (conjugated Gamma)
    phi2 <- gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt)
    phi2_prerun[r] <- phi2
    W2 <- 1/phi2
    
    
    # Sample theta1 (IRLS + Chan)
    irls_res <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                         phi1, y, tol, M_irls_max, chan_smoothing_theta1)
    
    theta1_tilde <- irls_res$theta1_tilde
    theta1  <- chan_sample_from_build(irls_res$theta1_build, Tt)
    
    # Sample theta2 (Chan)
    theta2_build <- chan_smoothing_theta2(theta1, phi1, phi2, theta_02)
    theta2 <- chan_sample_from_build(theta2_build, Tt)
  }
  
  return(list(phi1_prerun = phi1_prerun, 
              phi2_prerun = phi2_prerun,
              theta_01 = theta_01,
              theta_02 = theta_02,
              W1 = W1,
              W2 = W2,
              theta1 = theta1,
              theta2 = theta2,
              irls_res = irls_res))
}


sample_sir_collapsed <- function(y, N, burnin, R_prerun, 
                                 M_is, M_sir, M_irls_max, tol,
                                 mu_01, sigma2_01, mu_02, sigma2_02,
                                 nu_01, eta_01, nu_02, eta_02,
                                 W1, W2, theta_01, theta_02, theta1, theta2,
                                 theta1_tilde) {
  Tt <- length(y)
  
  # Prepare Chan static objects
  chan_smoothing_theta1 <- make_chan_theta1_smoother(Tt)
  chan_smoothing_theta2 <- make_chan_theta2_smoother(Tt)
  log_det_K0 <- chan_log_det_K0(Tt)
  
  # Auxiliary variables
  theta_01_hist <- numeric(N)
  theta_02_hist <- numeric(N)
  W1_hist <- numeric(N)
  W2_hist <- numeric(N)
  theta1_hist <- matrix(0, N, Tt)
  theta2_hist <- matrix(0, N, Tt)
  accepted1_hist <- logical(N)    # W1 MH acceptance
  accepted2_hist <- logical(N)    # W2 MH acceptance
  itr_irls_hist <- numeric(N)
  ess_is_hist  <- numeric(N)      # IS of W1's integrated likelihood
  ess_sir_hist <- numeric(N)      # SIR of theta1
  
  # Pre-run stage (to calibrate CE parameters)
  prerun_res <- collapsed_prerun(R_prerun, y, Tt,
                                 M_irls_max, tol,
                                 mu_01, sigma2_01, mu_02, sigma2_02,
                                 nu_01, eta_01, nu_02, eta_02,
                                 W1, W2, theta_01, theta_02, theta1, theta2, 
                                 theta1_tilde,
                                 chan_smoothing_theta1, chan_smoothing_theta2)
  
  # Recover pre-run parameters
  phi1_prerun  <- prerun_res$phi1_prerun
  phi2_prerun  <- prerun_res$phi2_prerun
  theta_01     <- prerun_res$theta_01
  theta_02     <- prerun_res$theta_02
  W1           <- prerun_res$W1
  W2           <- prerun_res$W2
  theta1       <- prerun_res$theta1
  theta2       <- prerun_res$theta2
  irls_cur     <- prerun_res$irls_res
  theta1_tilde <- prerun_res$irls_res$theta1_tilde
  phi1 <- 1/W1
  phi2 <- 1/W2
  
  # Cross entropy parameters calibration
  ce1_params <- calibrate_ce_gamma(phi1_prerun)
  ce2_params <- calibrate_ce_gamma(phi2_prerun)
  
  # Integrated likelihood for W1
  res_lik1 <- is_log_lik(irls_cur, y, phi1, theta2, theta_01, theta_02, M_is)
  log_lik1_cur <- res_lik1$log_lik
  
  # Integrated likelihood for W2
  res_lik2 <- log_marginal_lik_w2(theta1, phi1, phi2, theta_02, 
                                  chan_smoothing_theta2, log_det_K0)
  log_lik2_cur <- res_lik2$log_lik
  theta2_build_cur <- res_lik2$build
  
  
  # Gibbs sampling
  for (n in 1:N) {
    
    # Sample theta_02 (conjugated Normal)
    sigma2_02_bar <- (1 / sigma2_02 + 1 / W1 + 1 / W2)^(-1)
    mu_02_bar <- sigma2_02_bar * ((theta1[1] - theta_01) / W1 +
                                    theta2[1] / W2 + mu_02 / sigma2_02)
    theta_02 <- rnorm(1, mean = mu_02_bar, sd = sqrt(sigma2_02_bar))
    
    
    # Sample theta_01 (conjugated Normal)
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01 / sigma2_01 +
                                    (theta1[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))
    
    
    # Collapsed MH for phi1 (approximate integrated likelihood, marginalizes theta1)
    mh1_res <- mh_phi1_collapsed(phi1, log_lik1_cur, irls_cur,
                                 ce1_params, nu_01, eta_01,
                                 theta2, theta_01, theta_02,
                                 theta1_tilde, y, tol, M_irls_max, M_is,
                                 chan_smoothing_theta1)
    
    phi1 <- mh1_res$phi1
    W1 <- 1/phi1
    log_lik1_cur <- mh1_res$log_lik
    irls_cur <- mh1_res$irls_res
    
    
    # Collapsed MH for phi2 (EXACT integrated likelihood, marginalizes theta2)
    # Recompute log_lik2_cur first: theta_02 and phi1 both may have changed
    # since the last time it was updated (theta_01/theta_02 sampled above,
    # phi1 updated by the W1 MH step just now).
    res_lik2_cur <- log_marginal_lik_w2(theta1, phi1, phi2, theta_02, 
                                        chan_smoothing_theta2,
                                        log_det_K0)
    log_lik2_cur <- res_lik2_cur$log_lik
    theta2_build_cur <- res_lik2_cur$build
    
    mh2_res <- mh_phi2_collapsed(phi2, log_lik2_cur, theta2_build_cur,
                                 ce2_params, nu_02, eta_02,
                                 theta1, phi1, theta_02,
                                 chan_smoothing_theta2, log_det_K0)
    
    phi2 <- mh2_res$phi2
    W2 <- 1/phi2
    
    
    # Sample theta2 (Chan sampler)
    # must happen BEFORE SIR of theta1, so SIR sees the fresh theta2
    theta2_build <- chan_smoothing_theta2(theta1, phi1, phi2, theta_02)
    theta2 <- chan_sample_from_build(theta2_build, Tt)
    
    
    # Sample theta1* via SIR given accepted phi1 and the fresh theta2
    res_sir <- sir_theta1(irls_cur, y, phi1, theta2, theta_01, theta_02, M_sir)
    theta1 <- res_sir$theta1
    theta1_tilde <- irls_cur$theta1_tilde
   
    
    # Update irls_cur and log_lik1_cur for next iteration
    # (theta2, theta_01, theta_02 all changed in this iteration)
    irls_cur <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                         phi1, y, tol, M_irls_max, chan_smoothing_theta1)
    theta1_tilde <- irls_cur$theta1_tilde
    
    res_lik1 <- is_log_lik(irls_cur, y, phi1, theta2, theta_01, theta_02, M_is)
    log_lik1_cur <- res_lik1$log_lik
    
    
    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1
    theta2_hist[n, ] <- theta2
    itr_irls_hist[n] <- irls_cur$itr
    ess_is_hist[n] <- res_lik1$ess_is
    ess_sir_hist[n] <- res_sir$ess
    accepted1_hist[n] <- mh1_res$accepted
    accepted2_hist[n] <- mh2_res$accepted
  }
  
  return(list(
    theta_01_hist = theta_01_hist,
    theta_02_hist = theta_02_hist,
    W1_hist = W1_hist,
    W2_hist = W2_hist,
    theta1_hist = theta1_hist,
    theta2_hist = theta2_hist,
    itr_irls_hist = itr_irls_hist,
    ess_is_hist = ess_is_hist,
    ess_sir_hist = ess_sir_hist,
    accepted1_hist = accepted1_hist,
    accepted2_hist = accepted2_hist))
}
