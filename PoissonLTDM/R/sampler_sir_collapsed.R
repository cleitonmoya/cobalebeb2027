# Poisson - 2nd Order Polynomial Dynamic Model
#
# MCMC: Gibbs with SIR for theta1 and Collapsed Samplers for W1 and W2
#
# Strategy:
#   - Pre-run: standard Gibbs (R_prerun iterations) to collect phi1 and phi2 samples
#   - CE calibration: fit Gamma(c_hat, d_hat) proposals for both phi1 and phi2
#   - Gibbs sampling:
#       * (theta_01, theta1): EXTENDED-state joint block ((T+1)-dimensional),
#         sampled together via chan_smoothing_theta1/chan_sample_from_build -
#         still needs SIR (the Poisson likelihood on theta1[1..T] makes it
#         only a Laplace/IRLS approximation; theta_01 itself needs no
#         approximation, no likelihood on it).
#       * (theta_02, theta2): EXTENDED-state joint block ((T+1)-dimensional),
#         sampled together via chan_smoothing_theta2_ext/chan_sample_from_build -
#         exact, no SIR needed (no Poisson likelihood involved).
#       * W1: Collapsed MH with CE proposal marginalizing (theta_01, theta1)
#             JOINTLY; Integrated likelihod: IS with Laplace approx.
#       * W2: Collapsed MH with CE proposal marginalizing theta2, theta_02
#             held FIXED (ANCHORED, T-dimensional block - mathematically
#             distinct from the joint theta_02/theta2 block above, see
#             chan_smoothing_theta2_w2lik); Integrate likelihood: Gaussian
#
# Author: Cleiton Moya de Almeida


# IRLS: build the Laplace approximation for (theta_01, theta1) JOINTLY via the
# extended block. theta_01 itself needs no linearization (no likelihood) -
# only its exact Gaussian prior and the exact process link into theta1[1].
run_irls <- function(theta1_tilde, theta2, theta_02, phi1,
                     y, tol, M_irls_max, chan_smoothing_theta1,
                     mu_01, sigma2_01) {
    for (j in 1:M_irls_max) {
        f_t   <- exp(-theta1_tilde)
        phi_V <- 1 / f_t
        z_t   <- theta1_tilde + f_t * y - 1   # Poisson pseudo-observation

        theta1_build <- chan_smoothing_theta1(z_t, phi_V, phi1, mu_01, sigma2_01, theta_02, theta2)
        theta1_tilde_new <- theta1_build$theta1_hat[-1]   # drop node 0 (theta_01)
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


# IS estimator of log p(y | phi1, theta2) [integrated over (theta_01, theta1)
# JOINTLY]. Approximate: the Poisson likelihood makes p(theta1|y,phi1,theta2)
# non-Gaussian. theta_01 itself needs no approximation (no likelihood), but
# once the proposal q integrates it jointly (via the extended block), the
# target here must integrate it jointly too - holding it fixed would break
# the target/proposal cancellation, not just lose efficiency.
is_log_lik <- function(irls_res, y, phi1, theta2, theta_02, M_is,
                       Tt, Ttp1, mu_01, sigma2_01) {
    
    theta1_build <- irls_res$theta1_build
    eta_hat <- theta1_build$theta1_hat
    ch <- theta1_build$ch
    
    W1 <- 1/phi1

    log_det_H <- 2 * as.numeric(Matrix::determinant(ch, logarithm = TRUE)$modulus)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])
    log_norm_H <- -Ttp1 / 2 * log(2 * pi) + 0.5 * log_det_H

    log_w <- numeric(M_is)
    u_mat <- matrix(rnorm(M_is * Ttp1), nrow = M_is)
    d_ch  <- Matrix::diag(ch)

    for (i in 1:M_is) {
        u <- u_mat[i, ]
        w <- u / sqrt(d_ch)
        x <- as.vector(Matrix::solve(ch, w, system = "Lt"))
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

    log_lik <- logsumexp(log_w) - log(M_is)
    w_norm  <- exp(log_w - logsumexp(log_w))
    ess_is  <- 1 / sum(w_norm^2)

    return(list(log_lik = log_lik, ess_is = ess_is))
}



# ANCHORED, T-dimensional Chan smoother for theta2 given theta_02 held FIXED
# (single consumer: log_marginal_lik_w2 below, for the W2 marginal
# likelihood - mathematically distinct from the extended, jointly-sampled
# (theta_02, theta2) block used elsewhere in this sampler, so it is kept as
# its own dedicated, locally-scoped factory rather than in utils.R)
make_chan_theta2_smoother <- function(Tt) {
	
	res <- chan_build_static_objects(Tt)
	P2_matrix      <- res$K0
	Ch02_factor    <- res$Ch0_factor
	main_diag_base <- res$main_diag_base
	sub_diag_base  <- res$sub_diag_base
	idx_diag       <- res$idx_diag
	idx_sub        <- res$idx_sub
	
	chan_smoothing_theta2 <- function(theta1, phi1, phi2, theta_02) {
		z <- diff(theta1)   # z_t = theta1[t+1] - theta1[t], t=1,...,T-1
		
		diag_obs <- c(rep(phi1, Tt-1), 0)
		P2_matrix@x[idx_diag] <- (main_diag_base*phi2) + diag_obs
		P2_matrix@x[idx_sub]  <- -phi2
		
		Ch2_factor <- Matrix::update(Ch02_factor, P2_matrix)
		
		b <- numeric(Tt)
		b[1:(Tt-1)] <- z * phi1
		b[1] <- b[1] + theta_02 * phi2
		
		theta2_hat <- as.numeric(Matrix::solve(Ch2_factor, b, system="A"))
		list(theta2_hat = theta2_hat, ch = Ch2_factor, z = z)
	}
}


# EXACT integrated likelihood log p(z | phi2, phi1, theta_02) [integrated over
# theta2, theta_02 held FIXED]. Mathematically distinct from the joint
# (theta_02, theta2) extended block used elsewhere in this sampler - the W2
# collapsed MH needs the likelihood conditional on the CURRENT theta_02, so
# chan_smoothing_theta2 here must be the ANCHORED, T-dimensional smoother
# (make_chan_theta2_smoother, never the extended one). Exact: z|theta2 is
# Gaussian evaluated via the Chan sparse Cholesky factor at the conditional
# mode theta2_hat:
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
    Tt <- length(theta1)

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


# Collapsed MH step for W1 (independence chain with CE proposal, approximate
# IS likelihood, marginalizing (theta_01, theta1) JOINTLY)
mh_phi1_collapsed <- function(phi1_cur, log_lik1_cur, irls_cur,
                            ce1_params, nu_01, eta_01,
                            theta2, theta_02,
                            theta1_tilde, y, tol, M_irls_max, M_is,
                            chan_smoothing_theta1, mu_01, sigma2_01, Tt, Ttp1) {

    phi1_prop <- rgamma(1, shape = ce1_params$shape, rate = ce1_params$rate)

    irls_prop    <- run_irls(theta1_tilde, theta2, theta_02,
                             phi1_prop, y, tol, M_irls_max, chan_smoothing_theta1,
                             mu_01, sigma2_01)
    res_lik_prop <- is_log_lik(irls_prop, y, phi1_prop, theta2, theta_02, M_is,
                               Tt, Ttp1, mu_01, sigma2_01)
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


# Collapsed MH step for W2 (independence chain with CE proposal, EXACT
# likelihood, theta_02 held FIXED - see log_marginal_lik_w2)
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


# Sampling Importance Resampling for (theta_01, theta1) JOINTLY (extended,
# (T+1)-dimensional; corrects the Laplace/IRLS Gaussian approximation via
# importance resampling)
sir_theta1 <- function(irls_res, y, phi1, theta2, theta_02, M_sir,
                       mu_01, sigma2_01, Tt, Ttp1) {
    
    eta_hat <- irls_res$theta1_build$theta1_hat
    ch <- irls_res$theta1_build$ch
    W1 <- 1/phi1

    log_det_H <- 2*as.numeric(Matrix::determinant(ch, logarithm = TRUE)$modulus)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])
    log_norm_H <- -(Ttp1/2)*log(2*pi) + 0.5*log_det_H

    log_w <- numeric(M_sir)
    draws <- matrix(0, M_sir, Ttp1)
    d_ch <- Matrix::diag(ch)

    for (i in 1:M_sir) {
        u <- rnorm(Ttp1)
        w <- u / sqrt(d_ch)
        x <- as.vector(Matrix::solve(ch, w, system = "Lt"))
        draw_i <- eta_hat + x
        draws[i, ] <- draw_i
        theta_01_i <- draw_i[1]
        th         <- draw_i[-1]

        log_py <- sum(y * th - exp(th))
        log_prior_01 <- -0.5 * log(2 * pi * sigma2_01) - (theta_01_i - mu_01)^2 / (2 * sigma2_01)
        th_lag1 <- c(theta_01_i, th[-Tt])
        eps <- th - th_lag1 - th_lag2_fixed
        log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
        log_q <- log_norm_H - 0.5 * sum(u^2)

        log_w[i] <- log_py + log_prior_01 + log_prior_th - log_q
    }

    w_norm <- exp(log_w - logsumexp(log_w))
    ess_sir <- 1 / sum(w_norm^2)
    idx <- sample.int(M_sir, size = 1, prob = w_norm)

    list(theta_01 = draws[idx, 1], theta1 = draws[idx, -1], ess = ess_sir)
}

# PRE-RUN: standard Gibbs to collect phi1 and phi2 samples for CE calibration
collapsed_prerun <- function(R_prerun, y, Tt, Ttp1,
                             M_irls_max, tol, M_sir,
                             mu_01, sigma2_01, mu_02, sigma2_02,
                             nu_01, eta_01, nu_02, eta_02,
                             W1, W2, theta_01, theta_02, theta1, theta2, 
                             theta1_tilde,
                             chan_smoothing_theta1, chan_smoothing_theta2_ext) {
  
  phi1 <- 1/W1
  phi2 <- 1/W2
  
  phi1_prerun <- numeric(R_prerun)
  phi2_prerun <- numeric(R_prerun)
  
  for (r in 1:R_prerun) {
    
    # Sample phi1 (conjugated Gamma)
    phi1 <- gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1,
                              theta_02, theta2, Tt)
    phi1_prerun[r] <- phi1
    W1 <- 1/phi1
    
    
    # Sample phi2 (conjugated Gamma)
    phi2 <- gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt)
    phi2_prerun[r] <- phi2
    W2 <- 1/phi2
    
    
    # (theta_01, theta1) jointly (SIR, extended)
    irls_res <- run_irls(theta1_tilde, theta2, theta_02,
                         phi1, y, tol, M_irls_max, chan_smoothing_theta1,
                         mu_01, sigma2_01)
    theta1_tilde <- irls_res$theta1_tilde
    res_sir <- sir_theta1(irls_res, y, phi1, theta2, theta_02, M_sir,
                          mu_01, sigma2_01, Tt, Ttp1)
    theta_01 <- res_sir$theta_01
    theta1   <- res_sir$theta1
    
    # (theta_02, theta2) jointly via extended block
    build2 <- chan_smoothing_theta2_ext(theta1, phi1, phi2, mu_02, sigma2_02, theta_01)
    theta2_draw <- chan_sample_from_build(build2, Ttp1)
    theta_02 <- theta2_draw[1]
    theta2   <- theta2_draw[-1]
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
                                 theta1_tilde_scal) {
  Tt <- length(y)
  Ttp1 <- Tt + 1
  theta1_tilde <- rep(theta1_tilde_scal, Tt)
  
  # Prepare Chan static objects
  # (theta_01, theta1) and (theta_02, theta2): extended, (T+1)-dimensional
  # joint blocks
  chan_smoothing_theta1 <- make_chan_theta1_smoother_ext(Ttp1)
  chan_smoothing_theta2_ext <- make_chan_theta2_smoother_ext(Ttp1)
  # W2's integrated likelihood needs theta_02 held FIXED (not jointly
  # marginalized) - a mathematically distinct computation from the joint
  # block above, so it gets its own dedicated ANCHORED, T-dimensional smoother
  chan_smoothing_theta2_w2lik <- make_chan_theta2_smoother(Tt)
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
  prerun_res <- collapsed_prerun(R_prerun, y, Tt, Ttp1,
                                 M_irls_max, tol, M_sir,
                                 mu_01, sigma2_01, mu_02, sigma2_02,
                                 nu_01, eta_01, nu_02, eta_02,
                                 W1, W2, theta_01, theta_02, theta1, theta2, 
                                 theta1_tilde,
                                 chan_smoothing_theta1, chan_smoothing_theta2_ext)
  
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
  res_lik1 <- is_log_lik(irls_cur, y, phi1, theta2, theta_02, M_is,
                         Tt, Ttp1, mu_01, sigma2_01)
  log_lik1_cur <- res_lik1$log_lik
  
  # Integrated likelihood for W2
  res_lik2 <- log_marginal_lik_w2(theta1, phi1, phi2, theta_02, 
                                  chan_smoothing_theta2_w2lik, log_det_K0)
  log_lik2_cur <- res_lik2$log_lik
  theta2_build_cur <- res_lik2$build
  
  
  # Gibbs sampling
  for (n in 1:N) {
    
    # Collapsed MH for phi1 (approximate integrated likelihood, marginalizes
    # (theta_01, theta1) JOINTLY)
    mh1_res <- mh_phi1_collapsed(phi1, log_lik1_cur, irls_cur,
                                 ce1_params, nu_01, eta_01,
                                 theta2, theta_02,
                                 theta1_tilde, y, tol, M_irls_max, M_is,
                                 chan_smoothing_theta1, mu_01, sigma2_01, Tt, Ttp1)
    
    phi1 <- mh1_res$phi1
    W1 <- 1/phi1
    log_lik1_cur <- mh1_res$log_lik
    irls_cur <- mh1_res$irls_res
    
    
    # Collapsed MH for phi2 (EXACT integrated likelihood, marginalizes
    # theta2, theta_02 held FIXED)
    # Recompute log_lik2_cur first: phi1 may have changed since the last time
    # it was updated (phi1 updated by the W1 MH step just now).
    res_lik2_cur <- log_marginal_lik_w2(theta1, phi1, phi2, theta_02, 
                                        chan_smoothing_theta2_w2lik,
                                        log_det_K0)
    log_lik2_cur <- res_lik2_cur$log_lik
    theta2_build_cur <- res_lik2_cur$build
    
    mh2_res <- mh_phi2_collapsed(phi2, log_lik2_cur, theta2_build_cur,
                                 ce2_params, nu_02, eta_02,
                                 theta1, phi1, theta_02,
                                 chan_smoothing_theta2_w2lik, log_det_K0)
    
    phi2 <- mh2_res$phi2
    W2 <- 1/phi2
    
    
    # (theta_02, theta2) jointly via extended block given accepted phi2.
    # Sampled BEFORE (theta_01, theta1): theta1's SIR step needs the FRESH
    # theta2 to keep the phi2-theta1 coupling tight - doing it in the
    # opposite order (theta1 first) was found to degrade ESS(W2)
    # substantially (ESS(W2) ~4500 vs ~350).
    build2 <- chan_smoothing_theta2_ext(theta1, phi1, phi2, mu_02, sigma2_02, theta_01)
    draw2  <- chan_sample_from_build(build2, Ttp1)
    theta_02 <- draw2[1]
    theta2   <- draw2[-1]
    
    
    # (theta_01, theta1) jointly (SIR, extended, given accepted phi1 and the
    # fresh theta2)
    res_sir <- sir_theta1(irls_cur, y, phi1, theta2, theta_02, M_sir,
                          mu_01, sigma2_01, Tt, Ttp1)
    theta_01 <- res_sir$theta_01
    theta1   <- res_sir$theta1
    theta1_tilde <- irls_cur$theta1_tilde
   
    
    # Update irls_cur and log_lik1_cur for next iteration
    # (theta2, theta_02 both changed in this iteration)
    irls_cur <- run_irls(theta1_tilde, theta2, theta_02,
                         phi1, y, tol, M_irls_max, chan_smoothing_theta1,
                         mu_01, sigma2_01)
    theta1_tilde <- irls_cur$theta1_tilde
    
    res_lik1 <- is_log_lik(irls_cur, y, phi1, theta2, theta_02, M_is,
                           Tt, Ttp1, mu_01, sigma2_01)
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
