# Poisson - 2nd Order Polynomial Dynamic Model
#
# Strategy:
#  - theta1: Importance Sampling
#      Importance density: Laplace (normal) approximation
#  - theta2: Precision sampling (Chan)
#
# Author: Cleiton Moya de Almeida


sample_sir_laplace <- function(y, N, burnin, M_is, M_irls_max, tol,
                               mu_01, sigma2_01, mu_02, sigma2_02,
                               nu_01, eta_01, nu_02, eta_02,
                               W1, W2, theta_01, theta_02, theta1, theta2,
                               theta1_tilde_scal) {
  
  Tt <- length(y)
  Ttp1 <- Tt + 1
  theta1_tilde <- rep(theta1_tilde_scal, Tt)
  
  # Prepare Chan static objects
  chan_smoothing_theta1 <- make_chan_theta1_smoother_ext(Ttp1)
  chan_smoothing_theta2 <- make_chan_theta2_smoother_ext(Ttp1)
  
  # Auxiliary variables
  W1_hist <- numeric(N)
  W2_hist <- numeric(N)
  theta_01_hist <- numeric(N)
  theta_02_hist <- numeric(N)
  theta1_hist <- matrix(0, N, Tt)
  theta2_hist <- matrix(0, N, Tt)
  ess_is <- numeric(N)
  itr_irls <- numeric(N) # number of iterations of IRLS (for each gibbs step)
  
  # Gibbs sampling
  for (n in 1:N) {
    
    # Sample phi1 (conjugated Gamma)
    phi1 <- gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1,
                              theta_02, theta2, Tt)
    W1 <- 1/phi1
    sd_W1 <- sqrt(W1)
    
    # Sample phi2 (conjugated Gamma)
    phi2 <- gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt)
    W2 <- 1/phi2
    
    
    #
    # Importance Sampling for  theta_t1
    #
    
    # Local approximation (IRLS), (theta_01, theta1) jointly via extended block
    theta1_tilde_old <- theta1_tilde
    for (j in 1:M_irls_max) {
      
      f_t <- exp(-theta1_tilde)         # observational variance
      z_t <- theta1_tilde + f_t*y - 1   # pseudo-observation
      
      res <- chan_smoothing_theta1(z_t, 1/f_t, phi1, mu_01, sigma2_01, theta_02, theta2)
      theta1_tilde <- res$theta1_hat[-1]   # drop node 0 (theta_01) for the next linearization point
      
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
        # sample (theta_01, theta1) proposed jointly
        draw_prop <- chan_sample_from_build(res, Ttp1)
        theta_01_trajectories[i] <- draw_prop[1]
        theta1_prop <- draw_prop[-1]
        trajectories[i, ] <- theta1_prop
        
        # Log-weights: log p(y|theta1) - log g(y|theta1)
        # (theta_01's prior/process terms cancel exactly between target
        # and proposal, so the weight formula does not change)
        log_p <- sum(y * theta1_prop - exp(theta1_prop))
        log_g <- sum(-0.5 * log(2*pi*f_t) - 0.5 * (theta1_prop - z_t)^2 / f_t)
        log_w[i] <- log_p - log_g
      }
      
      log_w <- log_w - logsumexp(log_w)
      w <- exp(log_w)
      ess_is[n] <- 1 / sum(exp(2*log_w))
      
      idx <- sample(1:M_is, 1, prob=w)   # index for (theta_01*, theta1*)
      theta_01 <- theta_01_trajectories[idx]
      theta1   <- trajectories[idx, ] # theta1
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
    
    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1
    theta2_hist[n, ] <- theta2
    
  }
  
  return(list(
    theta_01_hist = theta_01_hist,
    theta_02_hist = theta_02_hist,
    W1_hist = W1_hist,
    W2_hist = W2_hist,
    theta1_hist = theta1_hist,
    theta2_hist = theta2_hist,
    ess_is = ess_is,
    itr_irls = itr_irls
  ))
  
}

