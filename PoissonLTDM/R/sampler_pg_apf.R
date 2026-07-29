# Poisson Local Trend Dynamic Model
#
# Model:
#  y_t ~ Poisson(exp{theta_t1})
#  theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1, omega_t1 ~ N(0, W1)
#  theta_t2 =                 theta_{t-1,1} + omega_t2, omega_t2 ~ N(0, W2)
#
# Priors:
#  theta_01 | D_0 ~ N(mu_01, sigma2_01)
#  theta_01 | D_0 ~ N(mu_02, sigma2_02)
#  1/W1 ~ gamma(shape=nu_01, rate=eta_01)
#  1/W2 ~ gamma(shape=nu_02, rate=eta_02)
#
# PMCMC:
#   Particle Gibbs (PG) with Backward Sampling
#   Strategy:
#    - theta1: Auxiliary Particle Filter (APF);
#    - theta2: Precision sampler (Chan Method)
#
# Reference: Andrieu, C., Doucet, A., & Holenstein, R. (2010).
#    Particle Markov Chain Monte Carlo Methods. Journal of the Royal
#    Statistical Society Series B: Statistical Methodology, 72(3), 269-342.
#
# Author: Cleiton Moya de Almeida


sample_pg_apf <- function(y, N, burnin, K,
                          mu_01, sigma2_01, mu_02, sigma2_02,
                          nu_01, eta_01, nu_02, eta_02,
                          W1, W2, theta_01, theta_02, theta1, theta2) {
  
  Tt <- length(y)
  Ttp1 <- Tt + 1
  
  # Prepare Chan static objects
  chan_smoothing_theta2 <- make_chan_theta2_smoother_ext(Ttp1)
  
  # Auxiliary variables
  W1_hist <- numeric(N)
  W2_hist <- numeric(N)
  theta_01_hist <- numeric(N)
  theta_02_hist <- numeric(N)
  theta1_hist <- matrix(0, N, Tt)
  theta2_hist <- matrix(0, N, Tt)
  ess_smc <- numeric(N)
  
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
    # Conditional SMC for theta1, AUGMENTED with a t=0 layer so that theta_01
    # is sampled JOINTLY with theta1 via backward sampling (instead of a
    # separate conjugate Normal draw)
    #
    theta_1_k   <- matrix(0, Tt, K)
    log_w_tilde <- matrix(0, Tt, K)
    
    # t = 0
    theta_0_k <- rnorm(K, mean = mu_01, sd = sqrt(sigma2_01))
    theta_0_k[K] <- theta_01   # reference path
    log_w_tilde_0 <- rep(-log(K), K)  # no likelihood at t = 0 (uniform weights)
    
    # t = 1
    # Predictor (auxiliary variable)
    theta_hat_11_k <- theta_0_k + theta_02
    
    # Auxiliary weights
    log_lambda_1_k <- y[1] * theta_hat_11_k - exp(theta_hat_11_k)
    
    # First stage resampling
    log_aux_1 <- log_w_tilde_0 + log_lambda_1_k
    A_1 <- sample(1:K, K, replace = TRUE, prob = exp(log_aux_1 - max(log_aux_1)))
    A_1[K] <- K   # reference path
    
    # Propagation
    theta_1_k[1, ] <- rnorm(K, mean = theta_0_k[A_1] + theta_02, sd = sd_W1)
    theta_1_k[1, K] <- theta1[1]
    
    # Updated weights
    log_w_1 <- log_p_yt(y[1], theta_1_k[1, ]) - log_lambda_1_k[A_1]
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
      b <- sample(1:K, 1, prob = bw)
      theta1[t] <- theta_1_k[t, b]
    }
    
    # backward sampling for theta_01
    # Backward weights: w_0^k * N(theta1[1] | theta_0^k + theta_02, W1)
    log_bw_0 <- log_w_tilde_0 +
      dnorm(theta1[1],
            mean = theta_0_k + theta_02,
            sd = sd_W1,
            log = TRUE
      )
    log_bw_0 <- log_bw_0 - max(log_bw_0)
    bw_0 <- exp(log_bw_0)
    bw_0 <- bw_0 / sum(bw_0)
    b_0 <- sample(1:K, 1, prob = bw_0)
    theta_01 <- theta_0_k[b_0]
    
    # (theta_02, theta2) jointly via extended block (Chan Method)
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
    ess_smc = ess_smc
  ))
  
}


