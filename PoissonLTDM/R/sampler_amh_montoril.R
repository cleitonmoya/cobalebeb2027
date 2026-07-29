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
# MCMC:
#  theta_t1: Component-Wise Metropolis within Gibbs (Randm Walking)
#  theta2:   Chan
#
# Author: Cleiton Moya de Almeida


sample_amh_montoril <- function(y, N, burnin, varsigma2_scal, ac_ref,
                         mu_01, sigma2_01, mu_02, sigma2_02,
                         nu_01, eta_01, nu_02, eta_02,
                         W1, W2, theta_01, theta_02, theta1, theta2){

    Tt <- length(y)
    Ttp1 <- Tt + 1
    varsigma2 <- rep(varsigma2_scal, Tt)

    # Prepare Chan static objects
    chan_smoothing_theta2 <- make_chan_theta2_smoother_ext(Ttp1)

    # Auxiliary vectors and matrix to store the results
    theta1_hist <- matrix(nrow=N, ncol=Tt)
    theta2_hist <- matrix(nrow=N, ncol=Tt)
    W1_hist <- numeric(N)
    W2_hist <- numeric(N)
    theta_01_hist <-numeric(N)
    theta_02_hist <-numeric(N)
    ac_hist <- matrix(0, nrow=N, ncol=Tt)

    # Gibbs sampling
    for (n in 1:N) {

        # Sample theta_01 (conjugated Normal)
        theta_01 <- gibbs_sample_theta01(mu_01, sigma2_01, theta1[1],
                                         theta_02, W1)

        # Sample phi1 (conjugated Gamma)
        phi1 <- gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1,
                                  theta_02, theta2, Tt)
        W1 <- 1/phi1

        # Sample phi2 (conjugated Gamma)
        phi2 <- gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt)
        W2 <- 1/phi2

        # Sample theta1 (Component-wise Metropolis)
        res_theta1 <- cwmh_sample_theta1(y, theta_01, theta_02,
                                     theta1, theta2, W1, varsigma2, Tt)
        theta1 <-res_theta1$theta1
        ac_hist[n,] <- res_theta1$ac

        # Adaptive stage of varsigma2
        delta <- min(0.01, 1/sqrt(n))
        ls <- log(varsigma2)/2 + delta*(ac_hist[n,] - ac_ref)
        varsigma2 <- exp(2*ls)
        
        # (theta_02, theta2) jointly via extended block (Chan Method)
        build2 <- chan_smoothing_theta2(theta1, phi1, phi2, mu_02, sigma2_02, theta_01)
        draw2  <- chan_sample_from_build(build2, Ttp1)
        theta_02 <- draw2[1]
        theta2   <- draw2[-1]

        # Store the sampled values
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
        ac_hist = ac_hist
    ))
}
