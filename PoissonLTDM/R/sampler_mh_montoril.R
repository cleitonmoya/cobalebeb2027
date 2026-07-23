# Poisson - 2nd Order Polynomial Dynamic Model
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


sample_mh_montoril <- function(y, N, burnin, varsigma2,
                         mu_01, sigma2_01, mu_02, sigma2_02,
                         nu_01, eta_01, nu_02, eta_02,
                         W1, W2, theta_01, theta_02, theta1, theta2){

    Tt <- length(y)

    # Prepare Chan static objects
    chan_sample_theta2 <- make_chan_theta2_sampler(Tt)

    # Auxiliary vectors and matrix to store the results
    theta1_hist <- matrix(nrow=N, ncol=Tt)
    theta2_hist <- matrix(nrow=N, ncol=Tt)
    W1_hist <- numeric(N)
    W2_hist <- numeric(N)
    theta_01_hist <-numeric(N)
    theta_02_hist <-numeric(N)
    ac_hist <- numeric(N)

    # Gibbs sampling
    for (n in 1:N) {

        # Sample theta_01 (conjugated Normal)
        theta_01 <- gibbs_sample_theta01(mu_01, sigma2_01, theta1[1],
                                         theta_02, W1)

        # Sample theta_02 (conjugated Normal)
        theta_02 <- gibbs_sample_theta02(mu_02, sigma2_02, theta_01, theta1[1],
                                         theta2[1], W1, W2)

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
        n_ac <- res_theta1$n_ac

        # Sample theta2 (Chan Method)
        theta2 <- chan_sample_theta2(theta1, phi1, phi2, theta_02, Tt)

        # Store the sampled values
        theta_01_hist[n] <- theta_01
        theta_02_hist[n] <- theta_02
        W1_hist[n] <- W1
        W2_hist[n] <- W2
        theta1_hist[n, ] <- theta1
        theta2_hist[n, ] <- theta2
        ac_hist[n] <- n_ac/Tt # mean acceptance ratio of theta_t1
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
