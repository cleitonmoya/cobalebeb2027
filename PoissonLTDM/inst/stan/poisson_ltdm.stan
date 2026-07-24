// Poisson Local Trend Model
//
// Model:
//   y_t ~ Poisson(exp(theta_t1))
//   theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1,  omega_t1 ~ N(0, W1)
//   theta_t2 =               + theta_{t-1,2} + omega_t2,  omega_t2 ~ N(0, W2)
//
// Priors:
//   theta_01 ~ N(mu_01, sigma2_01)
//   theta_02 ~ N(mu_02, sigma2_02)
//   phi1 = 1/W1 ~ Gamma(nu_01, eta_01)
//   phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
//
// Author: Cleiton Moya de Ameida

data {
  int<lower=1> Tt; // number of observations
  array[Tt] int<lower=0> y; // observed data y_1, ..., y_T

  // Prior hyperparameters
  real mu_01;
  real<lower=0> sigma2_01;
  
  real mu_02;
  real<lower=0> sigma2_02;
  
  real<lower=0> nu_01;
  real<lower=0> eta_01;

  real<lower=0> nu_02;
  real<lower=0> eta_02;
}

parameters {
  array[Tt] real theta1;
  array[Tt] real theta2;
  real theta_01;
  real theta_02;
  real<lower=0> phi1;
  real<lower=0> phi2;
}

transformed parameters {
  real<lower=0> W1;
  W1 = 1.0 / phi1;

  real<lower=0> W2;
  W2 = 1.0 / phi2;
}

model {
  // Priors
  theta_01 ~ normal(mu_01, sqrt(sigma2_01));
  theta_02 ~ normal(mu_02, sqrt(sigma2_02));

  phi1 ~ gamma(nu_01, eta_01);
  phi2 ~ gamma(nu_02, eta_02);

  // System equations
  theta1[1] ~ normal(theta_01 + theta_02, sqrt(W1));
  theta2[1] ~ normal(theta_02, sqrt(W2));
  for (t in 2:Tt) {
    theta1[t] ~ normal(theta1[t-1] + theta2[t-1], sqrt(W1));
    theta2[t] ~ normal(theta2[t-1], sqrt(W2));
  }

  // Observation equation
  // y_t | theta_t1 ~ Poisson(exp(theta_t1))
  for (t in 1:Tt) {
    y[t] ~ poisson_log(theta1[t]);
  }
}
