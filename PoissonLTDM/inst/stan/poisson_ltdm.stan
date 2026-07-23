// =========================================================================
// Poisson 2nd Order Polynomial Model
// y_t ~ Poisson(exp(theta_t1))
// theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1,  omega_t1 ~ N(0, W1)
// theta_t2 =               + theta_{t-1,2} + omega_t2,  omega_t2 ~ N(0, W2)
// theta_01 ~ N(mu_01, sigma2_01)
// theta_02 ~ N(mu_02, sigma2_02)
// phi1 = 1/W1 ~ Gamma(nu_01, eta_01)
// phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
// =========================================================================

data {
  int<lower=1> T;          // number of observations
  array[T] int<lower=0> y; // observed data y_1, ..., y_T

  // Prior hyperparameters
  real mu_01;               // média prior de theta_01
  real<lower=0> sigma2_01;  // variância prior de theta_01
  real<lower=0> nu_01;      // shape da Gamma para phi = 1/W1
  real<lower=0> eta_01;     // rate  da Gamma para phi = 1/W1

  real mu_02;               // média prior de theta_02
  real<lower=0> sigma2_02;  // variância prior de theta_02
  real<lower=0> nu_02;      // shape da Gamma para phi = 1/W2
  real<lower=0> eta_02;     // rate  da Gamma para phi = 1/W2
}

parameters {
  array[T] real theta1;     // estados latentes theta_11, ..., theta_T1
  real theta1_init;         // theta_01 (estado inicial)
  real<lower=0> phi1;       // phi = 1/W1 (precisão do sistema)

  array[T] real theta2;     // estados latentes theta_12, ..., theta_T2
  real theta2_init;         // theta_02
  real<lower=0> phi2;       // phi = 1/W2 (precisão do sistema)
}

transformed parameters {
  real<lower=0> W1;
  W1 = 1.0 / phi1;

  real<lower=0> W2;
  W2 = 1.0 / phi2;
}

model {
  // ---- Priors ----

  // theta_01 ~ N(mu_01, sigma2_01)
  // theta_02 ~ N(mu_02, sigma2_02)
  theta1_init ~ normal(mu_01, sqrt(sigma2_01));
  theta2_init ~ normal(mu_02, sqrt(sigma2_02));

  // phi1 = 1/W1 ~ Gamma(nu_01, eta_01)
  // phi2 = 1/W2 ~ Gamma(nu_02, eta_02)
  phi1 ~ gamma(nu_01, eta_01);
  phi2 ~ gamma(nu_02, eta_02);

  // ---- System equation ----
  // theta_t1 = theta_{t-1,1} + theta_{t-1,2} + omega_t1, omega_t1 ~ N(0,W1)
  // theta_t2 =                 theta_{t-1,2} + omega_t2, omega_t1 ~ N(0,W2)
  theta1[1] ~ normal(theta1_init + theta2_init, sqrt(W1));
  for (t in 2:T) {
    theta1[t] ~ normal(theta1[t-1] + theta2[t-1], sqrt(W1));
  }

  theta2[1] ~ normal(theta2_init, sqrt(W2));
  for (t in 2:T) {
    theta2[t] ~ normal(theta2[t-1], sqrt(W2));
  }

  // ---- Likelihood ----
  // y_t | theta_t1 ~ Poisson(exp(theta_t1))
  for (t in 1:T) {
    y[t] ~ poisson_log(theta1[t]);
  }
}

generated quantities {
  array[T] real lambda_hat;
  for (t in 1:T) {
    lambda_hat[t] = exp(theta1[t]);
  }
}
