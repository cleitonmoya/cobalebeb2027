# Functions to compute the simulation metrics

# Effective sample size
metrics_ess <- function(sample) {
	coda::effectiveSize(coda::mcmc(sample))
}


# Gewege convergence diagnostic
metrics_geweke <- function(sample, frac1=0.1, frac2=0.5) {
	unname(coda::geweke.diag(sample, frac1=frac1, frac2=frac2)[[1]])
}


# Log-likelihood for the Poisson model
metrics_loglik <- function(y, lambda_mean) {
	sum(dpois(y, lambda_mean, log=TRUE))
}


# Root mean square error (RMSE)
metrics_rmse <- function(estimated, true) {
	sqrt(mean((estimated - true)^2))
}


# Mean absolute error (MAE)
metrics_mae <- function(estimated, true) {
	mean(abs(estimated - true))
}


# Credible interval (estimated) for theta1 and theta2
metrics_theta_ci <- function(theta_samples, alpha) {
	ci_lower <- apply(theta_samples, 2, quantile, probs = alpha/2)
	ci_upper <- apply(theta_samples, 2, quantile, probs = 1 - alpha/2)
	return(list(ci_lower=ci_lower, ci_upper=ci_upper))
}


# Credible interval (estimated) W1 and W2
metrics_W_ci <- function(W_samples, alpha) {
	ci_lower <- quantile(W_samples, probs = alpha/2)
	ci_upper <- quantile(W_samples, probs = 1 - alpha/2)
	return(list(ci_lower=ci_lower, ci_upper=ci_upper))
} 