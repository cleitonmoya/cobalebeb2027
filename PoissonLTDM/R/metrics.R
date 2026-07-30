# Functions to compute the simulation metrics

# Effective sample size
compute_ess <- function(sample) {
	coda::effectiveSize(coda::mcmc(sample))
}


# Gewege convergence diagnostic
compute_geweke <- function(sample, frac1=0.1, frac2=0.5) {
	unname(coda::geweke.diag(sample, frac1=frac1, frac2=frac2)[[1]])
}


# Log-likelihood for the Poisson model
compute_loglik <- function(y, lambda_mean) {
	sum(dpois(y, lambda_mean, log=TRUE))
}


# Root mean square error (RMSE)
compute_rmse <- function(estimated, true) {
	sqrt(mean((estimated - true)^2))
}


# Mean absolute error (MAE)
compute_mae <- function(estimated, true) {
	mean(abs(estimated - true))
}


# Empirical IC coverage
compute_cov <- function(alpha) {
	# To do
	-1
}