# Functions to compute the simulation metrics

# ---- Convergence diagnostics (R_hat, ESS bulk/tail) ----
#
# Replaces the previous coda::effectiveSize()/coda::geweke.diag()-based
# metrics_ess()/metrics_geweke(). R_hat and ESS bulk/tail (Vehtari et al.
# 2021, rank-normalized split-R_hat) are computed via rhat_ess_fast(), a
# C++ (Rcpp + RcppArmadillo + FFTW3 + OpenMP) reimplementation of
# posterior::rhat()/ess_bulk()/ess_tail() -- see src/convergence.cpp.
# Matches posterior's output to floating-point precision (~1e-13 to 1e-16
# relative difference, validated against posterior::summarise_draws()
# across multiple cases including odd n_post and N_chains=2), but avoids
# the R-level per-variable loop overhead that made posterior's own
# summarise_draws() take 90-100+ seconds per metric at Tt=1600 variables.
#
# Two variants, matching the two ways convergence is assessed across chains:
#
# - metrics_convergence(): R_hat and ESS bulk/tail POOLED across chains.
#   R_hat answers "do the chains agree with each other" and is always
#   pooled by construction. Pooled ESS is appropriate for scalar
#   parameters (theta_01, theta_02, W1, W2) where a single bulk/tail ESS
#   per parameter is the usual summary.
#
# - metrics_convergence_by_chain(): ESS bulk/tail computed SEPARATELY PER
#   CHAIN (not pooled). Used for theta1_hist/theta2_hist: a chain can have
#   high individual ESS while still sitting in a different mode from the
#   others -- pooling would mask that (R_hat, computed separately, is what
#   is supposed to catch it). Does not return R_hat (use
#   metrics_convergence() for that on the same hist_name).

# Build the (n_post, N_chains, Tt) array rhat_ess_fast() expects, from a
# result_list of N_chains elements each holding hist_name as either a
# vector of length N (scalar parameter, e.g. "W1_hist") or an N x Tt
# matrix (time-indexed parameter, e.g. "theta1_hist"). Internal, not
# exported.
.build_convergence_array <- function(result_list, hist_name, burnin) {
	N_chains <- length(result_list)
	first <- result_list[[1]][[hist_name]]
	if (is.null(first)) {
		stop(sprintf("'%s' not found in result_list[[1]].", hist_name))
	}

	if (is.null(dim(first))) {
		# scalar parameter: vector of length N -> array (n_post, N_chains, 1)
		N <- length(first)
		post_idx <- (burnin + 1):N
		n_post <- length(post_idx)
		arr <- array(NA_real_, dim = c(n_post, N_chains, 1))
		for (c in 1:N_chains) {
			arr[, c, 1] <- result_list[[c]][[hist_name]][post_idx]
		}
	} else {
		# time-indexed parameter: N x Tt matrix -> array (n_post, N_chains, Tt)
		N <- nrow(first)
		Tt <- ncol(first)
		post_idx <- (burnin + 1):N
		n_post <- length(post_idx)
		arr <- array(NA_real_, dim = c(n_post, N_chains, Tt))
		for (c in 1:N_chains) {
			arr[, c, ] <- result_list[[c]][[hist_name]][post_idx, ]
		}
	}
	arr
}


# R_hat and ESS bulk/tail, POOLED across chains.
#
# result_list: list of N_chains elements, each containing hist_name (either
#   a length-N vector for scalar parameters, or an N x Tt matrix for
#   time-indexed parameters).
# hist_name: name of the field to extract from each element of result_list
#   (e.g. "W1_hist", "theta1_hist").
# burnin: number of initial iterations to discard.
#
# Returns a list with rhat, ess_bulk, ess_tail -- each a scalar for scalar
# parameters, or a length-Tt vector for time-indexed parameters.
metrics_convergence <- function(result_list, hist_name, burnin) {
	arr <- .build_convergence_array(result_list, hist_name, burnin)
	res <- rhat_ess_fast(arr)
	if (dim(arr)[3] == 1) {
		list(rhat = res$rhat[1], ess_bulk = res$ess_bulk[1], ess_tail = res$ess_tail[1])
	} else {
		list(rhat = res$rhat, ess_bulk = res$ess_bulk, ess_tail = res$ess_tail)
	}
}


# ESS bulk/tail computed SEPARATELY PER CHAIN (not pooled). Does not
# return R_hat -- use metrics_convergence() on the same hist_name for that.
#
# Same arguments as metrics_convergence(). Returns a list with ess_bulk and
# ess_tail, each a Tt x N_chains matrix (or a length-N_chains vector when
# hist_name is a scalar parameter).
metrics_convergence_by_chain <- function(result_list, hist_name, burnin) {
	N_chains <- length(result_list)
	arr <- .build_convergence_array(result_list, hist_name, burnin)
	Tt <- dim(arr)[3]

	ess_bulk_mat <- matrix(NA_real_, nrow = Tt, ncol = N_chains)
	ess_tail_mat <- matrix(NA_real_, nrow = Tt, ncol = N_chains)
	for (c in 1:N_chains) {
		arr_c <- array(arr[, c, ], dim = c(dim(arr)[1], 1, Tt))
		res_c <- rhat_ess_fast(arr_c)
		ess_bulk_mat[, c] <- res_c$ess_bulk
		ess_tail_mat[, c] <- res_c$ess_tail
	}

	if (Tt == 1) {
		list(ess_bulk = ess_bulk_mat[1, ], ess_tail = ess_tail_mat[1, ])
	} else {
		list(ess_bulk = ess_bulk_mat, ess_tail = ess_tail_mat)
	}
}


# ESS bulk/tail for a SINGLE chain (no R_hat -- R_hat requires comparing
# multiple chains and is not meaningful here). Used in the main simulation
# stage (simulation.R), where convergence itself (R_hat, Geweke) is only
# assessed once, during the pilot/calibration stage that decides N and
# burnin -- not recomputed per replica. For that pilot/calibration stage,
# use metrics_convergence()/metrics_convergence_by_chain() with the
# multi-chain result_list instead.
#
# samples: either a vector of length n_post (scalar parameter, e.g. W1
#   after burn-in removal) or an n_post x Tt matrix (time-indexed
#   parameter, e.g. theta1 after burn-in removal).
#
# Returns a list with ess_bulk, ess_tail -- each a scalar for a vector
# input, or a length-Tt vector for a matrix input.
metrics_ess_single_chain <- function(samples) {
	if (is.null(dim(samples))) {
		n_post <- length(samples)
		arr <- array(samples, dim = c(n_post, 1, 1))
	} else {
		n_post <- nrow(samples)
		Tt <- ncol(samples)
		arr <- array(samples, dim = c(n_post, 1, Tt))
	}
	res <- rhat_ess_fast(arr)
	list(ess_bulk = res$ess_bulk, ess_tail = res$ess_tail)
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
metrics_theta1_ci <- function(theta_samples, alpha) {
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
