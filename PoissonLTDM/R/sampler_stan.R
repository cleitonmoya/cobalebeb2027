# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: NUTS (Stan)
# Author: Cleiton Moya de Almeida


sample_stan <- function(model, y, N, burnin, seed, N_chains, n_cores, chain_inits,
						mu_01, sigma2_01, mu_02, sigma2_02,
						nu_01, eta_01, nu_02, eta_02){

	Tt <- length(y)

	# Mapping stan data
	stan_data <- list(
		Tt        = Tt,
		y         = y,
		mu_01     = mu_01,
		sigma2_01 = sigma2_01,
		mu_02     = mu_02,
		sigma2_02 = sigma2_02,
		nu_01     = nu_01,
		eta_01    = eta_01,
		nu_02     = nu_02,
		eta_02    = eta_02
	)

	# Initial values: one list per chain, taken from chain_inits (same
	# dispersed initialization scheme as the other samplers -- see
	# make_chain_inits() in test_stan.R)
	initial_values <- lapply(chain_inits, function(init) {
		list(
			theta_01 = init$theta_01,
			theta_02 = init$theta_02,
			theta1   = init$theta1,
			theta2   = init$theta2,
			phi1     = 1 / init$W1,
			phi2     = 1 / init$W2
		)
	})

	invisible(capture.output(
		fit <- suppressWarnings(
			rstan::sampling(
				object  = model,
				data    = stan_data,
				chains  = N_chains,
				cores   = n_cores,
				iter    = N,
				warmup  = burnin,
				thin    = 1,
				seed    = seed,
				init    = initial_values,
				refresh = 0
			)
		)
	))

	warmup_time <- rstan::get_elapsed_time(fit)[, 1]
	sample_time <- rstan::get_elapsed_time(fit)[, 2]
	elapsed_time <- warmup_time + sample_time

	# Extract the samples, per chain, in the same list-of-N_chains-results
	# format used by the *_cpp() samplers (test_cpp.R), so that
	# print_and_plot_diagnostics() works unchanged.
	samples <- rstan::extract(fit, permuted = FALSE, inc_warmup = TRUE)
	sampler_params <- rstan::get_sampler_params(fit, inc_warmup = TRUE)

	results <- lapply(1:N_chains, function(c) {
		list(
			theta_01_hist = samples[, c, "theta_01"],
			theta_02_hist = samples[, c, "theta_02"],
			W1_hist       = samples[, c, "W1"],
			W2_hist       = samples[, c, "W2"],
			theta1_hist   = samples[, c, paste0("theta1[", 1:Tt, "]")],
			theta2_hist   = samples[, c, paste0("theta2[", 1:Tt, "]")],
			ac_hist       = sampler_params[[c]][, "accept_stat__"],
			elapsed_time  = elapsed_time[c]
		)
	})

	return(list(results = results, fit = fit, elapsed_time = elapsed_time))
}
