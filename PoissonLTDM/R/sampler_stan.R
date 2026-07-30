# Poisson - 2nd Order Polynomial Dynamic Model
# MCMC: NUTS (Stan)
# Author: Cleiton Moya de Almeida


sample_stan <- function(model, y, N, burnin, seed,
						mu_01, sigma2_01, mu_02, sigma2_02,
						nu_01, eta_01, nu_02, eta_02,
						W1, W2, theta_01, theta_02, theta1, theta2){

	Tt <- length(y)
	phi1 <- 1/W1
	phi2 <- 1/W2
	
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

	# Initial values	
	initial_values <- list(
		list(
			theta_01 = theta_01,
			theta_02 = theta_02,
			theta1   = theta1,
			theta2   = theta2,
			phi1     = phi1,
			phi2     = phi2
		)
	)
	
	invisible(capture.output(
		fit <- suppressWarnings(
			rstan::sampling(
				object = model,
				data   = stan_data,
				chains = 1,
				iter   = N,
				warmup = burnin,
				thin   = 1,
				seed   = seed,
				init   = initial_values
			)
		)
	))
	
	warmup_time <- rstan::get_elapsed_time(fit)[1, 1]
	sample_time <- rstan::get_elapsed_time(fit)[1, 2]
	elapsed_time <- warmup_time + sample_time
	
	# Extract the samples
	samples <- rstan::extract(fit, permuted = FALSE, inc_warmup = TRUE)
	theta_01_hist <- samples[, 1, "theta_01"]
	theta_02_hist <- samples[, 1, "theta_02"]
	W1_hist <- samples[, 1, "W1"]
	W2_hist <- samples[, 1, "W2"]
	theta1_hist <- samples[, 1, paste0("theta1[", 1:Tt, "]")]
	theta2_hist <- samples[, 1, paste0("theta2[", 1:Tt, "]")]
	
	sampler_params <- rstan::get_sampler_params(fit, inc_warmup = TRUE)
	ac_hist <- sampler_params[[1]][, "accept_stat__"]

	return(list(
		theta_01_hist = theta_01_hist,
		theta_02_hist = theta_02_hist,
		W1_hist = W1_hist,
		W2_hist = W2_hist,
		theta1_hist = theta1_hist,
		theta2_hist = theta2_hist,
		ac_hist = ac_hist,
		elapsed_time = elapsed_time,
		fit = fit
	))
}
