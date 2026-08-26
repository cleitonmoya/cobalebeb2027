# application/real_data_run.R
#
# Application of the five samplers (amh_montoril, pg_apf, sir_laplace,
# sir_collapsed, stan) to a real count series (campy -- Ferland, Latour &
# Oraichi, 2006), using the same priors and per-method N/burnin/K already
# validated by the calibration study (R_config below, copied verbatim from
# simulation_run.R's own R_config).
#
# Unlike the simulation grid (1 chain per replica, since R=200 replicas
# already supply the convergence information R_hat would otherwise give),
# this is a SINGLE real dataset -- convergence must be assessed WITHIN this
# one application, via K_chains=3 dispersed-initialization chains per
# method, mirroring calibration_run.R's make_chain_inits() scheme (see
# metricas_theta1.tex, Section 5, for the full rationale and equations).
#
# IMPORTANT (per discussion): unlike run_task() in simulation_run.R, which
# only stores ESS BULK, every convergence/efficiency quantity here is
# computed and stored for BOTH bulk and tail -- ess_bulk_*, ess_tail_*,
# ess_sec_*_bulk, ess_sec_*_tail. This was flagged as a gap in the
# simulation stage; it is not retroactively fixed here (out of scope), but
# must not be repeated in this script.
#
# NO batchtools here, deliberately (per discussion): batchtools/PBS
# templates exist in simulation_run.R/calibration_run.R to manage SWEEPS of
# many jobs (registry, per-job checkpointing, selective resubmission,
# chunking). This script needs exactly ONE job -- 5 methods x K_chains=3 =
# 15 chain-units, all fitting in a single ncpus=15 exclusive node
# (place=excl reserves the whole node regardless of ncpus requested; see
# discussion for the memory estimate at Tt=140, ~235MB for montoril's
# heaviest chain -- no risk of the OOM issue that forced smaller chunks at
# Tt=1600 in the simulation grid). A sweep-management tool bringing a
# registry, job collections, and a PBS-script-generating template adds
# nothing here. This script simply runs all 15 chain-units in local
# parallel (run_chain_units(), one core per chain-unit) whenever it is
# invoked -- interactively, or inside a plain PBS job via
# real_data_pbs.sh (qsub real_data_pbs.sh), which merely calls
# `Rscript -e 'source("real_data_run.R")'` on a reserved node. No R-level
# "run_mode" branching is needed: it is the SAME code path either way.
#
# Seeds: simplified scheme (per discussion), NOT reusing calibration/
# production's structural seed formula (unnecessary complexity for a
# single one-off application) -- seed(m, k) = 900000 + m*1000 + k, m =
# method index (methods_grid_ref order), k = chain id (1..K_chains). Both
# calibration's and production's seed formulas start at method_idx*1e7
# (minimum 1e7); this scheme's maximum possible value (905003) is
# structurally, not incidentally, far below either grid's minimum -- zero
# collision risk by construction.

.defs_only_flag <- exists("REAL_DATA_DEFS_ONLY", inherits = FALSE) &&
					isTRUE(REAL_DATA_DEFS_ONLY)

rm(list = setdiff(ls(), ".defs_only_flag"))
options(error = function() traceback(2))

setwd(dirname(this.path::this.path()))

pkgload::load_all("../PoissonLTDM", debug = FALSE, quiet = TRUE)

printf <- function(...) cat(paste(sprintf(...), "\n"))

path_data           <- "../data/real"
path_results         <- "../results/real_data"
path_results_chains  <- sprintf("%s/%s", path_results, "chains")
dir.create(path_data, showWarnings = FALSE, recursive = TRUE)
dir.create(path_results_chains, showWarnings = FALSE, recursive = TRUE)

verbose <- TRUE

# ---- Prior hyperparameters (identical to simulation_run.R/calibration_run.R) ----
mu_01     <- 0        # theta_01 ~ N(mu_01, sigma2_01)
sigma2_01 <- 100
mu_02     <- 0        # theta_02 ~ N(mu_02, sigma2_02)
sigma2_02 <- 100
nu_01  <- 2           # phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
eta_01 <- 0.01
nu_02  <- 2           # phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
eta_02 <- 0.0001

# Adaptive Metropolis hyperparameters (montoril)
varsigma2_scal <- 0.02
ac_ref <- 0.44

# SIR Laplace and SIR Collapsed hyperparameters
M_is <- 3
M_irls_max <- 20
tol <- 1e-4
R_prerun <- 3000      # sir_collapsed only
M_sir <- 3            # sir_collapsed only
theta1_tilde_scal <- 0

# ---- Per-method N/burnin/K, calibration-validated (copied verbatim from
# simulation_run.R's R_config -- point 1 of the discussion) ----
R_config <- data.frame(
	method = c("montoril", "pg_apf",  "sir_laplace", "sir_collapsed", "stan"),
	N      = c(110000,     22000,     11000,         11000,           11000),
	burnin = c(10000,      2000,      1000,          1000,            1000),
	K      = c(NA,         200,       NA,            NA,              NA),
	stringsAsFactors = FALSE
)

methods_grid_ref <- R_config$method
K_chains <- 3

# ==========================================================================
# ---- Real data: campy (Ferland, Latour & Oraichi, 2006) ----
# ==========================================================================
#
# Campylobacterosis case counts, north of Quebec (Canada), 28-day
# intervals, January 1990 to end of October 2000 -- 140 observations, no
# covariates, no seasonal component (see discussion and metricas_theta1.tex
# for why this series was chosen: real Poisson counts, no seasonality,
# compatible with the level+trend-only PoissonLTDM). campy.rds (a plain
# numeric vector, length 140) must already exist under path_data -- it is
# NOT built here; place the file at ../data/real/campy.rds before running
# this script.

data_file <- sprintf("%s/campy.rds", path_data)
if (!file.exists(data_file)) {
	stop(sprintf("Real data file not found: %s -- place campy.rds there before running.", data_file))
}
y <- readRDS(data_file)
Tt <- length(y)

# ==========================================================================
# ---- Chain grid: one row per (method, chain_id) -- 5 x K_chains = 15 rows ----
# ==========================================================================

chain_grid <- do.call(rbind, lapply(seq_len(nrow(R_config)), function(i) {
	data.frame(method = R_config$method[i], N = R_config$N[i], burnin = R_config$burnin[i],
			   K = R_config$K[i], chain_id = seq_len(K_chains), stringsAsFactors = FALSE)
}))
chain_grid$seed <- 900000 + match(chain_grid$method, methods_grid_ref) * 1000 + chain_grid$chain_id

chain_result_filename <- function(chain) {
	sprintf("%s/%s_chain%d.rds", path_results_chains, chain$method, chain$chain_id)
}

# ---- Dispersed chain initialization (identical scheme to
# calibration_run.R's make_chain_inits(), applied to the real series
# instead of a simulated calibration series) ----
make_chain_inits <- function(chain_id, seed, y, Tt) {
	theta1_ref   <- log(y + 0.5)
	theta2_ref   <- c(diff(theta1_ref), 0)
	theta_01_ref <- theta1_ref[1]
	theta_02_ref <- theta2_ref[1]
	W1_ref <- var(diff(theta1_ref))
	W2_ref <- var(diff(theta2_ref))

	set.seed(seed)

	theta_01_init <- theta_01_ref + rnorm(1, 0, sd = 2 * sqrt(abs(theta_01_ref) + 1))
	theta_02_init <- theta_02_ref + rnorm(1, 0, sd = 2 * sqrt(abs(theta_02_ref) + 1))
	W1_init <- W1_ref * exp(rnorm(1, 0, sd = 1))
	W2_init <- W2_ref * exp(rnorm(1, 0, sd = 1))

	theta1_init <- theta1_ref + rnorm(Tt, 0, sd = 2 * sqrt(W1_ref))
	theta2_init <- theta2_ref + rnorm(Tt, 0, sd = 2 * sqrt(W2_ref))

	list(chain_id = chain_id, seed = seed,
		 theta_01 = theta_01_init, theta_02 = theta_02_init,
		 W1 = W1_init, W2 = W2_init,
		 theta1 = theta1_init, theta2 = theta2_init,
		 theta1_tilde = rep(0, Tt),
		 varsigma2 = rep(varsigma2_scal, Tt))
}

# ---- One-chain runner: unified across all 5 methods (mirrors
# calibration_run.R's run_one_chain(), and simulation_run.R's stan-model
# loading pattern) ----
run_one_chain <- function(chain) {
	chain_result_file <- chain_result_filename(chain)
	if (file.exists(chain_result_file)) {
		if (verbose) printf("\tChain already run (%s, chain %d), loading checkpoint", chain$method, chain$chain_id)
		return(readRDS(chain_result_file))
	}

	init <- make_chain_inits(chain$chain_id, chain$seed, y, Tt)
	method <- chain$method
	N <- chain$N
	burnin <- chain$burnin
	K <- chain$K

	if (verbose) printf("Running chain %d (seed=%d), method=%s", chain$chain_id, chain$seed, method)
	set.seed(chain$seed)

	if (method == "stan") {
		# Loaded/compiled here, guarded by file existence -- run_chain_units()
		# pre-builds this cache SERIALLY before spawning the parallel cluster,
		# so by the time this runs (possibly 3x in parallel, one per stan
		# chain), the cache already exists and every worker just reads it --
		# no race on ../cache/poisson_ltdm.rds.
		if (file.exists("../cache/poisson_ltdm.rds")) {
			stan_model <- readRDS("../cache/poisson_ltdm.rds")
		} else {
			if (verbose) printf("Building the Stan model")
			dir.create("../cache", showWarnings = FALSE, recursive = TRUE)
			stan_model <- rstan::stan_model(file = "../PoissonLTDM/inst/stan/poisson_ltdm.stan", model_name = "PoissonLTDM")
			saveRDS(stan_model, file = "../cache/poisson_ltdm.rds")
		}
	}

	execution_bench <- system.time({
		res <- switch(method,

			"montoril" = amh_montoril_cpp(y = y, N = N,
										   mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
										   nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
										   theta1 = init$theta1, theta2 = init$theta2,
										   theta_01 = init$theta_01, theta_02 = init$theta_02,
										   W1 = init$W1, W2 = init$W2, ac_ref = ac_ref, varsigma2 = init$varsigma2,
										   verbose = FALSE, print_every = 1000),

			"pg_apf" = pg_as_cpp(y = y, K = K, N = N,
								  mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
								  nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
								  theta1 = init$theta1, theta2 = init$theta2,
								  theta_01 = init$theta_01, theta_02 = init$theta_02,
								  W1 = init$W1, W2 = init$W2, verbose = FALSE, print_every = 1000),

			"sir_laplace" = sir_laplace_cpp(y = y, N = N,
											 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
											 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
											 theta1 = init$theta1, theta2 = init$theta2,
											 theta_01 = init$theta_01, theta_02 = init$theta_02,
											 W1 = init$W1, W2 = init$W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol,
											 verbose = FALSE, print_every = 1000),

			"sir_collapsed" = sir_collapsed_cpp(y = y, R_prerun = R_prerun, N = N,
												 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
												 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
												 theta1 = init$theta1, theta2 = init$theta2, theta1_tilde = init$theta1_tilde,
												 theta_01 = init$theta_01, theta_02 = init$theta_02, W1 = init$W1, W2 = init$W2,
												 M_is_lik = M_is, M_sir_theta1 = M_sir, M_irls_max = M_irls_max, tol = tol,
												 verbose = FALSE, print_every = 1000),

			"stan" = sample_stan(stan_model, y, N, burnin, chain$seed,
								  N_chains = 1, n_cores = 1,
								  chain_inits = list(list(theta_01 = init$theta_01, theta_02 = init$theta_02,
														   theta1 = init$theta1, theta2 = init$theta2,
														   W1 = init$W1, W2 = init$W2)),
								  mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
								  nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
								  verbose = FALSE)$results[[1]],

			stop(sprintf("Unknown method: %s", method))
		)
	})
	elapsed_time <- execution_bench[["user.self"]]
	if (method == "stan") elapsed_time <- res$elapsed_time

	chain_result <- list(
		method = method, chain_id = chain$chain_id, seed = chain$seed, elapsed_time = elapsed_time,
		theta_01_hist = res$theta_01_hist, theta_02_hist = res$theta_02_hist,
		W1_hist = res$W1_hist, W2_hist = res$W2_hist,
		theta1_hist = res$theta1_hist, theta2_hist = res$theta2_hist
	)
	saveRDS(chain_result, file = chain_result_file)
	chain_result
}


# ---- resolve_n_cores(): physical core count (identical to
# simulation_run.R's version -- not parallel::detectCores(logical=FALSE),
# confirmed unreliable on Euler) ----
resolve_n_cores <- function(n_cores, n_needed) {
	if (!is.null(n_cores)) return(min(n_cores, n_needed))

	phys <- NA_integer_
	if (Sys.info()[["sysname"]] == "Linux" && nzchar(Sys.which("lscpu"))) {
		out <- tryCatch(
			system("lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#'", intern = TRUE),
			error = function(e) character(0)
		)
		if (length(out) > 0) phys <- length(unique(out))
	}
	if (is.na(phys)) phys <- parallel::detectCores(logical = FALSE)
	if (is.na(phys)) {
		logi <- parallel::detectCores(logical = TRUE)
		phys <- max(1, logi %/% 2)
		printf("Physical core count unavailable from OS; falling back to floor(logical/2) = %d", phys)
	}
	min(phys, n_needed)
}

# ---- run_chain_units(): runs every row of `chains` in parallel, one core
# per chain-unit, inside a single local PSOCK cluster. Called identically
# whether this script runs interactively (login node / laptop) or inside a
# PBS job submitted via real_data_pbs.sh (qsub) -- ncpus=15 in that PBS
# script matches nrow(chain_grid), so every chain-unit gets its own core
# either way. No chunking: all 15 chain-units fit in one such call. ----
run_chain_units <- function(chains) {
	n <- nrow(chains)

	# Pre-build the stan cache ONCE, serially, before spawning workers -- see
	# run_one_chain()'s stan branch for why this avoids a build race when
	# multiple stan chains start in parallel.
	if (any(chains$method == "stan") && !file.exists("../cache/poisson_ltdm.rds")) {
		if (verbose) printf("Pre-building the Stan model (once, before spawning workers)")
		dir.create("../cache", showWarnings = FALSE, recursive = TRUE)
		stan_model <- rstan::stan_model(file = "../PoissonLTDM/inst/stan/poisson_ltdm.stan", model_name = "PoissonLTDM")
		saveRDS(stan_model, file = "../cache/poisson_ltdm.rds")
	}

	n_cores_used <- resolve_n_cores(NULL, n)
	if (verbose) printf("Running %d chain-unit(s) using %d core(s)", n, n_cores_used)

	current_wd <- getwd()
	cl <- parallel::makeCluster(n_cores_used)
	on.exit(parallel::stopCluster(cl), add = TRUE)
	parallel::clusterExport(cl, "current_wd", envir = environment())
	parallel::clusterExport(cl, ls(envir = environment(run_one_chain)), envir = environment(run_one_chain))
	parallel::clusterEvalQ(cl, {
		setwd(current_wd)
		pkgload::load_all("../PoissonLTDM", debug = FALSE, quiet = TRUE)
	})

	`%dorng%` <- doRNG::`%dorng%`
	doParallel::registerDoParallel(cl)
	res <- foreach::foreach(i = seq_len(n)) %dorng% {
		run_one_chain(chains[i, ])
	}
	res
}


# ==========================================================================
# ---- Aggregation: pool a method's K chains into R_hat/ESS (bulk AND
# tail), ESS/s (bulk AND tail), log-likelihood and log-CPO ----
# ==========================================================================
aggregate_method <- function(method_chains) {
	method <- method_chains[[1]]$method
	burnin <- R_config$burnin[R_config$method == method]

	# Pooled convergence diagnostics -- Eq. 30, metricas_theta1.tex.
	# metrics_convergence() (package function) implements the
	# rank-normalized, multi-chain R_hat/ESS bulk+tail via rhat_ess_fast().
	conv_theta_01 <- metrics_convergence(method_chains, "theta_01_hist", burnin)
	conv_theta_02 <- metrics_convergence(method_chains, "theta_02_hist", burnin)
	conv_W1       <- metrics_convergence(method_chains, "W1_hist", burnin)
	conv_W2       <- metrics_convergence(method_chains, "W2_hist", burnin)
	conv_theta1   <- metrics_convergence(method_chains, "theta1_hist", burnin)  # length-Tt vectors
	conv_theta2   <- metrics_convergence(method_chains, "theta2_hist", burnin)

	rhat_theta1_max <- max(conv_theta1$rhat)
	rhat_theta2_max <- max(conv_theta2$rhat)

	# Pooled posterior mean of W1, W2 (scalars, K chains stacked) -- W1/W2
	# don't vary with t, so their cross-method agreement (Delta_max) is a
	# single number per parameter, not a curve over t (unlike theta1/theta2).
	W1_pooled <- unlist(lapply(method_chains, function(ch) ch$W1_hist[-(1:burnin)]))
	W2_pooled <- unlist(lapply(method_chains, function(ch) ch$W2_hist[-(1:burnin)]))
	W1_mean <- mean(W1_pooled)
	W2_mean <- mean(W2_pooled)

	# Total CPU time (SUM across the K chains, not wall time) -- denominator
	# of Eq. 31.
	total_time <- sum(vapply(method_chains, function(ch) ch$elapsed_time, numeric(1)))

	# Per-chain CPU time -- K_chains independent replicate timings (distinct
	# seeds/inits, same N/burnin per method), giving mean and SE for the
	# cross-method CPU-time comparison against the tscount literature
	# benchmarks (see chat discussion). CAVEAT: with K_chains=3 (2 d.f.), the
	# Student-t multiplier for a 95% CI is t_{0.975,2}~4.30, not the usual
	# 1.96 -- report accordingly, do not treat this SE as if K were large.
	chain_times     <- vapply(method_chains, function(ch) ch$elapsed_time, numeric(1))
	total_time_mean <- mean(chain_times)
	total_time_se   <- sd(chain_times) / sqrt(length(chain_times))

	# ESS/s -- BOTH bulk and tail, for every parameter (the point flagged at
	# the top of this script: simulation_run.R only stored bulk).
	ess_sec <- function(conv) list(bulk = metrics_ess_per_sec(conv$ess_bulk, total_time),
									tail = metrics_ess_per_sec(conv$ess_tail, total_time))
	ess_sec_theta_01 <- ess_sec(conv_theta_01)
	ess_sec_theta_02 <- ess_sec(conv_theta_02)
	ess_sec_W1       <- ess_sec(conv_W1)
	ess_sec_W2       <- ess_sec(conv_W2)

	# Pooled post-burn-in draws (K chains stacked) for theta1_mean(t),
	# lambda_t = exp(theta1_t) samples (log-CPO), and log-likelihood at the
	# posterior mean -- Eq. 33-34.
	theta1_pooled <- do.call(rbind, lapply(method_chains, function(ch) {
		ch$theta1_hist[-(1:burnin), , drop = FALSE]
	}))
	theta1_mean_t  <- colMeans(theta1_pooled)
	theta2_pooled  <- do.call(rbind, lapply(method_chains, function(ch) {
		ch$theta2_hist[-(1:burnin), , drop = FALSE]
	}))
	theta2_mean_t  <- colMeans(theta2_pooled)
	lambda_mean_t  <- exp(theta1_mean_t)
	lambda_samples <- exp(theta1_pooled)

	log_lik <- metrics_loglik(y, lambda_mean_t)
	cpo <- metrics_log_cpo(y, lambda_samples)

	# 95% HPD band of theta1(t) and theta2(t), pooled across the K chains --
	# for the fit-overlay figure (application/figures_real_data.R), mirroring
	# by_t's band_lower/band_upper in the simulation study. metrics_theta_ci()
	# (package function, HDInterval::hdi() applied column-wise) already
	# exists for exactly this.
	ci_theta1 <- metrics_theta_ci(theta1_pooled, credMass = 0.95)
	ci_theta2 <- metrics_theta_ci(theta2_pooled, credMass = 0.95)

	list(
		method = method, total_time = total_time,
		total_time_mean = total_time_mean, total_time_se = total_time_se,

		rhat_theta_01 = conv_theta_01$rhat, rhat_theta_02 = conv_theta_02$rhat,
		rhat_W1 = conv_W1$rhat, rhat_W2 = conv_W2$rhat,
		rhat_theta1 = conv_theta1$rhat, rhat_theta2 = conv_theta2$rhat,
		rhat_theta1_max = rhat_theta1_max, rhat_theta2_max = rhat_theta2_max,

		ess_bulk_theta_01 = conv_theta_01$ess_bulk, ess_tail_theta_01 = conv_theta_01$ess_tail,
		ess_bulk_theta_02 = conv_theta_02$ess_bulk, ess_tail_theta_02 = conv_theta_02$ess_tail,
		ess_bulk_W1 = conv_W1$ess_bulk, ess_tail_W1 = conv_W1$ess_tail,
		ess_bulk_W2 = conv_W2$ess_bulk, ess_tail_W2 = conv_W2$ess_tail,
		ess_bulk_theta1 = conv_theta1$ess_bulk, ess_tail_theta1 = conv_theta1$ess_tail,
		ess_bulk_theta2 = conv_theta2$ess_bulk, ess_tail_theta2 = conv_theta2$ess_tail,
		ess_bulk_theta1_mean = mean(conv_theta1$ess_bulk), ess_bulk_theta1_min = min(conv_theta1$ess_bulk),
		ess_tail_theta1_mean = mean(conv_theta1$ess_tail), ess_tail_theta1_min = min(conv_theta1$ess_tail),
		ess_bulk_theta2_mean = mean(conv_theta2$ess_bulk), ess_bulk_theta2_min = min(conv_theta2$ess_bulk),
		ess_tail_theta2_mean = mean(conv_theta2$ess_tail), ess_tail_theta2_min = min(conv_theta2$ess_tail),

		ess_sec_theta_01_bulk = ess_sec_theta_01$bulk, ess_sec_theta_01_tail = ess_sec_theta_01$tail,
		ess_sec_theta_02_bulk = ess_sec_theta_02$bulk, ess_sec_theta_02_tail = ess_sec_theta_02$tail,
		ess_sec_W1_bulk = ess_sec_W1$bulk, ess_sec_W1_tail = ess_sec_W1$tail,
		ess_sec_W2_bulk = ess_sec_W2$bulk, ess_sec_W2_tail = ess_sec_W2$tail,
		ess_sec_theta1_bulk_mean = mean(metrics_ess_per_sec(conv_theta1$ess_bulk, total_time)),
		ess_sec_theta1_bulk_min  = min(metrics_ess_per_sec(conv_theta1$ess_bulk, total_time)),
		ess_sec_theta1_tail_mean = mean(metrics_ess_per_sec(conv_theta1$ess_tail, total_time)),
		ess_sec_theta1_tail_min  = min(metrics_ess_per_sec(conv_theta1$ess_tail, total_time)),
		ess_sec_theta2_bulk_mean = mean(metrics_ess_per_sec(conv_theta2$ess_bulk, total_time)),
		ess_sec_theta2_bulk_min  = min(metrics_ess_per_sec(conv_theta2$ess_bulk, total_time)),
		ess_sec_theta2_tail_mean = mean(metrics_ess_per_sec(conv_theta2$ess_tail, total_time)),
		ess_sec_theta2_tail_min  = min(metrics_ess_per_sec(conv_theta2$ess_tail, total_time)),

		theta1_mean = theta1_mean_t, theta2_mean = theta2_mean_t,
		theta1_ci_lower = ci_theta1$ci_lower, theta1_ci_upper = ci_theta1$ci_upper,
		theta2_ci_lower = ci_theta2$ci_lower, theta2_ci_upper = ci_theta2$ci_upper,
		W1_mean = W1_mean, W2_mean = W2_mean,
		log_lik = log_lik, log_cpo = cpo$log_cpo, cpo_t = cpo$cpo
	)
}


# ==========================================================================
# ---- Run all 15 chain-units, then aggregate per method + cross-method
# agreement. Single, unconditional code path -- no run_mode branching:
# invoking this script IS the "local" case, and invoking it inside
# real_data_pbs.sh's PBS job (see application/real_data_pbs.sh) is the
# only other way it ever runs, with the exact same code. ----
# ==========================================================================

if (!.defs_only_flag) {

	printf("Running all %d chain-unit(s).", nrow(chain_grid))
	invisible(run_chain_units(chain_grid))

	# ---- Aggregation: checkpoints make this idempotent -- re-running the
	# script after an interrupted run just reloads whichever chains are
	# already done and re-runs the rest. ----
	all_chain_files <- vapply(seq_len(nrow(chain_grid)), function(i) chain_result_filename(chain_grid[i, ]), character(1))

	if (all(file.exists(all_chain_files))) {
		if (verbose) printf("All %d chain-unit(s) done -- aggregating.", nrow(chain_grid))

		chain_results <- lapply(all_chain_files, readRDS)

		method_summaries <- lapply(methods_grid_ref, function(m) {
			method_chains <- Filter(function(ch) ch$method == m, chain_results)
			aggregate_method(method_chains)
		})
		names(method_summaries) <- methods_grid_ref

		# Cross-method agreement -- Eq. 32, metricas_theta1.tex. theta1/theta2
		# are curves over t; W1/W2 are single posterior-mean scalars per
		# method, so their Delta_max is a single number, not a length-Tt
		# vector (metrics_agreement_max() handles both the same way).
		theta1_means_by_method <- lapply(method_summaries, function(s) s$theta1_mean)
		theta2_means_by_method <- lapply(method_summaries, function(s) s$theta2_mean)
		W1_means_by_method <- lapply(method_summaries, function(s) s$W1_mean)
		W2_means_by_method <- lapply(method_summaries, function(s) s$W2_mean)
		delta_max_theta1 <- metrics_agreement_max(theta1_means_by_method)
		delta_max_theta2 <- metrics_agreement_max(theta2_means_by_method)
		delta_max_W1 <- metrics_agreement_max(W1_means_by_method)
		delta_max_W2 <- metrics_agreement_max(W2_means_by_method)

		summary_df <- do.call(rbind, lapply(method_summaries, function(s) {
			data.frame(
				method = s$method, total_time = s$total_time,
				total_time_mean = s$total_time_mean, total_time_se = s$total_time_se,

				rhat_theta1_max = s$rhat_theta1_max, rhat_theta2_max = s$rhat_theta2_max,

				ess_bulk_theta1_mean = s$ess_bulk_theta1_mean, ess_tail_theta1_mean = s$ess_tail_theta1_mean,
				ess_bulk_theta1_min  = s$ess_bulk_theta1_min,  ess_tail_theta1_min  = s$ess_tail_theta1_min,
				ess_bulk_theta2_mean = s$ess_bulk_theta2_mean, ess_tail_theta2_mean = s$ess_tail_theta2_mean,

				ess_sec_theta1_bulk_mean = s$ess_sec_theta1_bulk_mean, ess_sec_theta1_tail_mean = s$ess_sec_theta1_tail_mean,
				ess_sec_theta2_bulk_mean = s$ess_sec_theta2_bulk_mean, ess_sec_theta2_tail_mean = s$ess_sec_theta2_tail_mean,

				log_lik = s$log_lik, log_cpo = s$log_cpo,
				stringsAsFactors = FALSE
			)
		}))
		rownames(summary_df) <- NULL

		saveRDS(method_summaries, file = sprintf("%s/method_summaries.rds", path_results))
		saveRDS(list(theta1 = delta_max_theta1, theta2 = delta_max_theta2, W1 = delta_max_W1, W2 = delta_max_W2),
				file = sprintf("%s/delta_max.rds", path_results))
		write.csv(summary_df, file = sprintf("%s/summary.csv", path_results), row.names = FALSE)

		if (verbose) {
			printf("Summary saved to %s/summary.csv", path_results)
			print(summary_df)
		}
	} else {
		printf("%d of %d chain-unit(s) missing -- re-run this script to continue from checkpoints.",
			   sum(!file.exists(all_chain_files)), length(all_chain_files))
	}
}
