# application/predictive_loglik_pf.R
#
# One-step-ahead PREDICTIVE log-likelihood for each of the 5 methods --
# implements Eqs. pf_propagate_theta1/theta2, pf_weight, pf_contribution,
# pf_final of metricas_paper.tex (Section "Verossimilhanca preditiva
# um-passo-a-frente"). Each of a method's S^(m) post-burn-in posterior
# samples is a legitimate draw of (theta_01, theta_02, W1, W2); every
# particle carries its own (W1, W2) through every time step.
#
# Reads the RAW CHAIN CHECKPOINTS directly (../results/real_data/chains/),
# NOT method_summaries.rds. Does not depend on, modify, or re-run any of
# the 5 samplers themselves -- purely a post-hoc evaluation of their
# already-computed posterior samples.
#
# DESIGN (per discussion):
#  - Every method uses the SAME particle count N_COMMON (fair comparison --
#    using each method's raw S^(m) directly would confound "quality of the
#    estimated posterior" with "how large a particle swarm that method
#    happened to produce", since e.g. montoril's calibration-driven
#    N=110000 gives it ~10x more raw draws than sir_laplace/sir_collapsed/
#    stan).
#  - N_COMMON is set generously above the level that caused pathological
#    degeneracy (N=20000 gave an effective sample size of only ~57
#    particles at the outbreak, t=100) -- N_COMMON=200000 gives a ~10x
#    margin (effective size ~570 at that same point), comfortable without
#    being wasteful. The total compute cost to reach a target standard
#    error is roughly INVARIANT to how the budget is split between N and
#    replicates (SD scales ~1/sqrt(N), so N x replicates needed for a
#    fixed SE is ~constant) -- N only needs to clear the degeneracy floor;
#    beyond that, replicates do the rest of the precision work more
#    cheaply. Sampling WITH replacement when N_COMMON exceeds a method's
#    S^(m) is not a real limitation: every particle receives fresh,
#    independent process noise at each of the Tt=140 steps, so particles
#    that start as exact duplicates diverge immediately after the first
#    predict step -- the real diversity driver is the forward noise path,
#    not the starting draw alone.
#  - N_REPLICATES is also generous: SD (a property of a single filter
#    pass, controlled by N_COMMON) does not shrink by adding replicates --
#    only the standard ERROR of the mean does (SD/sqrt(replicates)).
#  - Runs in PARALLEL across the full (method x replicate) task grid, one
#    core per task, same pattern as run_chain_units() in real_data_run.R
#    (resolve_n_cores(), makeCluster, clusterExport, foreach %dorng%).
#
# VALIDATED (see chat/metricas_paper.tex) against KFAS's exact closed-form-
# adjacent computation for the nested level-only case (W2=0): bootstrap PF
# matches within Monte Carlo noise. Also checked for particle degeneracy at
# the outbreak (t=100) -- this is exactly why N_COMMON clears that floor
# with margin rather than being set minimally.

setwd(dirname(this.path::this.path()))

printf <- function(...) cat(paste(sprintf(...), "\n"))

path_data <- "../data/real"
path_results <- "../results/real_data"
path_results_chains <- sprintf("%s/chains", path_results)

y <- readRDS(sprintf("%s/campy.rds", path_data))
Tt <- length(y)

R_config <- data.frame(
	method = c("montoril", "pg_apf", "sir_laplace", "sir_collapsed", "stan"),
	burnin = c(10000, 2000, 1000, 1000, 1000),
	stringsAsFactors = FALSE
)
methods_grid_ref <- R_config$method
K_chains <- 3

# ---- Generous, cluster-scale settings -- N clears the degeneracy floor
# with a comfortable margin (see design note above); replicates then do
# the remaining precision work cheaply. ----
N_COMMON <- 200000
N_REPLICATES <- 100


# ---- Pool a method's K=3 chains' post-burn-in draws of theta_01, theta_02,
# W1, W2 -- these become the pool from which particles are drawn. Loaded
# ONCE per method (reads the ~100-250MB chain files a single time), then
# only the small pooled numeric vectors (a few MB) are exported to workers
# -- not re-read from disk per task. ----
pool_method_samples <- function(method, burnin) {
	chain_files <- sprintf("%s/%s_chain%d.rds", path_results_chains, method, seq_len(K_chains))
	if (!all(file.exists(chain_files))) {
		stop(sprintf("Missing chain checkpoint(s) for method '%s': %s", method,
					 paste(chain_files[!file.exists(chain_files)], collapse = ", ")))
	}
	chains <- lapply(chain_files, readRDS)
	list(
		theta_01 = unlist(lapply(chains, function(ch) ch$theta_01_hist[-(1:burnin)])),
		theta_02 = unlist(lapply(chains, function(ch) ch$theta_02_hist[-(1:burnin)])),
		W1       = unlist(lapply(chains, function(ch) ch$W1_hist[-(1:burnin)])),
		W2       = unlist(lapply(chains, function(ch) ch$W2_hist[-(1:burnin)]))
	)
}

printf("Pooling posterior samples for all 5 methods (reads chain checkpoints once)...")
all_samples <- lapply(methods_grid_ref, function(m) {
	burnin <- R_config$burnin[R_config$method == m]
	s <- pool_method_samples(m, burnin)
	printf("  %s: S_full = %d", m, length(s$theta_01))
	s
})
names(all_samples) <- methods_grid_ref


# ---- Bootstrap particle filter -- Eqs. pf_propagate_theta2/theta1,
# pf_weight, pf_contribution of metricas_paper.tex. Order confirmed against
# amh_montoril.cpp: theta1_t uses theta2_{t-1} (the value BEFORE this
# iteration's update), not the freshly-updated theta2_t. ----
subsample_particles <- function(samples, N, seed) {
	set.seed(seed)
	S_full <- length(samples$theta_01)
	idx <- sample.int(S_full, N, replace = (N > S_full))
	list(theta_01 = samples$theta_01[idx], theta_02 = samples$theta_02[idx],
		 W1 = samples$W1[idx], W2 = samples$W2[idx])
}

pf_predictive_loglik <- function(y, theta_01_samples, theta_02_samples, W1_samples, W2_samples, seed = 1) {
	set.seed(seed)
	Tt <- length(y)
	S <- length(theta_01_samples)

	loglik_t <- numeric(Tt)
	theta1 <- theta_01_samples
	theta2 <- theta_02_samples
	W1 <- W1_samples
	W2 <- W2_samples

	for (t in seq_len(Tt)) {
		# predict -- theta1_t uses theta2_{t-1}
		theta1 <- theta1 + theta2 + rnorm(S, 0, sqrt(W1))
		theta2 <- theta2 + rnorm(S, 0, sqrt(W2))

		# weight
		w <- dpois(y[t], exp(theta1))
		loglik_t[t] <- log(mean(w))

		# resample -- particles carry their own (W1[i], W2[i]) forward
		idx <- sample.int(S, S, replace = TRUE, prob = w)
		theta1 <- theta1[idx]
		theta2 <- theta2[idx]
		W1 <- W1[idx]
		W2 <- W2[idx]
	}
	sum(loglik_t)
}


# ---- resolve_n_cores(): physical core count (identical convention to
# real_data_run.R -- not parallel::detectCores(logical=FALSE), unreliable
# on Euler). ----
resolve_n_cores <- function(n_cores, n_needed) {
	if (!is.null(n_cores)) return(min(n_cores, n_needed))
	phys <- NA_integer_
	if (Sys.info()[["sysname"]] == "Linux" && nzchar(Sys.which("lscpu"))) {
		out <- tryCatch(system("lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#'", intern = TRUE),
						 error = function(e) character(0))
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


# ---- Full task grid: one task per (method, replicate). ----
task_grid <- do.call(rbind, lapply(methods_grid_ref, function(m) {
	data.frame(method = m, replicate_id = seq_len(N_REPLICATES), stringsAsFactors = FALSE)
}))
n_tasks <- nrow(task_grid)

n_cores_used <- resolve_n_cores(20, n_tasks)  # hard-capped to match #PBS -l select=1:ncpus=20 --
                                               # DO NOT pass NULL here: with n_tasks this large,
                                               # automatic detection is not capped by n_needed and
                                               # can report the whole node's physical cores rather
                                               # than the cgroup-restricted allocation, oversubscribing
                                               # and crashing the job. Update both this value and the
                                               # PBS script together if you change ncpus.
printf("Running %d tasks (%d methods x %d replicates) using %d core(s)", n_tasks, length(methods_grid_ref), N_REPLICATES, n_cores_used)

# ---- FORK, not PSOCK: unlike real_data_run.R's run_chain_units() (PSOCK,
# needed there for the compiled Rcpp samplers), this script is pure R with
# no compiled-code dependency, so there is no reason to pay for PSOCK's
# network-socket handshake between the master and each worker. The
# previous run failed with "invalid connection" / sendData.SOCKnode right
# after all 20 workers reported starting -- i.e. a socket handshake
# problem on this node, not an error in the task code itself (no worker
# ever printed an error, even with outfile=""). parallel::mclapply() uses
# fork() directly (Linux/Unix only, fine on this cluster) -- no sockets,
# no clusterExport needed (forked children inherit the parent's memory via
# copy-on-write, so all_samples/y/Tt/N_COMMON/the two functions are already
# visible), and no handshake step that can fail this way.
t0 <- Sys.time()
task_results_list <- parallel::mclapply(seq_len(n_tasks), function(i) {
	m <- task_grid$method[i]
	r <- task_grid$replicate_id[i]
	seed <- 10000 * match(m, c("montoril", "pg_apf", "sir_laplace", "sir_collapsed", "stan")) + r

	ll <- tryCatch({
		sub <- subsample_particles(all_samples[[m]], N_COMMON, seed = seed)
		pf_predictive_loglik(y, sub$theta_01, sub$theta_02, sub$W1, sub$W2, seed = seed + 5000000)
	}, error = function(e) {
		message(sprintf("Task %d (method=%s, replicate=%d) FAILED: %s", i, m, r, conditionMessage(e)))
		NA_real_
	})

	data.frame(method = m, replicate_id = r, predictive_loglik = ll)
}, mc.cores = n_cores_used, mc.set.seed = TRUE)
task_results <- do.call(rbind, task_results_list)
elapsed_total <- as.numeric(Sys.time() - t0, units = "secs")
printf("All %d tasks done in %.1fs (%.1f min)", n_tasks, elapsed_total, elapsed_total / 60)


# ---- Aggregate: mean, SD (per-replicate noise), SE of the mean, per method ----
summary_df <- do.call(rbind, lapply(methods_grid_ref, function(m) {
	vals <- task_results$predictive_loglik[task_results$method == m]
	data.frame(
		method = m, S_full = length(all_samples[[m]]$theta_01), N = N_COMMON, n_replicates = length(vals),
		predictive_loglik_mean = mean(vals), predictive_loglik_sd = sd(vals),
		predictive_loglik_se = sd(vals) / sqrt(length(vals)),
		stringsAsFactors = FALSE
	)
}))
rownames(summary_df) <- NULL

saveRDS(task_results, file = sprintf("%s/predictive_loglik_pf_replicates.rds", path_results))
saveRDS(summary_df, file = sprintf("%s/predictive_loglik_pf.rds", path_results))
write.csv(summary_df, file = sprintf("%s/predictive_loglik_pf.csv", path_results), row.names = FALSE)

printf("Done. Summary:")
print(summary_df)
