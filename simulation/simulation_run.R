# simulation/simulation_run.R
#
# Full simulation grid (Estrategia C):
#   - R=200 for all five methods (montoril, pg_apf, sir_laplace,
#     sir_collapsed, stan) -- stan's R was raised from an initial R=50
#     reference subsample to the full R=200 once compute budget allowed it;
#     R_config below reflects that final value, matching the committed
#     results/simulation/ output.
#   - 4 functions x 4 Tt values
#
# Architecture mirrors calibration_phase.R: single file, re-sourced by
# batchtools inside every job (source = "simulation_run.R"), with the same
# guard mechanism to prevent recursive dispatch. See that file's comments
# for the detailed rationale behind each guard; not re-explained here
# except where this script differs.
#
# Two independent choices, both in the "Run control" block below:
#
# 1. run_mode: "local" or "cluster".
#      "local"   -- runs every task of the (possibly subset) job grid in
#                   THIS session, all in one maximally-parallel batch
#                   using every physical core available on this machine
#                   (resolve_n_cores(), same lscpu-based detection as
#                   calibration_phase.R -- NOT parallel::detectCores(),
#                   which was confirmed unreliable on Euler during the
#                   ncpus=40 pilot: it reported 40 logical threads as 40
#                   "physical" cores on a 20-physical-core node).
#      "cluster" -- submits every chunk of the (possibly subset) job grid
#                   as ONE PBS job each, via batchtools, using
#                   simulation_pbs.tmpl (ncpus=20, place=excl). Each job
#                   runs run_tasks() for its own chunk's ~20 tasks, in
#                   parallel on that node's physical cores. submitJobs()
#                   returns as soon as jobs are dispatched.
#
# 2. grid_subset: NULL (run the full job grid) or a subset filter applied
#      to it (one row per CHUNK, not per task) -- e.g. only category ==
#      "pesado" (stan only), only Tt == 1600, or a single chunk_id, to run
#      a smaller batch without editing the grid definition. Applies
#      identically in both run_mode values.
#
# Chunking (cluster mode): k* = 20 tasks/job -- 1 core/task, 20 physical
# cores/node (2 sockets x 10 cores, HT on -- confirmed via lscpu + the
# ncpus=40 pilot; requesting ncpus=40 in the PBS resources counts logical
# threads, not physical cores). Chunks are homogeneous in (category, Tt),
# not just category: cost varies up to ~17x across Tt within stan alone,
# so mixing Tt within one chunk would leave workers idle waiting for the
# slowest task once the fast ones finish.
#
# NOTE: the previous multi-node "tarefas estrela" design (a second,
# heavier checkpoint level storing full chain history for a handful of
# (f, Tt) combinations) has been removed -- superseded by the full raw
# data already saved during calibration for Tt in {200, 1600}, both
# functions, all methods (see planning discussion). The SSH-based
# make_pbs_cluster() (multi-node, taskset-pinned workers) that supported
# spanning several nodes in one job has also been removed: chunking now
# allocates exactly 1 node per job (select=1), for which the plain local
# parallel::makeCluster() path (already proven by calibration_phase.R's
# non-stan branch) is simpler and avoids re-exposing the SSH-related bugs
# found and fixed during the ncpus=40 pilot (malformed rscript=, core
# dedup collapsing across sockets, do.call(c, ...) dropping the cluster
# class).

.defs_only_flag <- exists("SIMULATION_DEFS_ONLY", inherits = FALSE) &&
                    isTRUE(SIMULATION_DEFS_ONLY)

# Survives rm() below the same way .defs_only_flag does -- ensures the
# grid summary (printed further down) only shows once per R process, not
# once per internal re-source (loadRegistry()/makeRegistry() each
# re-source this whole file, so a single real invocation of this script
# prints it 3x otherwise).
.grid_already_printed <- exists(".grid_already_printed", inherits = FALSE) &&
                          isTRUE(.grid_already_printed)

rm(list = setdiff(ls(), c(".defs_only_flag", ".grid_already_printed")))
options(error = function() traceback(2))

setwd(dirname(this.path::this.path()))

# IMPORTANT: load_all() compiles C++ in DEBUG mode by default (-g -O0),
# silently overriding any -O3 in PoissonLTDM/src/Makevars, causing a
# ~4-5x slowdown of every sampler -- critical here since this is the main,
# long-running simulation.
pkgload::load_all("../PoissonLTDM", debug = FALSE, quiet = TRUE)

printf <- function(...) cat(paste(sprintf(...), "\n"))

# ---- archive_registry(): same as calibration_phase.R's version -- moves
# a finished registry into registry_archive/ with a timestamp instead of
# deleting it, preserving logs/history. Defined at top level (survives a
# loadRegistry()-triggered re-source of this file). ----
registry_archive_dir <- "registry_archive"
archive_registry <- function(reason, registry_dir) {
	dir.create(registry_archive_dir, showWarnings = FALSE, recursive = TRUE)
	archived_path <- file.path(registry_archive_dir,
		paste0("registry_simulation_", format(Sys.time(), "%Y-%m-%d_%H%M%S")))
	printf("%s -- archiving to: %s", reason, archived_path)
	file.rename(registry_dir, archived_path)
}

# ==========================================================================
# ---- Run control ----
# ==========================================================================
#
# run_mode and grid_subset are defined in simulation_grid_config.R (same directory),
# not inline here -- so changing what runs (which grid subset, local vs
# cluster) only requires editing that small file directly on the cluster,
# never re-sending simulation_run.R itself. See simulation_grid_config.R for the syntax
# and commented examples.

source("simulation_grid_config.R", local = FALSE)

# ==========================================================================

path_data <- "../data/simulated"
path_results <- "../results/simulation"
path_results_partial <- sprintf("%s/%s", path_results, "partial")
dir.create(path_results_partial, showWarnings = FALSE, recursive = TRUE)

verbose <- TRUE   # print information in the console

# N/burnin/K: per-method now (see R_config below), NOT global constants --
# this file previously had N<-10000/burnin<-1000/K<-50 as single global
# values used for every method's task, which never reflected the
# calibration study's per-method findings at all (amh_montoril needed
# N=110000/burnin=10000 to reach its efficiency plateau; pg_as needed
# N=22000/burnin=2000/K=200 to get rhat comfortably under 1.01; both were
# silently running at N=10000 here, 11x and 2.2x too small respectively --
# caught before launching production, not after).

# Adaptive Metropolis hyperparameters (montoril)
varsigma2_scal <- 0.02     # initial varsigma2
ac_ref <- 0.44             # acceptance ratio target

# SIR Laplace and SIR Collapsed hyperparameters
M_is <- 3             # Number of particles - IS for W1 integrated likelihood
M_irls_max <- 20
tol <- 1e-4

# Only for SIR Collapsed
R_prerun <- 3000      # pre-run iterations to calibrate CE proposals (phi1 and phi2)
M_sir <- 3            # Number of particles - SIR of theta1

# Prior hyperparameters
mu_01     <- 0        # theta_01 ~ N(mu_01, sigma2_01)
sigma2_01 <- 100
mu_02     <- 0        # theta_02 ~ N(mu_02, sigma2_02)
sigma2_02 <- 100
nu_01  <- 2           # phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
eta_01 <- 0.01
nu_02  <- 2           # phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
eta_02 <- 0.0001

# Initialization (other parameters is funciton of Tt, see the main loop)
theta1_tilde_scal <- 0 # sir_laplace and sir_collapsed


#####
# Build the task grid (Estrategia C) and chunk it (k* = 20)

# ---- Per-method category, replica count, and sampler hyperparameters ----
#
# N/burnin/K per method, not global -- each row here is the config the
# calibration study actually validated/decided for that method:
#   - montoril:      N=110000/burnin=10000, the efficiency-plateau pick
#     (peak ESS/s among the 5 configs tested; Tt=1600's rhat>1.01 is a
#     documented structural finding for this method, not something a
#     larger N resolves -- see the rhat-vs-N extrapolation from the
#     N=520000 calibration point).
#   - pg_apf:        N=22000/burnin=2000/K=200, prioritizing the lowest
#     rhat achieved (1.005) over the cheaper K=100 config (rhat=1.008 at
#     roughly half the cost) -- explicit choice, not the cost-optimal one.
#   - sir_laplace/sir_collapsed/stan: N=11000/burnin=1000 is the only
#     config calibration tested for these three, so it is the
#     calibration-validated config by construction, not a choice among
#     alternatives the way montoril/pg_apf's configs were.
R_config <- data.frame(
	method   = c("montoril", "pg_apf",  "sir_laplace", "sir_collapsed", "stan"),
	category = c("leve",     "medio",   "leve",        "leve",          "pesado"),
	R        = c(200,        200,       200,           200,             200),
	N        = c(110000,     22000,     11000,         11000,           11000),
	burnin   = c(10000,      2000,      1000,          1000,            1000),
	K        = c(NA,         200,       NA,            NA,              NA),
	stringsAsFactors = FALSE
)

# Full reference grids for the seed formula -- stay complete regardless of
# grid_subset, so the formula can place ANY (method, f, Tt) at its correct
# index, not just the ones in the current subset.
Tt_grid_ref        <- c(200, 400, 800, 1600)
functions_grid_ref <- c("constant", "linear", "quadratic", "sinusoidal")
methods_grid_ref   <- c("montoril", "pg_apf", "sir_laplace", "sir_collapsed", "stan")

task_grid <- do.call(rbind, lapply(seq_len(nrow(R_config)), function(i) {
	expand.grid(
		Tt = Tt_grid_ref, f = functions_grid_ref, replica = seq_len(R_config$R[i]),
		method = R_config$method[i], category = R_config$category[i],
		N = R_config$N[i], burnin = R_config$burnin[i], K = R_config$K[i],
		stringsAsFactors = FALSE
	)
}))

	# ---- Seed formula ----
	#
	#   seed(m, g, tau, r) = m*1e7 + g*1e6 + tau*1e5 + 10000 + r*10
	#
	# where m = match(method, methods_grid_ref), g = match(f,
	# functions_grid_ref), tau = match(Tt, Tt_grid_ref), r = replica.
	#
	# Deliberately matches calibration_run.R's task_seed_base() as closely
	# as possible:
	#
	#   seed_calib(m, g, tau, c, r) = m*1e7 + g*1e6 + tau*1e5 + c*1e3 + r*10
	#
	# Same leading multipliers for method/f/Tt (1e7/1e6/1e5) -- these were
	# kept byte-for-byte identical on purpose, so the leading digits of a
	# seed identify (method, f, Tt) the same way in both files, and the
	# same replica*10 trailing term, reused verbatim from calibration's
	# own formula rather than invented fresh.
	#
	# The one deliberate difference is the "+10000" in place of
	# calibration's "c*1e3" (c = config_idx, calibration's per-method
	# index into its own candidate N/burnin/K configs -- a dimension
	# production does not have, since each method here runs exactly ONE
	# config, the one calibration decided on). This is NOT a stylistic
	# choice or a "config_idx=0" placeholder -- it exists to rule out a
	# real bias, confirmed in discussion:
	#
	#   Naively using config_idx=0 (i.e. just m*1e7+g*1e6+tau*1e5+r*10,
	#   dropping the c*1e3 term entirely) does NOT prevent collisions,
	#   because calibration's own replica is FIXED AT 1 for the whole
	#   calibration phase (see calibration_run.R) -- so calibration's
	#   entire occupied sub-block below tau's digit is small, spanning
	#   only c*1000 + 1*10 + chain_offset(0..2) for c = 1..6 (pg_as has
	#   the most candidate configs, 6), i.e. values 1010-6012. Production's
	#   replica genuinely varies (1..200 here), so replica*10 ALONE spans
	#   10-2000, which overlaps that exact 1010-6012 band once replica
	#   exceeds ~100 -- confirmed as a real, reproducible collision (18
	#   exact seed matches across the full production grid against every
	#   calibration seed, all of them production replica=101 landing
	#   exactly on some method/f/Tt's calibration config_idx=1, chain=1
	#   seed). This is a genuine look-ahead/data-snooping risk, not merely
	#   cosmetic: for the 3 methods calibration tested only ONE config for
	#   (sir_laplace, sir_collapsed, stan), that single config IS the one
	#   production uses, so config_idx=1 there is not a hypothetical
	#   collision target -- it is exactly the calibration chain whose
	#   good convergence at that specific seed is *why* that config was
	#   chosen. Reusing that same seed for production's replica=1 would
	#   hand that one replica a selection advantage the other 199 never
	#   had.
	#
	#   The fix is not a bigger/smaller config_idx marker -- no single
	#   marker value avoids the overlap, since it is replica*10 itself
	#   (not a nonzero config_idx digit) that encroaches on calibration's
	#   occupied range. Adding a flat +10000 offset instead pushes
	#   production's ENTIRE sub-block (10010-12000 for replica=1..200)
	#   above calibration's maximum possible occupied value (config_idx
	#   up to 6, or generously up to 9, gives at most ~9012) -- a
	#   structural, not incidental, separation. Verified by brute-force
	#   intersection of the full sets: all 13600 production task seeds
	#   against all 204 calibration seeds (68 tasks x 3 chains) share ZERO
	#   values.
	#
	# Concrete examples (f=constant, Tt=1600):
	#   amh_montoril: calibration's chosen config (N=110000) is
	#     config_idx=3 there -> seed 11403010 (chain 1). Production
	#     replica=1 -> seed 11410010. Production replica=101 -> 11411010.
	#     Neither collides with 11403010, nor with any of the other 4
	#     candidate configs' seeds (11401010, 11402010, 11404010, 11405010).
	#   sir_laplace: calibration's ONLY config (N=11000, the one production
	#     also uses) is config_idx=1 -> seed 31401010 (chain 1). Production
	#     replica=1 -> seed 31410010 -- close in magnitude but NOT equal,
	#     by the +10000 construction above, not by chance.
	#
	# Confirmed this poses no RNG-range problem either: both phases only
	# ever call R's own set.seed(seed) (never pass a seed into the Rcpp
	# samplers directly), and the samplers themselves use RNGScope +
	# R::norm_rand() (confirmed directly in sir_laplace.cpp's source),
	# i.e. R's own centralized RNG state -- so the only limit that matters
	# is R's own .Machine$integer.max (2147483647). The largest seed this
	# formula can produce (method=5, f=4, Tt=4, replica=200) is 54410500,
	# about 2.5% of that limit.
	task_grid$seed <- match(task_grid$method, methods_grid_ref) * 1e7 +
						match(task_grid$f, functions_grid_ref) * 1e6 +
						match(task_grid$Tt, Tt_grid_ref) * 1e5 +
						10000 +
						task_grid$replica * 10

CHUNK_SIZE <- 20L

# ---- Per-(method, Tt) chunk size override ----
#
# CHUNK_SIZE=20 assumes 20 tasks can run concurrently on one 128GB-RAM
# node (see simulation_pbs.tmpl: place=excl reserves the WHOLE node
# regardless of the ncpus requested, so all 20 tasks in a chunk share
# that one node's RAM budget). This holds for every method EXCEPT
# montoril at Tt=1600: N=110000 there means each task's raw theta1_hist+
# theta2_hist alone is N*Tt*8*2 bytes = ~2.82GB, so 20 of them running at
# once is ~56.3GB of raw history BEFORE R's own overhead (temporary
# copies, the parallel diagnostics computation, gc lag) -- confirmed as
# the actual cause of a real production failure: 40 consecutive chunk
# jobs (batchtools job.id 297-336) all failed identically with "Error in
# unserialize(socklist[[n]]) : error reading from connection" right after
# "Running 20 task(s) using 20 core(s)", with zero task-level progress in
# any of them -- the exact signature of an OOM-killed parallel worker
# (its socket becomes unreadable mid-read). 40 is not a coincidence: it
# is EXACTLY montoril's Tt=1600 chunk count (200 replicas x 4 functions =
# 800 tasks / CHUNK_SIZE=20 = 40 chunks), and no other method/Tt
# combination showed any failures. pg_apf (N=22000) and sir_laplace/
# sir_collapsed (N=11000) at the same Tt=1600 are ~5-10x lighter per task
# (0.56GB and 0.28GB respectively) and were never at risk.
#
# The fix: halve concurrency (10 tasks/chunk instead of 20) specifically
# for montoril+Tt=1600, giving each task roughly double the RAM headroom
# on the same 128GB node (10 tasks x ~2.82GB raw = ~28.2GB, leaving ample
# room for R's overhead) -- at the cost of doubling wall-clock time and
# PBS job count for just this one (method, Tt) slice, not the whole grid.
chunk_size_for <- function(method, Tt) {
	ifelse(method == "montoril" & Tt == 1600, 10L, CHUNK_SIZE)
}

task_grid <- task_grid[order(task_grid$category, task_grid$Tt, task_grid$method,
							  task_grid$f, task_grid$replica), ]
rownames(task_grid) <- NULL

# ncpus per row: the chunk size for THIS task's (method, Tt) -- almost
# always CHUNK_SIZE, except montoril/Tt=1600 (see chunk_size_for()).
task_grid$ncpus <- chunk_size_for(task_grid$method, task_grid$Tt)

# Grouped by (category, Tt, method) now, not just (category, Tt) -- so a
# method with an overridden chunk_size_for() gets ITS OWN chunk sequence,
# sized correctly, rather than being mixed into chunks sized for whatever
# CHUNK_SIZE the rest of that category/Tt combination uses. (This also
# generalizes cleanly what already held true by coincidence before: since
# montoril's task count per Tt, 800, happens to be an exact multiple of
# the old global CHUNK_SIZE=20, its chunks were already "pure" -- never
# mixed with sir_laplace/sir_collapsed's "leve" tasks at the same Tt --
# but that was incidental to task ordering + divisibility, not something
# the grouping key itself guaranteed until now.)
task_grid$chunk_index <- ave(seq_len(nrow(task_grid)),
							  paste(task_grid$category, task_grid$Tt, task_grid$method),
							  FUN = function(idx) ceiling(seq_along(idx) / task_grid$ncpus[idx[1]]))

task_grid$chunk_id <- sprintf("%s_Tt%04d_%s_%03d", task_grid$category, task_grid$Tt, task_grid$method, task_grid$chunk_index)

grid_size <- nrow(task_grid)

# ---- Job-level grid (one row per chunk/PBS job) ----
job_grid <- unique(task_grid[, c("chunk_id", "category", "Tt", "ncpus")])
job_grid <- job_grid[order(job_grid$category, job_grid$Tt, job_grid$chunk_id), ]
rownames(job_grid) <- NULL

if (!is.null(grid_subset)) {
	job_grid <- grid_subset(job_grid)
	rownames(job_grid) <- NULL
}

n_chunks <- nrow(job_grid)
if (verbose && !.grid_already_printed) {
	printf("Job grid: %d chunk(s)%s", n_chunks, if (!is.null(grid_subset)) " (subset applied)" else "")
	for (cat in unique(job_grid$category)) {
		printf("  %s: %d chunk(s)", cat, sum(job_grid$category == cat))
	}
}
assign(".grid_already_printed", TRUE, envir = .GlobalEnv)

# ---- Walltime per category ----
#
# CRITICAL, confirmed by Euler support via email: walltime <= 6h routes to
# "fila paralela curta", which had only 2 nodes allocated at the time of
# that reply (max ~4 concurrent place=excl jobs, and in practice as few as
# 2 were observed, likely because the other node was already occupied by
# someone else). walltime > 6h routes to "fila paralela longa" instead (14
# free nodes at the time of that reply), where this Category-2 account's
# 160-ncpus quota is the actual binding constraint -- 160/20 = 8 concurrent
# place=excl jobs, which is the concurrency the whole CHUNK_SIZE-based
# chunking design here is built around. Every value below is deliberately
# > 6h for this reason -- NOT because any single chunk is expected to take
# that long (each is still a generous margin over the worst measured
# per-task time within that category; see planning discussion), but
# because requesting walltime <= 6h, even for the fastest chunks, would
# silently confine ALL of them to the short queue's ~2-4-node ceiling
# regardless of how quickly they actually finish. There is no downside to
# requesting more walltime than a job needs (it releases its node as soon
# as it's done either way) -- the downside only runs the other direction,
# requesting too little.
walltime_hours_for <- function(category) {
	switch(category,
		leve   = 6.5,
		medio  = 7,
		pesado = 8,
		stop(sprintf("Unknown category: %s", category))
	)
}


task_result_filename <- function(task) {
	sprintf("%s/%s_%s_%s_%s.rds", path_results_partial, task$method, task$f, task$Tt, task$replica)
}

data_filename <- function(task) {
	sprintf("%s/%s_%s_%s.rds", path_data, task$f, task$Tt, task$replica)
}


run_task <- function(task) {

	Tt <- task$Tt
	replica <- task$replica
	f <- task$f
	method <- task$method
	seed <- task$seed
	N <- task$N
	burnin <- task$burnin
	K <- task$K

	task_result_file <- task_result_filename(task)

	# Checkpoint: skip the task if it is already done
	if (file.exists(task_result_file)) {
		if (verbose) printf("\tTask already run, loading results")
		result <- readRDS(task_result_file)
	} else {

		# Load the data
		data_file <- data_filename(task)
		data <- readRDS(data_file)
		y <- data$y
		theta1_true <- data$theta

		# Set the initial values
		theta1 <- log(y + 0.5)
		theta2 <- c(diff(theta1), 0)
		theta_01 <- theta1[1]
		theta_02 <- theta2[1]
		W1 <- var(diff(theta1))
		W2 <- var(diff(theta2))
		
		# Run the task
		set.seed(seed)

		if (method == 'stan') {
			# Load or compile the model
			if (file.exists("../cache/poisson_ltdm.rds")) {
				model <- readRDS("../cache/poisson_ltdm.rds")
			} else {
				if (verbose) printf("Building the Stan model")
				file <- "../PoissonLTDM/inst/stan/poisson_ltdm.stan"
				model <- rstan::stan_model(file = file, model_name = "PoissonLTDM")
				saveRDS(model, file = "../cache/poisson_ltdm.rds")
			}
		}

		execution_bench <- system.time({
			res <- switch(method,

				# The four branches below call the Rcpp samplers directly
				# (<method>_cpp), matching calibration_phase.R's
				# proven-working signatures -- the sample_<method>()
				# wrapper names used previously did not exist in the
				# compiled package.

				"montoril" = amh_montoril_cpp(y = y, N = N,
											   mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
											   nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
											   theta1 = theta1, theta2 = theta2,
											   theta_01 = theta_01, theta_02 = theta_02,
											   W1 = W1, W2 = W2, ac_ref = ac_ref, varsigma2 = rep(varsigma2_scal, Tt),
											   verbose = FALSE, print_every = 1000),

				"pg_apf" = pg_as_cpp(y = y, K = K, N = N,
									  mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
									  nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
									  theta1 = theta1, theta2 = theta2,
									  theta_01 = theta_01, theta_02 = theta_02,
									  W1 = W1, W2 = W2, verbose = FALSE, print_every = 1000),

				"sir_laplace" = sir_laplace_cpp(y = y, N = N,
												 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
												 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
												 theta1 = theta1, theta2 = theta2,
												 theta_01 = theta_01, theta_02 = theta_02,
												 W1 = W1, W2 = W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol,
												 verbose = FALSE, print_every = 1000),

				"sir_collapsed" = sir_collapsed_cpp(y = y, R_prerun = R_prerun, N = N,
													 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
													 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
													 theta1 = theta1, theta2 = theta2, theta1_tilde = rep(theta1_tilde_scal, Tt),
													 theta_01 = theta_01, theta_02 = theta_02, W1 = W1, W2 = W2,
													 M_is_lik = M_is, M_sir_theta1 = M_sir, M_irls_max = M_irls_max, tol = tol,
													 verbose = FALSE, print_every = 1000),

				# sample_stan()'s real signature (sampler_stan.R) takes
				# N_chains/n_cores/chain_inits (a list of per-chain init
				# lists), not theta1/theta2/W1/W2/theta_01/theta_02 as
				# direct arguments -- the previous call didn't match this
				# at all ("unused arguments (theta1, theta2)", confirmed
				# via a real cluster job error). For 1 task = 1 chain, we
				# build a single-element chain_inits from the same
				# variables already in scope. The return value is also
				# nested (list(results=<list of N_chains>, fit,
				# elapsed_time)), not flat like the other 4 samplers --
				# $results[[1]] unwraps it to match what the rest of this
				# function expects (res$theta_01_hist, res$elapsed_time, etc.).
				"stan" = sample_stan(model, y, N, burnin, seed,
									 N_chains = 1, n_cores = 1,
									 chain_inits = list(list(theta_01 = theta_01, theta_02 = theta_02,
															  theta1 = theta1, theta2 = theta2,
															  W1 = W1, W2 = W2)),
									 mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
									 nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
									 verbose = FALSE)$results[[1]],

				stop(sprintf("Unknow method: %s", method))
			)
		})
		elapsed_time <- execution_bench[["user.self"]]
		if (verbose) printf("All tasks ran")

		# Common
		theta_01_hist <- res$theta_01_hist
		theta_02_hist <- res$theta_02_hist
		W1_hist <- res$W1_hist
		W2_hist <- res$W2_hist
		theta1_hist <- res$theta1_hist
		theta2_hist <- res$theta2_hist

		if (method == "stan") {
			# overwrite the previous elapsed time computed by system.time()
			elapsed_time <- res$elapsed_time
		}

		# Samples without burn-in
		theta_01_samples <- theta_01_hist[-(1:burnin)]
		theta_02_samples <- theta_02_hist[-(1:burnin)]
		W1_samples <- W1_hist[-(1:burnin)]
		W2_samples <- W2_hist[-(1:burnin)]
		theta1_samples <- theta1_hist[-(1:burnin), ]
		theta2_samples <- theta2_hist[-(1:burnin), ]

		# Effective Sample Size (bulk/tail; R_hat/Geweke not assessed per
		# replica -- only during the pilot/calibration stage)
		ess_theta_01 <- metrics_ess_single_chain(theta_01_samples)$ess_bulk
		ess_theta_02 <- metrics_ess_single_chain(theta_02_samples)$ess_bulk
		ess_W1 <- metrics_ess_single_chain(W1_samples)$ess_bulk
		ess_W2 <- metrics_ess_single_chain(W2_samples)$ess_bulk
		ess_theta1 <- metrics_ess_single_chain(theta1_samples)$ess_bulk
		ess_theta1_mean <- mean(ess_theta1)
		ess_theta1_min <- min(ess_theta1)

		ess_theta2 <- metrics_ess_single_chain(theta2_samples)$ess_bulk
		ess_theta2_mean <- mean(ess_theta2)
		ess_theta2_min <- min(ess_theta2)

		# Effective Sample Size per second
		ess_sec_theta_01 <- ess_theta_01 / elapsed_time
		ess_sec_theta_02 <- ess_theta_02 / elapsed_time
		ess_sec_W1 <- ess_W1 / elapsed_time
		ess_sec_W2 <- ess_W2 / elapsed_time
		ess_sec_theta1_mean <- ess_theta1_mean / elapsed_time
		ess_sec_theta2_mean <- ess_theta2_mean / elapsed_time

		# Fit metrics
		W1_mean <- mean(W1_samples)
		W1_median <- median(W1_samples)
		W1_var <- var(W1_samples)
		W1_ci <- metrics_W_ci(W1_samples, 0.95)
		W1_ci_lower <- W1_ci$ci_lower
		W1_ci_upper <- W1_ci$ci_upper

		W2_mean <- mean(W2_samples)
		W2_median <- median(W2_samples)
		W2_var <- var(W2_samples)
		W2_ci <- metrics_W_ci(W2_samples, 0.95)
		W2_ci_lower <- W2_ci$ci_lower
		W2_ci_upper <- W2_ci$ci_upper

		theta1_mean <- colMeans(theta1_samples)
		theta2_mean <- colMeans(theta2_samples)
		lambda_mean <- exp(theta1_mean)

		log_lik <- metrics_loglik(y, lambda_mean)
		rmse_theta1 <- metrics_rmse(theta1_mean, theta1_true)
		mae_theta1 <- metrics_mae(theta1_mean, theta1_true)

		theta1_ci <- metrics_theta_ci(theta1_samples, 0.95)
		theta1_ci_lower <- theta1_ci$ci_lower
		theta1_ci_upper <- theta1_ci$ci_upper

		# Coverage indicator C(t,r) (Eq. 7, metricas_theta1.tex), reduced to
		# a scalar rate per replica -- Bug 3's fix: keeping the raw
		# length-Tt ci_lower/ci_upper vectors directly as data.frame()
		# columns silently recycled every other (scalar) column to Tt rows.
		coverage_rate_theta1 <- mean(theta1_true >= theta1_ci_lower & theta1_true <= theta1_ci_upper)
		theta1_ci_width_mean <- mean(theta1_ci_upper - theta1_ci_lower)

		task_result <- data.frame(
			f = f, Tt = Tt, method = method, replica = replica,
			elapsed_time = elapsed_time,

			ess_theta_01 = ess_theta_01, ess_theta_02 = ess_theta_02,
			ess_W1 = ess_W1, ess_W2 = ess_W2,
			ess_theta1_mean = ess_theta1_mean, ess_theta1_min = ess_theta1_min,
			ess_theta2_mean = ess_theta2_mean, ess_theta2_min = ess_theta2_min,

			ess_sec_theta_01 = ess_sec_theta_01, ess_sec_theta_02 = ess_sec_theta_02,
			ess_sec_W1 = ess_sec_W1, ess_sec_W2 = ess_sec_W2,
			ess_sec_theta1_mean = ess_sec_theta1_mean, ess_sec_theta2_mean = ess_sec_theta2_mean,

			W1_mean = W1_mean, W1_median = W1_median, W1_var = W1_var,
			W1_ci_lower = W1_ci_lower, W1_ci_upper = W1_ci_upper,
			W2_mean = W2_mean, W2_median = W2_median, W2_var = W2_var,
			W2_ci_lower = W2_ci_lower, W2_ci_upper = W2_ci_upper,

			log_lik = log_lik, rmse_theta1 = rmse_theta1, mae_theta1 = mae_theta1,
			coverage_rate_theta1 = coverage_rate_theta1,
			theta1_ci_width_mean = theta1_ci_width_mean,

			stringsAsFactors = FALSE
		)

		# List-columns (length Tt each, one row per replica): kept for
		# simulation_aggregate.R -- Bias(t)/width(t) and the estimator
		# dispersion band (Sec. 3.7, metricas_theta1.tex) need these per-t
		# vectors; they are deliberately NOT part of the scalar columns
		# above (a CSV cannot hold vector cells cleanly).
		task_result$theta1_mean <- list(theta1_mean)
		task_result$theta1_ci_lower <- list(theta1_ci_lower)
		task_result$theta1_ci_upper <- list(theta1_ci_upper)

		saveRDS(task_result, file = task_result_file)

		result <- task_result
	}

	return(result)
}


#####
# Parallel execution

# ---- resolve_n_cores(): physical core count, same lscpu-based detection
# as calibration_phase.R -- NOT parallel::detectCores(logical=FALSE),
# confirmed unreliable on Euler (reported 40 instead of 20 physical cores
# during the ncpus=40 pilot). n_needed caps the result (never spin up more
# workers than there are tasks to run). ----
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

# ---- run_tasks(): runs every row of `tasks` in parallel, using a local
# PSOCK cluster (parallel::makeCluster()) -- the same mechanism already
# proven by calibration_phase.R's non-stan branch. n_cores = NULL (local
# mode) auto-detects and uses every physical core available, capped by the
# number of tasks; n_cores = CHUNK_SIZE (cluster mode, inside a PBS job)
# matches the node's ncpus allocation exactly. ----
run_tasks <- function(tasks, n_cores = NULL) {

	n <- nrow(tasks)
	n_cores_used <- resolve_n_cores(n_cores, n)
	if (verbose) printf("Running %d task(s) using %d core(s)", n, n_cores_used)

	current_wd <- getwd()
	`%dorng%` <- doRNG::`%dorng%`

	cl <- parallel::makeCluster(n_cores_used)
	on.exit(parallel::stopCluster(cl), add = TRUE)
	doParallel::registerDoParallel(cl)
	parallel::clusterExport(cl, "current_wd", envir = environment())
	# foreach's auto-export only inspects the immediate calling frame
	# (run_tasks()'s own local environment) -- it does not walk further up
	# the lexical chain to where run_task() and its dependencies (N,
	# mu_01, ..., task_result_filename, etc.) are actually defined. That
	# worked when this logic lived at the script's top level (same scope
	# as those globals); wrapping it in run_tasks() breaks the implicit
	# auto-export, so those names must be exported explicitly.
	parallel::clusterExport(cl, ls(envir = environment(run_task)), envir = environment(run_task))

	parallel::clusterEvalQ(cl, {
		setwd(current_wd)
		pkgload::load_all("../PoissonLTDM", debug = FALSE, quiet = TRUE)
	})

	foreach::foreach(i = seq_len(n), .combine = rbind) %dorng% {
		run_task(tasks[i, ])
	}
}


# ==========================================================================
# ---- Dispatch: run_mode == "local" or "cluster" ----
# ==========================================================================

running_inside_pbs_job <- Sys.getenv("PBS_NODEFILE") != ""
making_registry <- Sys.getenv("SIMULATION_MAKING_REGISTRY") != ""

if (!running_inside_pbs_job && !.defs_only_flag && !making_registry) {

if (run_mode == "local") {

	# Every task of the (filtered) job grid, in ONE maximally-parallel
	# batch using every physical core available on this machine -- chunk
	# boundaries only matter for PBS job packaging (cluster mode), not
	# here.
	task_grid_filtered <- task_grid[task_grid$chunk_id %in% job_grid$chunk_id, ]
	printf("Running %d task(s) locally, across %d chunk(s).", nrow(task_grid_filtered), n_chunks)
	df_results <- run_tasks(task_grid_filtered, n_cores = NULL)

} else if (run_mode == "cluster") {

	library(batchtools)

	Sys.setenv(SIMULATION_MAKING_REGISTRY = "TRUE")
	on.exit(Sys.unsetenv("SIMULATION_MAKING_REGISTRY"), add = TRUE)

	registry_dir <- "registry_simulation"

	if (dir.exists(registry_dir)) {
		old_reg <- tryCatch(loadRegistry(registry_dir, writeable = FALSE),
							 error = function(e) NULL)
		registry_dir <- "registry_simulation"  # re-assign after the re-source

		if (is.null(old_reg)) {
			archive_registry(sprintf("Unreadable leftover registry directory: %s", registry_dir),
							  registry_dir)
		} else {
			# findNotDone() = done is NA OR error is set -- this lumps
			# genuinely-still-running jobs together with ones that already
			# TERMINATED without ever confirming success (e.g. the
			# "Expired" misclassification seen during testing: batch.id
			# vanished from the scheduler's active list, but done stayed
			# NA due to what looks like an NFS result-visibility race,
			# even though the job's real output was already written
			# correctly). findOnSystem() checks the opposite condition
			# (done is NA AND batch.id IS STILL in the scheduler's active
			# list) -- the correct signal for "is something truly still
			# running right now", which is what this guard is actually
			# meant to protect against.
			n_pending <- nrow(batchtools::findOnSystem(reg = old_reg))
			if (n_pending > 0) {
				stop(sprintf(
					"Registry '%s' still has %d job(s) actually active on the scheduler. ",
					registry_dir, n_pending),
					"Refusing to archive it automatically -- wait for those jobs to finish, ",
					"or archive/remove the directory yourself if you mean to abandon them.")
			}
			archive_registry(sprintf("Previous registry '%s' has no active jobs", registry_dir),
							  registry_dir)
		}
	}

	reg <- makeRegistry(
		file.dir = registry_dir,
		source = "simulation_run.R",
		seed = 1
	)

	reg$cluster.functions <- makeClusterFunctionsTORQUE("simulation_pbs.tmpl")

	# ---- Skip chunks whose tasks are ALL already checkpointed ----
	chunk_done <- vapply(job_grid$chunk_id, function(cid) {
		chunk_tasks <- task_grid[task_grid$chunk_id == cid, ]
		all(file.exists(task_result_filename(chunk_tasks)))
	}, logical(1))

	n_skipped <- sum(chunk_done)
	if (n_skipped > 0) {
		printf("Skipping %d already-checkpointed chunk(s) (no job submitted for them).", n_skipped)
	}
	job_grid_to_submit <- job_grid[!chunk_done, ]

	if (nrow(job_grid_to_submit) == 0) {
		printf("Every chunk in the (filtered) grid is already checkpointed -- nothing to submit.")
	} else {

		ids <- batchMap(
			fun = function(chunk_id) {
				chunk_tasks <- task_grid[task_grid$chunk_id == chunk_id, ]
				run_tasks(chunk_tasks, n_cores = chunk_tasks$ncpus[1])
			},
			chunk_id = job_grid_to_submit$chunk_id,
			reg = reg
		)
		# batchMap()'s returned ids has only a job.id column (not the mapped
		# chunk_id) -- but job.id order matches the input vector's order
		# exactly (confirmed empirically), so category/ncpus can be attached
		# directly by position, no join needed.
		ids$category <- job_grid_to_submit$category
		ids$ncpus <- job_grid_to_submit$ncpus

		# max.concurrent.jobs = floor(160 / ncpus) varies now too, since
		# ncpus is no longer a single global CHUNK_SIZE for every
		# submission -- montoril/Tt=1600's chunks (ncpus=10) get a higher
		# ceiling (16) than everything else (ncpus=20, ceiling 8), correctly
		# reflecting that twice as many of the smaller jobs fit in the same
		# 160-ncpus quota. Submission is grouped by (category, ncpus) pairs
		# instead of category alone, since submitJobs() requires uniform
		# resources within one call, and ncpus now varies WITHIN the "leve"
		# category (montoril/Tt=1600 vs. everything else in "leve").
		for (cat in unique(ids$category)) {
			for (ncpus_val in unique(ids$ncpus[ids$category == cat])) {

				grp_ids <- ids[ids$category == cat & ids$ncpus == ncpus_val, "job.id", drop = FALSE]
				max_concurrent_grp <- 160 %/% ncpus_val

				submitJobs(grp_ids,
						   resources = list(ncpus = ncpus_val, walltime_hours = walltime_hours_for(cat),
											 max.concurrent.jobs = max_concurrent_grp),
						   reg = reg)

				printf("Submitted %d job(s) for category '%s' (ncpus=%d, max.concurrent=%d).",
					   nrow(grp_ids), cat, ncpus_val, max_concurrent_grp)
			}
		}

		printf("Use check_progress.R to monitor (see /areas/euler-cluster-simulation.md).")
		printf("NOTE: submitJobs() blocks (in THIS R session) until every requested job has")
		printf("      been dispatched -- not until they finish. max.concurrent.jobs is a real,")
		printf("      global throttle (checked against the live scheduler via getBatchIds(),")
		printf("      not a per-category counter) -- it varies by ncpus group now (160 %%/%% ncpus),")
		printf("      160/20=8 for most chunks, 160/10=16 for montoril/Tt=1600's smaller ones --")
		printf("      but the true per-group cap is respected regardless of how many (category,")
		printf("      ncpus) groups are submitted in one run. This call can still block for a")
		printf("      long time (hours, for the larger categories), so run this inside")
		printf("      tmux/screen or with nohup, never directly in an SSH session that might")
		printf("      disconnect.")
	}

} else {
	stop(sprintf("Unknown run_mode: '%s'. Use \"local\" or \"cluster\".", run_mode))
}

}
