# Poisson Local Trend Dynamic Model
# Pilot script: test ncpus=40 / place=excl packing on a single node,
# using 40 real sir_laplace/Tt=200 tasks (the cheapest tasks in the grid).
#
# This does NOT modify simulation.R. The helper functions below are
# copied verbatim from simulation.R (make_pbs_cluster, filename helpers,
# a trimmed run_task) so the pilot exercises the same real code path
# (checkpointing, NFS writes, worker pinning) without running the full
# production grid.
#
# Author: Cleiton Moya de Almeida

options(error = function() traceback(2))
rm(list = ls())

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

pkgload::load_all("../PoissonLTDM", debug = FALSE)

path_data <- "../data/simulated"
path_results_partial <- "results/pilot_ncpus40"
dir.create(path_results_partial, showWarnings = FALSE, recursive = TRUE)

verbose <- TRUE
printf <- function(...) cat(paste(sprintf(...), "\n"))

cluster <- Sys.getenv("PBS_NODEFILE") != ""
if (!cluster) stop("This pilot is meant to run inside a PBS job (PBS_NODEFILE not set).")

if (verbose) {
	printf("PBS_QUEUE: %s", Sys.getenv("PBS_QUEUE"))
	printf("PBS_JOBID: %s", Sys.getenv("PBS_JOBID"))
	printf("hostname: %s", Sys.info()[["nodename"]])
	# Quick, cheap check (no SSH needed, since select=1): does ncpus=40
	# requested from PBS correspond to physical cores or logical
	# (hyperthreaded) threads? parallel::detectCores() reads /proc/cpuinfo
	# directly on this same node.
	printf("detectCores(logical=TRUE):  %d", parallel::detectCores(logical = TRUE))
	printf("detectCores(logical=FALSE): %d", parallel::detectCores(logical = FALSE))
	printf("--- lscpu (authoritative, kernel-reported) ---")
	lscpu_out <- tryCatch(system2("lscpu", stdout = TRUE), error = function(e) sprintf("lscpu failed: %s", conditionMessage(e)))
	cat(paste(lscpu_out, collapse = "\n"), "\n")
	printf("-----------------------------------------------")
}

# Same production hyperparameters used for sir_laplace (calibration-decided)
N <- 11000
burnin <- 1000
M_is <- 3
M_irls_max <- 20
tol <- 1e-4

mu_01 <- 0; sigma2_01 <- 100
mu_02 <- 0; sigma2_02 <- 100
nu_01 <- 2; eta_01 <- 0.01
nu_02 <- 2; eta_02 <- 0.0001
W1 <- 0.01; W2 <- 0.01
theta_01 <- 0; theta_02 <- 0
theta1_tilde_scal <- 0

# ---- Pilot task grid ----
# Grid size is NOT fixed at 40 here: it is set dynamically below, after
# make_pbs_cluster() reports how many physical-core workers it actually
# created for this node (which may be < 40 if ncpus=40 counts logical
# hyperthreaded threads rather than physical cores).
method <- "sir_laplace"
f <- "constant"
Tt <- 200
N_replicas_max <- 40  # upper bound, limited by how many replicas of
                       # constant/Tt=200 are already available in path_data

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

	theta1 <- numeric(Tt)
	theta2 <- numeric(Tt)

	task_result_file <- task_result_filename(task)

	if (file.exists(task_result_file)) {
		if (verbose) printf("\tTask already run, loading results")
		result <- readRDS(task_result_file)
	} else {

		data_file <- data_filename(task)
		if (!file.exists(data_file)) {
			stop(sprintf("Missing simulated data file: %s -- generate it before running this pilot, or reduce N_replicas_pilot to the number of replicas already available.", data_file))
		}
		data <- readRDS(data_file)
		y <- data$y

		set.seed(seed)

		execution_bench <- system.time({
			res <- sir_laplace_cpp(y = y, N = N,
									mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
									nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
									theta1 = theta1, theta2 = theta2,
									theta_01 = theta_01, theta_02 = theta_02,
									W1 = W1, W2 = W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol,
									verbose = FALSE, print_every = 1000)
		})
		elapsed_time <- execution_bench[["user.self"]]

		result <- data.frame(
			f = f, Tt = Tt, method = method, replica = replica,
			elapsed_time = elapsed_time,
			worker_host = Sys.info()[["nodename"]],
			pid = Sys.getpid(),
			stringsAsFactors = FALSE
		)

		saveRDS(result, file = task_result_file)
	}

	return(result)
}

# ---- Cluster infrastructure (copied verbatim from simulation.R) ----
make_pbs_cluster <- function() {

	nodefile <- Sys.getenv("PBS_NODEFILE")
	if (nodefile == "") stop("PBS_NODEFILE is not set - is this running inside a PBS job?")
	hosts <- unique(readLines(nodefile))

	get_physical_cpu_ids <- function(host) {
		remote_cmd <- "grep -E '^processor|^core id|^physical id' /proc/cpuinfo"
		out <- system2("ssh", args = c(host, shQuote(remote_cmd)), stdout = TRUE)
		proc_ids <- as.integer(gsub(".*:\\s*", "", out[grepl("^processor", out)]))
		phys_ids <- as.integer(gsub(".*:\\s*", "", out[grepl("^physical id", out)]))
		core_ids <- as.integer(gsub(".*:\\s*", "", out[grepl("^core id", out)]))
		df <- data.frame(processor = proc_ids, physical = phys_ids, core = core_ids)
		df <- df[!duplicated(df[, c("physical", "core")]), ]
		df$processor
	}

	cpu_map <- lapply(hosts, get_physical_cpu_ids)
	names(cpu_map) <- hosts

	worker_specs <- do.call(rbind, lapply(hosts, function(h) {
		data.frame(host = h, cpu = cpu_map[[h]], stringsAsFactors = FALSE)
	}))

	if (verbose) printf("PBS cluster: %d nodes, %d workers (physical cores, no HT)",
						 length(hosts), nrow(worker_specs))

	worker_script_dir <- path.expand("~/cobalebeb2027/simulation/.worker_scripts")
	dir.create(worker_script_dir, showWarnings = FALSE, recursive = TRUE)

	# rscript= must be a single executable path with no embedded spaces --
	# parallel::makePSOCKcluster() quotes it as one token when building the
	# remote ssh command, so anything with arguments (e.g. "env VAR=x script"
	# or "bash -c '...'") breaks regardless of how it is escaped on our side.
	# Generating one small wrapper per cpu, with the core number baked into
	# its content instead of passed as an argument, avoids this entirely.
	make_pinned_worker <- function(host, cpu) {
		script_path <- sprintf("%s/worker_%s_%d.sh", worker_script_dir, Sys.getenv("PBS_JOBID"), cpu)
		writeLines(c(
			"#!/bin/bash",
			"source ~/setup_env.sh > /dev/null 2>&1",
			sprintf('exec taskset -c %d Rscript "$@"', cpu)
		), script_path)
		Sys.chmod(script_path, mode = "0755")
		parallel::makePSOCKcluster(host, rscript = script_path)
	}

	workers <- mapply(make_pinned_worker,
					   worker_specs$host, worker_specs$cpu,
					   SIMPLIFY = FALSE)

	# c() on classed list objects (SOCKcluster/cluster) drops the class --
	# there is no c.cluster S3 method -- so the combined object needs its
	# class reassigned explicitly, or checkCluster() rejects it downstream.
	cl <- do.call(c, workers)
	class(cl) <- c("SOCKcluster", "cluster")
	cl
}

# ---- Run the pilot ----
start_time <- proc.time()

current_wd <- getwd()
`%dorng%` <- doRNG::`%dorng%`

cl <- make_pbs_cluster()
n_workers <- length(cl)
printf("make_pbs_cluster() created %d workers (physical cores detected).", n_workers)
if (n_workers != 40) {
	printf("NOTE: worker count != 40 -- ncpus=40 likely counted logical (hyperthreaded) threads, not physical cores. k* for the production packing plan should be revisited accordingly.")
}

# Build the task grid now, sized to the number of physical-core workers
# actually available (capped by N_replicas_max), instead of assuming 40.
N_replicas_pilot <- min(n_workers, N_replicas_max)
task_grid <- data.frame(
	Tt = Tt,
	f = f,
	replica = 1:N_replicas_pilot,
	method = method,
	stringsAsFactors = FALSE
)
task_grid$seed <- task_grid$replica * 10  # pilot-only seed, no collision risk here
grid_size <- nrow(task_grid)
if (verbose) printf("Pilot grid size: %d", grid_size)

doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, c("current_wd", "printf"))

diag_results <- parallel::clusterEvalQ(cl, {
	setwd(current_wd)
	pkgload::load_all("../PoissonLTDM", debug = FALSE)
	list(
		node = Sys.info()[["nodename"]],
		pid = Sys.getpid(),
		exists_search = exists("sample_sir_laplace"),
		exists_ns = exists("sample_sir_laplace", envir = asNamespace("PoissonLTDM")),
		poissonltdm_attached = any(grepl("PoissonLTDM", search())),
		search_path = paste(search(), collapse = " | ")
	)
})

for (d in diag_results) {
	printf("[worker pid=%s @ %s] search=%s ns=%s attached=%s",
		   d$pid, d$node, d$exists_search, d$exists_ns, d$poissonltdm_attached)
}
printf("search() on worker 1: %s", diag_results[[1]]$search_path)

df_results <- foreach::foreach(
	i = 1:grid_size,
	.combine = rbind) %dorng% {
	run_task(task_grid[i, ])
}

parallel::stopCluster(cl)

end_time <- proc.time()
wall_time <- (end_time - start_time)[[3]]

# ---- Diagnostics ----
printf("=== Pilot summary ===")
printf("Wall time (s): %.2f", wall_time)
printf("Tasks completed: %d", nrow(df_results))
printf("Distinct worker PIDs used: %d", length(unique(df_results$pid)))
printf("Distinct hosts used: %d", length(unique(df_results$worker_host)))
printf("Mean task elapsed_time (s): %.2f", mean(df_results$elapsed_time))
printf("Max task elapsed_time (s): %.2f", max(df_results$elapsed_time))
saveRDS(df_results, file = sprintf("%s/pilot_summary.rds", path_results_partial))
