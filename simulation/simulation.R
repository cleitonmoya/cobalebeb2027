# Poisson Local Trend Dynamic Model
# Main simulation script
# Author: Cleiton Moya de Almeida

# Provide more informative traceback
options(error = function() traceback(2)) 
rm(list = ls())     # clear the environment

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# IMPORTANT: load_all() compiles C++ in DEBUG mode by default (-g -O0),
# silently overriding any -O3 in PoissonLTDM/src/Makevars, causing a
# ~4-5x slowdown of every sampler -- critical here since this is the main,
# long-running simulation. Using pkgload::load_all() directly (not
# devtools::load_all()) because devtools::load_all()'s `...` does not
# reliably forward debug= to pkgload in this environment (devtools 2.5.2 /
# pkgload 1.5.3), producing a "must be used" warning even though the
# argument is valid.
pkgload::load_all("../PoissonLTDM", debug = FALSE) # package with the samplers

path_data <- "../data/simulated"
path_results <- "results"
path_results_partial <- sprintf("%s/%s", path_results, "partial")
path_results_star <- sprintf("%s/%s", path_results, "star")

verbose <- TRUE   # print information in the console
parallel <- TRUE  # used only in local mode (cluster = FALSE)
n_cores <- 2      # used only in local mode (cluster = FALSE)


# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))

# ---- Execution environment auto-detection----
# cluster = TRUE runs in multi-node PBS mode (Euler), with 1 worker per
# physical core (no hyperthreading), pinned via taskset.
# cluster = FALSE runs locally via a simple PSOCK cluster (n_cores).
# Auto-detected from the presence of PBS_NODEFILE (this variable only
# exists inside the execution of a PBS job).
cluster <- Sys.getenv("PBS_NODEFILE") != ""
if (cluster) parallel <- TRUE  # cluster mode implies parallel execution

# General simulation parameters
N <- 10000          # number of iterations
burnin <- 1000

# Adaptive Metropolis hyperparameters
# (for montoril)
varsigma2_scal <- 0.02     # initial varsigma2
ac_ref <- 0.44             # acceptance ratio target

# Particle Gibbs (sir_pg) hyperparameter
K <- 50  # Number of particles

# SIR Laplace and SIR Collapsed  hyperparameters
M_is <- 3             # Number of particles - IS for W1 integrated likelihood
M_irls_max <- 20
tol <- 1e-4

# Only for SIR Collapsed
R_prerun <- 3000      # pre-run iterations to calibrate CE proposals (phi1 and phi2)
M_sir <- 3            # Number of particles - SIR of theta1

# Prior hyperparameters
# theta_01 ~ N(mu_01, sigma2_01)
mu_01     <- 0
sigma2_01 <- 100

# theta_02 ~ N(mu_2, sigma2_02)
mu_02     <- 0
sigma2_02 <- 100

# phi1 = W1^(-1) ~ Gamma(nu_01, eta_01)
nu_01  <- 2
eta_01 <- 0.01

# phi2 = W2^(-1) ~ Gamma(nu_02, eta_02)
nu_02  <- 2
eta_02 <- 0.0001

# Initialization
W2 <- 0.01
W1 <- 0.01
theta_01 <- 0
theta_02 <- 0
theta1_tilde_scal <- 0 # sir_laplace and sir_collapsed


#####
# Build the Task Grid 

# Used to build the grid and simulate
Tt_grid <- c(200, 400, 800, 1600)
functions_grid <- c("constant","linear", "quadratic", "sinusoidal")
methods_grid <- c("stan")
N_replicas <- 10

# Only for reference (used to compute the seed)
Tt_grid_ref <- c(200, 400, 800, 1600)
functions_grid_ref <- c("constant","linear", "quadratic", "sinusoidal")
methods_grid_ref <- c("montoril", "pg_apf", "sir_laplace", "sir_collapsed", "stan")

task_grid <- expand.grid(
	Tt = Tt_grid,
	f = functions_grid,
	replica = 1:N_replicas,
	method = methods_grid,
	stringsAsFactors = FALSE
)

grid_size <- nrow(task_grid)
if (verbose) printf("Total of tasks: %d", grid_size)

# Task stars: Store all simulated data for the (stars) selected functions and Tt 
f_star <- "quadratic"
Tt_star <- c(200, 1600)

task_grid$star <- with(task_grid,
	replica == 1 & Tt %in% Tt_star & f == f_star
)

if (verbose) printf("Total of star tasks: %d", sum(task_grid$star))


# Insert seed column
task_grid$seed <- match(task_grid$method, methods_grid_ref) * 1e5 +
					match(task_grid$f, functions_grid_ref) * 1e4 +
					match(task_grid$Tt, Tt_grid_ref) * 1e3 +
					task_grid$replica * 10


task_result_filename <- function(task) {
	file_name <- sprintf("%s/%s_%s_%s_%s.rds",
						 path_results_partial, task$method, task$f, task$Tt,
						 task$replica)
	return(file_name)
}

task_star_result_filename <- function(task) {
	file_name <- sprintf("%s/%s_%s_%s_%s.rds",
						 path_results_star, task$method, task$f, task$Tt,
						 task$replica)
	return(file_name)
}

data_filename <- function(task) {
	file_name <- sprintf("%s/%s_%s_%s.rds",
						 path_data, task$f, task$Tt, task$replica)
	return(file_name)
}


run_task <- function(task) {

	Tt <- task$Tt
	replica <- task$replica
	f <- task$f
	method <- task$method
	seed <- task$seed

	# Initial values for theta1, theta2 (need Tt) 
	theta1 <- numeric(Tt)
	theta2 <- numeric(Tt)
	
	
	task_result_file  <- task_result_filename(task)
	task_star_result_file <- task_star_result_filename(task)

	
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
				
				"montoril" = sample_amh_montoril(y, N, burnin, varsigma2_scal, ac_ref,
												    mu_01, sigma2_01, mu_02, sigma2_02,
												    nu_01, eta_01, nu_02, eta_02,
												    W1, W2, theta_01, theta_02, theta1, theta2),
								
				"pg_apf" = sample_pg_apf(y, N, burnin, K,
										 mu_01, sigma2_01, mu_02, sigma2_02,
										 nu_01, eta_01, nu_02, eta_02,
										 W1, W2, theta_01, theta_02, theta1, theta2),
				
				"sir_laplace" = sample_sir_laplace(y, N, burnin, M_is, M_irls_max, tol,
												   mu_01, sigma2_01, mu_02, sigma2_02,
												   nu_01, eta_01, nu_02, eta_02,
												   W1, W2, theta_01, theta_02, theta1, theta2,
												   theta1_tilde_scal),
								
				"sir_collapsed" = sample_sir_collapsed(y, N, burnin, R_prerun, 
													   M_is, M_sir, M_irls_max, tol,
													   mu_01, sigma2_01, mu_02, sigma2_02,
													   nu_01, eta_01, nu_02, eta_02,
													   W1, W2, theta_01, theta_02, theta1, theta2,
													   theta1_tilde_scal),
								
				"stan" = sample_stan(model, y, N, burnin, seed,
									 mu_01, sigma2_01, mu_02, sigma2_02,
									 nu_01, eta_01, nu_02, eta_02,
									 W1, W2, theta_01, theta_02, theta1, theta2),
				
				stop(sprintf("Unknow method: %s", method))
			)
		})
		elapsed_time <- execution_bench[["user.self"]]
		if (verbose) printf("All tasks ran")
		
		
		#
		# Extract  the results
		#
		
		# Common
		theta_01_hist = res$theta_01_hist
		theta_02_hist = res$theta_02_hist
		W1_hist = res$W1_hist
		W2_hist = res$W2_hist
		theta1_hist = res$theta1_hist
		theta2_hist = res$theta2_hist
		
		if (method == "stan") {
			# overwrite the previous elapsed time computed by system.time()
			elapsed_time <- res$elapsed_time
		}
		
		#
		# Compute the metrics of the replica
		#
		
		# Samples without burn-in
		theta_01_samples <- theta_01_hist[-(1:burnin)]
		theta_02_samples <- theta_02_hist[-(1:burnin)]
		W1_samples <- W1_hist[-(1:burnin)]
		W2_samples <- W2_hist[-(1:burnin)]
		theta1_samples <- theta1_hist[-(1:burnin), ]
		theta2_samples <- theta2_hist[-(1:burnin), ]
		
		
		# Effective Sample Size (bulk/tail; R_hat/Geweke not assessed per
		# replica -- only during the pilot/calibration stage that sets N and
		# burnin, using metrics_convergence() on a multi-chain run)
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
		ess_sec_theta_01 <- ess_theta_01/elapsed_time
		ess_sec_theta_02 <- ess_theta_02/elapsed_time
		ess_sec_W1 <- ess_W1/elapsed_time
		ess_sec_W2 <- ess_W2/elapsed_time
		ess_sec_theta1_mean <- ess_theta1_mean/elapsed_time
		ess_sec_theta2_mean <- ess_theta2_mean/elapsed_time
		
		# Fit metrics
		
		
		W1_mean <- mean(W1_samples)
		W1_median <- median(W1_samples)
		W1_var <- var(W1_samples)
		
		W2_mean <- mean(W2_samples)
		W2_median <- median(W2_samples)
		W2_var <- var(W2_samples)
		
		theta1_mean <- colMeans(theta1_samples)
		theta2_mean <- colMeans(theta2_samples)
		lambda_mean <- exp(theta1_mean)

		log_lik <- metrics_loglik(y, lambda_mean)
		rmse_theta1 <- metrics_rmse(theta1_mean, theta1_true)
		mae_theta1 <- metrics_mae(theta1_mean, theta1_true)
		
		theta1_ci <- metrics_theta_ci(theta1_samples, 0.05)
		theta1_ci_lower <- theta1_ci$ci_lower
		theta1_ci_upper <- theta1_ci$ci_upper
		
		theta2_ci <- metrics_theta_ci(theta2_samples, 0.05)
		theta2_ci_lower <- theta2_ci$ci_lower
		theta2_ci_upper <- theta2_ci$ci_upper
		
		
		# TO-DO: Coverage of the Empirical CI
		
		#
		# Save the results
		#
		task_result <- data.frame(
			f = f,
			Tt = Tt,
			method = method,
			replica = replica,
			elapsed_time = elapsed_time,
			
			ess_theta_01 = ess_theta_01,
			ess_theta_02 = ess_theta_02,
			ess_W1 = ess_W1,
			ess_W2 = ess_W2,
			ess_theta1_mean = ess_theta1_mean,
			ess_theta1_min = ess_theta1_min,
			ess_theta2_mean = ess_theta2_mean,
			ess_theta2_min = ess_theta2_min,
			
			ess_sec_theta_01 = ess_sec_theta_01,
			ess_sec_theta_02 = ess_sec_theta_02,
			ess_sec_W1 = ess_sec_W1,
			ess_sec_W2 = ess_sec_W2,
			ess_sec_theta1_mean = ess_sec_theta1_mean,
			ess_sec_theta2_mean = ess_sec_theta2_mean,
			
			W1_mean = W1_mean,
			W1_median = W1_median,
			W1_var = W1_var,
			W2_mean = W2_mean,  
			W2_median = W2_median,
			W2_var = W2_var,
			
			log_lik = log_lik,
			rmse_theta1 = rmse_theta1,
			mae_theta1 = mae_theta1,
			
			theta1_ci_lower = theta1_ci_lower,
			theta1_ci_upper = theta1_ci_upper,
			
			theta2_ci_lower = theta2_ci_lower,
			theta2_ci_upper = theta2_ci_upper,
			
			stringsAsFactors = FALSE
		)
		
		saveRDS(task_result, file = task_result_file)
		
		# Level two - Star tasks
		if (task$star) {
			
			result_star <- list(
				W1_hist = W1_hist,
				W2_hist = W2_hist,
				theta_01_hist = theta_01_hist,
				theta_02_hist = theta_02_hist,
				theta1_hist = theta1_hist,
				theta2_hist = theta2_hist)

			# Only for Montoril and Stan
			if (method %in% c("montoril", "stan")) {
				result_star$ac_hist = res$ac_hist
			}
			
			if (method == "pg_apf") {
				result_star$ess_smc = res$ess_smc
			}
			
			if (method %in% c("sir_laplace", "sir_collapsed")) {
				result_star$ess_is = res$ess_is
				result_star$itr_irls = res$itr_irls
			}
			
			if (method == "sir_collapsed") {
				result_star$accepted1_hist <- res$accepted1_hist
				result_star$accepted2_hist <- res$accepted2_hist
			}
			
			saveRDS(result_star, file = task_star_result_file)
		}
		
		result <- task_result
	}
	
	return(result)
}


#####
# Cluster infrastructure

# Local mode: simple PSOCK cluster, n_cores on this same machine
make_local_cluster <- function(n_cores) {
	parallel::makeCluster(n_cores)
}

# Cluster mode (PBS/Euler): 1 worker per physical core (no hyperthreading),
# pinned via taskset, spread across all nodes allocated to the job.
# Each remote worker needs to redo the environment setup (modules, SSL/
# toolchain variables), since the SSH session does not inherit the PBS
# job's shell.
make_pbs_cluster <- function() {
	
	nodefile <- Sys.getenv("PBS_NODEFILE")
	if (nodefile == "") stop("PBS_NODEFILE is not set - is this running inside a PBS job?")
	hosts <- unique(readLines(nodefile))

	get_physical_cpu_ids <- function(host) {
		remote_cmd <- "grep -E '^processor|^core id' /proc/cpuinfo"
		out <- system2("ssh", args = c(host, shQuote(remote_cmd)), stdout = TRUE)
		proc_ids <- as.integer(gsub(".*:\\s*", "", out[grepl("^processor", out)]))
		core_ids <- as.integer(gsub(".*:\\s*", "", out[grepl("^core id", out)]))
		df <- data.frame(processor = proc_ids, core = core_ids)
		df <- df[!duplicated(df$core), ]
		df$processor
	}

	cpu_map <- lapply(hosts, get_physical_cpu_ids)
	names(cpu_map) <- hosts

	worker_specs <- do.call(rbind, lapply(hosts, function(h) {
		data.frame(host = h, cpu = cpu_map[[h]], stringsAsFactors = FALSE)
	}))

	if (verbose) printf("PBS cluster: %d nodes, %d workers (physical cores, no HT)",
						 length(hosts), nrow(worker_specs))

	make_pinned_worker <- function(host, cpu) {
		parallel::makePSOCKcluster(host,
			rscript = sprintf("bash -c 'source ~/setup_env.sh > /dev/null 2>&1 && taskset -c %d Rscript'", cpu))
	}

	workers <- mapply(make_pinned_worker,
					   worker_specs$host, worker_specs$cpu,
					   SIMPLIFY = FALSE)
	do.call(c, workers)
}


#####
# Main exection

start_time = proc.time() # execution time
if (!parallel) {
	
	if (verbose) printf("Starting sequential execution")
	
	result_list <- vector("list", grid_size)
	for (i in 1:grid_size) {
		
		if (verbose) printf("Running task %d/%d: f=%s, Tt=%d, method=%s,replica=%d",
							i,
							grid_size,
							task_grid$f[i],
							task_grid$Tt[i],
							task_grid$method[i],
							task_grid$replica[i])
		
		result_list[[i]] <- run_task(task_grid[i, ])
	}
	
	# Stack the results (1-row data.frame) in a single data.frame
	df_results <- do.call(rbind, result_list)
	
} else {
	
	# Parallel execution
	if (verbose) printf("Starting parallel execution (cluster = %s)", cluster)
	current_wd <- getwd()
	
	`%dorng%` <- doRNG::`%dorng%`
	
	cl <- if (cluster) make_pbs_cluster() else make_local_cluster(n_cores)
	doParallel::registerDoParallel(cl)
	parallel::clusterExport(cl, "current_wd")
	
	parallel::clusterEvalQ(cl, {
		setwd(current_wd)
		# See comment near the top of this file: pkgload::load_all() with
		# debug=FALSE avoids the ~4-5x DEBUG-mode slowdown in every worker.
		pkgload::load_all("../PoissonLTDM", debug = FALSE)
	})
	
	df_results <- foreach::foreach(
		i = 1:grid_size,
		.combine = rbind) %dorng% {
		run_task(task_grid[i, ])
	}

	parallel::stopCluster(cl)
}
end_time <- proc.time()
execution_time <- (end_time - start_time)[[3]]

if (verbose) printf("Simulation complete in %.1f min", execution_time/60)

