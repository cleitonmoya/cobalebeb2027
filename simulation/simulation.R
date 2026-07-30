# Poisson Local Trend Dynamic Model
# Main simulation script
# Author: Cleiton Moya de Almeida

# library(rstan)
devtools::load_all("PoissonLTDM") # package with the samplers

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
path_data <- "../data/simulated"
path_results <- "results"
path_results_partial <- sprintf("%s/%s", path_results, "partial")
path_results_star <- sprintf("%s/%s", path_results, "star")

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))


# General simulation parameters
N <- 10000          # number of iterations
burnin <- 1000

# Adaptive Metropolis hyperparameters
# (for amh_cw and amg_montoril)
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

# General hyperparameters

#####
# Build the Task Grid 

Tt_grid <- c(200)
functions_grid <- c("quadratic")
methods_grid <- c("montoril", "pg_apf", "sir_collapsed", "sir_laplace", "stan")

# Tt_grid <- c(200, 400, 800, 2000)
# functions_grid <- c("constant", "linear", "quadratic", "sinusoidal")
# methods_grid <- c("montoril", "pg_apf", "sir_collapsed", "sir_laplace", "stan")
N_replicas <- 1

task_grid <- expand.grid(
	Tt = Tt_grid,
	f = functions_grid,
	replica = 1:N_replicas,
	method = methods_grid,
	stringsAsFactors = FALSE
)

grid_size <- nrow(task_grid)
printf("Total of tasks: %d\n", grid_size)

# Task stars: Store all simulated data for the (stars) selected functions and Tt 
f_star <- "quadratic"
Tt_star <- c(200, 2000)

task_grid$star <- with(task_grid,
	replica == 1 & Tt %in% Tt_star & f == f_star
)

printf("Total of star tasks: %d\n", sum(task_grid$star))


# Insert seed column
task_grid$seed <- 1:grid_size


task_result_filename <- function(task) {
	file_name <- sprintf("%s/%d_%s_%s_r%03d.rds",
						 path_results_partial, task$Tt, task$f,
						 task$method, task$replica)
	return(file_name)
}

task_star_result_filename <- function(task) {
	file_name <- sprintf("%s/Tt%d_%s_%s_r%03d.rds",
						 path_results_star, task$Tt, task$f,
						 task$method, task$replica)
	return(file_name)
}

data_filename <- function(task) {
	file_name <- sprintf("%s/%s_%d_%d.rds",
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
		result <- readRDS(task_result_file) #
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
				printf("Model loaded")
			} else {
				printf("Building the Stan model")
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
		elapsed_time <- execution_bench[["elapsed"]]
		
		
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
		
		# Effective Sample Size
		ess_theta_01 <- compute_ess(theta_01_hist[-(1:burnin)])
		ess_theta_02 <- compute_ess(theta_02_hist[-(1:burnin)])
		ess_W1 <- compute_ess(W1_hist[-(1:burnin)])
		ess_W2 <- compute_ess(W2_hist[-(1:burnin)])
		ess_theta1 <- compute_ess(theta1_hist[-(1:burnin), ])
		ess_theta2 <- compute_ess(theta2_hist[-(1:burnin), ])
		
		# Effective Sample Size per second
		ess_sec_theta_01 <- ess_theta_01/elapsed_time
		ess_sec_theta_02 <- ess_theta_02/elapsed_time
		ess_sec_W1 <- ess_W1/elapsed_time
		ess_sec_W2 <- ess_W2/elapsed_time
		ess_sec_theta1_mean <- mean(ess_theta1/elapsed_time)
		ess_sec_theta2_mean <- mean(ess_theta2/elapsed_time)
		
		# Geweke diagnostic
		z_W1 <- compute_geweke(W1_hist[-(1:burnin)])
		z_W2 <- compute_geweke(W2_hist[-(1:burnin)])
		z_theta_01 <- compute_geweke(theta_01_hist[-(1:burnin)])
		z_theta_02 <- compute_geweke(theta_02_hist[-(1:burnin)])
		
		z_theta1 <- compute_geweke(theta1_hist[-(1:burnin), ])
		z_theta2 <- compute_geweke(theta2_hist[-(1:burnin), ])
		z_theta1_out <- sum((z_theta1 < -1.96) | (z_theta1 > 1.96))/Tt
		z_theta2_out <- sum((z_theta2 < -1.96) | (z_theta2 > 1.96))/Tt
		
		# Fit metrics
		W1_mean <- mean(W1_hist[-(1:burnin)])
		W1_median <- median(W1_hist[-(1:burnin)])
		W2_mean <- mean(W2_hist[-(1:burnin)])
		W2_median <- median(W2_hist[-(1:burnin)])
		
		theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
		theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
		lambda_mean <- exp(theta1_mean)

		log_lik <- compute_loglik(y, lambda_mean)
		rmse_theta1 <- compute_rmse(theta1_mean, theta1_true)
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
			ess_theta2_mean = ess_theta2_mean,
			
			ess_sec_theta_01 = ess_sec_theta_01,
			ess_sec_theta_02 = ess_sec_theta_02,
			ess_sec_W1 = ess_sec_W1,
			ess_sec_W2 = ess_sec_W2,
			ess_sec_theta1_mean = ess_sec_theta1_mean,
			ess_sec_theta2_mean = ess_sec_theta2_mean,
			
			z_W1 = z_W1,
			z_W2 = z_W2,
			z_theta_01 = z_theta_01,
			z_theta_02 = z_theta_02,
			z_theta1_out = z_theta1_out,
			z_theta2_out = z_theta2_out,
			
			W1_mean = W1_mean,
			W1_median = W1_median,
			W2_mean = W2_mean,  
			W2_median = W2_median,
			
			log_lik = log_lik,
			rmse_theta1 = rmse_theta1,
			
			stringsAsFactors = FALSE
		)
		
		# Only for Montoril and Stan
		if (method %in% c("montoril", "stan")) {
			task_result$ac_hist = res$ac_hist
		}

		if (method == "pg_apf") {
			task_result$ess_smc = res$ess_smc
		}
		
		if (method %in% c("sir_laplace", "sir_collapsed")) {
			task_result$ess_is = res$ess_is
			task_result$itr_irls = res$itr_irls
		}
		
		if (method == "sir_collapsed") {
			task_result$accepted1_hist <- res$accepted1_hist
			task_result$accepted2_hist <- res$accepted2_hist
		}
		
		saveRDS(task_result, file = task_result_file)
		
		# Level two - Star tasks
		if (task$star) {
			saveRDS(list(
				W1_hist = W1_hist,
				W2_hist = W2_hist,
				theta_01_hist = theta_01_hist,
				theta_02_hist = theta_02_hist,
				theta1_hist = theta1_hist,
				theta2_hist = theta2_hist
			), file = task_star_result_file)
		}
		
		result <- task_result
	}
	
	return(result)
}


#####
# Sequential execution
#

result_list <- vector("list", grid_size)

for (i in 1:grid_size) {
	
	printf("Running task %d/%d: f=%s, Tt=%d, method=%s,replica=%d",
		   i,
		   grid_size,
		   task_grid$Tt[i], 
		   task_grid$f[i],
		   task_grid$method[i],
		   task_grid$replica[i])

	result_list[[i]] <- run_task(task_grid[i, ])
}

# Stack the results (1-row data.frame) in a single data.frame
df_results <- do.call(rbind, result_list)


#
# Aggregate the results
#
file_rds <- sprintf("%s/simulation_results.rds", path_results) 
saveRDS(df_results, file = file_rds)

file_csv <- sprintf("%s/simulation_results.csv", path_results) 
write.csv(df_results, file = file_csv, row.names = FALSE)

printf("Simulation complete!")
