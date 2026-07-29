# Poisson Local Trend Dynamic Model
# Main simulation script
# Author: Cleiton Moya de Almeida

library(Matrix)
library(coda)
library(invgamma)
# library(rstan)
devtools::load_all("PoissonLTDM") # package with the samplers

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

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

# Build the Task Grid 
Tt_grid <- c(200, 400, 800, 2000)
functions_grid <- c("constant", "linear", "quadratic", "sinusoidal")
methods_grid <- c("mh_cw", "mh_montoril", "pg_apf", "sir_collapsed", "sir_laplace", "stan")
N_replicas <- 50

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
f_star <- "piece_trend"
Tt_star <- c(200, 2000)

task_grid$star <- with(task_grid,
	replica == 1 & Tt %in% Tt_star & f == f_star
)

printf("Total of star tasks: %d\n", sum(task_grid$star))


# Insert seed column
task_grid$seed <- 1:grid_size


task_file_name <- function(task) {
	path_partial <- "results/partial"
	file_name <- sprintf("%s/%d_%s_%s_r%03d.rds",
						 path_partial, task$Tt, task$f,
						 task$method, task$replica)
	return(file_name)
}

task_star_file_name <- function(task) {
	path_star <- "results/star"
	file_name <- sprintf("%s/Tt%d_%s_%s_r%03d.rds",
						 path_star, task$Tt, task$f,
						 task$method, task$replica)
	return(file_name)
}

data_file_name <- function(task) {
	path_data <- "../data/simulated"
	file_name <- sprintf("%s/%s_%d_%d.rds",
						 path_data, task$f, task$Tt, task$replica)
	return(file_name)
}


run_task <- function(task) {

	Tt <- task$Tt
	f <- task$f
	method <- task$method
	replica <- task$replica
	seed <- task$seed

	# Initial values for theta1, theta2 (need Tt) 
	theta1 <- numeric(Tt)
	theta2 <- numeric(Tt)
	
	
	file_name  <- task_file_name(task)
	file_star_name <- task_star_file_name(task)

	
	# Checkpoint: skip the task if it is already done
	if (file.exists(file_name)) {
		result = readRDS(file_name) #
	} else {
		
		# Load the data
		data <- readRDS(data_file_name(task))
		y <- data$y
		
		# Run the task
		set.seed(seed)
		
		if (method == 'stan') {
			# Load or compile the model
			if (file.exists("../cache/poisson_ltdm.rds")) {
				model <- readRDS("../cache/poisson_ltdm.rds")
				printf("Model loaded")
			} else {
				printf("Building the model")
				file <- "../PoissonLTDM/inst/stan/poisson_ltdm.stan"
				model <- rstan::stan_model(file = file, model_name = "PoissonLTDM")
				saveRDS(model, file = "../cache/poisson_ltdm.rds")
			}
		}
		
		execution_bench <- system.time({
			res <- switch(method,
						  
				"amh_cw" = sample_amh_cw(y, N, burnin, varsigma2_scal, ac_ref,
										 mu_01, sigma2_01, mu_02, sigma2_02,
										 nu_01, eta_01, nu_02, eta_02,
										 W1, W2, theta_01, theta_02, theta1, theta2),
								
				"amh_montoril" = sample_amh_montoril(y, N, burnin, varsigma2_scal, ac_ref,
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
								
				"stan" = sample_stan(model, y, N, burnin,
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
		
		# Only for amh_cw, amd_montoril and stan
		if (method %in% c("amh_cw", "amg_montoril", "stan")) {
			ac_hist = res$ac_hist
		}
		
		if (method == "pg_apf") {
			ess_smc <- res$ess_smc
		}
		
		if (method %in% c("sir_laplace", "sir_collapsed")) {
			ess_is        <- res$ess_is
			itr_irls      <- res$itr_irls
		}
		
		if (method == "sir_collapsed") {
			accepted1_hist <- res$accepted1_hist
			accepted2_hist <- res$accepted2_hist
		}
		
		if (method == "stan") {
			stan_elapsed_time <- res$elapsed_time
		}
		
		#
		# Compute the metrics
		#
		
		# Effective Sample Size
		ess_theta_01 <- compute_ess(theta_01_hist[-(1:burnin)])
		ess_theta_01 <- compute_ess(theta_01_hist[-(1:burnin)])
		
		ess_W1 <- compute_ess(w1_hist[-(1:burnin)])
		ess_W2 <- compute_ess(w2_hist[-(1:burnin)])
		
		ess_theta1 <- compute_ess(theta1_hist[-(1:burnin), ])
		ess_theta2 <- compute_ess(theta1_hist[-(1:burnin), ])
		
		
		# Effective Sample Size per second
		
		
		# Geweke diagnostic
		
		
		
		
		# 
		rmse_theta1   <- calcular_rmse(restheta1_hist, theta1_true)
		cobertura     <- calcular_cobertura(resultado$theta1_hist, theta1_true)
		geweke_z      <- geweke_agregado(resultado)
		
		
		
		# --- 3.4 montar linha de resultado (data.frame de 1 linha) ---
		linha_resultado <- data.frame(
			Tt             = Tt,
			funcao        = funcao,
			metodo        = metodo,
			replica       = replica,
			tempo_s       = tempo_total,
			ess_theta1    = ess_theta1,
			ess_W1        = ess_W1,
			ess_W2        = ess_W2,
			ess_theta1_s  = ess_theta1 / tempo_total,
			rmse_theta1   = rmse_theta1,
			cobertura     = cobertura,
			geweke_z      = geweke_z,
			stringsAsFactors = FALSE
		)
		
		# --- 3.5 checkpoint: salva o resultado individual (resumo) em disco ---
		saveRDS(linha_resultado, file = arquivo_saida)
		
		# Level two for star tasks
		if (task$star) {
			saveRDS(list(
				theta1_hist = theta1_hist,
				theta2_hist = theta2_hist,
				W1_hist     = W1_hist,
				W2_hist     = W2_hist
			), file = file_name)
		}
		
		result=linha_resultado
	}
	return(result)
}


# 4. EXECUCAO (sequencial por enquanto) ########################################

# Nesta etapa: for() simples, sequencial, so para validar que a estrutura
# funciona ponta a ponta antes de paralelizar.
#
# Cada chamada a rodar_tarefa() ja salva seu proprio resultado em disco
# (bloco 3.5) e pula tarefas ja concluidas (bloco 3.0) -- entao interromper
# e reiniciar este for() a qualquer momento e seguro, sem perder progresso
# nem reprocessar tarefas ja feitas.
#
# Quando for paralelizar (ver bloco 4b comentado abaixo), a unica mudanca
# necessaria e trocar este for() por foreach(...) %dorng% { rodar_tarefa(...) },
# sem tocar em rodar_tarefa() nem no restante do script -- o checkpointing
# por arquivo individual ja e seguro em paralelo, pois cada worker escreve
# em seu proprio arquivo, nunca compartilhado com outro processo.

resultados_lista <- vector("list", nrow(task_grid))

for (i in seq_len(nrow(task_grid))) {
	cat(sprintf("Tarefa %d/%d: Tt=%d, funcao=%s, metodo=%s, replica=%d\n",
	            i, nrow(task_grid),
	            task_grid$Tt[i], task_grid$funcao[i],
	            task_grid$metodo[i], task_grid$replica[i]))

	resultados_lista[[i]] <- rodar_tarefa(task_grid[i, ])
}

resultados <- do.call(rbind, resultados_lista)



# 5. COMBINAR RESULTADOS PARCIAIS E SALVAR CONSOLIDADO #########################

# Le todos os arquivos individuais salvos em resultados/parcial/ e combina
# num unico data.frame. Roda independente de como as tarefas foram executadas
# (sequencial, interrompida e retomada, ou paralela) -- a fonte da verdade
# sao os arquivos em disco, nao o objeto `resultados` em memoria.

arquivos_parciais <- list.files(dir_parcial, pattern = "\\.rds$", full.names = TRUE)
cat(sprintf("Arquivos parciais encontrados: %d de %d tarefas esperadas\n",
            length(arquivos_parciais), nrow(task_grid)))

resultados <- do.call(rbind, lapply(arquivos_parciais, readRDS))

dir.create("resultados", showWarnings = FALSE)
saveRDS(resultados, file = "resultados/resultados_grid_principal.rds")
write.csv(resultados, file = "resultados/resultados_grid_principal.csv", row.names = FALSE)

cat("Simulacao concluida.\n")
