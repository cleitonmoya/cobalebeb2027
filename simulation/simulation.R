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

# Build the Task Grid 
Tt_grid <- c(200, 400, 800, 2000)
functions_grid <- c("block", "piecewise_linear", "heavisine", "piecewise_pol")
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
	path_partial <- "results/partial/"
	file_name <- sprintf("%s/Tt%d_%s_%s_r%03d.rds",
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


run_task <- function(task) {

	Tt <- task$Tt
	f <- task$f
	method <- task$method
	replica <- task$replica
	seed <- task$seed

	file_name  <- task_file_name(task)
	file_star_name <- task_star_file_name(task)

	# Checkpoint: skip the task if it is already done
	if (file.exists(file_name)) {
		result = readRDS(file_name) #
	} else {
		
		# Load the data
		f_teste <- obter_funcao_teste(funcao)
		dados   <- gerar_dados_poisson(f = f_teste, Tt = Tt, seed = seed_data)
		y            <- dados$y
		theta1_true  <- dados$theta1_true
		
		# hiperparametros de priori calibrados a partir de f(t) (ja decidido em conversas anteriores)
		priors <- calibrar_priors(f_teste, Tt)
		
		# --- 3.2 rodar o metodo correspondente ---
		set.seed(seed)
		
		tempo <- system.time({
			resultado <- switch(metodo,
								"mh_cw"          = sample_mh_cw(y, priors, n_iter = 5000, burnin = 1000),
								"mh_montoril"    = sample_mh_montoril(y, priors, n_iter = 5000, burnin = 1000),
								"pg_apf"         = sample_pg_apf(y, priors, n_iter = 5000, burnin = 1000, n_particulas = 200),
								"sir_collapsed"  = sample_sir_collapsed(y, priors, n_iter = 5000, burnin = 1000),
								"sir_laplace"    = sample_sir_laplace(y, priors, n_iter = 5000, burnin = 1000),
								"stan"           = stan_sample(modelo_stan, y, N = 5000, burnin = 1000, priors),
								stop(sprintf("Metodo desconhecido: %s", metodo))
			)
		})
		tempo_total <- tempo[["elapsed"]]
		
		# --- 3.3 calcular metricas ---
		ess_theta1    <- calcular_ess(resultado$theta1_hist)
		ess_W1        <- calcular_ess(resultado$W1_hist)
		ess_W2        <- calcular_ess(resultado$W2_hist)
		rmse_theta1   <- calcular_rmse(resultado$theta1_hist, theta1_true)
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
				theta1_hist = resultado$theta1_hist,
				theta2_hist = resultado$theta2_hist,
				W1_hist     = resultado$W1_hist,
				W2_hist     = resultado$W2_hist
			), file = arquivo_cadeia)
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
