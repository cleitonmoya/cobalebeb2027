# Poisson Local Trend Dynamic Model
# Main simulation script
# Author: Cleiton Moya de Almeida

library(Matrix)
library(coda)
library(invgamma)
# library(rstan)

devtools::load_all("PoissonLTDM")      # package with the samplers

source("simulacoes/funcoes_teste.R")   # f(t): near-constant, linear, inflexao, oscilatoria
source("simulacoes/gerar_dados.R")     # gerar_dados_poisson(f, Tt, seed)
source("simulacoes/metricas.R")        # calcular_ess, calcular_rmse, calcular_cobertura, geweke_agregado

# Build the Task Grid 
Tt_grid <- c(200, 400, 800, 2000)
functions_grid <- c("piece_trend", "linear", "inflexao", "sine")
methods_grid <- c("mh_cw", "mh_montoril", "pg_apf", "sir_collapsed", "sir_laplace", "stan")
N_replicas <- 50

task_grid <- expand.grid(
	Tt = Tt_grid,
	f = functions_grid,
	method = methods_grid,
	replica = 1:N_replicas,
	stringsAsFactors = FALSE
)

grid_size <- nrow(task_grid)
task_grid$seed <- 1:grid_size

# Task stars: Store all simulated data for the (stars) selected functions and Tt 
f_star <- "piece_trend"
Tt_star <- c(200, 2000)

task_grid$star_task <- with(task_grid,
	replica == 1 & Tt %in% Tt_star & f == f_star
)


printf("Total of tasks: %d\n", grid_size)
printf("Total of star tasks: %d\n", sum(task_grid$star_tesk))


# pasta onde cada tarefa salva seu proprio resultado individual (checkpointing)
dir_parcial <- "results/partial"
dir.create(dir_parcial, showWarnings = FALSE, recursive = TRUE)

# pasta onde as tarefas "representativas" (ver salvar_cadeia acima) salvam a
# cadeia bruta completa, separada dos resumos -- mantida a parte por ser
# muito maior em volume e nao fazer parte do checkpointing leve do grid inteiro
dir_cadeias <- "results/star"
dir.create(dir_cadeias, showWarnings = FALSE, recursive = TRUE)

# nome de arquivo unico e deterministico por tarefa -- usado tanto para salvar
# quanto para checar se a tarefa ja foi rodada (permite retomar apos interrupcao)
nome_arquivo_tarefa <- function(tarefa_row) {
	sprintf("%s/Tt%d_%s_%s_r%03d.rds",
	        dir_parcial, tarefa_row$Tt, tarefa_row$funcao,
	        tarefa_row$metodo, tarefa_row$replica)
}

# nome de arquivo para a cadeia bruta (mesmo padrao, pasta diferente)
nome_arquivo_cadeia <- function(tarefa_row) {
	sprintf("%s/Tt%d_%s_%s_r%03d_cadeia.rds",
	        dir_cadeias, tarefa_row$Tt, tarefa_row$funcao,
	        tarefa_row$metodo, tarefa_row$replica)
}


# 3. FUNCAO DE UMA TAREFA (unidade de paralelizacao) ###########################

# Esta funcao encapsula TUDO que uma tarefa precisa: gerar os dados (se ja nao
# tiver sido gerado nesta combinacao Tt x funcao x replica), rodar o metodo,
# calcular as metricas, e retornar uma linha de resultado.
#
# Por ser autocontida (sem depender de estado externo alem dos argumentos),
# esta funcao pode ser chamada tanto num for() sequencial quanto num
# foreach(...) %dopar% sem nenhuma mudanca de logica interna.

rodar_tarefa <- function(tarefa_row) {

	Tt <- tarefa_row$Tt
	funcao  <- tarefa_row$funcao
	metodo  <- tarefa_row$metodo
	replica <- tarefa_row$replica
	seed_data <- tarefa_row$seed_data
	seed_mcmc  <- tarefa_row$seed_mcmc

	arquivo_saida  <- nome_arquivo_tarefa(tarefa_row)
	arquivo_cadeia <- nome_arquivo_cadeia(tarefa_row)
	precisa_cadeia <- isTRUE(tarefa_row$salvar_cadeia)

	# --- 3.0 checkpoint: pula se esta tarefa ja foi rodada e salva ---
	# (considera tambem a cadeia bruta, se esta tarefa exige salva-la)
	resumo_pronto <- file.exists(arquivo_saida)
	cadeia_pronta <- !precisa_cadeia || file.exists(arquivo_cadeia)
	if (resumo_pronto && cadeia_pronta) {
		return(readRDS(arquivo_saida))
	}

	# --- 3.1 carregar dados (determinisitico dado seed_data) ---
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

	# --- 3.6 nivel 2: salva a cadeia bruta completa, se esta tarefa for representativa ---
	if (precisa_cadeia) {
		saveRDS(list(
			theta1_hist = resultado$theta1_hist,
			theta2_hist = resultado$theta2_hist,
			W1_hist     = resultado$W1_hist,
			W2_hist     = resultado$W2_hist
		), file = arquivo_cadeia)
	}

	linha_resultado
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
