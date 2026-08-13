# Resubmete, dentro do MESMO registro (sem arquivar/recriar nada), os jobs
# que o batchtools ainda considera pendentes -- ex.: chunks travados num
# nó com problema, deletados via qdel mas nunca reportaram conclusão pro
# batchtools.
#
# Rode de dentro de simulation/ no cluster.

library(batchtools)
CALIBRATION_PHASE_DEFS_ONLY <- TRUE
reg <- loadRegistry("registry_calibration", writeable = TRUE)

# CRÍTICO: loadRegistry() NÃO restaura reg$cluster.functions a partir do
# disco -- é um objeto de sessão (closure), não um dado serializável
# normal. No fluxo usual do calibration_phase.R, isso é re-estabelecido
# pelo próprio bloco de despacho (reg$cluster.functions <-
# makeClusterFunctionsTORQUE(...)) -- mas esse bloco só roda quando
# CALIBRATION_PHASE_DEFS_ONLY está DESARMADO, e aqui ele precisa estar
# ARMADO (pra loadRegistry() não re-disparar o despacho completo e tentar
# arquivar/recriar o registro). Sem esta linha, reg$cluster.functions cai
# no padrão do batchtools ("Interactive"), que executa o job ALI MESMO na
# sessão atual (login node) em vez de submeter pro PBS -- confirmado
# acontecendo na prática antes desta correção (evidência: a mensagem
# "using cluster functions 'Interactive'" e o erro "object 'pend_tab' not
# found" logo em seguida, causado pelo próprio job rodando inline e
# disparando o mesmo re-source que apaga .GlobalEnv).
reg$cluster.functions <- makeClusterFunctionsTORQUE("calibration_pbs.tmpl")
cat("cluster.functions confirmado:", reg$cluster.functions$name, "\n")
if (reg$cluster.functions$name != "TORQUE") {
    stop("cluster.functions não ficou como TORQUE -- NÃO prosseguir. Confira calibration_pbs.tmpl.")
}

# Todo o resto do estado (pend_tab, ids_cat, etc.) fica local a esta
# função -- mesma razão documentada em check_progress_calibration.R e
# check_progress.R: se algo aqui dentro disparar um re-source de
# calibration_phase.R, qualquer variável vivendo em .GlobalEnv seria
# apagada; variáveis locais de uma função sobrevivem a isso.
resubmit_pending <- function(reg) {
    # findNotDone(): submitted, never got a done/error signal (e.g. the
    # 'Expired' chunk -- died without reporting back).
    # findErrors(): submitted, DID run, but errored for real (e.g. the
    # chunk that accidentally ran inline before the cluster.functions fix
    # above -- it has a genuine recorded error, not just a missing
    # signal, so findNotDone() alone would silently skip it).
    # Resubmitting covers both: batchtools clears the previous
    # done/error/started markers for a job the moment it's resubmitted.
    to_retry <- unique(rbind(findNotDone(reg = reg), findErrors(reg = reg)))
    cat("Job(s) a resubmeter:", nrow(to_retry), "\n")

    if (nrow(to_retry) == 0) {
        cat("Nada pendente/com erro -- nada a resubmeter.\n")
        return(invisible(NULL))
    }

    tab <- getJobTable(reg = reg)
    pend_tab <- tab[tab$job.id %in% to_retry$job.id, ]
    chunk_ids <- vapply(pend_tab$job.pars, function(p) p$chunk_id, character(1))
    pend_tab$category <- job_grid$category[match(chunk_ids, job_grid$chunk_id)]

    print(pend_tab[, c("job.id", "job.hash")])
    cat("Categorias:", paste(pend_tab$category, collapse = ", "), "\n\n")

    for (cat_name in unique(pend_tab$category)) {
        ids_cat <- pend_tab$job.id[pend_tab$category == cat_name]
        cat(sprintf("Resubmetendo %d job(s) da categoria '%s' (walltime=%gh)...\n",
                    length(ids_cat), cat_name, walltime_hours_for(cat_name)))
        submitJobs(
            ids = data.frame(job.id = ids_cat),
            resources = list(ncpus = CHUNK_SIZE_CALIB,
                              walltime_hours = walltime_hours_for(cat_name),
                              max.concurrent.jobs = 8),
            reg = reg
        )
    }
    cat("\nResubmissão concluída. Acompanhe com check_progress_calibration.R.\n")
}

resubmit_pending(reg)
