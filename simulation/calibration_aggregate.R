# simulation/calibration_aggregate.R
#
# Builds/updates results/calibration/summary.csv from every .rds saved
# under results/calibration/ (one per task, written by
# run_calibration_task() whether the task ran locally or as a cluster
# job). Needed specifically for run_mode == "cluster" in
# calibration_phase.R: concurrent PBS jobs do not each write to the shared
# summary.csv (see run_and_save_task()'s header comment there for why), so
# this script is the step that produces it, run once all jobs have
# finished. Also needed for calibration_phase.R's checkpoint: a
# checkpointed (already-run) task is skipped ENTIRELY there -- no
# diagnostics recomputation, no summary.csv row -- so this script is what
# actually produces/fills that row.
#
# INCREMENTAL: reads any existing summary.csv first and keeps its rows for
# tasks that are still present there -- those are NOT recomputed. Only
# .rds files with no matching (method, f, Tt, replica) row in the existing
# summary.csv are processed. This matters because diagnostics
# recomputation, while much cheaper than re-running a sampler, is not free
# (it re-does the R_hat/ESS pass over every chain's full history) --
# there's no reason to redo it for tasks whose row is already known good.
# Delete summary.csv (or the specific row) by hand to force a task to be
# reprocessed.
#
# Does NOT re-run any sampler: each .rds already contains the full
# per-chain results (theta1_hist, theta2_hist, theta_01_hist, theta_02_hist,
# W1_hist, W2_hist, and any algorithm-specific fields), so R_hat/ESS are
# recomputed from that (via print_and_plot_diagnostics(..., plots = FALSE),
# which is cheap -- just the metrics_convergence()/metrics_convergence_by_chain()
# calls, no MCMC/SMC/HMC work), and build_summary_row() (defined in
# calibration_phase.R) turns the result into one summary.csv row per task.

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

# Reuses build_summary_row(), flatten_chains(), nz(), path_results, and
# print_and_plot_diagnostics() from calibration_phase.R -- sourcing it here
# pulls in those definitions (plus method_grid/Tt_grid/function_grid,
# unused here but harmless) without duplicating them. Setting
# CALIBRATION_PHASE_DEFS_ONLY <- TRUE first suppresses calibration_phase.R's
# own dispatch block (run_mode == "local"/"cluster") -- see the guard
# comment there -- so sourcing it here only defines things, it does not
# run a task or submit cluster jobs.
CALIBRATION_PHASE_DEFS_ONLY <- TRUE
source("calibration_phase.R")

printf("Aggregating calibration results from: %s", path_results)

rds_files <- list.files(path_results, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) {
    stop(sprintf("No .rds files found in %s -- has calibration_phase.R been run yet?", path_results))
}
printf("Found %d task .rds file(s).", length(rds_files))

# ---- Existing summary.csv: keep rows already present, only add new ones ----
summary_file <- file.path(path_results, "summary.csv")
if (file.exists(summary_file)) {
    summary_existing <- read.csv(summary_file, stringsAsFactors = FALSE)
    printf("Existing summary.csv has %d row(s); those tasks will be skipped.", nrow(summary_existing))
} else {
    summary_existing <- NULL
}

# task_key(): (method, f, Tt, replica) collapsed into one string per row,
# used to match a .rds file's task against summary_existing's rows without
# a 4-column merge.
task_key <- function(method, f, Tt, replica) sprintf("%s|%s|%s|%s", method, f, Tt, replica)

existing_keys <- if (is.null(summary_existing)) character(0) else {
    with(summary_existing, task_key(method, f, Tt, replica))
}

new_rows <- list()

for (i in seq_along(rds_files)) {
    rds_file <- rds_files[i]

    # Peek at just enough of the .rds to know its task key, before
    # deciding whether to load/process the (potentially large) results
    # field -- readRDS() itself still reads the whole file (R has no
    # partial-read for .rds), but this keeps the skip/process decision in
    # one place and the loop body symmetric with the non-incremental
    # version this replaced.
    task <- readRDS(rds_file)
    method  <- task$method
    f       <- task$f
    Tt      <- task$Tt
    replica <- task$replica

    key <- task_key(method, f, Tt, replica)
    if (key %in% existing_keys) {
        printf("[%d/%d] %s -- already in summary.csv, skipping.", i, length(rds_files), basename(rds_file))
        next
    }

    printf("[%d/%d] %s -- new, computing diagnostics.", i, length(rds_files), basename(rds_file))

    # Fields saved by run_calibration_task(): results, elapsed_time,
    # method, f, Tt, replica, N, burnin, K, seed_base.
    results      <- task$results
    N            <- task$N
    burnin       <- task$burnin
    K            <- task$K
    elapsed_time <- task$elapsed_time

    # Reload the data for this task, needed by print_and_plot_diagnostics()
    # for theta1_true / t_obs / changepoints (same lookup as
    # run_calibration_task() in calibration_phase.R).
    source_name <- sprintf("%s_%s_%s", f, Tt, replica)
    data <- readRDS(file.path("..", "data", "simulated", paste0(source_name, ".rds")))
    theta1_true  <- data$theta
    changepoints <- c(0.25, 0.5, 0.75) * Tt

    if (Tt == 200)  t_obs <- c(50, 100, 150, 175)
    if (Tt == 400)  t_obs <- c(75, 100, 200, 300)
    if (Tt == 800)  t_obs <- c(100, 300, 500, 700)
    if (Tt == 1600) t_obs <- c(400, 800, 1200, 1600)

    # plots = FALSE: this step only needs the returned R_hat/ESS values,
    # not a re-rendered PDF (each task's plot was already saved once, when
    # the task itself ran).
    diag <- print_and_plot_diagnostics(
        results, data$y, changepoints,
        theta1_true = theta1_true, theta2_true = NULL,
        t_obs = t_obs, burnin = burnin, elapsed_time = elapsed_time,
        plot_chain = plot_chain,
        nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
        ac_ref = ac_ref, plots = FALSE,
        compute_rhat = compute_rhat, compute_ess = compute_ess)

    new_rows[[length(new_rows) + 1]] <- build_summary_row(diag, method, f, Tt, replica, N, burnin, K, elapsed_time)
}

if (length(new_rows) == 0) {
    printf("No new tasks to add -- summary.csv is already up to date (%d row(s)).",
           if (is.null(summary_existing)) 0 else nrow(summary_existing))
} else {
    summary_new <- do.call(rbind, new_rows)
    summary_all <- if (is.null(summary_existing)) summary_new else rbind(summary_existing, summary_new)
    write.csv(summary_all, summary_file, row.names = FALSE)
    printf("Summary written to: %s (%d existing + %d new = %d task(s) total)",
           summary_file,
           if (is.null(summary_existing)) 0 else nrow(summary_existing),
           length(new_rows), nrow(summary_all))
}
