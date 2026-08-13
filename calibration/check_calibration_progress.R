# calibration/check_calibration_progress.R
#
# Read-only progress/health check for the calibration phase. Safe to run
# anytime, as often as needed -- NEVER submits, modifies, or archives
# anything (CALIBRATION_PHASE_DEFS_ONLY suppresses calibration_run.R's
# dispatch entirely, same guard used throughout this project).
#
# A separate script from check_progress.R (not a shared one) because the
# two phases differ in more than just names: calibration_run.R's
# flattened chunking (see its own header comments) means there are TWO
# independent levels of completion worth reporting here, not one --
#
#   1. Task-level completion (ground truth for the calibration run as a
#      whole). Compares the final per-task .rds files in
#      results/calibration/ against the full expected calibration_grid,
#      broken down by category and by grid cell (method x f x Tt x
#      config_idx). aggregate_task() itself only ever trusts file
#      existence, never batchtools' bookkeeping, so this number is always
#      correct regardless of anything below.
#
#   2. Chain-level completion (visibility INSIDE an in-progress chunk).
#      A chunk's job is not "half done" or "not done" from task-level
#      completion alone -- it may have 15 of its 18 chain-units finished
#      and be about to produce several tasks' final .rds in quick
#      succession, or it may be stuck on its first chain-unit. This checks
#      results/calibration/chains_tmp/ directly against chain_grid to show
#      that finer-grained progress, which is otherwise invisible until a
#      whole chunk's tasks all land at once.
#
#   3. batchtools job status (operational visibility only). Queue/
#      running/error state for whatever is in registry_calibration/ right
#      now. Same "Expired" NFS-staleness safeguard as check_progress.R --
#      re-checked once after a short wait before being flagged as a real
#      problem.

CALIBRATION_PHASE_DEFS_ONLY <- TRUE
rm(list = setdiff(ls(), "CALIBRATION_PHASE_DEFS_ONLY"))
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

source("calibration_run.R", local = FALSE)

printf("========================================")
printf("Calibration progress check: %s", format(Sys.time()))
printf("========================================")

# ---- 1. Task-level completion (ground truth) ----

task_names <- with(calibration_grid, mapply(task_name_for, method, f, Tt, N, burnin, K))
calibration_grid$done <- file.exists(file.path(path_results, paste0(task_names, ".rds")))

overall_done <- sum(calibration_grid$done)
overall_total <- nrow(calibration_grid)
printf("\n--- Task completion (ground truth: final .rds files present) ---")
printf("Overall: %d / %d tasks done (%.1f%%)", overall_done, overall_total, 100 * overall_done / overall_total)

for (cat in unique(calibration_grid$category)) {
    sub <- calibration_grid[calibration_grid$category == cat, ]
    printf("  %-8s: %5d / %5d (%.1f%%)", cat, sum(sub$done), nrow(sub), 100 * sum(sub$done) / nrow(sub))
}

cell_key <- with(calibration_grid, paste(method, f, Tt, N, burnin, K))
cell_done <- tapply(calibration_grid$done, cell_key, sum)
cell_total <- tapply(calibration_grid$done, cell_key, length)
incomplete <- names(cell_done)[cell_done < cell_total]

if (length(incomplete) > 0) {
    printf("\nGrid cells not yet fully done (%d of %d):", length(incomplete), length(cell_done))
    for (k in incomplete) {
        printf("  %-55s: %4d / %4d", k, cell_done[[k]], cell_total[[k]])
    }
} else {
    printf("\nAll grid cells fully done.")
}

# ---- 2. Chain-level completion (visibility inside in-progress chunks) ----

printf("\n--- Chain-level completion (results/calibration/chains_tmp/) ---")

# NOTE: this necessarily UNDERCOUNTS chain-level progress for tasks that
# are already task-done -- aggregate_task() deletes every chain checkpoint
# except plot_chain's once a task's final .rds exists (see its own header
# comment), so a fully-done task always shows N_chains-1 "missing" chain
# files here even though every chain genuinely ran. This section is only
# meant to show progress WITHIN chunks that are still in flight; use
# section 1 above for ground-truth task completion.
incomplete_grid <- calibration_grid[!calibration_grid$done, ]
if (nrow(incomplete_grid) == 0) {
    printf("No incomplete tasks -- nothing to show here (see section 1 above).")
} else {
    chain_files <- character(0)
    for (i in seq_len(nrow(incomplete_grid))) {
        t <- incomplete_grid[i, ]
        chain_files <- c(chain_files,
                          vapply(seq_len(N_chains),
                                 function(k) chain_result_filename(t$method, t$f, t$Tt, t$N, t$burnin, t$K, k),
                                 character(1)))
    }
    n_chain_done <- sum(file.exists(chain_files))
    printf("Among %d not-yet-aggregated task(s): %d / %d chain-unit checkpoint(s) present (%.1f%%)",
           nrow(incomplete_grid), n_chain_done, length(chain_files), 100 * n_chain_done / length(chain_files))
}

# ---- 3. batchtools job status (operational visibility) ----

printf("\n--- batchtools status ---")

if (!dir.exists("registry_calibration")) {

    printf("No active registry_calibration/ found (nothing currently submitted, or already archived after finishing).")

} else {

    library(batchtools)

    # ALL of this section's state (reg, expired_ids, still_expired,
    # err_ids) is kept local to check_batchtools_status()'s own call
    # frame, deliberately -- not because of any nesting preference, but
    # because loadRegistry() re-sources calibration_run.R internally
    # (same mechanism as the initial source() above), which runs
    # rm(list = setdiff(ls(), ".defs_only_flag")) at its own top. That rm()
    # wipes EVERYTHING in .GlobalEnv except that one flag -- confirmed in
    # practice TWICE now, on two different kinds of casualty:
    #   1. A named helper function (an earlier version wrapped the
    #      "assign flag + loadRegistry()" pattern in a `safe_load_registry()`
    #      function -- that function object, living in .GlobalEnv, was
    #      wiped by its own FIRST call's internal re-source, so the SECOND
    #      call site failed with "could not find function
    #      'safe_load_registry'").
    #   2. A plain variable (`expired_ids`, needed AFTER the second
    #      loadRegistry() call for the re-check comparison, was ALSO a
    #      .GlobalEnv variable and got wiped the same way, failing with
    #      "object 'expired_ids' not found" on the exact next line that
    #      used it).
    # Both failures share one root cause: anything living in .GlobalEnv
    # when loadRegistry() is called does not reliably survive the call.
    # Rather than keep discovering casualties one at a time and inlining
    # around each individually, this wraps the WHOLE section in a
    # function -- reg/expired_ids/still_expired/err_ids then live in this
    # function's own local execution frame (a child environment of
    # .GlobalEnv, created fresh when the function is called), which
    # calibration_run.R's rm(list = ..., i.e. operating on .GlobalEnv,
    # via sys.source(fn, envir = .GlobalEnv)) never touches. Verified
    # directly: a local variable inside a wrapping function survives an
    # rm(list=ls()) that targets .GlobalEnv, even though an identically-
    # named .GlobalEnv variable does not.
    check_batchtools_status <- function() {
        assign("CALIBRATION_PHASE_DEFS_ONLY", TRUE, envir = .GlobalEnv)
        reg <- loadRegistry("registry_calibration", writeable = FALSE)
        print(getStatus(reg = reg))

        expired_ids <- findExpired(reg = reg)
        if (nrow(expired_ids) > 0) {
            printf("\n%d job(s) show 'Expired' -- re-checking after a short wait", nrow(expired_ids))
            printf("(known NFS staleness false-positive)...")
            Sys.sleep(10)
            assign("CALIBRATION_PHASE_DEFS_ONLY", TRUE, envir = .GlobalEnv)
            reg <- loadRegistry("registry_calibration", writeable = FALSE)  # re-sync
            still_expired <- findExpired(ids = expired_ids, reg = reg)
            if (nrow(still_expired) > 0) {
                printf("Still %d 'Expired' after re-check -- worth a closer look (job.id: %s).",
                       nrow(still_expired), paste(still_expired$job.id, collapse = ", "))
            } else {
                printf("Resolved after re-check -- was staleness, not a real problem.")
            }
        }

        err_ids <- findErrors(reg = reg)
        if (nrow(err_ids) > 0) {
            printf("\n%d job(s) with a REAL error (not staleness):", nrow(err_ids))
            print(getErrorMessages(ids = err_ids, reg = reg))
        }
    }
    check_batchtools_status()
}

printf("\n========================================")
