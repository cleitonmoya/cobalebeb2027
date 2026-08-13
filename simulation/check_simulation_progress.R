# simulation/check_simulation_progress.R
#
# Read-only progress/health check for the production simulation run.
# Safe to run anytime, as often as needed -- NEVER submits, modifies, or
# archives anything (SIMULATION_DEFS_ONLY suppresses simulation.R's
# dispatch entirely, same guard used throughout this project).
#
# Two independent checks:
#
#   1. Data-level completeness (ground truth). Compares actual .rds
#      checkpoint files in results/partial/ against the full expected
#      task_grid, broken down by category and by grid cell
#      (method x f x Tt). This is what actually matters -- run_task()
#      itself only ever trusts file existence, never batchtools'
#      bookkeeping, so this number is always correct regardless of
#      anything below.
#
#   2. batchtools job status (operational visibility only). Queue/
#      running/error state for whatever is in registry_simulation/ right
#      now. Includes a specific safeguard against the known "Expired"
#      false positive (NFS directory-listing staleness in batchtools'
#      sync() -- see /areas/euler-cluster-simulation.md and the planning
#      discussion): any job reported as Expired is re-checked once after
#      a short wait before being flagged as a real problem, since the
#      underlying files are frequently already correct and simply
#      weren't visible yet at the first sync. Genuine errors (an actual
#      R error captured with a message) are NOT subject to this
#      staleness issue and are reported directly.

SIMULATION_DEFS_ONLY <- TRUE
rm(list = setdiff(ls(), "SIMULATION_DEFS_ONLY"))
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

source("simulation_run.R", local = FALSE)

printf("========================================")
printf("Progress check: %s", format(Sys.time()))
printf("========================================")

# ---- 1. Data-level completeness (ground truth) ----

task_grid$done <- file.exists(task_result_filename(task_grid))

overall_done <- sum(task_grid$done)
overall_total <- nrow(task_grid)
printf("\n--- Data completeness (ground truth: .rds files present) ---")
printf("Overall: %d / %d tasks done (%.1f%%)", overall_done, overall_total, 100 * overall_done / overall_total)

for (cat in unique(task_grid$category)) {
	sub <- task_grid[task_grid$category == cat, ]
	printf("  %-8s: %5d / %5d (%.1f%%)", cat, sum(sub$done), nrow(sub), 100 * sum(sub$done) / nrow(sub))
}

cell_key <- paste(task_grid$method, task_grid$f, task_grid$Tt)
cell_done <- tapply(task_grid$done, cell_key, sum)
cell_total <- tapply(task_grid$done, cell_key, length)
incomplete <- names(cell_done)[cell_done < cell_total]

if (length(incomplete) > 0) {
	printf("\nGrid cells not yet fully done (%d of %d):", length(incomplete), length(cell_done))
	for (k in incomplete) {
		printf("  %-40s: %4d / %4d", k, cell_done[[k]], cell_total[[k]])
	}
} else {
	printf("\nAll grid cells fully done.")
}

# ---- 2. batchtools job status (operational visibility) ----

printf("\n--- batchtools status ---")

if (!dir.exists("registry_simulation")) {

	printf("No active registry_simulation/ found (nothing currently submitted, or already archived after finishing).")

} else {

	library(batchtools)

	# ALL of this section's state (reg, expired_ids, still_expired,
	# err_ids) is kept local to check_batchtools_status()'s own call
	# frame, deliberately -- not because of any nesting preference, but
	# because loadRegistry() re-sources simulation.R internally (same
	# mechanism as the initial source() above), which runs
	# rm(list = setdiff(ls(), ".defs_only_flag")) at its own top. That rm()
	# wipes EVERYTHING in .GlobalEnv except that one flag -- confirmed in
	# practice on calibration_run.R's identical pattern (see
	# check_calibration_progress.R), on two different kinds of casualty:
	#   1. A named helper function (an earlier version wrapped the
	#      "assign flag + loadRegistry()" pattern in a `safe_load_registry()`
	#      function -- wiped by its own FIRST call's internal re-source,
	#      so the SECOND call site failed with "could not find function
	#      'safe_load_registry'").
	#   2. A plain variable (`expired_ids`, needed AFTER the second
	#      loadRegistry() call for the re-check comparison, wiped the same
	#      way -- "object 'expired_ids' not found" on the exact next line
	#      that used it, even after fixing (1) by inlining instead of a
	#      shared function).
	# Both failures share one root cause: anything living in .GlobalEnv
	# when loadRegistry() is called does not reliably survive the call.
	# Wrapping the WHOLE section in a function sidesteps this entirely --
	# reg/expired_ids/still_expired/err_ids then live in this function's
	# own local execution frame (a child environment of .GlobalEnv,
	# created fresh when the function is called), which the rm() (via
	# sys.source(fn, envir = .GlobalEnv)) never touches. Verified directly
	# (see check_calibration_progress.R's testing): a local variable
	# inside a wrapping function survives an rm(list=ls()) that targets
	# .GlobalEnv, even though an identically-named .GlobalEnv variable
	# does not.
	check_batchtools_status <- function() {
		assign("SIMULATION_DEFS_ONLY", TRUE, envir = .GlobalEnv)
		reg <- loadRegistry("registry_simulation", writeable = FALSE)
		print(getStatus(reg = reg))

		expired_ids <- findExpired(reg = reg)
		if (nrow(expired_ids) > 0) {
			printf("\n%d job(s) show 'Expired' -- re-checking after a short wait", nrow(expired_ids))
			printf("(known NFS staleness false-positive; see /areas/euler-cluster-simulation.md)...")
			Sys.sleep(10)
			assign("SIMULATION_DEFS_ONLY", TRUE, envir = .GlobalEnv)
			reg <- loadRegistry("registry_simulation", writeable = FALSE)  # re-sync
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
