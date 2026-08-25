# simulation/compute_pointwise_coverage.R
#
# Builds a pointwise-by-t coverage table (method, f, Tt, t, coverage) from
# the SAME per-replica .rds checkpoints simulation_aggregate.R already
# reads -- theta1_ci_lower/theta1_ci_upper (length-Tt list-columns) are
# already saved by run_task() (simulation_run.R), so no re-simulation is
# needed. This is a read-only, standalone companion to
# simulation_aggregate.R: it does not modify or overwrite
# summary_replicas.csv / summary_by_t.csv / summary_aggregated.csv.
#
# Motivation: summary_aggregated.csv's "coverage" column is the coverage
# rate averaged over BOTH replicas and t (Eq. cov_global). That average
# cannot distinguish "coverage is uniformly ~0.95 across t" from
# "coverage is ~1.0 almost everywhere but drops sharply near a few t's" --
# exactly the question of whether overcoverage on linear/quadratic/
# sinusoidal is a uniform model-level effect, versus the constant
# function's coverage being dragged down by localized breakpoint failures
# (see chat discussion). This script produces the by-t breakdown needed to
# tell those two stories apart.
#
# Output: summary_pointwise_coverage.csv, one row per (method, f, Tt, t):
#   coverage = fraction of replicas for which theta1_true(t) fell inside
#              that replica's HPD interval at that specific t.

SIMULATION_DEFS_ONLY <- TRUE
rm(list = setdiff(ls(), "SIMULATION_DEFS_ONLY"))
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

# Reuses path_results, path_results_partial, path_data, and printf --
# same sourcing pattern as simulation_aggregate.R.
source("simulation_run.R", local = FALSE)

printf("Computing pointwise coverage from: %s", path_results_partial)

rds_files <- list.files(path_results_partial, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) {
	stop(sprintf("No .rds files found in %s -- has simulation_run.R been run yet?", path_results_partial))
}
printf("Found %d replica .rds file(s).", length(rds_files))

replicas <- do.call(rbind, lapply(rds_files, readRDS))
rownames(replicas) <- NULL

if (is.null(replicas$theta1_ci_lower) || is.null(replicas$theta1_ci_upper)) {
	stop("theta1_ci_lower/theta1_ci_upper not found in the saved .rds -- ",
		 "these are required for pointwise coverage (see check_rds_structure.R).")
}

# theta1_true(t) per (f, Tt): deterministic, identical across replicas --
# read once per (f, Tt), from replica 1's data file. Same helper as
# simulation_aggregate.R.
theta1_true_for <- function(f, Tt) {
	d <- readRDS(sprintf("%s/%s_%s_1.rds", path_data, f, Tt))
	return(d$theta)
}

cells <- unique(replicas[, c("method", "f", "Tt")])
rownames(cells) <- NULL

pointwise_rows <- vector("list", nrow(cells))

for (i in seq_len(nrow(cells))) {

	cell <- cells[i, ]
	idx <- which(replicas$method == cell$method & replicas$f == cell$f & replicas$Tt == cell$Tt)
	R <- length(idx)
	Tt <- cell$Tt
	theta1_true <- theta1_true_for(cell$f, Tt)

	printf("[%d/%d] method=%s, f=%s, Tt=%d (R=%d)", i, nrow(cells), cell$method, cell$f, Tt, R)

	# ---- Stack per-t vectors across the R replicas of this cell: R x Tt matrices ----
	ci_lower_mat <- do.call(rbind, replicas$theta1_ci_lower[idx])
	ci_upper_mat <- do.call(rbind, replicas$theta1_ci_upper[idx])

	# Pointwise coverage indicator C(t, r), averaged over r (replicas) --
	# same indicator as coverage_rate_theta1 in run_task(), but kept per
	# t instead of collapsed with mean() over t.
	true_mat <- matrix(theta1_true, nrow = R, ncol = Tt, byrow = TRUE)
	covered_mat <- (true_mat >= ci_lower_mat) & (true_mat <= ci_upper_mat)
	coverage_t <- colMeans(covered_mat)

	pointwise_rows[[i]] <- data.frame(
		method = cell$method, f = cell$f, Tt = Tt, t = seq_len(Tt),
		coverage = coverage_t, R = R,
		stringsAsFactors = FALSE
	)
}

summary_pointwise_coverage <- do.call(rbind, pointwise_rows)

write.csv(summary_pointwise_coverage, file.path(path_results, "summary_pointwise_coverage.csv"), row.names = FALSE)

printf("Wrote summary_pointwise_coverage.csv (%d rows).", nrow(summary_pointwise_coverage))
