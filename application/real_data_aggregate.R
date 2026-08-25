# application/real_data_aggregate.R
#
# DISPOSABLE utility script -- NOT part of the permanent pipeline. Use only
# when the 15 chain-unit checkpoints (../results/real_data/chains/*.rds)
# already exist from a prior run, and aggregate_method() in real_data_run.R
# has changed since then (e.g., the theta1/theta2 95% HPD band added for
# application/figures_real_data.R) -- this re-derives method_summaries.rds,
# delta_max.rds and summary.csv from those checkpoints WITHOUT re-running
# any sampling and WITHOUT spawning a parallel cluster (sequential,
# lightweight -- safe to run directly on the login node, no qsub needed).
#
# Anyone replicating the application from scratch does NOT need this
# script: real_data_run.R alone runs the 15 chains AND aggregates them
# (with the current aggregate_method(), band included) in one pass.
#
# Uses real_data_run.R's REAL_DATA_DEFS_ONLY escape hatch to load every
# definition this needs (chain_grid, chain_result_filename(),
# aggregate_method(), methods_grid_ref, R_config, y, Tt, and the package's
# metrics_* functions) without triggering its run_chain_units() /
# aggregation block.

REAL_DATA_DEFS_ONLY <- TRUE
source("real_data_run.R")

printf("Re-aggregating from %d existing chain-unit checkpoint(s) (no resampling, no parallel cluster).", nrow(chain_grid))

all_chain_files <- vapply(seq_len(nrow(chain_grid)), function(i) chain_result_filename(chain_grid[i, ]), character(1))
missing <- all_chain_files[!file.exists(all_chain_files)]
if (length(missing) > 0) {
	stop(sprintf("Missing %d chain checkpoint(s), cannot aggregate:\n%s",
				 length(missing), paste(missing, collapse = "\n")))
}

chain_results <- lapply(all_chain_files, readRDS)

method_summaries <- lapply(methods_grid_ref, function(m) {
	method_chains <- Filter(function(ch) ch$method == m, chain_results)
	aggregate_method(method_chains)
})
names(method_summaries) <- methods_grid_ref

# Cross-method agreement -- Eq. 32, metricas_theta1.tex. theta1/theta2 are
# curves over t; W1/W2 are single posterior-mean scalars per method, so
# their Delta_max is a single number, not a length-Tt vector
# (metrics_agreement_max() handles both the same way).
theta1_means_by_method <- lapply(method_summaries, function(s) s$theta1_mean)
theta2_means_by_method <- lapply(method_summaries, function(s) s$theta2_mean)
W1_means_by_method <- lapply(method_summaries, function(s) s$W1_mean)
W2_means_by_method <- lapply(method_summaries, function(s) s$W2_mean)
delta_max_theta1 <- metrics_agreement_max(theta1_means_by_method)
delta_max_theta2 <- metrics_agreement_max(theta2_means_by_method)
delta_max_W1 <- metrics_agreement_max(W1_means_by_method)
delta_max_W2 <- metrics_agreement_max(W2_means_by_method)

summary_df <- do.call(rbind, lapply(method_summaries, function(s) {
	data.frame(
		method = s$method, total_time = s$total_time,

		rhat_theta1_max = s$rhat_theta1_max, rhat_theta2_max = s$rhat_theta2_max,

		ess_bulk_theta1_mean = s$ess_bulk_theta1_mean, ess_tail_theta1_mean = s$ess_tail_theta1_mean,
		ess_bulk_theta1_min  = s$ess_bulk_theta1_min,  ess_tail_theta1_min  = s$ess_tail_theta1_min,
		ess_bulk_theta2_mean = s$ess_bulk_theta2_mean, ess_tail_theta2_mean = s$ess_tail_theta2_mean,

		ess_sec_theta1_bulk_mean = s$ess_sec_theta1_bulk_mean, ess_sec_theta1_tail_mean = s$ess_sec_theta1_tail_mean,
		ess_sec_theta2_bulk_mean = s$ess_sec_theta2_bulk_mean, ess_sec_theta2_tail_mean = s$ess_sec_theta2_tail_mean,

		log_lik = s$log_lik, log_cpo = s$log_cpo,
		stringsAsFactors = FALSE
	)
}))
rownames(summary_df) <- NULL

saveRDS(method_summaries, file = sprintf("%s/method_summaries.rds", path_results))
saveRDS(list(theta1 = delta_max_theta1, theta2 = delta_max_theta2, W1 = delta_max_W1, W2 = delta_max_W2),
		file = sprintf("%s/delta_max.rds", path_results))
write.csv(summary_df, file = sprintf("%s/summary.csv", path_results), row.names = FALSE)

printf("Re-aggregation done. Summary saved to %s/summary.csv", path_results)
print(summary_df)
