# simulation/simulation_aggregate.R
#
# Builds the three consolidated output files from every per-replica .rds
# checkpoint under results/partial/ (written by run_task() in
# simulation.R, whether the task ran locally or as a cluster job):
#
#   - summary_replicas.csv   -- 1 row per (method, f, Tt, replica); every
#                                SCALAR metric already computed by
#                                run_task() (coverage_rate_theta1,
#                                theta1_ci_width_mean, RMSE, MAE, ESS/s,
#                                etc.)
#   - summary_by_t.csv       -- 1 row per (method, f, Tt, t); Bias(t) with
#                                its Monte Carlo CI, width(t), and the
#                                estimator-dispersion band (mean +/- HDI
#                                across replicas) -- see metricas_theta1.tex
#                                Sections 2.2 and 3.6-3.7 for the formulas.
#   - summary_aggregated.csv -- 1 row per (method, f, Tt); RMSE/MAE
#                                mean+SD, global coverage + SE + Wilson CI,
#                                global width, ESS/s mean+SD, mean
#                                elapsed_time -- one row per grid cell,
#                                ready for the article/thesis tables.
#
# Unlike calibration_aggregate.R, this script does NOT recompute any
# diagnostics from raw chains -- run_task() already computes everything
# eagerly and saves only the reduced per-replica summary (keeping full
# chain histories for ~13,600 tasks was deliberately ruled out -- see
# planning discussion). Aggregation here is a plain recompute-from-scratch
# over the already-reduced per-replica values, cheap enough (a few
# thousand small .rds reads + array math) that there is no incremental/
# skip logic like calibration_aggregate.R's -- every run rebuilds all
# three files from every .rds present.

SIMULATION_DEFS_ONLY <- TRUE
rm(list = setdiff(ls(), "SIMULATION_DEFS_ONLY"))
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

# Reuses path_results, path_results_partial, path_data, printf, and
# metrics_theta_ci() (via the PoissonLTDM load_all() at the top of
# simulation.R) -- SIMULATION_DEFS_ONLY suppresses simulation.R's own
# dispatch block (run_mode == "local"/"cluster"), same guard mechanism as
# submit_simulation.R used earlier.
source("simulation_run.R", local = FALSE)
# No rm(SIMULATION_DEFS_ONLY, ...) needed here: simulation.R's own
# rm(list = setdiff(ls(), ".defs_only_flag")) at its top already destroys
# it as a side effect of being sourced.

printf("Aggregating simulation results from: %s", path_results_partial)

rds_files <- list.files(path_results_partial, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) {
	stop(sprintf("No .rds files found in %s -- has simulation_run.R been run yet?", path_results_partial))
}
printf("Found %d replica .rds file(s).", length(rds_files))

replicas <- do.call(rbind, lapply(rds_files, readRDS))
rownames(replicas) <- NULL

# ---- File 1: summary_replicas.csv (scalar columns only -- the 3
# list-columns below cannot be written to CSV) ----
list_cols <- c("theta1_mean", "theta1_ci_lower", "theta1_ci_upper")
scalar_cols <- setdiff(names(replicas), list_cols)
summary_replicas <- replicas[, scalar_cols]
write.csv(summary_replicas, file.path(path_results, "summary_replicas.csv"), row.names = FALSE)
printf("Wrote summary_replicas.csv (%d rows).", nrow(summary_replicas))

# ---- theta1_true(t) per (f, Tt): deterministic, identical across
# replicas -- read once per (f, Tt), from replica 1's data file. ----
theta1_true_for <- function(f, Tt) {
	d <- readRDS(sprintf("%s/%s_%s_1.rds", path_data, f, Tt))
	d$theta
}

cells <- unique(replicas[, c("method", "f", "Tt")])
rownames(cells) <- NULL

by_t_rows <- vector("list", nrow(cells))
agg_rows <- vector("list", nrow(cells))

z_wilson <- qnorm(1 - 0.05 / 2)  # gamma = 0.05: 95% CI on the coverage estimate (Eq. wilson_final, metricas_theta1.tex)
credmass_band <- 0.95            # alpha for the HDI band reused on replica-level means (Eq. banda_dispersao)

for (i in seq_len(nrow(cells))) {

	cell <- cells[i, ]
	idx <- which(replicas$method == cell$method & replicas$f == cell$f & replicas$Tt == cell$Tt)
	cell_rows <- replicas[idx, ]
	R <- nrow(cell_rows)
	Tt <- cell$Tt
	theta1_true <- theta1_true_for(cell$f, Tt)

	printf("[%d/%d] method=%s, f=%s, Tt=%d (R=%d)", i, nrow(cells), cell$method, cell$f, Tt, R)

	# ---- Stack per-t vectors across the R replicas of this cell: R x Tt matrices ----
	theta1_mean_mat <- do.call(rbind, replicas$theta1_mean[idx])
	ci_lower_mat    <- do.call(rbind, replicas$theta1_ci_lower[idx])
	ci_upper_mat    <- do.call(rbind, replicas$theta1_ci_upper[idx])

	# ---- summary_by_t.csv rows ----

	# Bias(t) (Eq. bias) + Monte Carlo SE/CI (Eq. bias_mcse)
	post_mean_agg <- colMeans(theta1_mean_mat)
	bias_t        <- post_mean_agg - theta1_true
	sd_t          <- apply(theta1_mean_mat, 2, sd)
	mcse_bias_t   <- sd_t / sqrt(R)
	t_crit        <- qt(1 - 0.05 / 2, df = R - 1)
	bias_ci_lower <- bias_t - t_crit * mcse_bias_t
	bias_ci_upper <- bias_t + t_crit * mcse_bias_t

	# width(t) (Eq. width_agg, marginal-by-time direction)
	width_t <- colMeans(ci_upper_mat - ci_lower_mat)

	# Estimator dispersion band (Eq. banda_dispersao): HDI across
	# replicas' theta1_mean, per t -- reuses metrics_theta_ci() exactly as
	# documented, with the R replica-level means as the input "sample"
	# instead of one replica's posterior draws.
	band <- metrics_theta_ci(theta1_mean_mat, credmass_band)

	by_t_rows[[i]] <- data.frame(
		method = cell$method, f = cell$f, Tt = Tt, t = seq_len(Tt),
		bias = bias_t, bias_mcse = mcse_bias_t,
		bias_ci_lower = bias_ci_lower, bias_ci_upper = bias_ci_upper,
		width_mean = width_t,
		post_mean_agg = post_mean_agg,
		band_lower = band$ci_lower, band_upper = band$ci_upper,
		stringsAsFactors = FALSE
	)

	# ---- summary_aggregated.csv row ----

	rmse_mean <- mean(cell_rows$rmse_theta1)
	rmse_sd   <- sd(cell_rows$rmse_theta1)
	mae_mean  <- mean(cell_rows$mae_theta1)
	mae_sd    <- sd(cell_rows$mae_theta1)

	# Global coverage (Eq. cov_global) + SE (Eq. se_cobertura) + Wilson CI
	# (Eqs. wilson_final/wilson_center)
	coverage    <- mean(cell_rows$coverage_rate_theta1)
	coverage_se <- sqrt(coverage * (1 - coverage) / R)
	wilson_center <- (R * coverage + z_wilson^2 / 2) / (R + z_wilson^2)
	wilson_margin <- (z_wilson / (R + z_wilson^2)) * sqrt(R * coverage * (1 - coverage) + z_wilson^2 / 4)
	coverage_wilson_lower <- wilson_center - wilson_margin
	coverage_wilson_upper <- wilson_center + wilson_margin

	width_mean_agg <- mean(cell_rows$theta1_ci_width_mean)

	agg_rows[[i]] <- data.frame(
		method = cell$method, f = cell$f, Tt = Tt, R = R,
		rmse_mean = rmse_mean, rmse_sd = rmse_sd,
		mae_mean = mae_mean, mae_sd = mae_sd,
		coverage = coverage, coverage_se = coverage_se,
		coverage_wilson_lower = coverage_wilson_lower,
		coverage_wilson_upper = coverage_wilson_upper,
		width_mean = width_mean_agg,
		ess_sec_theta1_mean = mean(cell_rows$ess_sec_theta1_mean),
		ess_sec_theta1_sd = sd(cell_rows$ess_sec_theta1_mean),
		ess_sec_theta_01_mean = mean(cell_rows$ess_sec_theta_01),
		elapsed_time_mean = mean(cell_rows$elapsed_time),
		stringsAsFactors = FALSE
	)
}

summary_by_t <- do.call(rbind, by_t_rows)
summary_aggregated <- do.call(rbind, agg_rows)

write.csv(summary_by_t, file.path(path_results, "summary_by_t.csv"), row.names = FALSE)
write.csv(summary_aggregated, file.path(path_results, "summary_aggregated.csv"), row.names = FALSE)

printf("Wrote summary_by_t.csv (%d rows) and summary_aggregated.csv (%d rows).",
	   nrow(summary_by_t), nrow(summary_aggregated))
