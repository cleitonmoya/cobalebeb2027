# calibration/calibration_run.R
#
# Official, reproducible run of the calibration-phase grid: a cross-join
# of method_configs (one row per (method, N, burnin, K) configuration --
# a method may appear more than once, to compare configurations) with
# (f, Tt) (2 functions x 2 Tt), replica = 1, N_chains = 3. Documents the
# N/burnin/K configuration(s) decided for each method during the
# pre-calibration discussion. Reuses, unchanged, the same
# chain-initialization scheme, hyperparameters, and diagnostic pipeline as
# tests/test_cpp.R and tests/test_stan.R -- this script is not a
# reimplementation, it is those two scripts made parametrizable over the
# calibration grid, plus disk output (raw results + plots + a one-row
# summary) for later aggregation.
#
# Two independent choices, both set in the "Run control" block below:
#
# 1. run_mode: "local" or "cluster".
#      "local"   -- runs ONE task (task_id below) in this R session, exactly
#                   as before: same process, same console output, updates
#                   results/calibration/summary.csv directly. For quick
#                   testing/debugging of a single grid row.
#      "cluster" -- submits every row of the (possibly subset) calibration
#                   grid as ONE PBS job each, via batchtools, using
#                   calibration_pbs.tmpl (ncpus=3, place=excl -- see that
#                   file). Each job runs run_and_save_task() for its own
#                   grid row, in its own R process, saving .rds + .pdf as
#                   usual. submitJobs() returns immediately after
#                   dispatching -- it does NOT wait for jobs to finish and
#                   does NOT write summary.csv. Aggregation is a separate
#                   step: run calibration_aggregate.R once the jobs are
#                   done (see that file; check progress with
#                   batchtools::getStatus(reg)).
#
# 2. grid_subset: NULL (run the full calibration_grid) or a subset filter
#      applied to it (see "Grid subset" below) -- e.g. only the "stan" rows,
#      or only method=="amh_montoril" & N==110000, to submit/run a smaller
#      batch without editing the grid definition itself. Applies
#      identically in both run_mode values: in "local" mode task_id indexes
#      into the FILTERED grid; in "cluster" mode every row of the FILTERED
#      grid becomes one job.

# Capture CALIBRATION_PHASE_DEFS_ONLY (if the caller -- e.g.
# calibration_aggregate.R -- set it before source()ing this file) BEFORE
# rm(list = ls()) below wipes it out along with everything else in this
# environment. See the "Guard against recursive dispatch" block further
# down for how .defs_only_flag is used.
.defs_only_flag <- exists("CALIBRATION_PHASE_DEFS_ONLY", inherits = FALSE) &&
                    isTRUE(CALIBRATION_PHASE_DEFS_ONLY)

rm(list = setdiff(ls(), ".defs_only_flag"))
options(error = function() traceback(2))

# NOTE (untested in a real batchtools/PBS job): this.path::this.path() is
# designed for direct Rscript/source() invocation. Whether it still
# resolves correctly when batchtools re-sources this file from inside
# doJobCollection() (a different call stack) has not been verified here.
# If a cluster-mode job's log shows a setwd()/this.path error, the fix is
# to hardcode the known absolute path instead, e.g.:
#   setwd("/mnt/nfs/home/cleiton/poisson_ltdm/simulation")
# rather than relying on this.path() inside job execution.
setwd(dirname(this.path::this.path()))

# See test_cpp.R for why debug=FALSE matters here (elapsed_time feeds
# directly into ESS/second diagnostics inside print_and_plot_diagnostics()).
pkgload::load_all("../PoissonLTDM", debug = FALSE)

source("../PoissonLTDM/R/sampler_stan.R")
source("../tests/plot_diagnostics.R")

printf <- function(...) cat(paste(sprintf(...), "\n"))

# ---- archive_registry(): move a registry directory into registry_archive/
# with a timestamp, preserving its logs and full job history -- used by
# run_mode == "cluster"'s auto-clean step (see there for the full
# rationale). Defined here, at top level (like printf above), rather than
# inside the "cluster" dispatch branch, because loadRegistry() -- called
# from within that auto-clean step -- re-sources this whole file via
# loadRegistryDependencies(), which resets .GlobalEnv; anything defined
# only inside the (now-guarded) dispatch block would NOT survive that
# reset and be gone by the time archive_registry() is actually called
# afterwards. This was caught during testing as an "object
# 'archive_registry' not found"-style failure, the same root cause as the
# "object 'registry_dir' not found" error fixed alongside it.
registry_archive_dir <- "registry_archive"
archive_registry <- function(reason, registry_dir) {
    dir.create(registry_archive_dir, showWarnings = FALSE, recursive = TRUE)
    archived_path <- file.path(registry_archive_dir,
        paste0("registry_calibration_", format(Sys.time(), "%Y-%m-%d_%H%M%S")))
    printf("%s -- archiving to: %s", reason, archived_path)
    file.rename(registry_dir, archived_path)
}

# ---- task_name_for(): canonical task/file-name for one (method, f, Tt, N,
# burnin, K) combination ----
#
# <method>_<f>_<Tt>_N<N>_b<burnin>[_K<K>] -- same convention already
# applied by hand to the pre-existing calibration .rds/.pdf files via
# rename_calibration_results.sh (K suffix present only when not NA, e.g.
# pg_as). Centralized here (used for task_rds, plot_file, and by
# calibration_aggregate.R) so the naming rule is defined in exactly one
# place. Defined at top level, like printf/archive_registry above, for the
# same reason: it must survive a loadRegistry()/makeRegistry()-triggered
# re-source of this file.
task_name_for <- function(method, f, Tt, N, burnin, K) {
    base <- sprintf("%s_%s_%s_N%s_b%s", method, f, Tt, N, burnin)
    if (!is.na(K)) base <- paste0(base, "_K", K)
    base
}
# NOTE: this name is built from (method, f, Tt, N, burnin, K) only, not
# config_idx -- relies on every row of method_configs having a distinct
# (N, burnin, K) combination within the same method (true for all 13 rows
# above). If a future edit ever adds two config rows for the same method
# with identical (N, burnin, K) -- e.g. to vary something not reflected in
# the file name -- their .rds/.pdf paths (and checkpoint-skip logic) would
# collide silently.

# ---- walltime_hours_for(): per-category walltime bucket, same three-tier
# design as simulation.R's walltime_hours_for() ----
#
# CRITICAL, confirmed by Euler support via email: walltime <= 6h routes to
# "fila paralela curta", which had only 2 nodes allocated at the time of
# that reply (max ~4 concurrent place=excl jobs, and in practice as few as
# 2 were observed, likely because the other node was already occupied by
# someone else). walltime > 6h routes to "fila paralela longa" instead (14
# free nodes at the time of that reply), where this Category-2 account's
# 160-ncpus quota is the actual binding constraint -- 160/20 = 8 concurrent
# place=excl jobs, which is the concurrency the whole chunking/CHUNK_SIZE_
# CALIB design here is built around. Every value below is deliberately
# > 6h for this reason -- NOT because any single chunk is expected to take
# that long (see the per-category cost breakdown further below, largely
# unchanged from before this was discovered), but because requesting
# walltime <= 6h, even for the fastest chunks, would silently confine ALL
# of them to the short queue's ~2-4-node ceiling regardless of how quickly
# they actually finish. There is no downside to requesting more walltime
# than a job needs (it releases its node as soon as it's done either way)
# -- the downside only runs the other direction, requesting too little.
#
# Per-category cost model (unchanged from before, just re-based above the
# 6h threshold instead of tightly around actual cost): "leve" covers
# amh_montoril's N<=110000 configs and the sir_*/pg_as-K<=50 configs (well
# under a minute even at Tt=1600, calibration's N_chains=3 run in parallel
# so wall time tracks a single chain's time, not 3x it); "medio" covers
# amh_montoril's N=220000 config and pg_as's K=100/200 or N=22000 configs
# (estimated single-digit minutes, extrapolating linearly from the
# K=100/N=11000 production timings); "pesado" is stan alone (~2.1h
# measured at Tt=1600). Same as calibration_run.R's method_configs
# category comment: re-bucket a row by hand if a probe run shows it
# landing in the wrong tier before submitting the full grid.
walltime_hours_for <- function(category) {
    switch(category,
        leve   = 6.5,
        medio  = 7,
        pesado = 8,
        stop(sprintf("Unknown category: %s", category))
    )
}

# ==========================================================================
# ---- Run control ----
# ==========================================================================

run_mode <- "cluster"   # "local" or "cluster"

# Grid subset: NULL runs everything; otherwise a function(grid) -> grid
# (row filter) applied to calibration_grid right after it is built below.
# Example: only stan, only Tt=1600, only sir_collapsed + quadratic, etc.
grid_subset <- NULL
# grid_subset <- function(grid) subset(grid, method == "stan")

# Only used when run_mode == "local": which row of the (filtered) grid to
# run. A single integer runs just that row. "all" runs every row of the
# (filtered) grid, one after another, in this same R session -- useful for
# running a small grid_subset locally without going through batchtools/PBS
# at all (e.g. a handful of cheap tasks), or for a first end-to-end smoke
# test of the whole grid before committing to run_mode == "cluster".
task_id <- "all"

# ==========================================================================

# ---- Output paths ----
path_results <- "../results/calibration"
path_plots   <- file.path(path_results, "plots")
dir.create(path_results, showWarnings = FALSE, recursive = TRUE)
dir.create(path_plots,   showWarnings = FALSE, recursive = TRUE)


# ---- Per-method configuration(s), as decided during calibration ----
#
# A TABLE, not a single config per method: a method can appear more than
# once, with different (N, burnin, K), when more than one configuration is
# being compared for it -- e.g. amh_montoril here has three rows, to
# document the R_hat plateau found when increasing N/burnin (55000/5000 ->
# 75000/5000 -> 110000/10000) at Tt=1600, where R_hat(theta1_max) stayed
# structurally elevated (Fisher information argument -- see calibration
# discussion) rather than genuinely converging with more iterations.
#
# N / burnin / K fixed across the ENTIRE Tt grid for a given config row --
# deliberately NOT re-tuned per Tt, so that ESS/time comparisons across Tt
# stay apples-to-apples (amh_montoril and pg_as both show degradation at
# Tt=1600 that is treated as a genuine finding, not "fixed" by inflating N
# further within a single config row).
#
# category (leve/medio/pesado) drives which walltime bucket a task's PBS
# job is submitted under (see walltime_hours_for() and the "cluster"
# dispatch branch below) -- same three-tier design as simulation.R's
# R_config, adapted here per CONFIG ROW rather than per method, since a
# single method's wall-clock cost now varies substantially across its own
# rows (e.g. amh_montoril's N=55000 vs N=220000 configs, or pg_as's K=50
# vs K=200 configs) -- a per-method category would force the lightest and
# heaviest config of the same method into the same walltime bucket.
#
# Thresholds below are rough, deliberately conservative buckets based on
# the production/calibration timings already measured for these samplers
# at Tt=1600 (the worst case in the grid): amh_montoril and the sir_*
# methods are cheap regardless of N in the ranges used here (well under a
# minute even at the largest N); pg_as's cost scales with K*N and crosses
# into noticeably-longer territory at K=200; stan is in a class of its own
# (~2h at Tt=1600 already measured). Re-bucket a row by hand if a probe run
# (see "Sanity check" in the walltime discussion) shows it landing in the
# wrong tier.
method_configs <- rbind(
    data.frame(method = "amh_montoril",  N = 55000,  burnin = 5000,  K = NA, category = "leve"),
    data.frame(method = "amh_montoril",  N = 75000,  burnin = 5000,  K = NA, category = "leve"),
    data.frame(method = "amh_montoril",  N = 110000, burnin = 10000, K = NA, category = "leve"),
    data.frame(method = "amh_montoril",  N = 220000, burnin = 20000, K = NA, category = "medio"),
    data.frame(method = "pg_as",         N = 11000,  burnin = 1000,  K = 50,  category = "leve"),
    data.frame(method = "pg_as",         N = 11000,  burnin = 1000,  K = 100, category = "medio"),
    data.frame(method = "pg_as",         N = 11000,  burnin = 1000,  K = 200, category = "medio"),
    data.frame(method = "pg_as",         N = 22000,  burnin = 2000,  K = 50,  category = "medio"),
    data.frame(method = "pg_as",         N = 22000,  burnin = 2000,  K = 100, category = "medio"),
    data.frame(method = "pg_as",         N = 22000,  burnin = 2000,  K = 200, category = "pesado"),
    data.frame(method = "sir_laplace",   N = 11000,  burnin = 1000,  K = NA, category = "leve"),
    data.frame(method = "sir_collapsed", N = 11000,  burnin = 1000,  K = NA, category = "leve"),
    data.frame(method = "stan",          N = 11000,  burnin = 1000,  K = NA, category = "pesado")
)

# config_idx: 1st, 2nd, 3rd... config row for a given method, in the order
# written above -- used by the seed formula below to distinguish rows that
# share (method, f, Tt) but differ in (N, burnin, K), so that e.g.
# amh_montoril's N=55000 and N=220000 configs at the same (f, Tt) get
# different (and reproducible) seeds instead of colliding.
method_configs$config_idx <- ave(seq_len(nrow(method_configs)),
                                  method_configs$method, FUN = seq_along)


# ---- Calibration grid ----
#
# Full calibration phase: cross-join of method_configs (13 rows above)
# with (f, Tt) (2 x 2 = 4 combinations), giving 52 tasks, plus 12 more
# from the pg_as x sinusoidal extension below (6 pg_as configs x 2 Tt
# values) -- 64 tasks total. method_grid/Tt_grid/
# function_grid below are the FULL reference grids (used for the
# deterministic seed formula) and stay complete regardless of grid_subset
# -- the seed formula must be able to place ANY (method, f, Tt) at its
# correct index, not just the ones in the current subset.
method_grid   <- c("amh_montoril", "pg_as", "sir_laplace", "sir_collapsed", "stan")
Tt_grid       <- c(200, 400, 800, 1600)          # full reference grid (for seed encoding)
function_grid <- c("constant", "linear", "quadratic", "sinusoidal") # full reference grid

f_Tt_grid <- expand.grid(
    f  = c("constant", "quadratic"),
    Tt = c(200, 1600),
    stringsAsFactors = FALSE
)

calibration_grid <- merge(method_configs, f_Tt_grid, by = NULL)

# pg_as-only extension: sinusoidal has no changepoints, unlike the
# constant/quadratic pair above -- the cleanest scenario for separating
# "K too small" (degeneracy accumulates regardless of local dynamics) from
# "degeneracy near changepoints" (the original calibration hypothesis) as
# the explanation for pg_as's collapse at large Tt documented in the
# simulation results. Tt in {200, 1600} only (not the full 4-value grid),
# matching the same "worst/best case" spirit as f_Tt_grid above.
pg_as_sinusoidal_grid <- merge(
    method_configs[method_configs$method == "pg_as", ],
    expand.grid(f = "sinusoidal", Tt = c(200, 1600), stringsAsFactors = FALSE),
    by = NULL
)
calibration_grid <- rbind(calibration_grid, pg_as_sinusoidal_grid)

calibration_grid <- calibration_grid[order(calibration_grid$method, calibration_grid$f,
                                            calibration_grid$Tt, calibration_grid$N), ]
rownames(calibration_grid) <- NULL

replica <- 1  # fixed for the whole calibration phase

# ---- Grid subset (see "Run control" above) ----
if (!is.null(grid_subset)) {
    calibration_grid <- grid_subset(calibration_grid)
    rownames(calibration_grid) <- NULL
}

n_tasks <- nrow(calibration_grid)
printf("Calibration grid: %d tasks%s", n_tasks,
       if (!is.null(grid_subset)) " (subset applied)" else "")

# task_id: stable row identifier, used below to map flattened chain-units
# back to the task they belong to (calibration_grid's own row order is not
# safe to rely on for this once grid_subset has possibly reordered/dropped
# rows).
calibration_grid$task_id <- seq_len(n_tasks)

# replica: same value for the whole grid ("fixed for the whole calibration
# phase", per `replica <- 1` above) -- carried here as an explicit column
# (not read from the global `replica` inside cluster workers) so
# chain_grid/units below are fully self-contained data, not dependent on a
# clusterExport()'d snapshot staying in sync with it.
calibration_grid$replica <- replica


N_chains   <- 3
plot_chain <- 1

# ---- Chunking (cluster mode): flatten calibration_grid into ONE ROW PER
# (task, chain_id), then group into chunks of up to CHUNK_SIZE_CALIB
# chain-units each, for the flattened/local-cluster dispatch below (see
# run_chain_units()'s header comment for why chain-level, not task-level,
# is the unit that actually fixes node under-utilization). ----
#
# CHUNK_SIZE_CALIB=18 (not 20): leaves 2 of the node's 20 physical cores
# free for OS/NFS/scheduler overhead on a place=excl allocation, rather
# than saturating every core -- the same margin already used informally
# elsewhere in the project. TASKS_PER_CHUNK = 18 %/% 3 = 6 tasks/chunk when
# N_chains=3, matching the "6 tasks, 18 cores" figure directly -- but
# expressed as a chain-count division so it stays correct if N_chains ever
# changes without anyone needing to re-derive this by hand.
#
# Chunks are assigned WITHIN (category, Tt) groups, same reasoning as
# simulation.R's chunking: cost varies enormously across category (leve/
# medio/pesado) and across Tt within a method (e.g. stan's cost ratio
# between Tt=200 and Tt=1600), so mixing them in one chunk would leave that
# chunk's workers idle waiting for its slowest chain-unit once the fast
# ones finish -- exactly the sync-barrier problem the flattening is meant
# to avoid, just moved up one level (across tasks in a chunk) instead of
# solved. A task's own N_chains chain-units are NEVER split across chunks
# (chunk_index is computed on TASKS within a (category, Tt) group, not on
# raw chain-units), so aggregate_task() can always assume every chain file
# it needs was produced by the same chunk's job.
CHUNK_SIZE_CALIB <- 18L
TASKS_PER_CHUNK  <- CHUNK_SIZE_CALIB %/% N_chains

calibration_grid$chunk_index <- ave(seq_len(n_tasks),
                                     paste(calibration_grid$category, calibration_grid$Tt),
                                     FUN = function(idx) ceiling(seq_along(idx) / TASKS_PER_CHUNK))
calibration_grid$chunk_id <- sprintf("%s_Tt%04d_%03d",
                                      calibration_grid$category, calibration_grid$Tt,
                                      calibration_grid$chunk_index)

# ---- chain_grid: one row per (task, chain_id) -- the actual flat unit of
# parallel work dispatched by run_chain_units() below. ----
chain_grid <- do.call(rbind, lapply(1:N_chains, function(k) {
    g <- calibration_grid
    g$chain_id <- k
    g
}))
chain_grid <- chain_grid[order(chain_grid$chunk_id, chain_grid$task_id, chain_grid$chain_id), ]
rownames(chain_grid) <- NULL

# ---- Job-level grid (one row per chunk/PBS job) ----
job_grid <- unique(calibration_grid[, c("chunk_id", "category", "Tt")])
job_grid <- job_grid[order(job_grid$category, job_grid$Tt, job_grid$chunk_id), ]
rownames(job_grid) <- NULL
n_chunks <- nrow(job_grid)
printf("Calibration grid: %d chunk(s) (%d task(s)/chunk max, %d chain-unit(s)/chunk max)",
       n_chunks, TASKS_PER_CHUNK, CHUNK_SIZE_CALIB)
for (cat in unique(job_grid$category)) {
    printf("  %s: %d chunk(s)", cat, sum(job_grid$category == cat))
}

# Prior hyperparameters (identical across all methods and the whole grid)
mu_01     <- 0
sigma2_01 <- 100
mu_02     <- 0
sigma2_02 <- 100
nu_01  <- 2
eta_01 <- 0.01
nu_02  <- 2
eta_02 <- 0.0001

# amh_montoril adaptive-MH target acceptance
ac_ref <- 0.44

# sir_laplace / sir_collapsed
M_is       <- 3
M_irls_max <- 20
tol        <- 1e-4
R_prerun   <- 3000  # sir_collapsed only
M_sir      <- 3     # sir_collapsed only

verbose      <- FALSE
print_every  <- 1000
plots        <- TRUE
compute_rhat <- TRUE
compute_ess  <- TRUE


# ---- Physical core count (same logic as simulation.R's resolve_n_cores(),
# now genuinely needed here too: run_chain_units() calls this with an
# explicit n_cores=CHUNK_SIZE_CALIB from cluster mode, unlike the old code
# which only ever called it with n_cores=NULL -- capping by n_needed even
# when n_cores is explicit avoids opening more workers than there are
# chain-units to process for an incomplete trailing chunk.) ----
resolve_n_cores <- function(n_cores, n_needed) {
    if (!is.null(n_cores)) return(min(n_cores, n_needed))

    phys <- NA_integer_
    if (Sys.info()[["sysname"]] == "Linux" && nzchar(Sys.which("lscpu"))) {
        out <- tryCatch(
            system("lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#'", intern = TRUE),
            error = function(e) character(0)
        )
        if (length(out) > 0) phys <- length(unique(out))
    }
    if (is.na(phys)) phys <- parallel::detectCores(logical = FALSE)
    if (is.na(phys)) {
        logi <- parallel::detectCores(logical = TRUE)
        phys <- max(1, logi %/% 2)
        printf("Physical core count unavailable from OS; falling back to floor(logical/2) = %d", phys)
    }
    min(phys, n_needed)
}


# ---- build_summary_row(): one-row summary from a diag object ----
#
# Every field comes directly from print_and_plot_diagnostics()'s return
# value (diag) -- no recomputation. Per-chain vectors (length N_chains) are
# flattened into one column per chain (e.g. ess_w1_chain1, ess_w1_chain2,
# ess_w1_chain3) via flatten_chains(), so the whole summary stays a single
# data.frame row, matching every line print_and_plot_diagnostics() prints
# to the console for a task.
#
# Factored out of run_calibration_task() so calibration_aggregate.R can
# call it too, on a diag recomputed from each task's saved .rds (via
# print_and_plot_diagnostics(..., plots = FALSE)), without duplicating this
# construction.
flatten_chains <- function(x, name) {
    if (is.null(x)) {
        out <- list(NA_real_)
        names(out) <- name
        return(out)
    }
    out <- as.list(x)
    names(out) <- paste0(name, "_chain", seq_along(x))
    out
}

# nz(): NULL -> NA (scalar). data.frame() turns a NULL list element into a
# 0-row column, not NA, which breaks row-binding whenever a field is NULL
# for the current method (e.g. ac_ratio_mean is NULL for anything other
# than amh_montoril). Every scalar diag$* field that can be NULL must be
# wrapped in nz() before entering the data.frame() call below.
nz <- function(x) if (is.null(x)) NA_real_ else x

build_summary_row <- function(diag, method, f, Tt, replica, N, burnin, K, elapsed_time) {
    data.frame(
        c(
            list(
                method = method, f = f, Tt = Tt, replica = replica,
                N = N, burnin = burnin, K = K,
                elapsed_time_s = elapsed_time,

                W1_mean = nz(diag$W1_mean), W1_median = nz(diag$W1_median),
                W2_mean = nz(diag$W2_mean), W2_median = nz(diag$W2_median),
                loglik = nz(diag$loglik),

                rhat_theta01 = nz(diag$rhat_theta01),
                rhat_theta02 = nz(diag$rhat_theta02),
                rhat_w1 = nz(diag$rhat_w1),
                rhat_w2 = nz(diag$rhat_w2),
                rhat_theta1_mean = nz(diag$rhat_theta1_mean),
                rhat_theta1_max  = nz(diag$rhat_theta1_max),
                rhat_theta2_mean = nz(diag$rhat_theta2_mean),
                rhat_theta2_max  = nz(diag$rhat_theta2_max),
                rhat_max = nz(diag$rhat_max)
            ),
            flatten_chains(diag$ess_theta01, "ess_theta01"),
            flatten_chains(diag$ess_theta01_tail, "ess_theta01_tail"),
            list(ess_theta01_cv = nz(diag$ess_theta01_cv)),
            flatten_chains(diag$ess_theta02, "ess_theta02"),
            flatten_chains(diag$ess_theta02_tail, "ess_theta02_tail"),
            list(ess_theta02_cv = nz(diag$ess_theta02_cv)),
            flatten_chains(diag$ess_w1, "ess_w1"),
            flatten_chains(diag$ess_w1_tail, "ess_w1_tail"),
            list(ess_w1_cv = nz(diag$ess_w1_cv)),
            flatten_chains(diag$ess_w2, "ess_w2"),
            flatten_chains(diag$ess_w2_tail, "ess_w2_tail"),
            list(ess_w2_cv = nz(diag$ess_w2_cv)),

            flatten_chains(diag$ess_theta1_mean_over_t, "ess_theta1_mean_over_t"),
            flatten_chains(diag$ess_theta1_mean_over_t_tail, "ess_theta1_mean_over_t_tail"),
            list(ess_theta1_mean_over_t_cv = nz(diag$ess_theta1_mean_over_t_cv)),
            flatten_chains(diag$ess_theta1_min_over_t, "ess_theta1_min_over_t"),
            flatten_chains(diag$ess_theta1_min_over_t_tail, "ess_theta1_min_over_t_tail"),
            list(ess_theta1_min_over_t_cv = nz(diag$ess_theta1_min_over_t_cv)),

            flatten_chains(diag$ess_theta2_mean_over_t, "ess_theta2_mean_over_t"),
            flatten_chains(diag$ess_theta2_mean_over_t_tail, "ess_theta2_mean_over_t_tail"),
            list(ess_theta2_mean_over_t_cv = nz(diag$ess_theta2_mean_over_t_cv)),
            flatten_chains(diag$ess_theta2_min_over_t, "ess_theta2_min_over_t"),
            flatten_chains(diag$ess_theta2_min_over_t_tail, "ess_theta2_min_over_t_tail"),
            list(ess_theta2_min_over_t_cv = nz(diag$ess_theta2_min_over_t_cv)),

            list(
                ess_bulk_min = nz(diag$ess_bulk_min),

                ess_sec_w1_bulk = nz(diag$ess_sec_w1_bulk),
                ess_sec_w1_tail = nz(diag$ess_sec_w1_tail),
                ess_sec_w2_bulk = nz(diag$ess_sec_w2_bulk),
                ess_sec_w2_tail = nz(diag$ess_sec_w2_tail),
                ess_sec_theta1_mean_bulk = nz(diag$ess_sec_theta1_mean_bulk),
                ess_sec_theta1_mean_tail = nz(diag$ess_sec_theta1_mean_tail),
                ess_sec_theta1_min_bulk  = nz(diag$ess_sec_theta1_min_bulk),
                ess_sec_theta1_min_tail  = nz(diag$ess_sec_theta1_min_tail),
                ess_sec_theta2_mean_bulk = nz(diag$ess_sec_theta2_mean_bulk),
                ess_sec_theta2_mean_tail = nz(diag$ess_sec_theta2_mean_tail),
                ess_sec_theta2_min_bulk  = nz(diag$ess_sec_theta2_min_bulk),
                ess_sec_theta2_min_tail  = nz(diag$ess_sec_theta2_min_tail),

                ac_ratio_mean = nz(diag$ac_ratio_mean),
                ac_ratio_at_changepoints = nz(diag$ac_ratio_at_changepoints),
                w1_mh_acceptance_rate = nz(diag$w1_mh_acceptance_rate),
                w2_mh_acceptance_rate = nz(diag$w2_mh_acceptance_rate),
                ce1_shape = nz(diag$ce1_shape), ce1_rate = nz(diag$ce1_rate),
                ce2_shape = nz(diag$ce2_shape), ce2_rate = nz(diag$ce2_rate)
            )
        ),
        stringsAsFactors = FALSE
    )
}


# ---- Chain initialization (identical scheme to test_cpp.R / test_stan.R) ----
make_chain_inits <- function(chain_id, seed, y, Tt, method) {
    theta1_ref   <- log(y + 0.5)
    theta2_ref   <- c(diff(theta1_ref), 0)
    theta_01_ref <- theta1_ref[1]
    theta_02_ref <- theta2_ref[1]
    W1_ref <- var(diff(theta1_ref))
    W2_ref <- var(diff(theta2_ref))

    set.seed(seed)

    theta_01_init <- theta_01_ref + rnorm(1, 0, sd = 2 * sqrt(abs(theta_01_ref) + 1))
    theta_02_init <- theta_02_ref + rnorm(1, 0, sd = 2 * sqrt(abs(theta_02_ref) + 1))
    W1_init <- W1_ref * exp(rnorm(1, 0, sd = 1))
    W2_init <- W2_ref * exp(rnorm(1, 0, sd = 1))

    theta1_init <- theta1_ref + rnorm(Tt, 0, sd = 2 * sqrt(W1_ref))
    theta2_init <- theta2_ref + rnorm(Tt, 0, sd = 2 * sqrt(W2_ref))

    list(chain_id = chain_id, seed = seed,
         theta_01 = theta_01_init, theta_02 = theta_02_init,
         W1 = W1_init, W2 = W2_init,
         theta1 = theta1_init, theta2 = theta2_init,
         theta1_tilde = rep(0, Tt),
         varsigma2 = rep(0.02, Tt))
}


# ---- One-chain runner: unified across ALL 5 methods, including stan ----
#
# Previously stan was a special case inside run_calibration_task() (its own
# N_chains-way rstan::sampling(cores=) call, run OUTSIDE the
# parallel::makeCluster() block used by the other 4 methods) -- this made
# stan the one method that could not be decomposed into independent,
# individually-checkpointable (task, chain_id) units, which is exactly what
# the flattened chunking below needs from every method uniformly. Calling
# sample_stan() with N_chains=1, n_cores=1 (one call per chain, like
# simulation.R already does for production) removes that asymmetry: stan's
# per-chain cost is now paid inside a single-core worker exactly like
# amh_montoril/pg_as/sir_laplace/sir_collapsed, and its 3 calibration
# chains are combined afterwards by aggregate_task() exactly the same way
# regardless of method. The stan model object is loaded/compiled by the
# CALLER (run_one_chain_unit(), once per worker process) and passed in here
# -- loading/compiling it fresh inside every chain call would be wasteful
# and (worse) would race on ../cache/poisson_ltdm.rds if two stan chain
# units happened to run concurrently on the same node with no cache file
# yet.
run_one_chain <- function(init, method, y, N, burnin, K, stan_model = NULL) {
    printf("Chain %d (seed=%d), method=%s", init$chain_id, init$seed, method)
    set.seed(init$seed)

    if (method == "stan") {
        # burnin passed through UNCHANGED to sample_stan(), which maps it
        # straight to rstan::sampling(warmup = burnin) -- this is REAL HMC
        # warmup (step-size + mass-matrix adaptation), not a value stan
        # discards afterwards. sampler_stan.R calls rstan::extract(...,
        # inc_warmup = TRUE), so the returned theta1_hist/theta2_hist/etc.
        # already include those burnin draws, in N total rows -- exactly
        # the same shape the 4 Rcpp samplers return (N total draws, no
        # internal burnin trim of their own). aggregate_task()'s call to
        # print_and_plot_diagnostics(..., burnin = burnin) is what trims
        # the first `burnin` rows for every method uniformly, stan
        # included -- burnin must NOT be zeroed here, or rstan performs no
        # adaptation at all for the chain.
        sample_stan(
            model = stan_model, y = y, N = N, burnin = burnin, seed = init$seed,
            N_chains = 1, n_cores = 1,
            chain_inits = list(init),
            mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
            nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
            verbose = verbose)$results[[1]]

    } else if (method == "pg_as") {
        pg_as_cpp(y = y, K = K, N = N,
                  mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                  nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                  theta1 = init$theta1, theta2 = init$theta2,
                  theta_01 = init$theta_01, theta_02 = init$theta_02,
                  W1 = init$W1, W2 = init$W2, verbose = verbose, print_every = print_every)
    } else if (method == "amh_montoril") {
        amh_montoril_cpp(y = y, N = N,
                          mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                          nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                          theta1 = init$theta1, theta2 = init$theta2,
                          theta_01 = init$theta_01, theta_02 = init$theta_02,
                          W1 = init$W1, W2 = init$W2, ac_ref = ac_ref, varsigma2 = init$varsigma2,
                          verbose = verbose, print_every = print_every)
    } else if (method == "sir_laplace") {
        sir_laplace_cpp(y = y, N = N,
                         mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                         nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                         theta1 = init$theta1, theta2 = init$theta2,
                         theta_01 = init$theta_01, theta_02 = init$theta_02,
                         W1 = init$W1, W2 = init$W2, M_irls_max = M_irls_max, M_is = M_is, tol = tol,
                         verbose = verbose, print_every = print_every)
    } else if (method == "sir_collapsed") {
        sir_collapsed_cpp(y = y, R_prerun = R_prerun, N = N,
                           mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
                           nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
                           theta1 = init$theta1, theta2 = init$theta2, theta1_tilde = init$theta1_tilde,
                           theta_01 = init$theta_01, theta_02 = init$theta_02, W1 = init$W1, W2 = init$W2,
                           M_is_lik = M_is, M_sir_theta1 = M_sir, M_irls_max = M_irls_max, tol = tol,
                           verbose = verbose, print_every = print_every)
    } else {
        stop(sprintf("Unknown method: '%s'.", method))
    }
}


# ---- Chain-level seed: shared by run_one_chain_unit() and any code that
# needs to know a chain's seed without actually running it (e.g. tests) ----
#
# Same formula run_calibration_task() used before, now factored out since
# both run_one_chain_unit() (runs a chain) and aggregate_task() (loads a
# task's already-run chains back from disk, no need to recompute the seed,
# but needs seed_base for the final task .rds's seed_base field) need it.
task_seed_base <- function(method, f, Tt, config_idx, replica) {
    method_idx <- match(method, method_grid)
    Tt_idx     <- match(Tt, Tt_grid)
    f_idx      <- match(f, function_grid)
    method_idx * 1e7 + f_idx * 1e6 + Tt_idx * 1e5 + config_idx * 1e3 + replica * 10
}


# ---- chain_result_filename(): path for one (task, chain_id) checkpoint ----
#
# Lives under path_chains (../results/calibration/chains_tmp/), separate from
# path_results (the final per-TASK .rds, unchanged) -- these are two
# different kinds of checkpoint at two different granularities (chain vs.
# task), and aggregate_task() deletes most of the chain-level ones once the
# task-level one exists (see there), so keeping them in their own
# subdirectory makes that cleanup a simple glob instead of having to
# pattern-match filenames apart from the task-level .rds files living in
# the same directory.
path_chains <- file.path(path_results, "chains_tmp")
dir.create(path_chains, showWarnings = FALSE, recursive = TRUE)

chain_result_filename <- function(method, f, Tt, N, burnin, K, chain_id) {
    task_name <- task_name_for(method, f, Tt, N, burnin, K)
    file.path(path_chains, sprintf("%s_chain%d.rds", task_name, chain_id))
}


# ---- run_one_chain_unit(): the flat, individually-checkpointed unit of
# parallel work -- one (task, chain_id) pair ----
#
# This is what gets distributed across the ~18 workers of a chunk's local
# cluster (see run_chain_units() below) -- NOT run_calibration_task() /
# whole tasks, which is what made the OLD design block a full node (20
# cores) down to only N_chains=3 busy workers. Checkpointed exactly like
# run_calibration_task() used to be at the task level, just one level
# finer-grained: if this chain's .rds already exists, load and return it
# instead of re-running -- lets a chunk resume after a job dies partway
# through without redoing chains that already finished, and lets
# aggregate_task() (below) simply assume every chain file it needs is
# already on disk once run_chain_units() returns.
#
# stan_model is NULL for every method except stan; it is loaded/compiled
# ONCE per worker process by run_chain_units()'s clusterEvalQ() (or, in
# "local" single-task mode, once by run_calibration_task() itself), not
# once per chain -- see run_one_chain()'s header comment for why.
run_one_chain_unit <- function(method, f, Tt, N, burnin, K, config_idx, replica, chain_id, stan_model = NULL) {
    out_file <- chain_result_filename(method, f, Tt, N, burnin, K, chain_id)

    if (file.exists(out_file)) {
        return(readRDS(out_file))
    }

    source_name <- sprintf("%s_%s_%s", f, Tt, replica)
    data <- readRDS(file.path("..", "data", "simulated", paste0(source_name, ".rds")))
    y <- data$y

    seed_base <- task_seed_base(method, f, Tt, config_idx, replica)
    seed <- seed_base + (chain_id - 1)

    init <- make_chain_inits(chain_id, seed, y, Tt, method)

    start_time <- Sys.time()
    hist <- run_one_chain(init, method, y, N, burnin, K, stan_model = stan_model)
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")

    chain_result <- list(chain_id = chain_id, seed = seed, elapsed_time = elapsed_time, hist = hist)
    saveRDS(chain_result, file = out_file)
    chain_result
}


# ---- aggregate_task(): combines a task's N_chains already-run chain
# checkpoints into the final task-level .rds + diagnostic .pdf + summary
# row -- unchanged output contract from the old run_calibration_task(),
# just fed from chain_result_filename() checkpoints instead of an in-memory
# `results` list built by a single N_chains-way parallel block ----
#
# Callable only once every one of this task's N_chains chain files exists
# on disk (run_chain_units()/run_one_chain_unit() above are what create
# them) -- callers are responsible for that ordering (both
# run_calibration_task() below and the "cluster" dispatch's per-chunk
# aggregation loop guarantee it).
aggregate_task <- function(method, f, Tt, N, burnin, K, config_idx, replica = 1) {

    task_name <- task_name_for(method, f, Tt, N, burnin, K)
    task_rds  <- file.path(path_results, paste0(task_name, ".rds"))
    plot_file <- file.path(path_plots, paste0(task_name, ".pdf"))
    printf("==== Aggregating task: %s ====", task_name)

    # Same checkpoint spirit as before: a task already aggregated is
    # skipped entirely -- see run_calibration_task()'s old header comment
    # (unchanged rationale, just moved here since aggregation, not chain
    # execution, is now the step that produces the task-level artifacts).
    if (file.exists(task_rds)) {
        printf("Task already aggregated (%s exists) -- skipping entirely.", task_rds)
        return(NULL)
    }

    source_name <- sprintf("%s_%s_%s", f, Tt, replica)
    data <- readRDS(file.path("..", "data", "simulated", paste0(source_name, ".rds")))
    y <- data$y
    theta1_true <- data$theta
    changepoints <- c(0.25, 0.5, 0.75) * Tt

    if (Tt == 200)  t_obs <- c(50, 100, 150, 175)
    if (Tt == 400)  t_obs <- c(75, 100, 200, 300)
    if (Tt == 800)  t_obs <- c(100, 300, 500, 700)
    if (Tt == 1600) t_obs <- c(400, 800, 1200, 1600)

    seed_base <- task_seed_base(method, f, Tt, config_idx, replica)

    chain_files <- vapply(1:N_chains, function(k) chain_result_filename(method, f, Tt, N, burnin, K, k),
                           character(1))
    missing <- !file.exists(chain_files)
    if (any(missing)) {
        stop(sprintf(
            "aggregate_task(%s): %d/%d chain checkpoint(s) missing -- run_chain_units() ",
            task_name, sum(missing), N_chains),
            "must complete every chain of a task before aggregate_task() is called for it.")
    }
    chain_results <- lapply(chain_files, readRDS)

    # `results`: same shape print_and_plot_diagnostics()/build_summary_row()
    # already expect -- a plain list of N_chains history objects, in chain
    # order, exactly what the old N_chains-way parallel block produced
    # in-memory.
    results <- lapply(chain_results, `[[`, "hist")

    # elapsed_time: previously the WALL-CLOCK time of the whole N_chains-way
    # parallel block (chains genuinely running side by side, so this was
    # close to max(per-chain time), not their sum). Chains are now flat
    # units that may run at different points in time, possibly across
    # different chunk jobs on a resumed run -- max() over their individually
    # measured elapsed_time is the closest equivalent available and is what
    # ESS/second (a downstream diagnostic) is computed against, same as
    # before. NOT sum(): summing would understate ESS/second by conflating
    # "3 chains ran on 3 different cores" with "3 chains ran back-to-back on
    # 1 core", which is not what happened in either the old or new design.
    elapsed_time <- max(vapply(chain_results, `[[`, numeric(1), "elapsed_time"))

    if (plots) {
        pdf(plot_file, width = 8, height = 6)
        on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
    }

    diag <- print_and_plot_diagnostics(
        results, y, changepoints,
        theta1_true = theta1_true, theta2_true = NULL,
        t_obs = t_obs, burnin = burnin, elapsed_time = elapsed_time,
        plot_chain = plot_chain,
        nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
        ac_ref = ac_ref, plots = plots,
        compute_rhat = compute_rhat, compute_ess = compute_ess)

    if (plots) {
        dev.off()
        printf("Plots saved to: %s", plot_file)
    }

    saveRDS(list(results = results, elapsed_time = elapsed_time,
                 method = method, f = f, Tt = Tt, replica = replica,
                 N = N, burnin = burnin, K = K, seed_base = seed_base),
            file = task_rds)

    # ---- Clean up per-chain checkpoints, keeping only plot_chain's ----
    #
    # Every chain's raw history is already inside task_rds (results list
    # above) -- the per-chain files under path_chains/ were only ever a
    # resumability checkpoint for the flattened parallel step, not a
    # separate source of truth. Deleting all but plot_chain's avoids
    # doubling NFS usage for the heaviest configs (e.g. amh_montoril's
    # N=220000/Tt=1600 row: ~5.6GB per chain history -- keeping all 3
    # intermediate files alongside the already-heavy task_rds would nearly
    # double peak disk use for that row for no benefit). plot_chain's file
    # is kept deliberately: it is the one already re-read for the
    # diagnostic plot above, so keeping it lets that specific chain's raw
    # history be inspected/re-plotted later without reloading task_rds's
    # (larger, all-chains) results list.
    for (k in seq_len(N_chains)) {
        if (k != plot_chain) {
            f_chain <- chain_result_filename(method, f, Tt, N, burnin, K, k)
            if (file.exists(f_chain)) file.remove(f_chain)
        }
    }

    build_summary_row(diag, method, f, Tt, replica, N, burnin, K, elapsed_time)
}


# ---- run_chain_units(): runs a whole data.frame of (task, chain_id) rows
# in ONE flat, maximally-parallel local cluster ----
#
# This is the piece that actually delivers the "6 tasks per job, 18/20
# cores busy instead of 3/20" goal: instead of one parallel::makeCluster()
# PER TASK (N_chains=3 workers, opened and torn down once per task -- the
# old design, and the reason a whole place=excl node sat at 3/20 cores
# busy), ONE cluster is opened for the WHOLE chunk (up to CHUNK_SIZE_CALIB
# chain-units, e.g. 6 tasks x 3 chains = 18), and every chain-unit --
# regardless of which task it belongs to -- is pulled from one shared
# foreach queue. A worker that finishes task A's chain 2 early does not
# wait for task A's chains 1/3 to finish; it immediately picks up the next
# pending chain-unit in the queue, whatever task it belongs to. This is
# what actually removes the per-task synchronization barrier that a naive
# "pack 6 fixed tasks per job" (task-level, not chain-level, packing) would
# still have if those 6 tasks differ in cost -- mirrors simulation.R's
# run_tasks(), just with "one chain of one task" as the flat unit instead
# of "one whole task" (which is already atomic there, since production
# tasks are single-chain).
#
# stan_model is loaded/compiled ONCE here (not per chain, not per task) and
# exported to every worker via clusterExport() -- workers that never handle
# a stan chain-unit simply never touch it. is_stan: TRUE only if `units`
# contains at least one stan row, so the (potentially slow, first-time)
# model compilation is skipped entirely for chunks that never need it.
run_chain_units <- function(units, n_cores = NULL) {

    n <- nrow(units)
    n_cores_used <- resolve_n_cores(n_cores, n)
    printf("Running %d chain-unit(s) using %d core(s)", n, n_cores_used)

    is_stan <- any(units$method == "stan")
    stan_model <- NULL
    if (is_stan) {
        rstan::rstan_options(auto_write = FALSE)
        if (file.exists("../cache/poisson_ltdm.rds")) {
            stan_model <- readRDS("../cache/poisson_ltdm.rds")
        } else {
            printf("Building the Stan model")
            stan_model <- rstan::stan_model(
                file = "../PoissonLTDM/inst/stan/poisson_ltdm.stan",
                model_name = "PoissonLTDM")
            saveRDS(stan_model, file = "../cache/poisson_ltdm.rds")
        }
    }

    cl <- parallel::makeCluster(n_cores_used)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    `%dorng%` <- doRNG::`%dorng%`

    current_wd <- getwd()
    parallel::clusterExport(cl, "current_wd", envir = environment())
    if (is_stan) parallel::clusterExport(cl, "stan_model", envir = environment())

    # Same rationale as the old code's clusterExport() call: foreach's
    # auto-export only walks the LEXICAL environment of the %dorng% call,
    # not the globals run_one_chain()/run_one_chain_unit() actually use --
    # those must be exported explicitly from .GlobalEnv.
    parallel::clusterExport(cl,
        c("run_one_chain", "run_one_chain_unit", "chain_result_filename",
          "task_seed_base", "task_name_for", "make_chain_inits", "printf",
          "mu_01", "sigma2_01", "mu_02", "sigma2_02",
          "nu_01", "eta_01", "nu_02", "eta_02",
          "ac_ref", "M_is", "M_irls_max", "tol", "R_prerun", "M_sir",
          "verbose", "print_every", "N_chains", "path_chains",
          "method_grid", "Tt_grid", "function_grid"),
        envir = .GlobalEnv)

    parallel::clusterEvalQ(cl, {
        setwd(current_wd)
        pkgload::load_all("../PoissonLTDM", debug = FALSE)
        source("../PoissonLTDM/R/sampler_stan.R")
    })

    foreach::foreach(i = seq_len(n)) %dorng% {
        u <- units[i, ]
        run_one_chain_unit(u$method, u$f, u$Tt, u$N, u$burnin, u$K, u$config_idx, u$replica, u$chain_id,
                            stan_model = if (u$method == "stan") stan_model else NULL)
    }

    invisible(NULL)
}


# ---- Main entry point: runs ONE row of calibration_grid (all N_chains
# chains + aggregation) -- unchanged EXTERNAL call contract from before
# (same arguments, same return value), just internally a thin wrapper over
# run_chain_units() + aggregate_task() now. Used by "local" mode's
# task_id-based single-task runs (see dispatch block); the "cluster"
# dispatch below calls run_chain_units()/aggregate_task() directly instead,
# at chunk granularity, for the actual flattening benefit -- see
# run_chain_units()'s header comment for why that distinction matters. ----
run_calibration_task <- function(method, f, Tt, N, burnin, K, config_idx, replica = 1) {
    task_name <- task_name_for(method, f, Tt, N, burnin, K)
    task_rds  <- file.path(path_results, paste0(task_name, ".rds"))
    if (file.exists(task_rds)) {
        printf("Task already run (%s exists) -- skipping entirely.", task_name)
        return(NULL)
    }

    units <- data.frame(method = method, f = f, Tt = Tt, N = N, burnin = burnin, K = K,
                         config_idx = config_idx, replica = replica, chain_id = 1:N_chains,
                         stringsAsFactors = FALSE)
    run_chain_units(units, n_cores = NULL)

    aggregate_task(method, f, Tt, N, burnin, K, config_idx, replica)
}


# ---- run_and_save_task(): runs one grid row AND updates summary.csv ----
#
# ---- save_summary_row(): the summary.csv read-modify-write logic, on its
# own -- factored out of run_and_save_task() so callers that already have a
# summary_row in hand (e.g. the "local all" dispatch below, which calls
# aggregate_task() directly to avoid run_and_save_task()'s redundant
# run_chain_units() call -- see run_and_save_task()'s header comment) don't
# have to go through run_calibration_task() again just to reach this
# logic. Same behavior as before, unchanged. ----
save_summary_row <- function(summary_row, update_summary_csv = TRUE) {
    if (is.null(summary_row)) return(invisible(NULL))

    printf("Summary row: rhat_max=%.4f, ess_bulk_min=%.1f, elapsed_time=%.2fs",
           summary_row$rhat_max, summary_row$ess_bulk_min, summary_row$elapsed_time_s)

    if (update_summary_csv) {
        summary_file <- file.path(path_results, "summary.csv")
        if (file.exists(summary_file)) {
            summary_all <- read.csv(summary_file, stringsAsFactors = FALSE)
            # K may be NA (most methods) -- compare via a string key instead
            # of == directly, since NA == NA is NA in R, not TRUE, which
            # would silently fail to match/drop the old row for methods
            # without a K.
            existing_key <- with(summary_all, paste(method, f, Tt, N, burnin, K, replica))
            this_key <- paste(summary_row$method, summary_row$f, summary_row$Tt,
                               summary_row$N, summary_row$burnin, summary_row$K, summary_row$replica)
            keep <- existing_key != this_key
            summary_all <- rbind(summary_all[keep, , drop = FALSE], summary_row)
        } else {
            summary_all <- summary_row
        }
        write.csv(summary_all, summary_file, row.names = FALSE)
        printf("Summary written to: %s (%d task(s) total)", summary_file, nrow(summary_all))
    }

    summary_row
}


# ---- run_and_save_task(): runs one grid row's chains (via
# run_calibration_task(), which internally opens its own single-task
# cluster -- fine for LOCAL SINGLE-TASK use, but see the "local all"
# dispatch below for why it is NOT reused there) AND updates summary.csv.
# Used directly by run_mode == "local"'s single-task_id path. NOT used by
# run_mode == "cluster": each PBS job there calls run_chain_units()/
# aggregate_task() directly (see the chunk-level batchMap() body below)
# instead, since concurrent jobs writing to the same summary.csv at once
# would race/clobber each other, so cluster-mode jobs only produce their
# own .rds -- calibration_aggregate.R builds summary.csv from those .rds
# files afterwards, once all jobs are done.
#
# run_calibration_task() returns NULL for an already-checkpointed task
# (see its own header comment) -- there is nothing new to summarize in
# that case, so update_summary_csv is skipped and this function returns
# NULL too. That task's summary.csv row (if any) is left exactly as it
# was; calibration_aggregate.R is what fills in rows for checkpointed
# tasks, not this function.
run_and_save_task <- function(method, f, Tt, N, burnin, K, config_idx, replica = 1, update_summary_csv = TRUE) {
    summary_row <- run_calibration_task(method, f, Tt, N, burnin, K, config_idx, replica)
    save_summary_row(summary_row, update_summary_csv)
}


# ==========================================================================
# ---- Dispatch: run_mode == "local" or "cluster" (see "Run control") ----
# ==========================================================================

# ---- Guard against recursive dispatch inside a batchtools/PBS job, or
# when this file is sourced only for its definitions (e.g. by
# calibration_aggregate.R) ----
#
# makeRegistry(source = "calibration_run.R") makes every job re-source
# this WHOLE file (functions, config, AND the dispatch block below) before
# batchtools calls the batchMap() target function with that job's
# arguments. Without a guard, a job with run_mode == "cluster" would hit
# this same dispatch block and call submitJobs() again from inside itself.
# Sys.getenv("PBS_NODEFILE") is only non-empty when actually executing
# inside a running PBS job (same auto-detection already used in
# simulation.R) -- never true for the user's own interactive/Rscript
# invocation of this file, so it reliably distinguishes "being re-sourced
# inside a job" from "being run directly to submit/run a task".
#
# Separately, calibration_aggregate.R sources this file purely to reuse
# its function/config definitions (build_summary_row(), path_results,
# print_and_plot_diagnostics(), etc.) -- it must NOT trigger a real task
# run or a batchtools submission either, even though it runs outside any
# PBS job. It sets CALIBRATION_PHASE_DEFS_ONLY <- TRUE before sourcing this
# file for exactly that purpose; captured above (before rm(list=ls())) as
# .defs_only_flag, since that variable does not exist for a normal direct
# run of this script.
#
# A third, distinct re-entrancy case: makeRegistry() itself calls
# loadRegistryDependencies() -- which re-sources this file -- from INSIDE
# its own construction, before returning (see the comment right before the
# makeRegistry() call below). That happens in the very process that is
# submitting jobs, not inside a PBS job and not via
# calibration_aggregate.R, so neither of the two guards above catches it
# on its own. CALIBRATION_PHASE_MAKING_REGISTRY is set (via Sys.setenv(),
# so it survives sys.source()'s environment reset) immediately before that
# makeRegistry() call, for exactly this case.
running_inside_pbs_job <- Sys.getenv("PBS_NODEFILE") != ""
making_registry <- Sys.getenv("CALIBRATION_PHASE_MAKING_REGISTRY") != ""

if (!running_inside_pbs_job && !.defs_only_flag && !making_registry) {

if (run_mode == "local") {

    # ---- LOCAL: run task_id (a single row) or "all" rows of the (filtered) grid ----
    if (identical(task_id, "all")) {
        # Flatten across every chain-unit of the (filtered) grid and run
        # them ALL in one maximally-parallel local batch first -- mirrors
        # simulation.R's local mode ("every task... in ONE maximally-
        # parallel batch using every physical core"), adapted here to
        # chain-units (see run_chain_units()'s header comment). Aggregation
        # (per task, sequential, cheap) happens afterwards -- the
        # run_and_save_task() calls below will find every chain file
        # already checkpointed and skip straight to it.
        chain_grid_filtered <- chain_grid[chain_grid$task_id %in% calibration_grid$task_id, ]
        printf("Running %d chain-unit(s) locally (%d task(s)), in one parallel batch.",
               nrow(chain_grid_filtered), n_tasks)
        run_chain_units(chain_grid_filtered, n_cores = NULL)

        printf("Aggregating %d task(s), sequentially.", n_tasks)
        for (i in seq_len(n_tasks)) {
            task <- calibration_grid[i, ]
            printf("[%d/%d] method=%s, f=%s, Tt=%d, N=%d, burnin=%d, K=%s",
                   i, n_tasks, task$method, task$f, task$Tt, task$N, task$burnin,
                   if (is.na(task$K)) "NA" else task$K)
            # aggregate_task() directly (NOT run_and_save_task(), which
            # would call run_calibration_task() -> run_chain_units() again
            # -- pointless here since every chain-unit was already run and
            # checkpointed by the single batch call above; going through
            # run_and_save_task() would still work correctly (chain
            # checkpoints make it a fast no-op), but pays a full
            # makeCluster()/stopCluster() round-trip PER TASK just to reach
            # that no-op -- confirmed during testing to dominate this
            # loop's wall-clock time for a grid this size).
            summary_row <- aggregate_task(task$method, task$f, task$Tt, task$N, task$burnin, task$K,
                                           task$config_idx, task$replica)
            save_summary_row(summary_row)
        }
    } else {
        task <- calibration_grid[task_id, ]
        printf("Selected task %d/%d: method=%s, f=%s, Tt=%d, N=%d, burnin=%d, K=%s",
               task_id, n_tasks, task$method, task$f, task$Tt, task$N, task$burnin,
               if (is.na(task$K)) "NA" else task$K)

        summary_row <- run_and_save_task(task$method, task$f, task$Tt, task$N, task$burnin, task$K, task$config_idx, task$replica)
    }

} else if (run_mode == "cluster") {

    # ---- CLUSTER: submit every row of the (filtered) grid as one PBS job
    # each, via batchtools. Returns as soon as jobs are dispatched --
    # does NOT wait and does NOT write summary.csv (see
    # calibration_aggregate.R for that, run separately once jobs finish).
    library(batchtools)

    # ---- Guard against recursive dispatch, set BEFORE any
    # loadRegistry()/makeRegistry() call in this block ----
    #
    # Both loadRegistry() (used below to inspect a previous registry) and
    # makeRegistry() (used later to create the new one) call
    # loadRegistryDependencies() internally, which re-sources this whole
    # file via sys.source(fn, envir = .GlobalEnv). Without this guard set
    # BEFORE the very first such call, that re-source hits this same
    # dispatch block, which calls loadRegistry() again, which re-sources
    # again, ... -- infinite recursion (this was caught in testing: it
    # required Ctrl+C to stop, with a stack of hundreds of nested
    # loadRegistry() -> loadRegistryDependencies() -> sys.source() frames).
    # Setting the flag here, before the auto-clean section's
    # loadRegistry() call further down, closes that gap -- setting it only
    # right before the final makeRegistry() call (as in an earlier version
    # of this script) was NOT early enough.
    #
    # A Sys.setenv() var (not a plain R variable) is used because
    # sys.source(fn, envir = .GlobalEnv) resets the R environment on every
    # re-source but process-level env vars survive it -- same mechanism
    # already relied on for PBS_NODEFILE. on.exit() ensures the flag is
    # cleared even if this block errors out partway through (e.g. the
    # "jobs still queued/running" stop() below) -- otherwise it would stay
    # TRUE for the rest of the R session, silently suppressing this
    # dispatch block on any later run_and_save_task()-adjacent call that
    # happens to trigger a sys.source() of this file.
    Sys.setenv(CALIBRATION_PHASE_MAKING_REGISTRY = "TRUE")
    on.exit(Sys.unsetenv("CALIBRATION_PHASE_MAKING_REGISTRY"), add = TRUE)

    registry_dir <- "registry_calibration"  # fixed, hardcoded value -- see below for why it must be re-set after loadRegistry()

    # ---- Auto-clean a PREVIOUS, FINISHED registry ----
    #
    # makeRegistry() refuses to create a registry where one already exists
    # (the "File at path already exists" error seen during earlier
    # testing). Re-running this script after a previous cluster submission
    # has fully finished (every job either done or errored -- none still
    # queued/running) is the common case and should not require manually
    # rm -rf-ing the directory every time.
    #
    # This is intentionally NOT automatic when jobs are still
    # queued/running: archiving registry_calibration out from under jobs
    # that are still executing on the PBS side would not kill them, but it
    # would orphan them -- they keep running, but batchtools loses track
    # of them (no more getStatus(), no results collected). In that case,
    # stop with a clear message instead of silently touching anything, so
    # the person can decide (wait for the jobs, or explicitly archive/remove
    # the directory themselves if they really mean to abandon those jobs).
    #
    # A previous registry is ARCHIVED (moved into registry_archive/ with a
    # timestamp, via archive_registry() defined at the top of this file --
    # see there for why it must live there and not here), not deleted --
    # this preserves its logs and full job history (results/calibration/*.rds
    # and summary.csv are separate and were never affected either way)
    # while still freeing up the registry_calibration path for a fresh
    # makeRegistry() call.

    if (dir.exists(registry_dir)) {
        # NOTE: loadRegistry() below also re-sources this whole file (same
        # loadRegistryDependencies() mechanism as makeRegistry() -- see the
        # guard comment above), which resets .GlobalEnv via
        # rm(list = setdiff(ls(), ".defs_only_flag")) at the top of this
        # script. registry_dir (a plain R variable, set just above) does
        # NOT survive that reset the way env vars do -- so it must be
        # re-assigned immediately after any loadRegistry()/loadRegistryDependencies()
        # call, before it is used again. This was the cause of an "object
        # 'registry_dir' not found" error caught during testing.
        old_reg <- tryCatch(loadRegistry(registry_dir, writeable = FALSE),
                             error = function(e) NULL)
        registry_dir <- "registry_calibration"  # re-assign: see NOTE above

        if (is.null(old_reg)) {
            # Directory exists but isn't a valid/readable registry (e.g.
            # left over from an interrupted makeRegistry() call) -- safe
            # to archive, there is nothing coherent to lose track of.
            archive_registry(sprintf("Unreadable leftover registry directory: %s", registry_dir),
                              registry_dir)
        } else {
            n_pending <- nrow(batchtools::findNotDone(reg = old_reg))
            if (n_pending > 0) {
                stop(sprintf(
                    "Registry '%s' already exists with %d job(s) still queued/running. ",
                    registry_dir, n_pending),
                    "Refusing to archive it automatically -- wait for those jobs to finish ",
                    "(batchtools::getStatus(batchtools::loadRegistry(\"", registry_dir, "\"))), ",
                    "or archive/remove the directory yourself if you mean to abandon them.")
            }
            archive_registry(sprintf("Previous registry '%s' has no pending jobs", registry_dir),
                              registry_dir)
        }
    }

    reg <- makeRegistry(
        file.dir = registry_dir,
        source = "calibration_run.R",  # each job re-sources this whole
                                          # file (run_mode is irrelevant
                                          # inside the job -- batchMap()'s
                                          # target function is what runs,
                                          # not the dispatch block below)
        seed = 1
    )

    reg$cluster.functions <- makeClusterFunctionsTORQUE("calibration_pbs.tmpl")

    # ---- Skip chunks whose tasks are ALL already checkpointed ----
    #
    # aggregate_task() itself already skips a checkpointed task (see its
    # header comment) -- but only after a PBS job has been queued/started
    # for its whole chunk, i.e. after paying for a full place=excl node
    # allocation just to find nothing to do. Filtering here avoids creating
    # those jobs in the first place -- same spirit as the old task-level
    # filter, just checked per CHUNK (every task in the chunk must already
    # be done for the whole chunk to be skipped; a chunk with even one
    # unfinished task still needs its job, though run_chain_units() /
    # aggregate_task()'s own per-chain and per-task checkpoints mean that
    # job will skip straight past whatever already finished).
    task_names <- with(calibration_grid, mapply(task_name_for, method, f, Tt, N, burnin, K))
    task_rds_exists <- file.exists(file.path(path_results, paste0(task_names, ".rds")))
    chunk_done <- ave(task_rds_exists, calibration_grid$chunk_id, FUN = all)
    n_skipped_chunks <- sum(!duplicated(calibration_grid$chunk_id) & chunk_done)
    if (n_skipped_chunks > 0) {
        printf("Skipping %d already-fully-checkpointed chunk(s) (no job submitted for them).",
               n_skipped_chunks)
    }
    job_grid_to_submit <- job_grid[!(job_grid$chunk_id %in%
                                      unique(calibration_grid$chunk_id[chunk_done == TRUE])), ]

    if (nrow(job_grid_to_submit) == 0) {
        printf("Every chunk in the (filtered) grid is already checkpointed -- nothing to submit.")
        printf("Run calibration_aggregate.R to build/update summary.csv from the existing .rds files.")
    } else {

    # One job per CHUNK (not per task): each job runs run_chain_units() on
    # every chain-unit of its chunk in one flat, up-to-CHUNK_SIZE_CALIB-way
    # local cluster -- this is the actual fix for the old "3/20 cores busy"
    # problem, see run_chain_units()'s header comment -- then aggregates
    # every task in the chunk sequentially (cheap: just rhat/ESS + a plot
    # per task, no sampling). summary.csv is still NOT written here (same
    # reasoning as before: concurrent jobs racing on one file) --
    # calibration_aggregate.R remains the step that builds it from the
    # .rds files once every job is done.
    ids <- batchMap(
        fun = function(chunk_id) {
            chunk_chain_units <- chain_grid[chain_grid$chunk_id == chunk_id, ]
            run_chain_units(chunk_chain_units, n_cores = CHUNK_SIZE_CALIB)

            chunk_tasks <- calibration_grid[calibration_grid$chunk_id == chunk_id, ]
            for (i in seq_len(nrow(chunk_tasks))) {
                t <- chunk_tasks[i, ]
                aggregate_task(t$method, t$f, t$Tt, t$N, t$burnin, t$K, t$config_idx, t$replica)
            }
        },
        chunk_id = job_grid_to_submit$chunk_id,
        reg = reg
    )
    # batchMap()'s returned ids has only a job.id column (not the mapped
    # chunk_id/category) -- but job.id order matches the input vector's
    # order exactly (same behavior relied on in simulation.R), so category
    # can be attached directly by position, no join needed.
    ids$category <- job_grid_to_submit$category

    # Per-category walltime (leve/medio/pesado), mirroring simulation.R's
    # per-category submitJobs loop. ncpus=CHUNK_SIZE_CALIB (18), not
    # N_chains (3) -- this is the actual change that fixes node
    # under-utilization: the OLD code requested ncpus=N_chains=3 per job
    # while place=excl still reserved the whole 20-core node regardless,
    # leaving 17 cores idle for the job's entire walltime. Requesting 18
    # here doesn't change place=excl's node-level reservation either, but
    # it now matches what run_chain_units() actually uses inside the job,
    # which is the number that was wasted before.
    #
    # max.concurrent.jobs = floor(160 / CHUNK_SIZE_CALIB) = 8, NOT 10:
    # 10 * 18 = 180 ncpus, which exceeds this Category-2 account's 160-ncpus
    # quota -- letting batchtools push a 9th or 10th job to the scheduler
    # while 8 are already active can leave that job sitting in PBS's own
    # 'H' (Held) state (observed in practice: 7 Running + 2 Queued = 9
    # active chunks x 18 ncpus = 162, just over the 160 cap, right when an
    # 'H' job appeared in qstat). 8 is also the concurrency
    # walltime_hours_for()'s own header comment already assumes (the whole
    # chunking design is built around 8 concurrent place=excl jobs) --
    # max.concurrent.jobs=10 was never actually reachable in practice, it
    # was just a stale value copied from simulation.R's ORIGINAL (also
    # slightly-too-high) setting rather than derived from this file's own
    # CHUNK_SIZE_CALIB.
    max_concurrent_calib <- 160 %/% CHUNK_SIZE_CALIB
    for (cat in unique(ids$category)) {
        cat_ids <- ids[ids$category == cat, "job.id", drop = FALSE]
        submitJobs(cat_ids,
                   resources = list(ncpus = CHUNK_SIZE_CALIB, walltime_hours = walltime_hours_for(cat),
                                     max.concurrent.jobs = max_concurrent_calib),
                   reg = reg)
        printf("Submitted %d job(s) for category '%s' (walltime=%gh).",
               nrow(cat_ids), cat, walltime_hours_for(cat))
    }

    printf("Submitted %d job(s) total to the cluster (registry: %s).", nrow(ids), reg$file.dir)
    printf("Use check_progress_calibration.R to monitor.")
    printf("NOTE: submitJobs() blocks (in THIS R session) until every requested job has")
    printf("      been dispatched -- not until they finish. max.concurrent.jobs=%d is a", max_concurrent_calib)
    printf("      real, global throttle (checked against the live scheduler via")
    printf("      getBatchIds(), not a per-category counter), so submitting multiple")
    printf("      categories in one run is safe -- the true %d-job cap is respected", max_concurrent_calib)
    printf("      across all of them. But because this call can block for a long time")
    printf("      (hours, for the larger categories), run this inside tmux/screen or")
    printf("      with nohup, never directly in an SSH session that might disconnect.")
    printf("Once all jobs are done, run calibration_aggregate.R to build summary.csv.")

    }

} else {
    stop(sprintf("Unknown run_mode: '%s'. Use \"local\" or \"cluster\".", run_mode))
}

}  # end !running_inside_pbs_job && !.defs_only_flag && !making_registry guard
