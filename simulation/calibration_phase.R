# simulation/calibration_phase.R
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
task_id <- 1

# ==========================================================================

# ---- Output paths ----
path_results <- "results/calibration"
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
method_configs <- rbind(
    data.frame(method = "amh_montoril",  N = 55000,  burnin = 5000,  K = NA),
    data.frame(method = "amh_montoril",  N = 75000,  burnin = 5000,  K = NA),
    data.frame(method = "amh_montoril",  N = 110000, burnin = 10000, K = NA),
    data.frame(method = "pg_as",         N = 11000,  burnin = 1000,  K = 100),
    data.frame(method = "sir_laplace",   N = 11000,  burnin = 1000,  K = NA),
    data.frame(method = "sir_collapsed", N = 11000,  burnin = 1000,  K = NA),
    data.frame(method = "stan",          N = 11000,  burnin = 1000,  K = NA)
)


# ---- Calibration grid ----
#
# Full calibration phase: cross-join of method_configs (7 rows above) with
# (f, Tt) (2 x 2 = 4 combinations), giving 28 tasks. method_grid/Tt_grid/
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


N_chains   <- 3
plot_chain <- 1

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





# ---- Physical core count (identical logic to test_cpp.R / test_stan.R) ----
resolve_n_cores <- function(n_cores, N_chains) {
    if (!is.null(n_cores)) return(n_cores)

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
    min(phys, N_chains)
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
#
# theta1_dispersed distinguishes samplers where theta1 is genuine Markov
# chain state (amh_montoril, pg_as -- dispersion matters for R_hat) from
# samplers where theta1 is redrawn via SIR/IS or HMC every iteration
# (sir_laplace, sir_collapsed, stan -- zero-init is sufficient).
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

    theta1_dispersed <- method %in% c("amh_montoril", "pg_as")
    if (theta1_dispersed) {
        theta1_init <- theta1_ref + rnorm(Tt, 0, sd = 2 * sqrt(W1_ref))
    } else {
        theta1_init <- numeric(Tt)
    }
    theta2_init <- numeric(Tt)

    list(chain_id = chain_id, seed = seed,
         theta_01 = theta_01_init, theta_02 = theta_02_init,
         W1 = W1_init, W2 = W2_init,
         theta1 = theta1_init, theta2 = theta2_init,
         theta1_tilde = rep(0, Tt),
         varsigma2 = rep(0.02, Tt))
}


# ---- One-chain runner for the Rcpp samplers (mirrors test_cpp.R exactly) ----
run_one_chain_cpp <- function(init, method, y, N, K) {
    printf("Chain %d/%d (seed=%d)", init$chain_id, N_chains, init$seed)
    set.seed(init$seed)

    if (method == "pg_as") {
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
        stop(sprintf("Unknown Rcpp method: '%s'.", method))
    }
}


# ---- Main entry point: runs ONE row of calibration_grid ----
#
# This is the function that becomes the batchMap() target in the next
# step -- everything above it is shared setup, everything it does is
# self-contained (load data, run N_chains chains, compute diagnostics,
# save raw + plots + summary, return the summary row).
#
# N/burnin/K are now explicit arguments (one row of calibration_grid),
# not looked up from a per-method table -- a method can be run with more
# than one (N, burnin, K) configuration (see method_configs above), so the
# caller must say which one this particular task is.
run_calibration_task <- function(method, f, Tt, N, burnin, K, replica = 1) {

    task_name    <- task_name_for(method, f, Tt, N, burnin, K)
    task_rds     <- file.path(path_results, paste0(task_name, ".rds"))
    plot_file    <- file.path(path_plots, paste0(task_name, ".pdf"))
    printf("==== Task: %s ====", task_name)

    # ---- Checkpoint: skip entirely if this task's .rds already exists ----
    #
    # Same spirit as simulation.R's run_task(), but stronger: a task
    # already run is skipped COMPLETELY -- no data reload, no diagnostics
    # recomputation, no replotting -- not just "sampler skipped, everything
    # else redone". Replotting in particular is not cheap (Tt-point
    # traceplots/histograms), and is pure waste for a task whose .rds,
    # .pdf, and summary.csv row are all already on disk from a previous
    # run. Returns NULL; run_and_save_task() treats that as "nothing new
    # to report" and does not touch summary.csv for this task -- filling
    # summary.csv for already-checkpointed tasks is calibration_aggregate.R's
    # job (it reads existing rows straight from a prior summary.csv/the
    # .rds files, without recomputing anything either), not this script's.
    #
    # If N/burnin/K actually changed in method_config, this stale checkpoint
    # would silently keep reporting the OLD run's diagnostics -- delete the
    # .rds by hand to force a genuine re-run when that happens.
    if (file.exists(task_rds)) {
        printf("Task already run (%s exists) -- skipping entirely.", task_rds)
        return(NULL)
    }

    # ---- Load data (needed for diagnostics either way) ----
    source_name <- sprintf("%s_%s_%s", f, Tt, replica)
    data <- readRDS(file.path("..", "data", "simulated", paste0(source_name, ".rds")))
    y <- data$y
    theta1_true <- data$theta
    changepoints <- c(0.25, 0.5, 0.75) * Tt

    if (Tt == 200)  t_obs <- c(50, 100, 150, 175)
    if (Tt == 400)  t_obs <- c(75, 100, 200, 300)
    if (Tt == 800)  t_obs <- c(100, 300, 500, 700)
    if (Tt == 1600) t_obs <- c(400, 800, 1200, 1600)

    # ---- Seed (same deterministic formula used throughout the project) ----
    method_idx <- match(method, method_grid)
    Tt_idx     <- match(Tt, Tt_grid)
    f_idx      <- match(f, function_grid)
    seed_base  <- method_idx * 1e5 + f_idx * 1e4 + Tt_idx * 1e3 + replica * 10

    printf("Running %s for %s, seed_base=%d, N=%d, burnin=%d%s",
           method, source_name, seed_base, N, burnin,
           if (!is.na(K)) sprintf(", K=%d", K) else "")

    # ---- Chain initializations ----
    chain_inits <- lapply(1:N_chains, function(k) {
        make_chain_inits(k, seed_base + (k - 1), y, Tt, method)
    })

    # ---- Run ----
    start_time <- Sys.time()

    if (method == "stan") {
        options(mc.cores = resolve_n_cores(NULL, N_chains))
        rstan::rstan_options(auto_write = FALSE)

        if (file.exists("../cache/poisson_ltdm.rds")) {
            model <- readRDS("../cache/poisson_ltdm.rds")
        } else {
            printf("Building the Stan model")
            model <- rstan::stan_model(
                file = "../PoissonLTDM/inst/stan/poisson_ltdm.stan",
                model_name = "PoissonLTDM")
            saveRDS(model, file = "../cache/poisson_ltdm.rds")
        }

        stan_out <- sample_stan(
            model = model, y = y, N = N, burnin = burnin, seed = seed_base,
            N_chains = N_chains, n_cores = resolve_n_cores(NULL, N_chains),
            chain_inits = chain_inits,
            mu_01 = mu_01, sigma2_01 = sigma2_01, mu_02 = mu_02, sigma2_02 = sigma2_02,
            nu_01 = nu_01, eta_01 = eta_01, nu_02 = nu_02, eta_02 = eta_02,
            verbose = verbose)
        results <- stan_out$results

    } else {
        n_cores_used <- resolve_n_cores(NULL, N_chains)
        printf("Running %d chains in parallel on %d cores", N_chains, n_cores_used)

        cl <- parallel::makeCluster(n_cores_used)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        doParallel::registerDoParallel(cl)
        `%dorng%` <- doRNG::`%dorng%`

        current_wd <- getwd()
        parallel::clusterExport(cl, "current_wd", envir = environment())

        # run_one_chain_cpp() and the hyperparameter/config globals it uses
        # (mu_01, sigma2_01, ..., ac_ref, M_irls_max, M_is, tol, R_prerun,
        # M_sir, verbose, print_every) are defined at the script's top level
        # (.GlobalEnv), not inside run_calibration_task() -- foreach's
        # automatic variable detection only walks the LEXICAL environment of
        # the %dorng% call (i.e. run_calibration_task()'s own environment),
        # so it does not find them. They must be exported explicitly, from
        # .GlobalEnv, for the worker processes to see them.
        parallel::clusterExport(cl,
            c("run_one_chain_cpp", "printf",
              "mu_01", "sigma2_01", "mu_02", "sigma2_02",
              "nu_01", "eta_01", "nu_02", "eta_02",
              "ac_ref", "M_is", "M_irls_max", "tol", "R_prerun", "M_sir",
              "verbose", "print_every", "N_chains"),
            envir = .GlobalEnv)

        parallel::clusterEvalQ(cl, {
            setwd(current_wd)
            pkgload::load_all("../PoissonLTDM", debug = FALSE)
        })

        results <- foreach::foreach(init = chain_inits) %dorng% {
            run_one_chain_cpp(init, method, y, N, K)
        }
    }

    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    printf("Total wall-clock time (%d chains): %.2f s", N_chains, elapsed_time)

    gc(full = TRUE)

    # ---- Diagnostics: print (console) + plot (redirected to a PDF file) ----
    #
    # print_and_plot_diagnostics() prints to the console, plots directly to
    # the active graphics device, and returns the R_hat/ESS values it
    # computed (rhat_max, ess_bulk_min, and the per-parameter breakdowns --
    # see plot_diagnostics.R). Wrapping the call in pdf()/dev.off() is the
    # standard R way to redirect that plotting to a file without touching
    # plot_diagnostics.R.
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

    # ---- Save raw results ----
    saveRDS(list(results = results, elapsed_time = elapsed_time,
                 method = method, f = f, Tt = Tt, replica = replica,
                 N = N, burnin = burnin, K = K, seed_base = seed_base),
            file = task_rds)

    # ---- One-row summary for later aggregation across all 20 tasks ----
    summary_row <- build_summary_row(diag, method, f, Tt, replica, N, burnin, K, elapsed_time)

    return(summary_row)
}


# ---- run_and_save_task(): runs one grid row AND updates summary.csv ----
#
# Wraps run_calibration_task() with the summary.csv read-modify-write
# logic. Used directly by run_mode == "local" below. NOT used by
# run_mode == "cluster": each PBS job there calls this too (it is what
# batchMap() dispatches), but concurrent jobs writing to the same
# summary.csv at once would race/clobber each other, so cluster-mode jobs
# only produce their own .rds -- calibration_aggregate.R builds
# summary.csv from those .rds files afterwards, once all jobs are done.
#
# Any pre-existing summary.csv row for the same (method, f, Tt, replica)
# is dropped before appending the new one, so re-running a task after a
# config change (N, burnin, K, ...) updates that task's row instead of
# duplicating it. Column set may grow over time (e.g. a method-specific
# field appearing for the first time) -- rbind() with a mismatched column
# set errors out rather than silently misaligning columns, which is
# intentional here.
#
# run_calibration_task() returns NULL for an already-checkpointed task
# (see its own header comment) -- there is nothing new to summarize in
# that case, so update_summary_csv is skipped and this function returns
# NULL too. That task's summary.csv row (if any) is left exactly as it
# was; calibration_aggregate.R is what fills in rows for checkpointed
# tasks, not this function.
run_and_save_task <- function(method, f, Tt, N, burnin, K, replica = 1, update_summary_csv = TRUE) {
    summary_row <- run_calibration_task(method, f, Tt, N, burnin, K, replica)

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


# ==========================================================================
# ---- Dispatch: run_mode == "local" or "cluster" (see "Run control") ----
# ==========================================================================

# ---- Guard against recursive dispatch inside a batchtools/PBS job, or
# when this file is sourced only for its definitions (e.g. by
# calibration_aggregate.R) ----
#
# makeRegistry(source = "calibration_phase.R") makes every job re-source
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
        printf("Running all %d task(s) of the (filtered) grid, sequentially.", n_tasks)
        for (i in seq_len(n_tasks)) {
            task <- calibration_grid[i, ]
            printf("[%d/%d] method=%s, f=%s, Tt=%d, N=%d, burnin=%d, K=%s",
                   i, n_tasks, task$method, task$f, task$Tt, task$N, task$burnin,
                   if (is.na(task$K)) "NA" else task$K)
            run_and_save_task(task$method, task$f, task$Tt, task$N, task$burnin, task$K, replica)
        }
    } else {
        task <- calibration_grid[task_id, ]
        printf("Selected task %d/%d: method=%s, f=%s, Tt=%d, N=%d, burnin=%d, K=%s",
               task_id, n_tasks, task$method, task$f, task$Tt, task$N, task$burnin,
               if (is.na(task$K)) "NA" else task$K)

        summary_row <- run_and_save_task(task$method, task$f, task$Tt, task$N, task$burnin, task$K, replica)
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
        source = "calibration_phase.R",  # each job re-sources this whole
                                          # file (run_mode is irrelevant
                                          # inside the job -- batchMap()'s
                                          # target function is what runs,
                                          # not the dispatch block below)
        seed = 1
    )

    reg$cluster.functions <- makeClusterFunctionsTORQUE("calibration_pbs.tmpl")

    # ---- Skip rows whose .rds is already checkpointed ----
    #
    # run_calibration_task() itself already skips a checkpointed task
    # (see its header comment) -- but only after a PBS job has been
    # queued/started for it, i.e. after paying for a full place=excl node
    # allocation just to check one file and exit. Filtering here avoids
    # creating those jobs in the first place. Uses task_name_for(), the
    # same naming convention run_calibration_task() uses for its own
    # task_rds path.
    task_names <- with(calibration_grid, mapply(task_name_for, method, f, Tt, N, burnin, K))
    task_rds_exists <- file.exists(file.path(path_results, paste0(task_names, ".rds")))
    n_skipped <- sum(task_rds_exists)
    if (n_skipped > 0) {
        printf("Skipping %d already-checkpointed task(s) (no job submitted for them):", n_skipped)
        for (name in task_names[task_rds_exists]) printf("  %s", name)
    }
    calibration_grid_to_submit <- calibration_grid[!task_rds_exists, ]

    if (nrow(calibration_grid_to_submit) == 0) {
        printf("Every task in the (filtered) grid is already checkpointed -- nothing to submit.")
        printf("Run calibration_aggregate.R to build/update summary.csv from the existing .rds files.")
    } else {

    # update_summary_csv = FALSE: see run_and_save_task()'s header comment
    # -- concurrent jobs must not race on the same summary.csv. N/burnin/K
    # are mapped per-row (like method/f/Tt), not passed via more.args,
    # since they can now differ between rows of the same method (see
    # method_configs above) -- only replica is genuinely constant across
    # the whole grid.
    ids <- batchMap(
        fun = function(method, f, Tt, N, burnin, K, replica) {
            run_and_save_task(method, f, Tt, N, burnin, K, replica, update_summary_csv = FALSE)
        },
        method = calibration_grid_to_submit$method, f = calibration_grid_to_submit$f,
        Tt = calibration_grid_to_submit$Tt, N = calibration_grid_to_submit$N,
        burnin = calibration_grid_to_submit$burnin, K = calibration_grid_to_submit$K,
        more.args = list(replica = replica),
        reg = reg
    )

    # walltime_hours = 3 covers every method/Tt combination with generous
    # margin -- the most expensive case measured during calibration (stan,
    # Tt=1600) took ~2.1h for 3 chains; amh_montoril/pg_as/sir_laplace/
    # sir_collapsed are all well under an hour even at Tt=1600. Reduced
    # from an initial 6h: schedulers typically use requested walltime for
    # backfilling decisions, and a long requested walltime can make a job
    # harder to slot in even when physical nodes are free -- observed only
    # 2 jobs running concurrently out of 18 submitted despite place=excl
    # nodes apparently being available, which requesting less walltime
    # (while still safely covering the measured worst case) may help with.
    submitJobs(ids, resources = list(ncpus = N_chains, walltime_hours = 3), reg = reg)

    printf("Submitted %d job(s) to the cluster (registry: %s).", nrow(ids), reg$file.dir)
    printf("Check progress with batchtools::getStatus(loadRegistry(\"registry_calibration\")).")
    printf("Once all jobs are done, run calibration_aggregate.R to build summary.csv.")

    }

} else {
    stop(sprintf("Unknown run_mode: '%s'. Use \"local\" or \"cluster\".", run_mode))
}

}  # end !running_inside_pbs_job && !.defs_only_flag && !making_registry guard
