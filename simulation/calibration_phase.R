# simulation/calibration_phase.R
#
# Official, reproducible run of the calibration-phase grid (5 methods x 2
# functions x 2 Tt, replica = 1, N_chains = 3), documenting the N/burnin/K
# configuration decided for each method during the pre-calibration
# discussion. Reuses, unchanged, the same chain-initialization scheme,
# hyperparameters, and diagnostic pipeline as tests/test_cpp.R and
# tests/test_stan.R -- this script is not a reimplementation, it is those
# two scripts made parametrizable over the calibration grid, plus disk
# output (raw results + plots + a one-row summary) for later aggregation.
#
# THIS VERSION: local/manual use. Run it directly (Rscript calibration_phase.R
# or source() it in an R session) to execute ONE row of the grid, selected
# via task_id below -- exactly like test_cpp.R/test_stan.R let you pick one
# (method, f, Tt) combination at the top of the file. Once this works and
# is validated locally, task_id becomes the batchtools chunk/job index and
# run_calibration_task() becomes the function passed to batchMap() -- no
# other change needed.

rm(list = ls())
options(error = function() traceback(2))
setwd(dirname(this.path::this.path()))

# See test_cpp.R for why debug=FALSE matters here (elapsed_time feeds
# directly into ESS/second diagnostics inside print_and_plot_diagnostics()).
pkgload::load_all("../PoissonLTDM", debug = FALSE)

source("../PoissonLTDM/R/sampler_stan.R")
source("../tests/plot_diagnostics.R")

printf <- function(...) cat(paste(sprintf(...), "\n"))

# ---- Output paths ----
path_results <- "results/calibration"
path_plots   <- file.path(path_results, "plots")
dir.create(path_results, showWarnings = FALSE, recursive = TRUE)
dir.create(path_plots,   showWarnings = FALSE, recursive = TRUE)


# ---- Calibration grid ----
#
# Full calibration phase: 20 combinations (5 methods x 2 functions x 2 Tt),
# replica fixed at 1. method_grid/Tt_grid/function_grid below are the FULL
# reference grids (used for the deterministic seed formula) and are kept
# complete regardless of which subset calibration_grid actually runs.
#
# LOCAL TEST SUBSET (current): 8 tasks (sir_laplace, sir_collapsed) x 2
# functions x 2 Tt -- restricted here to validate the pipeline cheaply
# before running the full 5-method grid. Restore method = method_grid to
# run all 20.
method_grid   <- c("amh_montoril", "pg_as", "sir_laplace", "sir_collapsed", "stan")
Tt_grid       <- c(200, 400, 800, 1600)          # full reference grid (for seed encoding)
function_grid <- c("constant", "linear", "quadratic", "sinusoidal") # full reference grid

calibration_grid <- expand.grid(
    method = method_grid,
    f      = c("constant", "quadratic"),
    Tt     = c(200, 1600),
    stringsAsFactors = FALSE
)

replica <- 1  # fixed for the whole calibration phase

n_tasks <- nrow(calibration_grid)
printf("Calibration grid: %d tasks", n_tasks)


# ---- Per-method configuration, as decided during calibration ----
#
# N / burnin (and, where applicable, K = number of SMC particles) fixed per
# method across the ENTIRE Tt grid -- deliberately NOT re-tuned per Tt, so
# that ESS/time comparisons across Tt stay apples-to-apples (see
# discussion: amh_montoril and pg_as both showed degradation at Tt=1600
# that is treated as a genuine finding, not "fixed" by inflating N further).
method_config <- list(
    amh_montoril  = list(N = 55000, burnin = 5000),
    pg_as         = list(N = 11000, burnin = 1000, K = 100),
    sir_laplace   = list(N = 11000, burnin = 1000),
    sir_collapsed = list(N = 11000, burnin = 1000),
    stan          = list(N = 11000, burnin = 1000)
)

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
run_calibration_task <- function(method, f, Tt, replica = 1) {

    cfg <- method_config[[method]]
    if (is.null(cfg)) stop(sprintf("No configuration registered for method '%s'.", method))
    N      <- cfg$N
    burnin <- cfg$burnin
    K      <- if (!is.null(cfg$K)) cfg$K else NA_integer_

    task_name <- sprintf("%s_%s_%s", method, f, Tt)
    printf("==== Task: %s ====", task_name)

    # ---- Load data ----
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
    plot_file <- file.path(path_plots, paste0(task_name, ".pdf"))
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
            file = file.path(path_results, paste0(task_name, ".rds")))

    # ---- One-row summary for later aggregation across all 20 tasks ----
    #
    # Every field comes directly from print_and_plot_diagnostics()'s return
    # value (diag) -- no recomputation. Per-chain vectors (length N_chains)
    # are flattened into one column per chain (e.g. ess_w1_chain1,
    # ess_w1_chain2, ess_w1_chain3) via flatten_chains(), so the whole
    # summary stays a single data.frame row, matching every line printed to
    # the console for this task.
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

    # nz(): NULL -> NA (scalar). data.frame() turns a NULL list element into
    # a 0-row column, not NA, which breaks row-binding whenever a field is
    # NULL for the current method (e.g. ac_ratio_mean is NULL for anything
    # other than amh_montoril). Every scalar diag$* field that can be NULL
    # must be wrapped in nz() before entering the data.frame() call below.
    nz <- function(x) if (is.null(x)) NA_real_ else x

    summary_row <- data.frame(
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

    return(summary_row)
}


# =====================================================================
# LOCAL TEST DRIVER -- run ONE task at a time, exactly like test_cpp.R /
# test_stan.R let you pick one (method, f, Tt) at the top of the file.
# Change task_id below and re-run to test a different combination.
# =====================================================================

# ---- Choose which row of the calibration grid to run here ----
task_id <- 5
# ----------------------------------------------------------------------

task <- calibration_grid[task_id, ]
printf("Selected task %d/%d: method=%s, f=%s, Tt=%d",
       task_id, n_tasks, task$method, task$f, task$Tt)

summary_row <- run_calibration_task(task$method, task$f, task$Tt, replica)

printf("Summary row: rhat_max=%.4f, ess_bulk_min=%.1f, elapsed_time=%.2fs",
       summary_row$rhat_max, summary_row$ess_bulk_min, summary_row$elapsed_time_s)

# ---- Append/update the consolidated summary.csv ----
#
# Any pre-existing row for the same (method, f, Tt, replica) is dropped
# before appending the new one, so re-running a task_id after a config
# change (N, burnin, K, ...) updates that task's row instead of duplicating
# it. Column set may grow over time (e.g. a method-specific field appearing
# for the first time) -- rbind() with a mismatched column set errors out
# rather than silently misaligning columns, which is intentional here.
summary_file <- file.path(path_results, "summary.csv")
if (file.exists(summary_file)) {
    summary_all <- read.csv(summary_file, stringsAsFactors = FALSE)
    keep <- !(summary_all$method == summary_row$method &
              summary_all$f == summary_row$f &
              summary_all$Tt == summary_row$Tt &
              summary_all$replica == summary_row$replica)
    summary_all <- rbind(summary_all[keep, , drop = FALSE], summary_row)
} else {
    summary_all <- summary_row
}
write.csv(summary_all, summary_file, row.names = FALSE)
printf("Summary written to: %s (%d task(s) total)", summary_file, nrow(summary_all))
