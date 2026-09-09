# tests/test_rhat_ess_fast.R
#
# Validates PoissonLTDM::rhat_ess_fast() (src/convergence.cpp) against
# posterior::rhat()/ess_bulk()/ess_tail() -- confirms the C++
# reimplementation (RcppArmadillo + FFTW3 + OpenMP) is mathematically
# faithful to the Vehtari et al. (2021) algorithm used by the posterior
# package, within floating-point precision.
#
# Also benchmarks PoissonLTDM::metrics_convergence() /
# metrics_convergence_by_chain() (which call rhat_ess_fast() under the
# hood) against posterior::summarise_draws() at the real problem scale
# (n_post=9000, N_chains=3, Tt=1600).
#
# Run from inside cobalebeb2027/tests/:
#   Rscript test_rhat_ess_fast.R

# IMPORTANT: load_all() compiles C++ in DEBUG mode by default (-g -O0),
# silently overriding any -O3 in src/Makevars. This causes a ~4-5x
# slowdown vs an optimized build and would make Part 3's speed benchmark
# meaningless (compares an unoptimized rhat_ess_fast() against posterior,
# understating the real speedup by roughly that factor). Using
# pkgload::load_all() directly (not devtools::load_all()) because
# devtools::load_all()'s `...` does not reliably forward debug= to pkgload
# in this environment (devtools 2.5.2 / pkgload 1.5.3), producing a "must
# be used" warning even though the argument is valid. debug=FALSE forces
# an optimized rebuild; recompile is left at its default (FALSE) since
# load_all() already detects changed .cpp/.h files via timestamp and
# recompiles only what's needed.
pkgload::load_all("../PoissonLTDM", debug = FALSE)
library(posterior)

set.seed(123)

# ============================================================
# Part 1: correctness validation, calling rhat_ess_fast() directly
# ============================================================

compare_case <- function(n_post, N_chains, Tt, phi = 0.9, label = "") {
    gen_chain <- function(n, phi, sd = 1) {
        e <- rnorm(n, sd = sd)
        as.numeric(stats::filter(e, filter = phi, method = "recursive"))
    }
    arr <- array(NA_real_, dim = c(n_post, N_chains, Tt))
    for (tt in 1:Tt) for (cc in 1:N_chains) arr[, cc, tt] <- gen_chain(n_post, phi)

    draws <- posterior::as_draws_array(arr)
    ref <- posterior::summarise_draws(draws, rhat, ess_bulk, ess_tail)
    pkg_res <- PoissonLTDM::rhat_ess_fast(arr)

    cat(sprintf("--- %s (n_post=%d, N_chains=%d, Tt=%d) ---\n", label, n_post, N_chains, Tt))
    cat(sprintf("  max|rhat diff|:      %.2e\n", max(abs(ref$rhat - pkg_res$rhat))))
    cat(sprintf("  max|ess_bulk diff| (rel): %.2e\n", max(abs(ref$ess_bulk - pkg_res$ess_bulk) / ref$ess_bulk)))
    cat(sprintf("  max|ess_tail diff| (rel): %.2e\n\n", max(abs(ref$ess_tail - pkg_res$ess_tail) / ref$ess_tail)))
}

cat("========== Part 1: correctness validation (rhat_ess_fast vs posterior) ==========\n\n")

compare_case(200, 3, 5, 0.5, "Small, low autocorrelation")
compare_case(200, 3, 5, 0.98, "Small, high autocorrelation")
compare_case(9000, 3, 10, 0.9, "Real scale n_post=9000")
compare_case(201, 3, 3, 0.7, "ODD n_post (edge case)")
compare_case(500, 2, 5, 0.8, "N_chains=2")
compare_case(9000, 1, 5, 0.9, "N_chains=1 (used by metrics_ess_single_chain)")


# ============================================================
# Part 2: validation of metrics_convergence() / metrics_convergence_by_chain()
# (the layer plot_diagnostics.R actually calls)
# ============================================================

cat("========== Part 2: metrics_convergence() / metrics_convergence_by_chain() validation ==========\n\n")

n_post <- 300
burnin <- 50
N <- n_post + burnin
N_chains <- 3
Tt <- 6

gen_chain <- function(n, phi = 0.9, sd = 1) {
    e <- rnorm(n, sd = sd)
    as.numeric(stats::filter(e, filter = phi, method = "recursive"))
}

result_list <- lapply(1:N_chains, function(c) {
    list(
        W1_hist = gen_chain(N, phi = 0.7),
        theta1_hist = sapply(1:Tt, function(t) gen_chain(N, phi = 0.95))
    )
})

# --- metrics_convergence, scalar ---
res_w1 <- PoissonLTDM::metrics_convergence(result_list, "W1_hist", burnin)
post_idx <- (burnin + 1):N
w1_arr <- array(NA_real_, dim = c(length(post_idx), N_chains, 1))
for (c in 1:N_chains) w1_arr[, c, 1] <- result_list[[c]]$W1_hist[post_idx]
ref_w1 <- posterior::summarise_draws(posterior::as_draws_array(w1_arr), rhat, ess_bulk, ess_tail)

cat("--- metrics_convergence(), scalar (W1_hist) ---\n")
cat(sprintf("  rhat diff:     %.2e\n", abs(res_w1$rhat - ref_w1$rhat)))
cat(sprintf("  ess_bulk diff: %.2e\n", abs(res_w1$ess_bulk - ref_w1$ess_bulk)))
cat(sprintf("  ess_tail diff: %.2e\n\n", abs(res_w1$ess_tail - ref_w1$ess_tail)))

# --- metrics_convergence, matrix ---
res_t1 <- PoissonLTDM::metrics_convergence(result_list, "theta1_hist", burnin)
t1_arr <- array(NA_real_, dim = c(length(post_idx), N_chains, Tt))
for (c in 1:N_chains) t1_arr[, c, ] <- result_list[[c]]$theta1_hist[post_idx, ]
ref_t1 <- posterior::summarise_draws(posterior::as_draws_array(t1_arr), rhat, ess_bulk, ess_tail)

cat("--- metrics_convergence(), matrix (theta1_hist) ---\n")
cat(sprintf("  max|rhat diff|:     %.2e\n", max(abs(res_t1$rhat - ref_t1$rhat))))
cat(sprintf("  max|ess_bulk diff|: %.2e\n", max(abs(res_t1$ess_bulk - ref_t1$ess_bulk))))
cat(sprintf("  max|ess_tail diff|: %.2e\n\n", max(abs(res_t1$ess_tail - ref_t1$ess_tail))))

# --- metrics_convergence_by_chain, matrix (per-chain ESS, not pooled) ---
res_bychain <- PoissonLTDM::metrics_convergence_by_chain(result_list, "theta1_hist", burnin)
ess_bulk_ref <- matrix(NA_real_, nrow = Tt, ncol = N_chains)
ess_tail_ref <- matrix(NA_real_, nrow = Tt, ncol = N_chains)
for (c in 1:N_chains) {
    arr_c <- array(t1_arr[, c, ], dim = c(length(post_idx), 1, Tt))
    summ_c <- posterior::summarise_draws(posterior::as_draws_array(arr_c), ess_bulk, ess_tail)
    ess_bulk_ref[, c] <- summ_c$ess_bulk
    ess_tail_ref[, c] <- summ_c$ess_tail
}

cat("--- metrics_convergence_by_chain(), matrix (theta1_hist) ---\n")
cat(sprintf("  max|ess_bulk diff|: %.2e\n", max(abs(res_bychain$ess_bulk - ess_bulk_ref))))
cat(sprintf("  max|ess_tail diff|: %.2e\n\n", max(abs(res_bychain$ess_tail - ess_tail_ref))))

# --- metrics_ess_single_chain (used in simulation.R) ---
w1_single <- result_list[[1]]$W1_hist[post_idx]
res_single <- PoissonLTDM::metrics_ess_single_chain(w1_single)
single_arr <- array(w1_single, dim = c(length(w1_single), 1, 1))
ref_single <- posterior::summarise_draws(posterior::as_draws_array(single_arr), ess_bulk, ess_tail)

cat("--- metrics_ess_single_chain(), scalar (1 chain) ---\n")
cat(sprintf("  ess_bulk diff: %.2e\n", abs(res_single$ess_bulk - ref_single$ess_bulk)))
cat(sprintf("  ess_tail diff: %.2e\n\n", abs(res_single$ess_tail - ref_single$ess_tail)))


# ============================================================
# Part 3: speed benchmark at real scale
# ============================================================

cat("========== Part 3: speed benchmark (real scale) ==========\n\n")

n_post <- 9000
N_chains <- 3
Tt <- 1600

gen_chain_fast <- function(n, phi = 0.98, sd = 1) {
    e <- rnorm(n, sd = sd)
    as.numeric(stats::filter(e, filter = phi, method = "recursive"))
}

cat(sprintf("Generating data: n_post=%d, N_chains=%d, Tt=%d\n", n_post, N_chains, Tt))
arr <- array(NA_real_, dim = c(n_post, N_chains, Tt))
for (tt in 1:Tt) for (cc in 1:N_chains) arr[, cc, tt] <- gen_chain_fast(n_post)
cat("Data generated.\n\n")

t0 <- Sys.time()
pkg_res <- PoissonLTDM::rhat_ess_fast(arr)
t_pkg <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("PoissonLTDM::rhat_ess_fast(): %.2f s\n", t_pkg))

draws <- posterior::as_draws_array(arr)
t0 <- Sys.time()
ref <- posterior::summarise_draws(draws, rhat, ess_bulk, ess_tail, .cores = 1)
t_r_seq <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("posterior::summarise_draws (.cores=1): %.1f s\n", t_r_seq))
cat(sprintf("Speedup: %.1fx\n\n", t_r_seq / t_pkg))

cat("=== Correctness validation at real scale ===\n")
cat(sprintf("max|rhat diff|:      %.2e\n", max(abs(ref$rhat - pkg_res$rhat))))
cat(sprintf("max|ess_bulk diff| (rel): %.2e\n", max(abs(ref$ess_bulk - pkg_res$ess_bulk) / ref$ess_bulk)))
cat(sprintf("max|ess_tail diff| (rel): %.2e\n", max(abs(ref$ess_tail - pkg_res$ess_tail) / ref$ess_tail)))

cat("\n========== Done ==========\n")
