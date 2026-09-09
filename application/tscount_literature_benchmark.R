# application/tscount_literature_benchmark.R
#
# Reproduces the two observation-driven literature benchmarks reported in
# Table (predictive log-likelihood, campy application): INGARCH (identity
# link) and log-linear (log link) models, both under a Poisson observation
# distribution. All via tscount::tsglm() (Liboschik, Fokianos & Fried,
# 2017).
#
# IMPORTANT PROVENANCE NOTE (see chat discussion -- do not drop when citing
# these numbers):
#   - INGARCH Poisson: refits the MODEL SPECIFICATION of Ferland, Latour &
#     Oraichi (2006) -- same past_obs/past_mean structure -- but tsglm's
#     optimizer (constrOptim: barrier + BFGS, Liboschik et al. 2017, Eq. 8)
#     converges to different parameter estimates than their original Excel
#     Solver (GRG2) fit (e.g. c0_hat=1.4965 here vs. their reported 1.9280).
#     Their paper does not print a log-likelihood value at all -- the
#     -434.384 figure is entirely our own refit, not a number from the
#     paper.
#   - Log-linear Poisson: not a reproduction of a published fit to the
#     campy series. It follows the general model class described in
#     Liboschik et al. (2017) (log link = Fokianos & Tjostheim 2011) but is
#     a fit WE performed, not a result from a specific paper.
#
# Estimation method (both): conditional maximum likelihood via constrOptim
# (barrier method + BFGS using the analytical score, Liboschik et al. 2017,
# Sec. 3.1).
#
# logLik() on a fitted tsglm object returns the TRUE (not quasi) log-
# likelihood including all constant terms (Liboschik et al. 2017, Sec. 5)
# -- i.e. it already includes the -log(y_t!) normalizing term, same
# convention as metrics_loglik()/the particle-filter predictive log-lik
# used for the PoissonLTDM methods elsewhere in this project. It is a
# ONE-STEP-AHEAD PREDICTIVE quasi-likelihood (lambda_t computed forward
# from y_1,...,y_{t-1} only) -- comparable to the particle-filter
# predictive log-likelihood, NOT to log_lik/log_cpo (which use the
# smoothed posterior mean, conditioning on the whole series).

setwd(dirname(this.path::this.path()))
path_results <- "../results/application"
dir.create(path_results, showWarnings = FALSE, recursive = TRUE)

library(tscount)
data(campy)

printf <- function(...) cat(paste(sprintf(...), "\n"))

# ---- CPU-time measurement convention: same as elsewhere in the project
# (e.g. application/real_data_run.R's run_one_chain()) -- system.time(),
# keeping the "user.self" component (CPU time consumed by the process
# itself), not "elapsed" (wall-clock, which can include I/O wait /
# scheduling noise unrelated to the algorithm). Note this is NOT an
# ESS/s-comparable quantity (see chat discussion): tsglm() is a single
# deterministic optimizer run (constrOptim/BFGS to the conditional-ML
# point estimate), not a sampler producing a correlated draw sequence --
# there is no autocorrelation structure from which an effective sample
# size could be defined. cpu_time_mean_sec/cpu_time_se_sec below are
# reported purely as a computational-cost figure, parallel to (but not
# directly combinable with) the ESS/s columns used for the five
# PoissonLTDM methods elsewhere in the project.
#
# Repeated-timing scheme (see chat discussion): since tsglm() is
# deterministic, ALL variation in a single system.time() reading is
# measurement noise (scheduler jitter, cold/warm cache, one-off R
# bytecode-compilation cost on first call to a function) -- not
# algorithmic variation. Two consequences, both handled below:
#   (i) a warm-up call is issued once per model and discarded, so the
#       one-off bytecode-compilation cost is not folded into the timed
#       blocks;
#  (ii) proc.time()'s tick resolution on Linux (incl. CentOS 7/Euler) is
#       fixed at delta=10ms (USER_HZ=100, a kernel ABI constant,
#       independent of CONFIG_HZ) regardless of machine -- to keep the
#       per-call quantization error negligible without needing a pilot
#       test, B calls are grouped inside one system.time() and divided by
#       B, repeated for R independent blocks, giving R i.i.d. estimates of
#       per-call CPU time from which mean and SE are computed.
B_timing <- 50   # calls per block -- keeps block time >> delta=10ms regardless of per-call cost
R_timing <- 30   # independent blocks -- SE(mean) = sd(block_times)/sqrt(R_timing)

time_model_cpu <- function(fit_call) {
	invisible(fit_call())  # warm-up: pays one-off bytecode-compilation cost, discarded
	block_times <- numeric(R_timing)
	for (r in seq_len(R_timing)) {
		block_times[r] <- system.time(for (b in seq_len(B_timing)) invisible(fit_call()))[["user.self"]] / B_timing
	}
	list(mean = mean(block_times), se = sd(block_times) / sqrt(R_timing), block_times = block_times)
}

# ---- 1. INGARCH, Poisson (identity link) -- Ferland, Latour & Oraichi
# (2006) model specification, refit via tsglm ----
fit_ingarch_poisson <- tsglm(
	campy,
	model = list(past_obs = 1, past_mean = c(7, 13)),
	link = "identity",
	distr = "poisson"
)
timing_ingarch_poisson <- time_model_cpu(function() tsglm(
	campy, model = list(past_obs = 1, past_mean = c(7, 13)), link = "identity", distr = "poisson"))

# ---- 2. Log-linear model, Poisson (log link) -- Fokianos & Tjostheim
# (2011) model class ----
fit_loglinear_poisson <- tsglm(
	campy,
	model = list(past_obs = 1, past_mean = 13),
	link = "log",
	distr = "poisson"
)
timing_loglinear_poisson <- time_model_cpu(function() tsglm(
	campy, model = list(past_obs = 1, past_mean = 13), link = "log", distr = "poisson"))

# ---- Summary table ----
results <- data.frame(
	model = c("INGARCH, Poisson", "Log-linear, Poisson"),
	link = c("identity", "log"),
	distr = c("poisson", "poisson"),
	predictive_loglik = c(logLik(fit_ingarch_poisson), logLik(fit_loglinear_poisson)),
	AIC = c(AIC(fit_ingarch_poisson), AIC(fit_loglinear_poisson)),
	cpu_time_mean_sec = c(timing_ingarch_poisson$mean, timing_loglinear_poisson$mean),
	cpu_time_se_sec = c(timing_ingarch_poisson$se, timing_loglinear_poisson$se),
	stringsAsFactors = FALSE
)

printf("=== Predictive log-likelihood, 2 literature benchmarks (campy) ===")
print(results, row.names = FALSE)

results_file <- sprintf("%s/tscount_literature_benchmark.csv", path_results)
write.csv(results, file = results_file, row.names = FALSE)
printf("\nSaved to %s", results_file)
