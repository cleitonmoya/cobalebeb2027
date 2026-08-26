# application/tscount_literature_benchmark.R
#
# Reproduces the four observation-driven literature benchmarks reported in
# Table (predictive log-likelihood, campy application): INGARCH (identity
# link) and log-linear (log link) models, each under Poisson and negative
# binomial observation distributions. All via tscount::tsglm() (Liboschik,
# Fokianos & Fried, 2017).
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
#   - INGARCH NegBin, log-linear Poisson, log-linear NegBin: none of these
#     three are reproductions of a published fit to the campy series. They
#     follow the general model classes described in Liboschik et al.
#     (2017) (identity link = INGARCH per Ferland et al. 2006; log link =
#     Fokianos & Tjostheim 2011) but are fits WE performed, not results
#     from a specific paper.
#
# Estimation method (all four): conditional maximum likelihood via
# constrOptim (barrier method + BFGS using the analytical score, Liboschik
# et al. 2017, Sec. 3.1). For the negative binomial cases, per Christou &
# Fokianos (2014): the regression coefficients are estimated by POISSON
# quasi-likelihood (provably independent of the dispersion parameter), and
# the dispersion parameter sigma^2 = 1/phi is estimated in a SEPARATE
# second step via a Pearson chi^2 moment equation (Liboschik et al. 2017,
# Eq. 10) -- not joint likelihood maximization. This is why the regression
# coefficients below come out numerically IDENTICAL between the Poisson
# and NegBin fit within each link family (confirmed below).
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

# ---- 2. INGARCH, negative binomial (identity link) ----
fit_ingarch_nbin <- tsglm(
	campy,
	model = list(past_obs = 1, past_mean = c(7, 13)),
	link = "identity",
	distr = "nbinom"
)
timing_ingarch_nbin <- time_model_cpu(function() tsglm(
	campy, model = list(past_obs = 1, past_mean = c(7, 13)), link = "identity", distr = "nbinom"))

# ---- 3. Log-linear model, Poisson (log link) -- Fokianos & Tjostheim
# (2011) model class ----
fit_loglinear_poisson <- tsglm(
	campy,
	model = list(past_obs = 1, past_mean = 13),
	link = "log",
	distr = "poisson"
)
timing_loglinear_poisson <- time_model_cpu(function() tsglm(
	campy, model = list(past_obs = 1, past_mean = 13), link = "log", distr = "poisson"))

# ---- 4. Log-linear model, negative binomial (log link) ----
fit_loglinear_nbin <- tsglm(
	campy,
	model = list(past_obs = 1, past_mean = 13),
	link = "log",
	distr = "nbinom"
)
timing_loglinear_nbin <- time_model_cpu(function() tsglm(
	campy, model = list(past_obs = 1, past_mean = 13), link = "log", distr = "nbinom"))

# ---- Summary table ----
results <- data.frame(
	model = c("INGARCH, Poisson", "INGARCH, NegBin", "Log-linear, Poisson", "Log-linear, NegBin"),
	link = c("identity", "identity", "log", "log"),
	distr = c("poisson", "nbinom", "poisson", "nbinom"),
	predictive_loglik = c(logLik(fit_ingarch_poisson), logLik(fit_ingarch_nbin),
						   logLik(fit_loglinear_poisson), logLik(fit_loglinear_nbin)),
	AIC = c(AIC(fit_ingarch_poisson), AIC(fit_ingarch_nbin),
			AIC(fit_loglinear_poisson), AIC(fit_loglinear_nbin)),
	cpu_time_mean_sec = c(timing_ingarch_poisson$mean, timing_ingarch_nbin$mean,
						   timing_loglinear_poisson$mean, timing_loglinear_nbin$mean),
	cpu_time_se_sec = c(timing_ingarch_poisson$se, timing_ingarch_nbin$se,
						 timing_loglinear_poisson$se, timing_loglinear_nbin$se),
	stringsAsFactors = FALSE
)

printf("=== Predictive log-likelihood, 4 literature benchmarks (campy) ===")
print(results, row.names = FALSE)

# ---- Confirms the regression-coefficient invariance noted above: Poisson
# vs. NegBin coefficients are identical within each link family (only the
# dispersion parameter sigmasq differs, estimated in a separate step). ----
printf("\n=== Confirming Christou & Fokianos (2014): regression coefficients")
printf("    identical between Poisson and NegBin fits within each link family ===")
printf("INGARCH (identity link):")
print(rbind(Poisson = coef(fit_ingarch_poisson)[1:4], NegBin = coef(fit_ingarch_nbin)[1:4]))
printf("\nLog-linear (log link):")
print(rbind(Poisson = coef(fit_loglinear_poisson), NegBin = coef(fit_loglinear_nbin)[1:3]))

write.csv(results, file = "tscount_literature_benchmark.csv", row.names = FALSE)
printf("\nSaved to tscount_literature_benchmark.csv")
