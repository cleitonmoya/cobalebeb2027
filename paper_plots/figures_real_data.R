# ============================================================================
# figures_real_data.R
#
# Generates all real-data-application figures (fit + HPD band, R_hat(t),
# raw ESS, ESS/s, cross-method agreement Delta_max(t), total CPU time,
# log_lik/log_cpo) for the Cobal/EBEB 2027 manuscript, from
# method_summaries.rds, delta_max.rds, summary.csv (application/real_data_run.R's
# output) and the real data series itself (campy.rds).
#
# SEPARATE from ../simulation/figures.R (see discussion): the simulation
# study has R=200 replicas x 4 functions x 4 T -- distributions worth
# boxplotting. The real-data application has ONE dataset, ONE T, and a
# single point estimate per (method, parameter) -- there is no distribution
# to box; every "ESS/s" or "R_hat" figure here is a bar or a line, not a
# boxplot. The style block below (theme_paper, METHOD_*, PARAM_*,
# save_figure) is DUPLICATED verbatim from figures.R, not sourced, so this
# script stays self-contained (same convention already used for R_config
# across simulation_run.R / calibration_run.R / real_data_run.R) -- if the
# shared style changes in figures.R, mirror the change here too.
#
# Figures produced (see bottom of file):
#   real_fit.pdf              - y_t vs lambda_hat + 95% HPD band, facet by method
#   real_rhat.pdf              - R_hat(t), facet by parameter (theta_t1, theta_t2)
#   real_ess.pdf                - raw ESS (bulk + tail), facet by parameter (6)
#   real_ess_sec.pdf           - ESS/s (bulk + tail), facet by parameter (6)
#   real_delta_max.pdf         - cross-method agreement Delta_max(t)
#   real_total_time.pdf        - total CPU time by method
#   real_fit_metrics.pdf       - log_lik and log_cpo, dot plot
#   article_real_data_fit.pdf  - (article) representative fit + Delta_max(t)
#   article_real_data_ess.pdf  - (article) raw ESS + ESS/s, theta_t1/theta_t2/W1/W2
# ============================================================================

setwd(dirname(this.path::this.path()))

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(gridExtra)
library(grid)


# ---- Shared configuration (edit here to affect all figures) --------------

DATA_DIR <- "../data/real"
RESULTS_DIR <- "../results/real_data"
OUTPUT_DIR <- "plots/real_data"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

PATH_DATA <- file.path(DATA_DIR, "campy.rds")
PATH_SUMMARIES <- file.path(RESULTS_DIR, "method_summaries.rds")
PATH_DELTA_MAX <- file.path(RESULTS_DIR, "delta_max.rds")
PATH_SUMMARY_CSV <- file.path(RESULTS_DIR, "summary.csv")

# Method identifiers, plotting order, and English display labels -- DUPLICATED
# verbatim from ../simulation/figures.R (see header note).
METHOD_LEVELS <- c("montoril", "pg_apf", "sir_laplace", "sir_collapsed", "stan")

METHOD_LABELS <- c(
	montoril = "AM",
	pg_apf = "PG-AS",
	sir_laplace = "SIR",
	sir_collapsed = "SIR-Collapsed",
	stan = "Stan"
)

METHOD_COLORS <- c(
	montoril = "#D55E00",
	pg_apf = "#CC79A7",
	sir_laplace = "#E69F00",
	sir_collapsed = "#0072B2",
	stan = "#009E73"
)

METHOD_LINETYPES <- c(
	montoril = "solid",
	pg_apf = "twodash",
	sir_laplace = "dotted",
	sir_collapsed = "dotdash",
	stan = "longdash"
)

METHOD_SHAPES <- c(
	montoril = 16,
	pg_apf = 17,
	sir_laplace = 15,
	sir_collapsed = 18,
	stan = 4
)

# Parameter identifiers, order, and expression labels -- DUPLICATED verbatim
# from ../simulation/figures.R.
PARAM_LEVELS <- c("theta_t1", "theta_t2", "theta_01", "theta_02", "W1", "W2")

PARAM_LABELS <- c(
	theta_t1 = "theta[t1]",
	theta_t2 = "theta[t2]",
	theta_01 = 'theta["01"]',
	theta_02 = 'theta["02"]',
	W1 = "W[1]",
	W2 = "W[2]"
)

FULL_WIDTH <- 6
BASE_FONT_SIZE <- 9

PDF_DEVICE <- cairo_pdf


# ---- Shared ggplot theme (duplicated verbatim from figures.R) -------------

theme_paper <- function(base_size = BASE_FONT_SIZE) {
	theme_minimal(base_size = base_size) +
		theme(
			panel.grid.minor = element_blank(),
			panel.grid.major = element_line(color = "grey70", linewidth = 0.15),
			panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
			strip.text = element_text(face = "bold", size = rel(0.95)),
			strip.background = element_blank(),
			axis.title = element_text(size = rel(1.0)),
			axis.text = element_text(size = rel(0.85), color = "black"),
			legend.position = "bottom",
			legend.title = element_blank(),
			legend.text = element_text(size = rel(0.85)),
			legend.key.width = unit(1.3, "lines"),
			plot.title = element_blank(),
			plot.margin = margin(4, 6, 4, 4)
		)
}

# pch 18 (solid diamond, used for sir_collapsed) renders visually smaller
# than the other pch shapes at the same "size" value -- bump it up so all
# five method markers read as similar visual weight (mirrors figures_simulation.R).
METHOD_POINT_SIZES <- c(
	montoril = 1.6,
	pg_apf = 1.6,
	sir_laplace = 1.6,
	sir_collapsed = 2.2,
	stan = 1.6
)

# Theoretical ceiling on raw (bulk) ESS: S_total = S_per_chain * K_chains.
# The real-data application runs K = 3 chains per method (unlike the
# simulation study's K = 1 -- see ESS_CEILING_SIM in figures_simulation.R),
# per the calibration table:
#   amh_montoril: S = 100,000 * K = 3 -> 300,000
#   pg_apf: S = 20,000 * K = 3 -> 60,000
#   sir_laplace / sir_collapsed / stan: S = 10,000 * K = 3 -> 30,000
ESS_CEILING_REAL <- c(
	montoril = 300000,
	pg_apf = 60000,
	sir_laplace = 30000,
	sir_collapsed = 30000,
	stan = 30000
)

# Proper log-scale minor gridlines (2,3,4,...,9 within each decade -- the
# MATLAB-style log grid), instead of ggplot's default minor_breaks for a log
# scale, which is just the arithmetic midpoint (in log space) between two
# major breaks -- a single, visually meaningless line splitting each decade
# in half (duplicated verbatim from figures_simulation.R -- see header note).
log_minor_breaks <- function(x) {
	lo <- floor(log10(min(x)))
	hi <- ceiling(log10(max(x)))
	breaks <- as.vector(outer(2:9, 10^(lo:hi)))
	return(breaks[breaks >= min(x) & breaks <= max(x)])
}

scale_color_method <- function() {
	return(scale_color_manual(values = METHOD_COLORS, labels = METHOD_LABELS, breaks = METHOD_LEVELS))
}

scale_fill_method <- function() {
	return(scale_fill_manual(values = METHOD_COLORS, labels = METHOD_LABELS, breaks = METHOD_LEVELS))
}

scale_linetype_method <- function() {
	return(scale_linetype_manual(values = METHOD_LINETYPES, labels = METHOD_LABELS, breaks = METHOD_LEVELS))
}

scale_shape_method <- function() {
	return(scale_shape_manual(values = METHOD_SHAPES, labels = METHOD_LABELS, breaks = METHOD_LEVELS))
}

save_figure <- function(plot, filename, width = FULL_WIDTH, height) {
	ggsave(
		filename = file.path(OUTPUT_DIR, filename),
		plot = plot,
		width = width,
		height = height,
		units = "in",
		device = PDF_DEVICE
	)
	return(invisible(NULL))
}


# ---- Data loading -----------------------------------------------------------
#
# method_summaries.rds is a list keyed by method, each element holding
# per-t vectors (length Tt: rhat_theta1, ess_bulk_theta1, theta1_mean,
# theta1_ci_lower/upper, cpo_t, ...) and per-method scalars (rhat_theta1_max,
# ess_sec_theta1_bulk_mean, log_lik, log_cpo, ...). Reshaped here into tidy
# long-format data frames, one per figure family, so the make_fig_*()
# functions below read like ordinary ggplot pipelines.

load_data <- function() {
	y <- readRDS(PATH_DATA)
	Tt <- length(y)
	ms <- readRDS(PATH_SUMMARIES)
	dm <- readRDS(PATH_DELTA_MAX)

	# Real calendar dates for the x-axis, reconstructed from campy's actual
	# ts time base (tsp: start = 1990.000, frequency = 13 four-week periods
	# per year -- Ferland, Latour & Oraichi 2006 / R package tscount).
	# Standard fractional-year -> calendar-date conversion; verified against
	# the original .RData's time() values (identical), giving 28/29-day
	# spacing consistent with the dataset's known 4-week reporting cycle.
	campy_start_year <- 1990
	campy_frequency <- 13
	frac_year <- campy_start_year + (seq_len(Tt) - 1) / campy_frequency
	campy_dates <- as.Date(sprintf("%d-01-01", floor(frac_year))) + round((frac_year - floor(frac_year)) * 365.25)

	# ---- Fit: t, date, method, y, theta1_mean/ci, lambda_mean/ci ----
	fit <- bind_rows(lapply(METHOD_LEVELS, function(m) {
		s <- ms[[m]]
		data.frame(
			t = seq_len(Tt), date = campy_dates, method = m, y = y,
			theta1_mean = s$theta1_mean, theta1_ci_lower = s$theta1_ci_lower, theta1_ci_upper = s$theta1_ci_upper
		)
	})) %>%
		mutate(
			method = factor(method, levels = METHOD_LEVELS),
			lambda_mean = exp(theta1_mean),
			lambda_ci_lower = exp(theta1_ci_lower),
			lambda_ci_upper = exp(theta1_ci_upper)
		)

	# ---- R_hat(t): t, method, parameter (theta_t1/theta_t2), rhat ----
	rhat_t <- bind_rows(lapply(METHOD_LEVELS, function(m) {
		s <- ms[[m]]
		bind_rows(
			data.frame(t = seq_len(Tt), method = m, parameter = "theta_t1", rhat = s$rhat_theta1),
			data.frame(t = seq_len(Tt), method = m, parameter = "theta_t2", rhat = s$rhat_theta2)
		)
	})) %>%
		mutate(
			method = factor(method, levels = METHOD_LEVELS),
			parameter = factor(parameter, levels = PARAM_LEVELS)
		)

	# ---- Raw ESS and ESS/s: method x parameter (6 levels) x {bulk, tail} ----
	# theta_t1/theta_t2 use the already-summarised mean-over-t scalar (the
	# same convention simulation_run.R/figures.R use for a per-parameter
	# scalar summary); theta_01/theta_02/W1/W2 are themselves scalars.
	ess_raw <- bind_rows(lapply(METHOD_LEVELS, function(m) {
		s <- ms[[m]]
		data.frame(
			method = m,
			parameter = PARAM_LEVELS,
			bulk = c(s$ess_bulk_theta1_mean, s$ess_bulk_theta2_mean, s$ess_bulk_theta_01, s$ess_bulk_theta_02, s$ess_bulk_W1, s$ess_bulk_W2),
			tail = c(s$ess_tail_theta1_mean, s$ess_tail_theta2_mean, s$ess_tail_theta_01, s$ess_tail_theta_02, s$ess_tail_W1, s$ess_tail_W2)
		)
	})) %>%
		mutate(method = factor(method, levels = METHOD_LEVELS), parameter = factor(parameter, levels = PARAM_LEVELS)) %>%
		pivot_longer(cols = c(bulk, tail), names_to = "estimator", values_to = "value") %>%
		mutate(estimator = factor(estimator, levels = c("bulk", "tail")))

	ess_sec <- bind_rows(lapply(METHOD_LEVELS, function(m) {
		s <- ms[[m]]
		data.frame(
			method = m,
			parameter = PARAM_LEVELS,
			bulk = c(s$ess_sec_theta1_bulk_mean, s$ess_sec_theta2_bulk_mean, s$ess_sec_theta_01_bulk, s$ess_sec_theta_02_bulk, s$ess_sec_W1_bulk, s$ess_sec_W2_bulk),
			tail = c(s$ess_sec_theta1_tail_mean, s$ess_sec_theta2_tail_mean, s$ess_sec_theta_01_tail, s$ess_sec_theta_02_tail, s$ess_sec_W1_tail, s$ess_sec_W2_tail)
		)
	})) %>%
		mutate(method = factor(method, levels = METHOD_LEVELS), parameter = factor(parameter, levels = PARAM_LEVELS)) %>%
		pivot_longer(cols = c(bulk, tail), names_to = "estimator", values_to = "value") %>%
		mutate(estimator = factor(estimator, levels = c("bulk", "tail")))

	# ---- Delta_max(t): t, parameter (theta_t1/theta_t2), value ----
	delta_max <- bind_rows(
		data.frame(t = seq_len(Tt), parameter = "theta_t1", value = as.numeric(dm$theta1)),
		data.frame(t = seq_len(Tt), parameter = "theta_t2", value = as.numeric(dm$theta2))
	) %>%
		mutate(parameter = factor(parameter, levels = PARAM_LEVELS))

	# W1/W2 don't vary with t -- delta_max.rds$W1/$W2 are already single
	# scalars (see real_data_run.R's aggregate_method()).
	delta_max_W1 <- as.numeric(dm$W1)
	delta_max_W2 <- as.numeric(dm$W2)

	# ---- Total CPU time and fit metrics: one row per method (from summary.csv) ----
	summary_csv <- read.csv(PATH_SUMMARY_CSV, stringsAsFactors = FALSE) %>%
		mutate(method = factor(method, levels = METHOD_LEVELS))

	return(list(y = y, Tt = Tt, fit = fit, rhat_t = rhat_t, ess_raw = ess_raw, ess_sec = ess_sec,
				delta_max = delta_max, delta_max_W1 = delta_max_W1, delta_max_W2 = delta_max_W2,
				summary = summary_csv))
}


# ---- Figure 1: fit -- y_t vs lambda_hat + 95% HPD band, facet by method ----

make_fig_real_fit <- function(fit, methods = METHOD_LEVELS, facet_nrow = 2, facet_ncol = 3) {
	plot_data <- fit %>% filter(method %in% methods)
	obs_data <- plot_data %>% filter(method == methods[1]) %>% select(date, y)

	p <- ggplot() +
		geom_col(
			data = obs_data,
			aes(x = date, y = y),
			fill = "grey80", width = 28, alpha = 0.7
		) +
		geom_ribbon(
			data = plot_data,
			aes(x = date, ymin = lambda_ci_lower, ymax = lambda_ci_upper, fill = method),
			alpha = 0.25
		) +
		geom_line(
			data = plot_data,
			aes(x = date, y = lambda_mean, color = method, linetype = method),
			linewidth = 0.4
		) +
		facet_wrap(~method, nrow = facet_nrow, ncol = facet_ncol, labeller = as_labeller(METHOD_LABELS), scales = "free_x") +
		scale_color_method() +
		scale_fill_method() +
		scale_linetype_method() +
		scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
		labs(x = NULL, y = expression(y[t]~"and"~hat(lambda)[t])) +
		theme_paper() +
		theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

	return(p)
}


# ---- Figure 1b (article): single representative method + Delta_max(t) ----
#
# sir_collapsed is the paper's main methodological contribution and has the
# best log_cpo (Table, summary.csv) -- a principled, non-arbitrary choice
# of "representative" method, not cherry-picking the prettiest curve.
ARTICLE_REPRESENTATIVE_METHOD <- "sir_collapsed"

make_fig_article_real_fit <- function(fit, representative = ARTICLE_REPRESENTATIVE_METHOD) {
	plot_data <- fit %>% filter(method == representative)

	p_fit <- ggplot() +
		geom_col(
			data = plot_data,
			aes(x = date, y = y, fill = "Observed"),
			width = 28, alpha = 0.7
		) +
		geom_ribbon(
			data = plot_data,
			aes(x = date, ymin = lambda_ci_lower, ymax = lambda_ci_upper, fill = "95% HDI"),
			alpha = 0.25
		) +
		geom_line(
			data = plot_data,
			aes(x = date, y = lambda_mean, color = "Posterior mean"),
			linewidth = 0.45
		) +
		scale_fill_manual(
			values = c("Observed" = "grey80", "95% HDI" = METHOD_COLORS[[representative]]),
			breaks = c("Observed", "95% HDI")
		) +
		scale_color_manual(
			values = c("Posterior mean" = METHOD_COLORS[[representative]]),
			breaks = "Posterior mean"
		) +
		scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
		labs(x = NULL, y = expression(y[t]~"and"~hat(lambda)[t])) +
		theme_paper() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		guides(
			fill = guide_legend(order = 1, override.aes = list(alpha = c(0.7, 0.35))),
			color = guide_legend(order = 2)
		)

	return(p_fit)
}


# ---- Figure 1c (article, alternative): all 5 methods overlaid (ribbon +
# line per method), mirroring make_fig_overlay_hpd() from figures.R -- test
# variant against the single-representative-method version above. ----

make_fig_article_real_fit_overlay <- function(fit, methods = METHOD_LEVELS, font_size = 7) {
	plot_data <- fit %>% filter(method %in% methods)
	obs_data <- plot_data %>% filter(method == methods[1]) %>% select(date, y)

	p <- ggplot() +
		geom_col(
			data = obs_data,
			aes(x = date, y = y),
			fill = "grey55", width = 28, alpha = 0.7
		) +
		geom_line(
			data = plot_data,
			aes(x = date, y = lambda_mean, color = method, linetype = method),
			linewidth = 0.4
		) +
		scale_color_method() +
		scale_linetype_method() +
		scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
		scale_y_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60), expand = expansion(mult = c(0, 0))) +
		labs(x = NULL, y = expression(y[t]~"and"~hat(lambda)[t])) +
		theme_paper() +
		theme(
			text = element_text(size = font_size),
			axis.text.x = element_text(angle = 0, hjust = 0.5),
			legend.margin = margin(0, 0, 0, 0),
			legend.box.margin = margin(0, 0, 0, 0),
			legend.box.spacing = unit(2, "pt"),
			legend.spacing.x = unit(4, "pt"),
			legend.key.size = unit(0.7, "lines")
		) +
		guides(
			color = guide_legend(nrow = 1, byrow = TRUE),
			linetype = guide_legend(nrow = 1, byrow = TRUE)
		)

	return(p)
}


# ---- Figure 2: R_hat(t), facet by parameter (theta_t1, theta_t2) ----------

make_fig_real_rhat <- function(rhat_t) {
	p <- ggplot(rhat_t, aes(x = t, y = rhat, color = method, linetype = method)) +
		geom_hline(yintercept = 1.01, linetype = "dashed", color = "grey40", linewidth = 0.4) +
		geom_line(linewidth = 0.35) +
		facet_wrap(~parameter, nrow = 1, ncol = 2, labeller = as_labeller(PARAM_LABELS, label_parsed)) +
		scale_color_method() +
		scale_linetype_method() +
		labs(x = "t", y = expression(hat(R))) +
		theme_paper() +
		guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))

	return(p)
}


# ---- Shared helper: ESS-style bar chart (bulk + tail), facet by parameter -
#
# Bar, not boxplot: unlike the simulation study (R=200 replicas per cell),
# the real-data application has exactly one value per (method, parameter,
# estimator) -- there is no distribution to box.

make_fig_real_ess_style_bar <- function(plot_data, y_label, breaks_pow10 = -2:5, facet_nrow = 3, facet_ncol = 2) {
	p <- ggplot(plot_data, aes(x = method, y = value, fill = method, alpha = estimator)) +
		geom_col(position = position_dodge(width = 0.75), width = 0.7, color = "grey20", linewidth = 0.1) +
		facet_wrap(~parameter, nrow = facet_nrow, ncol = facet_ncol, labeller = as_labeller(PARAM_LABELS, label_parsed), scales = "free_x") +
		scale_fill_method() +
		scale_alpha_manual(values = c(bulk = 1, tail = 0.5), labels = c(bulk = "Bulk", tail = "Tail")) +
		scale_x_discrete(labels = METHOD_LABELS) +
		scale_y_log10(
			breaks = 10^breaks_pow10,
			labels = label_number(big.mark = ",", drop0trailing = TRUE)
		) +
		labs(x = NULL, y = y_label) +
		theme_paper() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
			legend.position = "bottom"
		) +
		guides(fill = "none", alpha = guide_legend(title = NULL, override.aes = list(fill = "grey30")))

	return(p)
}


# ---- Figure 3: raw ESS bar chart, facet by parameter -----------------------

make_fig_real_ess <- function(ess_raw, facet_nrow = 3, facet_ncol = 2) {
	return(make_fig_real_ess_style_bar(ess_raw, y_label = "ESS", facet_nrow = facet_nrow, facet_ncol = facet_ncol))
}


# ---- Figure 4: ESS/s bar chart, facet by parameter -------------------------

make_fig_real_ess_sec <- function(ess_sec, facet_nrow = 3, facet_ncol = 2) {
	return(make_fig_real_ess_style_bar(ess_sec, y_label = "ESS/s", facet_nrow = facet_nrow, facet_ncol = facet_ncol))
}


# ---- Figure 5: cross-method agreement Delta_max(t) -------------------------

make_fig_real_delta_max <- function(delta_max) {
	p <- ggplot(delta_max, aes(x = t, y = value)) +
		geom_line(color = "grey20", linewidth = 0.4) +
		facet_wrap(~parameter, nrow = 1, ncol = 2, labeller = as_labeller(PARAM_LABELS, label_parsed), scales = "free_y") +
		labs(x = "t", y = expression(Delta[max](t))) +
		theme_paper()

	return(p)
}


# ---- Figure 6: total CPU time by method ------------------------------------

make_fig_real_total_time <- function(summary_df) {
	p <- ggplot(summary_df, aes(x = method, y = total_time, fill = method)) +
		geom_col(width = 0.6, color = "grey20", linewidth = 0.15) +
		scale_fill_method() +
		scale_x_discrete(labels = METHOD_LABELS) +
		scale_y_log10(
			labels = label_number(big.mark = ",", drop0trailing = TRUE),
			breaks = 10^(-1:4)
		) +
		labs(x = NULL, y = "Total CPU time (s)") +
		theme_paper() +
		theme(legend.position = "none")

	return(p)
}


# ---- Figure 7: log_lik and log_cpo, bar chart with an explicit per-facet
# floor (not zero) -----------------------------------------------------------
#
# Cross-method differences here are small (0.5-2.2 nats) relative to the
# distance to zero (~-333 for log_lik, ~-388 for log_cpo) -- a bar from a
# 0 baseline would make every bar look identical. Same fix as
# make_fig_total_time() in figures.R (geom_rect() with an explicit floor
# instead of geom_col()'s implicit baseline), applied per facet since
# log_lik and log_cpo sit on very different scales.

make_fig_real_fit_metrics <- function(summary_df, floor_pad = 0.15) {
	plot_data <- summary_df %>%
		select(method, log_lik, log_cpo) %>%
		pivot_longer(cols = c(log_lik, log_cpo), names_to = "metric", values_to = "value") %>%
		mutate(
			metric = factor(metric, levels = c("log_lik", "log_cpo"),
							 labels = c("log-likelihood~(post.~mean)", "log-CPO~(leave-one-out)")),
			method_index = as.numeric(factor(method, levels = METHOD_LEVELS))
		) %>%
		group_by(metric) %>%
		mutate(floor = min(value) - floor_pad * diff(range(value))) %>%
		ungroup()

	p <- ggplot(plot_data) +
		geom_rect(
			aes(xmin = method_index - 0.35, xmax = method_index + 0.35, ymin = floor, ymax = value, fill = method),
			color = "grey20", linewidth = 0.15
		) +
		facet_wrap(~metric, nrow = 1, ncol = 2, labeller = label_parsed, scales = "free") +
		scale_fill_method() +
		scale_x_continuous(breaks = seq_along(METHOD_LEVELS), labels = METHOD_LABELS[METHOD_LEVELS]) +
		labs(x = NULL, y = NULL) +
		theme_paper() +
		theme(legend.position = "none")

	return(p)
}


# ---- Shared helper: pull the legend grob out of a ggplot, and stack several
# plots into ONE gtable with truly equal panel sizes -------------------------
#
# gridExtra::arrangeGrob(heights = ...) treats each plot as an opaque box and
# splits total height by the given weights -- since rows with less axis/strip
# overhead (e.g. a blanked x-axis) end up with a visibly BIGGER panel than
# rows with more overhead for the same weight, matching panel sizes required
# fragile trial-and-error weight ratios. gtable::gtable_rbind() instead merges
# the plots' gtables at the grid level, so each plot's "panel" row keeps its
# default 1-null sizing -- when the combined gtable is finally drawn, all
# panel rows compete equally for the same leftover space and end up exactly
# the same height, regardless of how much fixed-size axis/strip content
# surrounds them.

extract_legend <- function(p) {
	g <- ggplotGrob(p)
	idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
	return(g$grobs[[idx]])
}

stack_plots_equal_panels <- function(plots, legend = NULL, legend_gap_pt = 0) {
	grobs <- lapply(plots, ggplotGrob)

	max_widths <- grobs[[1]]$widths
	for (g in grobs[-1]) {
		max_widths <- grid::unit.pmax(max_widths, g$widths)
	}
	for (i in seq_along(grobs)) {
		grobs[[i]]$widths <- max_widths
	}

	combined <- grobs[[1]]
	for (g in grobs[-1]) {
		combined <- rbind(combined, g)
	}

	if (!is.null(legend)) {
		combined <- gtable::gtable_add_rows(combined, heights = grobHeight(legend) + unit(legend_gap_pt, "pt"))
		combined <- gtable::gtable_add_grob(combined, legend, t = nrow(combined), l = 1, r = ncol(combined))
	}

	return(combined)
}

# ---- Shared helper: ESS-style point chart (bulk + tail), facet by parameter
#
# Point, not bar: a single value per (method, parameter, estimator) doesn't
# need geom_col's visual weight (fill + border) -- a colored/shaped point
# per method reads just as clearly and far less cluttered, especially once
# bulk and tail are both shown (distinguished by alpha, dodged apart).

make_fig_real_ess_style_point <- function(plot_data, y_label, breaks_pow10, facet_nrow = 1, facet_ncol = 4, free_y = FALSE, limits = c(NA, NA)) {
	p <- ggplot(plot_data, aes(x = method, y = value, color = method, shape = method, size = method, alpha = estimator, group = estimator)) +
		geom_point(position = position_dodge(width = 0.5)) +
		facet_wrap(
			~parameter,
			nrow = facet_nrow,
			ncol = facet_ncol,
			scales = if (free_y) "free_y" else "fixed",
			labeller = as_labeller(PARAM_LABELS, label_parsed)
		) +
		scale_color_method() +
		scale_shape_method() +
		scale_size_manual(values = METHOD_POINT_SIZES, breaks = METHOD_LEVELS, guide = "none") +
		scale_alpha_manual(values = c(bulk = 1, tail = 0.45), labels = c(bulk = "Bulk", tail = "Tail")) +
		scale_x_discrete(labels = METHOD_LABELS) +
		scale_y_log10(
			breaks = 10^breaks_pow10,
			minor_breaks = log_minor_breaks,
			labels = label_number(big.mark = ",", drop0trailing = TRUE),
			limits = limits
		) +
		labs(x = NULL, y = y_label) +
		theme_paper() +
		theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1)) +
		guides(
			color = guide_legend(nrow = 1, override.aes = list(size = unname(METHOD_POINT_SIZES[METHOD_LEVELS]))),
			shape = guide_legend(nrow = 1, override.aes = list(size = unname(METHOD_POINT_SIZES[METHOD_LEVELS]))),
			alpha = guide_legend(title = NULL, override.aes = list(color = "grey30", shape = 16, size = 1.6))
		)

	return(p)
}


# ---- Figure 9 (article): raw ESS + ESS/s, theta_t1/theta_t2/W1/W2 ---------

make_fig_article_real_ess <- function(
	ess_raw, ess_sec,
	params = c("theta_t1", "theta_t2", "W1", "W2"),
	panel_spacing_pt = 2,
	row_margin = margin(1, 6, 1, 4),
	legend_gap_pt = 0,
	font_size = 7
) {
	facet_ncol <- length(params)
	shrink_text <- theme(text = element_text(size = font_size))
	minor_grid <- theme(panel.grid.minor = element_line(color = "grey80", linewidth = 0.12, linetype = "dashed"))

	plot_data_raw <- ess_raw %>% filter(parameter %in% params) %>% mutate(parameter = factor(parameter, levels = params))
	plot_data_sec <- ess_sec %>% filter(parameter %in% params) %>% mutate(parameter = factor(parameter, levels = params))

	# Theoretical ceiling on raw ESS (S_total = S_per_chain * K_chains,
	# K = 3 here), one horizontal dashed segment per method, centered on
	# that method's x category (spans most of its width; not dodged by
	# estimator -- the ceiling applies to bulk and tail alike).
	ceiling_data <- data.frame(method = factor(METHOD_LEVELS, levels = METHOD_LEVELS)) %>%
		mutate(method_index = as.numeric(method), ceiling = ESS_CEILING_REAL[as.character(method)])

	# Row 1 (raw ESS) and row 2 (ESS/s) share the same x variable (method),
	# so row 1's x-axis text/ticks/title are dropped -- row 2's are enough.
	# Both rows use one shared (non-free) y-axis across all 4 parameters, so
	# ggplot only draws y-axis text on the leftmost column.
	p_raw <- make_fig_real_ess_style_point(plot_data_raw, y_label = "ESS", breaks_pow10 = 2:6, facet_nrow = 1, facet_ncol = facet_ncol) +
		geom_segment(
			data = ceiling_data,
			aes(x = method_index - 0.35, xend = method_index + 0.35, y = ceiling, yend = ceiling),
			inherit.aes = FALSE,
			linetype = "dashed",
			color = "grey30",
			linewidth = 0.4
		) +
		shrink_text +
		minor_grid +
		theme(
			legend.position = "none",
			panel.spacing = unit(panel_spacing_pt, "pt"),
			plot.margin = row_margin,
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			axis.title.x = element_blank()
		)

	# Row 2 skips its own facet strip titles -- row 1's, directly above,
	# already name each parameter.
	p_sec_full <- make_fig_real_ess_style_point(plot_data_sec, y_label = "ESS/s", breaks_pow10 = 0:4, facet_nrow = 1, facet_ncol = facet_ncol, limits = c(1, NA)) +
		shrink_text +
		minor_grid +
		theme(
			panel.spacing = unit(panel_spacing_pt, "pt"),
			plot.margin = row_margin,
			strip.text = element_blank(),
			strip.background = element_blank(),
			axis.title.x = element_blank(),
			legend.margin = margin(0, 0, 0, 0),
			legend.box.margin = margin(0, 0, 0, 0),
			legend.box.spacing = unit(2, "pt"),
			legend.spacing.x = unit(4, "pt"),
			legend.key.size = unit(0.7, "lines")
		)
	p_sec <- p_sec_full + theme(legend.position = "none")

	# Shared bottom legend: built from synthetic dummy data (decoupled from
	# the real panels) so a "Theoretical max" key -- same grey dashed style
	# as the ceiling segments in row 1 -- can be appended after the method
	# keys. Methods get a "blank" linetype (they have no line in the real
	# plot, only points); the Bulk/Tail alpha legend is reproduced
	# separately since it isn't tied to method identity.
	legend_data_method <- data.frame(x = 1, y = 1, method = factor(c(METHOD_LEVELS, "max_ess"), levels = c(METHOD_LEVELS, "max_ess")))
	legend_data_estimator <- data.frame(x = 1, y = 1, estimator = factor(c("bulk", "tail"), levels = c("bulk", "tail")))

	legend_plot <- ggplot() +
		geom_line(data = legend_data_method, aes(x = x, y = y, color = method, linetype = method)) +
		geom_point(data = legend_data_method, aes(x = x, y = y, color = method, shape = method, size = method)) +
		geom_point(data = legend_data_estimator, aes(x = x, y = y, alpha = estimator), color = "grey30", shape = 16, size = 1.6) +
		scale_color_manual(values = c(METHOD_COLORS, max_ess = "grey30"), labels = c(METHOD_LABELS, max_ess = "Theoretical max"), breaks = c(METHOD_LEVELS, "max_ess")) +
		scale_shape_manual(values = c(METHOD_SHAPES, max_ess = NA), labels = c(METHOD_LABELS, max_ess = "Theoretical max"), breaks = c(METHOD_LEVELS, "max_ess")) +
		scale_linetype_manual(
			values = c(setNames(rep("blank", length(METHOD_LEVELS)), METHOD_LEVELS), max_ess = "dashed"),
			labels = c(METHOD_LABELS, max_ess = "Theoretical max"),
			breaks = c(METHOD_LEVELS, "max_ess")
		) +
		scale_size_manual(values = c(METHOD_POINT_SIZES, max_ess = 1.6), guide = "none") +
		scale_alpha_manual(values = c(bulk = 1, tail = 0.45), labels = c(bulk = "Bulk", tail = "Tail")) +
		theme_paper() +
		shrink_text +
		theme(
			legend.margin = margin(0, 0, 0, 0),
			legend.box.margin = margin(0, 0, 0, 0),
			legend.box.spacing = unit(2, "pt"),
			legend.spacing.x = unit(4, "pt"),
			legend.key.size = unit(0.7, "lines")
		) +
		guides(
			color = guide_legend(nrow = 1, override.aes = list(size = unname(c(METHOD_POINT_SIZES[METHOD_LEVELS], max_ess = 1.6)))),
			shape = guide_legend(nrow = 1, override.aes = list(size = unname(c(METHOD_POINT_SIZES[METHOD_LEVELS], max_ess = 1.6)))),
			linetype = guide_legend(nrow = 1),
			alpha = guide_legend(title = NULL, override.aes = list(color = "grey30", shape = 16, size = 1.6))
		)
	legend <- extract_legend(legend_plot)

	combined <- stack_plots_equal_panels(list(p_raw, p_sec), legend = legend, legend_gap_pt = legend_gap_pt)

	return(combined)
}


# ---- Article table: log_lik/log_cpo per method + Delta_max summary -------
#
# Replaces a figure with a \input{}-able LaTeX fragment (booktabs style,
# matching metricas_theta1.tex's convention) -- generated straight from
# summary.csv/delta_max.rds so the numbers in the manuscript can never
# drift from the actual results via manual transcription.
#
# Deliberately cautious wording is left to the manuscript text, not baked
# in here: with a single real dataset (no replicas), differences this small
# (log_lik span 0.51 nats, log_cpo span 2.23 nats) should be reported
# as-is, not framed as a statistically established ranking.

write_real_data_results_table <- function(summary_df, delta_max_W1, delta_max_W2, delta_max, filename = "real_data_results_table.tex") {
	tbl <- summary_df %>%
		mutate(method_label = METHOD_LABELS[as.character(method)]) %>%
		arrange(match(method, METHOD_LEVELS))

	rows <- sprintf("%s & %.2f & %.2f \\\\", tbl$method_label, tbl$log_lik, tbl$log_cpo)

	max_delta_theta1 <- max(delta_max$value[delta_max$parameter == "theta_t1"])
	max_delta_theta2 <- max(delta_max$value[delta_max$parameter == "theta_t2"])

	# Plain %.4f rounds anything below 5e-5 to "0.0000", hiding the value
	# (Delta_max(W2) is ~1e-6) -- switch to scientific notation for small
	# magnitudes so nothing silently disappears.
	fmt_delta <- function(x) if (abs(x) < 1e-3) formatC(x, format = "e", digits = 2) else sprintf("%.4f", x)

	lines <- c(
		"% Auto-generated by figures_real_data.R -- do not edit by hand.",
		"\\begin{table}[h]",
		"\\centering",
		"\\begin{tabular}{lrr}",
		"\\toprule",
		"Method & log-likelihood & log-CPO \\\\",
		"\\midrule",
		rows,
		"\\bottomrule",
		"\\end{tabular}",
		sprintf("\\caption{Log-likelihood at the posterior mean and leave-one-out log-CPO (Eq.~34) for the \\texttt{campy} application. Both are close across methods (spans of %.2f and %.2f nats, respectively); with a single real dataset, this ranking is suggestive rather than statistically established.}",
				diff(range(tbl$log_lik)), diff(range(tbl$log_cpo))),
		"\\end{table}",
		"",
		"% Cross-method agreement summary (Eq. 32) -- W1/W2 have no t dimension.",
		sprintf("Cross-method agreement was tight throughout: $\\max_t \\Delta_{\\max}(\\theta_{t1}) = %s$, $\\max_t \\Delta_{\\max}(\\theta_{t2}) = %s$, $\\Delta_{\\max}(W_1) = %s$, $\\Delta_{\\max}(W_2) = %s$.",
				fmt_delta(max_delta_theta1), fmt_delta(max_delta_theta2), fmt_delta(delta_max_W1), fmt_delta(delta_max_W2))
	)

	writeLines(lines, file.path(OUTPUT_DIR, filename))
	return(invisible(NULL))
}


# ---- Run everything ---------------------------------------------------------

data <- load_data()

save_figure(
	plot = make_fig_real_fit(data$fit),
	filename = "real_fit.pdf",
	height = 4.4
)

save_figure(
	plot = make_fig_real_rhat(data$rhat_t),
	filename = "real_rhat.pdf",
	height = 3.0
)

save_figure(
	plot = make_fig_real_ess(data$ess_raw),
	filename = "real_ess.pdf",
	height = 7.0
)

save_figure(
	plot = make_fig_real_ess_sec(data$ess_sec),
	filename = "real_ess_sec.pdf",
	height = 7.0
)

save_figure(
	plot = make_fig_real_delta_max(data$delta_max),
	filename = "real_delta_max.pdf",
	height = 2.8
)

save_figure(
	plot = make_fig_real_total_time(data$summary),
	filename = "real_total_time.pdf",
	height = 3.0
)

save_figure(
	plot = make_fig_real_fit_metrics(data$summary),
	filename = "real_fit_metrics.pdf",
	height = 3.2
)

save_figure(
	plot = make_fig_article_real_fit_overlay(data$fit),
	filename = "article_real_data_fit.pdf",
	height = 1.8
)

save_figure(
	plot = make_fig_article_real_ess(data$ess_raw, data$ess_sec),
	filename = "article_real_data_ess.pdf",
	height = 2.8
)

write_real_data_results_table(data$summary, data$delta_max_W1, data$delta_max_W2, data$delta_max)

cat("All figures written to", OUTPUT_DIR, "\n")
