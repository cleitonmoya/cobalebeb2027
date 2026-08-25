# ============================================================================
# figures_simulation.R
#
# Generates all simulation figures (RMSE, coverage, HPD overlay, total execution
# time, ESS/s) for the Cobal/EBEB 2027 manuscript from the three production
# simulation summary CSVs. Single shared style (colors, linetypes, shapes,
# method labels, ggplot theme) so all figures are visually consistent.
#
# Figures produced (see FIGURES_TO_RUN at the bottom):
#   rmse_boxplot.pdf   - RMSE per replica, boxplot, facet by function
#   coverage.pdf        - empirical 95% coverage vs T, facet by function
#   overlay_hpd.pdf     - all methods overlaid, mean + 95% HPD, facet by function
#   total_time_sum.pdf - total summed CPU time by method and T
#   boxplot_ess.pdf     - ESS/s boxplot by method, facet by parameter
# ============================================================================

setwd(dirname(this.path::this.path()))

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(gridExtra)
library(grid)


# ---- Shared configuration (edit here to affect all figures) --------------

DATA_DIR <- "../data"
RESULTS_DIR <- "../results"
OUTPUT_DIR <- "plots/simulation"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

PATH_AGGREGATED <- file.path(RESULTS_DIR, "summary_aggregated.csv")
PATH_BY_T <- file.path(RESULTS_DIR, "summary_by_t.csv")
PATH_REPLICAS <- file.path(RESULTS_DIR, "summary_replicas.csv")
PATH_POINTWISE <- file.path(RESULTS_DIR, "summary_pointwise_coverage.csv")

# Method identifiers, plotting order, and English display labels.
METHOD_LEVELS <- c("montoril", "pg_apf", "sir_laplace", "sir_collapsed", "stan")

METHOD_LABELS <- c(
	montoril = "AM",
	pg_apf = "PG-AS",
	sir_laplace = "SIR",
	sir_collapsed = "SIR-Collapsed",
	stan = "Stan"
)

# Colors (Okabe-Ito colorblind-safe subset) + linetypes + point shapes, one
# per method. The linetype/shape redundancy is what keeps the figures
# legible if printed in grayscale -- color alone is never load-bearing.
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

# Function (latent trend) identifiers and English display labels.
FUNCTION_LEVELS <- c("constant", "linear", "quadratic", "sinusoidal")

FUNCTION_LABELS <- c(
	constant = "Steps",
	linear = "p-Linear",
	quadratic = "p-Quadratic",
	sinusoidal = "Sinusoidal"
)

# Figure sizing (inches) -- standardized across the article. FULL_WIDTH
# matches the manuscript's text width (no two-column layout).
FULL_WIDTH <- 6
BASE_FONT_SIZE <- 9

PDF_DEVICE <- cairo_pdf


# ---- Shared ggplot theme ---------------------------------------------------

theme_paper <- function(base_size = BASE_FONT_SIZE) {
	theme_minimal(base_size = base_size) +
		theme(
			panel.grid.minor = element_blank(),
			panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
			panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
			strip.text = element_text(face = "bold", size = rel(0.95)),
			strip.background = element_rect(fill = "grey93", color = NA),
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

# Common scale helpers so every figure that maps `method` uses the exact
# same colors / linetypes / shapes / legend labels.
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

# Helper to save a plot with the standardized figure size.
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

load_data <- function() {
	aggregated <- read.csv(PATH_AGGREGATED, stringsAsFactors = FALSE) %>%
		mutate(
			method = factor(method, levels = METHOD_LEVELS),
			f = factor(f, levels = FUNCTION_LEVELS)
		)

	by_t <- read.csv(PATH_BY_T, stringsAsFactors= FALSE) %>%
		mutate(
			method = factor(method, levels = METHOD_LEVELS),
			f = factor(f, levels = FUNCTION_LEVELS)
		)

	replicas <- read.csv(PATH_REPLICAS, stringsAsFactors = FALSE) %>%
		mutate(
			method = factor(method, levels = METHOD_LEVELS),
			f = factor(f, levels = FUNCTION_LEVELS)
		)

	pointwise <- read.csv(PATH_POINTWISE, stringsAsFactors = FALSE) %>%
		mutate(
			method = factor(method, levels = METHOD_LEVELS),
			f = factor(f, levels = FUNCTION_LEVELS)
		)

	return(list(aggregated = aggregated, by_t = by_t, replicas = replicas, pointwise = pointwise))
}


# ---- Figure 1: RMSE per replica, boxplot, facet by function --------------

make_fig_rmse_boxplot <- function(replicas, functions = FUNCTION_LEVELS, facet_nrow = 2, facet_ncol = 2) {
	plot_data <- replicas %>%
		filter(f %in% functions) %>%
		mutate(Tt_factor = factor(Tt))

	p <- ggplot(plot_data, aes(x = Tt_factor, y = rmse_theta1, fill = method)) +
		geom_boxplot(
			notch = FALSE,
			outlier.shape = 21,
			outlier.fill = NA,
			outlier.stroke = 0.1,
			outlier.size = 0.5,
			outlier.alpha = 1,
			linewidth = 0.1,
			position = position_dodge2(padding = 0.15)
		) +
		facet_wrap(~f, nrow = facet_nrow, ncol = facet_ncol, labeller = as_labeller(FUNCTION_LABELS)) +
		scale_fill_method() +
		labs(x = "T", y = expression(RMSE~of~theta[t1])) +
		theme_paper() +
		guides(fill = guide_legend(nrow = 1))

	return(p)
}


# ---- Figure 2: empirical 95% coverage vs T, facet by function -------------
make_fig_coverage <- function(aggregated, functions = FUNCTION_LEVELS, facet_nrow = 2, facet_ncol = 2) {
	plot_data <- aggregated %>%
		filter(f %in% functions) %>%
		mutate(Tt_factor = factor(Tt))

	p <- ggplot(plot_data, aes(x = Tt_factor, y = coverage, color = method, group = method)) +
		geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40", linewidth = 0.4) +
		geom_errorbar(
			aes(ymin = coverage_wilson_lower, ymax = coverage_wilson_upper),
			position = position_dodge(width = 0.6),
			width = 0.35,
			linewidth = 0.35
		) +
		geom_point(
			aes(shape = method),
			position = position_dodge(width = 0.6),
			size = 1.6
		) +
		#geom_line(
		#	aes(x = as.integer(Tt_factor), linetype = method),
	    #		position = position_dodge(width = 0.6),
	#		linewidth = 0.35
	#	) +
		facet_wrap(~f, nrow = facet_nrow, ncol = facet_ncol, labeller = as_labeller(FUNCTION_LABELS)) +
		scale_color_method() +
		scale_shape_method() +
		scale_linetype_method() +
		labs(x = "T", y = expression(Empirical~95*"%"~coverage~of~theta[1])) +
		theme_paper() +
		guides(
			color = guide_legend(nrow = 1),
			shape = guide_legend(nrow = 1),
			linetype = guide_legend(nrow = 1)
		)

	return(p)
}


# ---- Figure 3: all methods overlaid, mean + 95% HPD, facet by function ----

load_true_curve <- function(Tt_selected, replica = 1) {
	true_curve <- lapply(FUNCTION_LEVELS, function(fn) {
		path <- file.path(DATA_DIR, "simulated", sprintf("%s_%d_%d.rds", fn, Tt_selected, replica))
		sim <- readRDS(path)
		data.frame(f = fn, t = seq_along(sim$theta1), true_theta1 = sim$theta1)
	}) %>%
		bind_rows() %>%
		mutate(f = factor(f, levels = FUNCTION_LEVELS))

	return(true_curve)
}

make_fig_overlay_hpd <- function(by_t, Tt_selected = 1600, facet_nrow = 2, facet_ncol = 2) {
	plot_data <- by_t %>%
		filter(Tt == Tt_selected)

	true_curve <- load_true_curve(Tt_selected)

	p <- ggplot() +
		geom_ribbon(
			data = plot_data,
			aes(x = t, ymin = band_lower, ymax = band_upper, fill = method),
			alpha = 0.2
		) +
		geom_ribbon(
			data = true_curve,
			aes(x = t, ymin = true_theta1, ymax = true_theta1, fill = "true"),
			alpha = 0.2
		) +
		geom_line(
			data = plot_data,
			aes(x = t, y = post_mean_agg, color = method, linetype = method),
			linewidth = 0.35
		) +
		geom_line(
			data = true_curve,
			aes(x = t, y = true_theta1, color = "true", linetype = "true"),
			linewidth = 0.35
		) +
		facet_wrap(~f, nrow = facet_nrow, ncol = facet_ncol, labeller = as_labeller(FUNCTION_LABELS)) +
		scale_color_manual(
			values = c(METHOD_COLORS, true = "black"),
			labels = c(METHOD_LABELS, true = "True"),
			breaks = c(METHOD_LEVELS, "true")
		) +
		scale_fill_manual(
			values = c(METHOD_COLORS, true = "white"),
			labels = c(METHOD_LABELS, true = "True"),
			breaks = c(METHOD_LEVELS, "true")
		) +
		scale_linetype_manual(
			values = c(METHOD_LINETYPES, true = "solid"),
			labels = c(METHOD_LABELS, true = "True"),
			breaks = c(METHOD_LEVELS, "true")
		) +
		scale_x_continuous(labels = function(x) prefix_first_label(as.character(x), "(t)")) +
		labs(x = NULL, y = expression(theta[t1])) +
		theme_paper() +
		theme(axis.title.x = element_blank()) +
		guides(
			color = guide_legend(nrow = 1, byrow = TRUE),
			fill = guide_legend(nrow = 1, byrow = TRUE),
			linetype = guide_legend(nrow = 1, byrow = TRUE)
		)

	return(p)
}


# ---- Figure 4: T CPU Time ---------------------

make_fig_total_time <- function(replicas, floor_hours = 0.1) {
	n_methods <- length(METHOD_LEVELS)
	group_width <- 0.8
	bar_width <- group_width / n_methods

	plot_data <- replicas %>%
		group_by(method, Tt) %>%
		summarise(total_hours = sum(elapsed_time) / 3600, .groups = "drop") %>%
		mutate(
			Tt_index = match(Tt, sort(unique(Tt))),
			method_index = match(method, METHOD_LEVELS),
			slot_center = Tt_index + (method_index - (n_methods + 1) / 2) * bar_width,
			xmin = slot_center - bar_width * 0.45,
			xmax = slot_center + bar_width * 0.45,
			ymin = floor_hours,
			ymax = total_hours
		)

	Tt_breaks <- sort(unique(plot_data$Tt))

	# NOTE: geom_col()'s default 0-baseline is expressed in *transformed*
	# coordinates once scale_y_log10() is applied, which silently becomes
	# the original-scale value 1 (log10(1) = 0) instead of true zero -- bars
	# just above 1 collapse to a sliver. geom_rect() with an explicit,
	# log-safe floor (ymin = floor_hours) avoids that pitfall entirely.
	p <- ggplot(plot_data) +
		geom_rect(
			aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = method),
			color = "grey20",
			linewidth = 0.15
		) +
		scale_fill_method() +
		scale_x_continuous(breaks = seq_along(Tt_breaks), labels = as.character(Tt_breaks)) +
		scale_y_log10(
			labels = label_number(big.mark = ",", drop0trailing = TRUE),
			breaks = 10^(-1:4),
			limits = c(floor_hours, 20000),
			expand = expansion(mult = c(0, 0.05))
		) +
		labs(x = "T", y = "Total CPU time (hours)") +
		theme_paper() +
		theme(
			legend.key.size = unit(0.6, "lines")
		) +
		guides(fill = guide_legend(nrow = 1))

	return(p)
}


# ---- Shared helper: ESS-style boxplot by method, facet by parameter -------

PARAM_LEVELS <- c("theta_t1", "theta_t2", "theta_01", "theta_02", "W1", "W2")

PARAM_LABELS <- c(
	theta_t1 = "theta[t1]",
	theta_t2 = "theta[t2]",
	theta_01 = 'theta["01"]',
	theta_02 = 'theta["02"]',
	W1 = "W[1]",
	W2 = "W[2]"
)

make_fig_ess_style_boxplot <- function(plot_data, y_label, breaks_pow10 = -2:3, facet_nrow = 3, facet_ncol = 2) {
	p <- ggplot(plot_data, aes(x = method, y = value, fill = method)) +
		geom_boxplot(
			notch = FALSE,
			outlier.shape = 21,
			outlier.fill = NA,
			outlier.stroke = 0.1,
			outlier.size = 0.4,
			outlier.alpha = 1,
			linewidth = 0.1
		) +
		facet_wrap(
			~parameter,
			nrow = facet_nrow,
			ncol = facet_ncol,
			labeller = as_labeller(PARAM_LABELS, label_parsed)
		) +
		scale_fill_method() +
		scale_x_discrete(labels = METHOD_LABELS) +
		scale_y_log10(
			breaks = 10^breaks_pow10,
			labels = function(x) {
				ifelse(x == 1, "1.00", label_number(big.mark = ",", drop0trailing = TRUE)(x))
			}
		) +
		labs(x = NULL, y = y_label) +
		theme_paper() +
		theme(
			axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
			legend.position = "none"
		)

	return(p)
}


# ---- Figure 5: ESS/s boxplot by method and T, facet by parameter ----------
# Broken down by T (rather than pooling all T together) to see whether ESS/s
# varies with T -- same x = T / dodge = method / facet = parameter layout as
# the RMSE boxplot (Figure 1), covering all 6 parameters.

make_fig_ess_boxplot <- function(replicas) {
	plot_data <- replicas %>%
		select(
			method, Tt,
			theta_t1 = ess_sec_theta1_mean,
			theta_t2 = ess_sec_theta2_mean,
			theta_01 = ess_sec_theta_01,
			theta_02 = ess_sec_theta_02,
			W1 = ess_sec_W1,
			W2 = ess_sec_W2
		) %>%
		pivot_longer(
			cols = -c(method, Tt),
			names_to = "parameter",
			values_to = "value"
		) %>%
		mutate(
			parameter = factor(parameter, levels = PARAM_LEVELS),
			Tt_factor = factor(Tt)
		)

	p <- ggplot(plot_data, aes(x = Tt_factor, y = value, fill = method)) +
		geom_boxplot(
			notch = FALSE,
			outlier.shape = 21,
			outlier.fill = NA,
			outlier.stroke = 0.1,
			outlier.size = 0.4,
			outlier.alpha = 1,
			linewidth = 0.1,
			position = position_dodge2(padding = 0.15)
		) +
		facet_wrap(
			~parameter,
			nrow = 3,
			ncol = 2,
			labeller = as_labeller(PARAM_LABELS, label_parsed)
		) +
		scale_fill_method() +
		scale_y_log10(
			breaks = 10^(-2:3),
			labels = function(x) {
				ifelse(x == 1, "1.00", label_number(big.mark = ",", drop0trailing = TRUE)(x))
			}
		) +
		labs(x = "T", y = "ESS/s") +
		theme_paper() +
		guides(fill = guide_legend(nrow = 1))

	return(p)
}


# ---- Figure 6: raw ESS boxplot by method, facet by parameter (Tt = 1600) --

make_fig_ess_raw_boxplot <- function(replicas) {
	plot_data <- replicas %>%
		select(
			method, Tt,
			theta_t1 = ess_theta1_mean,
			theta_t2 = ess_theta2_mean,
			theta_01 = ess_theta_01,
			theta_02 = ess_theta_02,
			W1 = ess_W1,
			W2 = ess_W2
		) %>%
		pivot_longer(
			cols = -c(method, Tt),
			names_to = "parameter",
			values_to = "value"
		) %>%
		mutate(
			parameter = factor(parameter, levels = PARAM_LEVELS),
			Tt_factor = factor(Tt)
		)

	p <- ggplot(plot_data, aes(x = Tt_factor, y = value, fill = method)) +
		geom_boxplot(
			notch = FALSE,
			outlier.shape = 21,
			outlier.fill = NA,
			outlier.stroke = 0.1,
			outlier.size = 0.4,
			outlier.alpha = 1,
			linewidth = 0.1,
			position = position_dodge2(padding = 0.15)
		) +
		facet_wrap(
			~parameter,
			nrow = 3,
			ncol = 2,
			scales = "free_y",
			labeller = as_labeller(PARAM_LABELS, label_parsed)
		) +
		scale_fill_method() +
		scale_y_log10(
			breaks = 10^(1:5),
			labels = label_number(big.mark = ",", drop0trailing = TRUE)
		) +
		labs(x = "T", y = "ESS") +
		theme_paper() +
		guides(fill = guide_legend(nrow = 1))

	return(p)
}


# ---- Figure 7: pointwise HPD coverage by t, facet by function (fixed T) ---

# Rolling-mean smoother for coverage(t) -- the raw per-t series is a mean
# over only R replicas (binomial noise), so a light rolling average makes
# the systematic pattern (dip at breakpoints or not) easier to see. window
# is in units of t (number of points), not seconds/hours.
rolling_mean <- function(x, window) {
	n <- length(x)
	half <- window %/% 2
	return(sapply(seq_len(n), function(i) {
		lo <- max(1, i - half)
		hi <- min(n, i + half)
		return(mean(x[lo:hi]))
	}))
}

# Breakpoints are defined as fractions of Tt (matching the data generator's
# 1/4, 1/2, 3/4 split), not fixed absolute t values, so they line up
# correctly regardless of which Tt is selected. Only the piecewise
# functions (constant, linear, quadratic) have breakpoints; sinusoidal
# has none.
BREAKPOINT_FRACTIONS <- c(0.25, 0.5, 0.75)

make_fig_pointwise_coverage <- function(pointwise, Tt_selected = 1600, window = 21, functions = FUNCTION_LEVELS, facet_nrow = 2, facet_ncol = 2, show_Tt_in_title = FALSE, show_breakpoint_legend = TRUE) {
	plot_data <- pointwise %>%
		filter(Tt == Tt_selected, f %in% functions) %>%
		group_by(method, f) %>%
		arrange(t) %>%
		mutate(coverage_smooth = rolling_mean(coverage, window = window)) %>%
		ungroup()

	breakpoint_data <- expand.grid(
		f = factor(intersect(c("constant", "linear", "quadratic"), functions), levels = FUNCTION_LEVELS),
		xintercept = BREAKPOINT_FRACTIONS * Tt_selected
	)

	facet_labels <- FUNCTION_LABELS
	if (show_Tt_in_title) {
		facet_labels <- setNames(paste0(FUNCTION_LABELS, sprintf(" (T = %d)", Tt_selected)), names(FUNCTION_LABELS))
	}

	p <- ggplot(plot_data, aes(x = t, y = coverage_smooth, color = method, linetype = method)) +
		geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40", linewidth = 0.4) +
		geom_vline(
			data = breakpoint_data,
			aes(xintercept = xintercept),
			inherit.aes = FALSE,
			linetype = "dotted",
			color = "red",
			linewidth = 0.4
		) +
		geom_line(linewidth = 0.35) +
		facet_wrap(~f, nrow = facet_nrow, ncol = facet_ncol, labeller = as_labeller(facet_labels)) +
		scale_color_method() +
		scale_linetype_method() +
		coord_cartesian(ylim = c(0, 1)) +
		labs(x = "t", y = expression(Pointwise~95*"%"~coverage~of~theta[t1])) +
		theme_paper() +
		guides(
			color = guide_legend(nrow = 1),
			linetype = guide_legend(nrow = 1)
		)

	# In-panel legend for the breakpoint vlines -- a plain ggplot legend can't
	# show them (they're not mapped to color/linetype = method), so a small
	# dotted-segment + label is drawn directly inside each panel that has
	# breakpoints, top-right. A white backing rect keeps it legible where it
	# crosses the near-1 coverage curves.
	if (nrow(breakpoint_data) > 0 && show_breakpoint_legend) {
		legend_inset <- data.frame(f = unique(breakpoint_data$f))

		p <- p +
			geom_rect(
				data = legend_inset,
				aes(xmin = 0.66 * Tt_selected, xmax = 0.99 * Tt_selected, ymin = 0.85, ymax = 0.95),
				inherit.aes = FALSE,
				fill = "white",
				alpha = 0.85,
				color = NA
			) +
			geom_segment(
				data = legend_inset,
				aes(x = 0.67 * Tt_selected, xend = 0.72 * Tt_selected, y = 0.90, yend = 0.90),
				inherit.aes = FALSE,
				linetype = "dotted",
				color = "red",
				linewidth = 0.4
			) +
			geom_text(
				data = legend_inset,
				aes(x = 0.735 * Tt_selected, y = 0.90, label = "Breakpoint"),
				inherit.aes = FALSE,
				hjust = 0,
				vjust = 0.5,
				size = 2.8,
				color = "grey30"
			)
	}

	return(p)
}


# ---- Shared helper: pull the legend grob out of a ggplot, for figures that --
# ---- stack several plots and show only one shared legend ------------------

extract_legend <- function(p) {
	g <- ggplotGrob(p)
	idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
	return(g$grobs[[idx]])
}

# Prepends "(T)"/"(t)" to the first axis break instead of a separate,
# dedicated axis-title line -- saves a full text row per plot. When facets
# share one (non-free) x scale, the prefixed break is repeated under every
# panel, not just the leftmost one.
prefix_first_label <- function(labels, prefix) {
	labels[1] <- paste(prefix, labels[1])
	return(labels)
}


# ---- Figure 8 (article): RMSE + coverage (vs T) + pointwise coverage
# ---- (vs t), 3x2 grid, constant & quadratic --------------------------------

make_fig_article_rmse_coverage <- function(
	replicas, aggregated, pointwise,
	functions = FUNCTION_LEVELS,
	Tt_selected = 1600,
	panel_spacing_pt = 2,
	row_margin = margin(1, 6, 1, 4),
	legend_gap_pt = 0,
	font_size = 7
) {
	facet_ncol <- length(functions)
	shrink_text <- theme(text = element_text(size = font_size))

	# Row 1 (RMSE) and row 2 (coverage vs T) share the same x variable (T),
	# so row 1's x-axis text/ticks/title are dropped -- row 2's are enough.
	p_rmse <- make_fig_rmse_boxplot(replicas, functions = functions, facet_nrow = 1, facet_ncol = facet_ncol) +
		shrink_text +
		theme(
			legend.position = "none",
			panel.spacing = unit(panel_spacing_pt, "pt"),
			plot.margin = row_margin,
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank(),
			axis.title.x = element_blank()
		)

	# Row 2 skips its own facet strip titles -- row 1's, directly above,
	# already name each function.
	p_coverage <- make_fig_coverage(aggregated, functions = functions, facet_nrow = 1, facet_ncol = facet_ncol) +
		scale_x_discrete(labels = function(x) prefix_first_label(x, "(T)")) +
		labs(y = expression(Emp.~95*"%"~Cov.~of~theta[1])) +
		shrink_text +
		theme(
			legend.position = "none",
			panel.spacing = unit(panel_spacing_pt, "pt"),
			plot.margin = row_margin,
			strip.text = element_blank(),
			strip.background = element_blank(),
			axis.title.x = element_blank()
		)

	p_pointwise <- make_fig_pointwise_coverage(
		pointwise,
		Tt_selected = Tt_selected,
		functions = functions,
		facet_nrow = 1,
		facet_ncol = facet_ncol,
		show_Tt_in_title = TRUE,
		show_breakpoint_legend = FALSE
	) +
		scale_x_continuous(labels = function(x) prefix_first_label(as.character(x), "(t)")) +
		labs(y = expression("Pt.wise"~95*"%"~Cov.~of~theta[t1])) +
		shrink_text +
		theme(legend.position = "none", panel.spacing = unit(panel_spacing_pt, "pt"), plot.margin = row_margin, axis.title.x = element_blank())

	# Shared bottom legend: built from synthetic dummy data (decoupled from
	# the real panels) so a "Breakpoint" key -- same red dotted linetype as
	# the vlines in row 3 -- can be appended after the method keys, instead
	# of a separate in-panel/caption legend.
	legend_data <- data.frame(
		x = 1,
		y = 1,
		method = factor(c(METHOD_LEVELS, "breakpoint"), levels = c(METHOD_LEVELS, "breakpoint"))
	)
	legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = method, linetype = method, shape = method)) +
		geom_line() +
		geom_point(size = 1.6) +
		scale_color_manual(values = c(METHOD_COLORS, breakpoint = "red"), labels = c(METHOD_LABELS, breakpoint = "Breakpoint"), breaks = c(METHOD_LEVELS, "breakpoint")) +
		scale_linetype_manual(values = c(METHOD_LINETYPES, breakpoint = "dotted"), labels = c(METHOD_LABELS, breakpoint = "Breakpoint"), breaks = c(METHOD_LEVELS, "breakpoint")) +
		scale_shape_manual(values = c(METHOD_SHAPES, breakpoint = NA), labels = c(METHOD_LABELS, breakpoint = "Breakpoint"), breaks = c(METHOD_LEVELS, "breakpoint")) +
		theme_paper() +
		shrink_text +
		guides(
			color = guide_legend(nrow = 1),
			linetype = guide_legend(nrow = 1),
			shape = guide_legend(nrow = 1)
		)
	legend <- extract_legend(legend_plot)

	# Rows 1-2 have less non-panel content than row 3 (no strip text on row 2,
	# no x-axis on row 1) so an equal "null" split renders their actual panels
	# taller than row 3's. Weighting them down brings the three panels to
	# roughly the same height.
	combined <- arrangeGrob(
		p_rmse,
		p_coverage,
		p_pointwise,
		legend,
		ncol = 1,
		heights = unit.c(unit(1, "null"), unit(1, "null"), unit(1.08, "null"), grobHeight(legend) + unit(legend_gap_pt, "pt"))
	)

	return(combined)
}


# ---- Figure 9 (article): raw ESS + ESS/s, 2x2 grid, theta_t1 & W1 ---------

make_fig_article_ess <- function(replicas, params = c("theta_t1", "W1"), Tt_selected = 1600) {
	plot_data_raw <- replicas %>%
		filter(Tt == Tt_selected) %>%
		select(
			method,
			theta_t1 = ess_theta1_mean,
			theta_t2 = ess_theta2_mean,
			theta_01 = ess_theta_01,
			theta_02 = ess_theta_02,
			W1 = ess_W1,
			W2 = ess_W2
		) %>%
		pivot_longer(cols = -method, names_to = "parameter", values_to = "value") %>%
		filter(parameter %in% params) %>%
		mutate(parameter = factor(parameter, levels = params))

	plot_data_sec <- replicas %>%
		select(
			method,
			theta_t1 = ess_sec_theta1_mean,
			theta_t2 = ess_sec_theta2_mean,
			theta_01 = ess_sec_theta_01,
			theta_02 = ess_sec_theta_02,
			W1 = ess_sec_W1,
			W2 = ess_sec_W2
		) %>%
		pivot_longer(cols = -method, names_to = "parameter", values_to = "value") %>%
		filter(parameter %in% params) %>%
		mutate(parameter = factor(parameter, levels = params))

	p_raw <- make_fig_ess_style_boxplot(plot_data_raw, y_label = "bulk-ESS", breaks_pow10 = 1:4, facet_nrow = 1, facet_ncol = 2)
	p_sec <- make_fig_ess_style_boxplot(plot_data_sec, y_label = "bulk-ESS/s", facet_nrow = 1, facet_ncol = 2)

	combined <- arrangeGrob(p_raw, p_sec, ncol = 1)

	return(combined)
}


# ---- Run everything ---------------------------------------------------------

data <- load_data()

save_figure(
	plot = make_fig_rmse_boxplot(data$replicas),
	filename = "rmse.pdf",
	height = 5.2
)

save_figure(
	plot = make_fig_coverage(data$aggregated),
	filename = "coverage.pdf",
	height = 5.2
)

save_figure(
	plot = make_fig_overlay_hpd(data$by_t, Tt_selected = 1600, facet_nrow = 1, facet_ncol = 4),
	filename = "article_hpd.pdf",
	height = 1.8
)

save_figure(
	plot = make_fig_total_time(data$replicas),
	filename = "total_cpu_time.pdf",
	height = 3.2
)

save_figure(
	plot = make_fig_ess_boxplot(data$replicas),
	filename = "ess_sec.pdf",
	height = 7.0
)

save_figure(
	plot = make_fig_ess_raw_boxplot(data$replicas),
	filename = "ess.pdf",
	height = 7.0
)

save_figure(
	plot = make_fig_pointwise_coverage(data$pointwise, Tt_selected = 1600),
	filename = "pointwise_coverage.pdf",
	height = 5.2
)

save_figure(
	plot = make_fig_article_rmse_coverage(data$replicas, data$aggregated, data$pointwise),
	filename = "article_rmse_coverage.pdf",
	height = 4.45
)

save_figure(
	plot = make_fig_article_ess(data$replicas),
	filename = "article_ess.pdf",
	height = 5.5
)

cat("All figures written to", OUTPUT_DIR, "\n")
