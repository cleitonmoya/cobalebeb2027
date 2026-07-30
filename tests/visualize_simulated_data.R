# =============================================================================
# plot_simulated_series.R
#
# Complementary script to visualize simulated series saved as .rds files
# (produced by simulate_state_functions.R).
#
# For each function class, displays one figure with a 2x2 grid (base R):
#   rows    -> replicas (2 selected replicas)
#   columns -> selected values of Tt (e.g. 200 and 2000)
#
# Each panel shows the observed counts y_t (points) together with the
# underlying mean exp(theta_t1) (line), matching y_t ~ Poisson(exp(theta_t1)).
# Assumes series are already saved on disk; does not regenerate anything.
# Plots are displayed only (RStudio graphics device), not saved to disk.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

input_dir    <- "../data/simulated/"
func_names   <- c("constant", "linear", "quadratic", "sinusoidal")
Tt_selected  <- c(200, 1600)
replicas_sel <- c(1, 2)

# -----------------------------------------------------------------------------
# 2. Helper: load one replica file
# -----------------------------------------------------------------------------

load_replica <- function(func_name, Tt, replica, dir = input_dir) {
	file_name <- sprintf("%s_%d_%d.rds", func_name, Tt, replica)
	file_path <- file.path(dir, file_name)
	
	if (!file.exists(file_path)) {
		stop(sprintf("File not found: %s", file_path))
	}
	
	readRDS(file_path)
}

# -----------------------------------------------------------------------------
# 3. Helper: draw a single panel for one replica
# -----------------------------------------------------------------------------

draw_panel <- function(replica_data, type = c("y", "theta"), row_label = NULL, col_label = NULL) {
	
	type <- match.arg(type)
	
	Tt     <- replica_data$Tt
	theta1 <- replica_data$theta1
	
	if (type == "theta") {
		plot(1:Tt, theta1,
			 type = "l", col = "red", lwd = 2,
			 xlab = if (!is.null(col_label)) col_label else "",
			 ylab = if (!is.null(row_label)) row_label else "",
			 font.lab = 2)
	} else {
		y      <- replica_data$y
		mean_t <- exp(theta1)
		
		plot(1:Tt, y,
			 type = "p", pch = 16, cex = 0.5,
			 xlab = if (!is.null(col_label)) col_label else "",
			 ylab = if (!is.null(row_label)) row_label else "",
			 font.lab = 2)
		
		lines(1:Tt, mean_t, col = "red", lwd = 2)
	}
}


# -----------------------------------------------------------------------------
# 4. Build 2x3 grid figure for one function class
# -----------------------------------------------------------------------------

plot_function_grid <- function(func_name, Tt_vals = Tt_selected, reps = replicas_sel) {
	
	old_par <- par(no.readonly = TRUE)
	on.exit(par(old_par))
	
	par(mfrow = c(length(reps), length(Tt_vals) + 1),
		oma = c(0, 0, 3, 0),
		mar = c(4, 4, 1, 1))
	
	for (i in seq_along(reps)) {
		
		row_label <- sprintf("Replica %d", reps[i])
		
		# First column: theta_t1, with Tt = Tt_vals[i] (row-specific)
		replica_data_theta <- load_replica(func_name, Tt_vals[i], reps[i])
		col_label <- if (i == length(reps)) "theta_t1" else ""
		draw_panel(replica_data_theta, type = "theta", row_label = row_label, col_label = col_label)
		
		# Remaining columns: y_t vs exp(theta_t1) for each selected Tt
		for (j in seq_along(Tt_vals)) {
			
			replica_data <- load_replica(func_name, Tt_vals[j], reps[i])
			col_label <- if (i == length(reps)) sprintf("T = %d", Tt_vals[j]) else ""
			draw_panel(replica_data, type = "y", row_label = "", col_label = col_label)
		}
	}
	
	mtext(sprintf("Function: %s", func_name), outer = TRUE, cex = 1.3, font = 2)
}

# -----------------------------------------------------------------------------
# 5. Display one figure per function class
# -----------------------------------------------------------------------------

for (func_name in func_names) {
	plot_function_grid(func_name, Tt_vals = Tt_selected, reps = replicas_sel)
}