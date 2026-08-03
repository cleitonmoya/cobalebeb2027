# Data Generation (fixed-parameters version)
#
# Generates Poisson observations y_t ~ Poisson(exp(theta_{t1}))
#
# Four function classes for theta_{t1} = f(t), all with FIXED parameters
# (no random draws of breakpoints, levels, curvature, or sinusoidal
# parameters). Only y varies across replicas, via the Poisson draw.
#
#   1) piecewise constant
#   2) piecewise linear
#   3) piecewise quadratic (polynomial, degree 2 per segment)
#   4) sinusoidal
#
# For each class and each Tt in Tt_grid, M replicas are generated with the
# SAME theta1 (fixed parameters); only y changes across replicas. Each
# replica is saved as an individual .rds file following the pattern:
# <function>_<Tt>_<replica>.rds
#
# Author: Cleiton Moya de Almeida

# Set directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# -----------------------------------------------------------------------------
# 1. Deterministic function generators (theta_{t1} = f(t))
# -----------------------------------------------------------------------------

piecewise_constant <- function(Tt, props, values) {
    breakpoints <- round(props * Tt)
    y <- numeric(Tt)
    K <- length(values)
    
    for (i in seq_len(K)) {
        t0 <- breakpoints[i]
        t1 <- breakpoints[i + 1]
        tt <- (t0 + 1):t1
        y[tt] <- values[i]
    }
    
    return(y)
}


piecewise_linear <- function(Tt, props, values) {
    breakpoints <- round(props * Tt)
    y <- numeric(Tt)
    y[1] <- values[1]
    K <- length(values) - 1
    
    for (i in seq_len(K)) {
        t0 <- breakpoints[i]
        t1 <- breakpoints[i + 1]
        v0 <- values[i]
        v1 <- values[i + 1]
        
        tt <- (t0 + 1):t1
        s  <- (tt - t0) / (t1 - t0)
        
        y[tt] <- v0 + (v1 - v0) * s
    }
    
    return(y)
}


piecewise_pol <- function(Tt, props, values, bulges) {
    breakpoints <- round(props * Tt)
    y <- numeric(Tt)
    y[1] <- values[1]
    
    for (i in seq_along(bulges)) {
        t0 <- breakpoints[i]
        t1 <- breakpoints[i + 1]
        v0 <- values[i]
        v1 <- values[i + 1]
        b  <- bulges[i]
        
        tt <- (t0 + 1):t1
        s  <- (tt - t0) / (t1 - t0)
        
        y[tt] <- v0 + (v1 - v0) * s + b * s * (1 - s)
    }
    
    return(y)
}


sinusoidal <- function(Tt, mean_level, amplitude, n_cycles, phase) {
    tt <- 1:Tt
    y  <- mean_level + amplitude * sin(2 * pi * n_cycles * tt / Tt + phase)
    return(y)
}

# -----------------------------------------------------------------------------
# 2. Fixed parameters
# -----------------------------------------------------------------------------

fixed_props <- c(0, 0.25, 0.5, 0.75, 1)  # K = 4 segments, shared by
                                          # constant, linear and quadratic

function_registry <- list(
    constant = list(
        params = list(props = fixed_props, values = c(1, 2, 0.5, 1.5)),
        gen_series = function(Tt, params) piecewise_constant(Tt, params$props, params$values)
    ),
    linear = list(
        params = list(props = fixed_props, values = c(1, 2, 0.5, 1.5, 1.25)),
        gen_series = function(Tt, params) piecewise_linear(Tt, params$props, params$values)
    ),
    quadratic = list(
        params = list(props = fixed_props, values = c(1, 2, 0.5, 1.5, 1.25),
                      bulges = c(1, -1, 1, -2)),
        gen_series = function(Tt, params) piecewise_pol(Tt, params$props, params$values, params$bulges)
    ),
    sinusoidal = list(
        params = list(mean_level = 1.25, amplitude = 0.75, n_cycles = 1, phase = 0),
        gen_series = function(Tt, params) sinusoidal(Tt, params$mean_level, params$amplitude,
                                                     params$n_cycles, params$phase)
    )
)

# -----------------------------------------------------------------------------
# 3. Driver: generate M replicas per function class and save to .rds
# -----------------------------------------------------------------------------

simulate_and_save_fixed <- function(Tt_grid, M, output_dir, seed_base = 1) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    for (func_name in names(function_registry)) {
        entry  <- function_registry[[func_name]]
        params <- entry$params
        
        for (Tt in Tt_grid) {
            
            # Reset the RNG for each (func_name, Tt) block, so increasing M
            # only appends new replicas within a block, without disturbing
            # the other blocks (which do not depend on M anymore).
            set.seed(seed_base + Tt)
            
            theta1 <- entry$gen_series(Tt, params)
            
            for (r in 1:M) {
                y <- rpois(Tt, lambda = exp(theta1))
                
                replica_data <- list(
                    function_name = func_name,
                    Tt            = Tt,
                    replica       = r,
                    params        = params,
                    theta1        = theta1,
                    y             = y
                )
                
                file_name <- sprintf("%s_%d_%d.rds", func_name, Tt, r)
                saveRDS(replica_data, file = file.path(output_dir, file_name))
            }
        }
    }
    
    invisible(NULL)
}

# -----------------------------------------------------------------------------
# 4. Execution: 4 functions x M = 100 replicas x 4 values of Tt
# -----------------------------------------------------------------------------

Tt_grid <- c(200, 400, 800, 1600)
M <- 200
output_dir <- "../data/simulated_fixed/"
seed_base <- 42

simulate_and_save_fixed(Tt_grid = Tt_grid, M = M, output_dir = output_dir, seed_base = seed_base)
