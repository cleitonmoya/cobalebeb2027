# Data Generation
#
# Generates Poisson observations y_t ~ Poisson(exp(theta_{t1}))
#
# Four function classes for theta_{t1} = f(t):
#   1) piecewise constant
#   2) piecewise linear
#   3) piecewise quadratic (polynomial, degree 2 per segment)
#   4) sinusoidal
#
# For each class and each Tt in Tt_grid, K replicas are generated with
# randomly or fixed drawn parameters (breakpoints via Dirichlet gaps, 
# levels/amplitudes # via uniform/normal draws). Each replica is saved as an 
# individual .rds file # following the pattern: <function>_<Tt>_<replica>.rds
#
# Author: Cleiton Moya de Almeida

set.seed(42)  # single global seed

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# -----------------------------------------------------------------------------
# 1. Parameters generators
# -----------------------------------------------------------------------------

sample_props <- function(K, alpha_conc, min_gap_frac = NULL) {
    # K segments -> K+1 breakpoints in [0, 1], strictly increasing, props[1] = 0,
    # props[K+1] = 1. Gaps g_i ~ Dirichlet(alpha_conc, ..., alpha_conc).
    repeat {
        g <- rgamma(K, shape = alpha_conc, rate = 1)
        g <- g / sum(g)
        if (is.null(min_gap_frac) || min(g) >= min_gap_frac) break
    }
    props <- c(0, cumsum(g))
    props[K + 1] <- 1  # correct floating point drift
    return(props)
}


generate_constant_params <- function(K = 5,
                                     alpha_conc = 8,
                                     value_range = c(0.5, 3.5),
                                     min_gap_frac = NULL) {
    props  <- sample_props(K, alpha_conc, min_gap_frac)
    values <- runif(K, min = value_range[1], max = value_range[2])
    
    return(list(props = props, values = values))
}


generate_linear_params <- function(K = 5,
                                   alpha_conc = 8,
                                   value_range = c(0.5, 3.5),
                                   min_gap_frac = NULL) {
    props  <- sample_props(K, alpha_conc, min_gap_frac)
    values <- runif(K + 1, min = value_range[1], max = value_range[2])
    
    return(list(props = props, values = values))
}


generate_piecewise_pol_params <- function(K = 5,
                                          alpha_conc = 8,
                                          value_range = c(0.5, 3.5),
                                          bulge_sd = 0.5,
                                          min_gap_frac = NULL) {
    props  <- sample_props(K, alpha_conc, min_gap_frac)
    values <- runif(K + 1, min = value_range[1], max = value_range[2])
    bulges <- rnorm(K, mean = 0, sd = bulge_sd)
    
    return(list(props = props, values = values, bulges = bulges))
}


generate_sinusoidal_params <- function(mean_range   = c(1.5, 3.0),
                                       amp_range    = c(0.3, 1.0),
                                       n_cycles_set = 2:2) {
    mean_level <- runif(1, min = mean_range[1], max = mean_range[2])
    # Keep amplitude below mean_level to avoid theta_t1 too close to/below 0
    amp_upper  <- min(amp_range[2], mean_level - 0.3)
    amplitude  <- runif(1, min = amp_range[1], max = max(amp_range[1] + 1e-6, amp_upper))
    n_cycles   <- sample(n_cycles_set, 1)
    phase      <- runif(1, min = 0, max = 2 * pi)
    
    return(list(mean_level = mean_level, amplitude = amplitude,
                n_cycles = n_cycles, phase = phase))
}

# -----------------------------------------------------------------------------
# 2. Deterministic function generators (theta_{t1} = f(t))
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
# 3. Registry mapping function name -> (parameter generator, series generator)
# -----------------------------------------------------------------------------

function_registry <- list(
    constant = list(
        gen_params = function() generate_constant_params(),
        gen_series = function(Tt, params) piecewise_constant(Tt, params$props, params$values)
    ),
    linear = list(
        gen_params = function() generate_linear_params(),
        gen_series = function(Tt, params) piecewise_linear(Tt, params$props, params$values)
    ),
    quadratic = list(
        gen_params = function() generate_piecewise_pol_params(),
        gen_series = function(Tt, params) piecewise_pol(Tt, params$props, params$values, params$bulges)
    ),
    sinusoidal = list(
        gen_params = function() generate_sinusoidal_params(),
        gen_series = function(Tt, params) sinusoidal(Tt, params$mean_level, params$amplitude,
                                                     params$n_cycles, params$phase)
    )
)

# -----------------------------------------------------------------------------
# 4. Driver: generate K replicas per function class and save to .rds
# -----------------------------------------------------------------------------

simulate_and_save <- function(Tt_grid, K, output_dir, fixed_parameters = FALSE) {
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    func_names <- names(function_registry)
    
    for (func_name in func_names) {
        entry <- function_registry[[func_name]]
        
        # If fixed_parameters = TRUE, parameters are drawn ONCE per function
        # (fixed across all K replicas and all Tt in the grid) -> same theta1
        # (true level) for every replica; only y varies, via the Poisson
        # draw seed. If FALSE, parameters are drawn once per (function,
        # replica) -> same waveform shape across all Tt in the grid, only
        # the discretization changes.
        if (fixed_parameters) params <- entry$gen_params()
        
        for (r in 1:K) {
            
            if (!fixed_parameters) params <- entry$gen_params()
            
            for (Tt in Tt_grid) {
                theta1 <- entry$gen_series(Tt, params)
                y      <- rpois(Tt, lambda = exp(theta1))
                
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
# 5. Execution: 4 functions x K = 50 replicas x 4 values of Tt
# -----------------------------------------------------------------------------

Tt_grid <- c(200, 400, 800, 1600)
K <- 100
output_dir = "../data/simulated/"
fixed_parameters <- TRUE  # same theta1 (across all K replicas per function

simulate_and_save(Tt_grid = Tt_grid, K = K, output_dir = output_dir,
                   fixed_parameters = fixed_parameters)
