# simulation/simulation_grid_config.R
#
# Edit this file directly on the cluster to control what simulation_run.R
# does. No need to re-send simulation_run.R for this -- it's sourced by
# simulation_run.R every time (direct run, batchtools job re-source, or
# check_simulation_progress.R/simulation_aggregate.R via SIMULATION_DEFS_ONLY).
#
# Just two variables:

run_mode <- "cluster"   # "local" or "cluster"

# Grid subset: NULL runs everything; otherwise a function(job_grid) ->
# job_grid (row filter), one row per chunk. job_grid has columns
# chunk_id, category ("leve"/"medio"/"pesado"), Tt.
# grid_subset <- function(g) rbind(head(subset(g, category == "medio"), 1), head(subset(g, category == "pesado" & Tt == 200), 1))   # [TESTE RELAMPAGO: 1 chunk medio + 1 chunk pesado/Tt=200]
grid_subset <- NULL                                            # everything
# grid_subset <- function(g) subset(g, category == "pesado")     # only stan
# grid_subset <- function(g) subset(g, Tt == 1600)                # only Tt=1600, all categories
# grid_subset <- function(g) g[1, ]                                # exactly one chunk
