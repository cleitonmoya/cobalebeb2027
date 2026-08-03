# Poisson Local Trend Dynamic Model
# Results aggregation script
# Author: Cleiton Moya de Almeida

rm(list = ls())     # clear the environment

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

path_data <- "../data/simulated"
path_results <- "results"
path_results_partial <- sprintf("%s/%s", path_results, "partial")
path_results_star <- sprintf("%s/%s", path_results, "star")

verbose <- TRUE

# Print auxiliary function
printf <- function(...) cat(paste(sprintf(...), "\n"))


#####
# Read and combine all partial results

partial_files <- list.files(path_results_partial, pattern = "\\.rds$", full.names = TRUE)

if (verbose) printf("Partial files found: %d", length(partial_files))

if (length(partial_files) == 0) {
	stop("No files found in ", path_results_partial)
}

df_results <- do.call(rbind, lapply(partial_files, readRDS))

if (verbose) printf("Total combined rows: %d", nrow(df_results))


#####
# Save the aggregated results

file_rds <- sprintf("%s/simulation_results.rds", path_results)
saveRDS(df_results, file = file_rds)

file_csv <- sprintf("%s/simulation_results.csv", path_results)
write.csv(df_results, file = file_csv, row.names = FALSE)

if (verbose) printf("Aggregated results saved to %s and %s", file_rds, file_csv)
