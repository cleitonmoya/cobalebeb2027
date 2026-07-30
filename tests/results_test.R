
# Provide more informative traceback
options(error = function() traceback(2)) 

# Change de directory to the same of the current file
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

path_results <- "../simulation/results"



