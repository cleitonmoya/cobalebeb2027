

# Effective sample size
compute_ess <- function(sample) {
	return(coda::effectiveSize(coda::mcmc(sample))
}

compute_ess_sec <- function()
