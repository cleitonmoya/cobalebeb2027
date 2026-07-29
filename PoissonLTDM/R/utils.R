# Common functions for the Poisson LTDM samplers

# Used in: pg_apf, sir_laplace, sir_collapsed 
logsumexp <- function(x) {
	cc <- max(x)
	return(cc + log(sum(exp(x - cc))))
}


# Used in: sampler_pg_apf, sir_laplace
# Log-likelihood
log_p_yt <- function(yt, theta_t1) {
	res <- yt * theta_t1 - exp(theta_t1)
	return(res)
}


# Used in: mh_montoril
gibbs_sample_theta01 <- function(mu_01, sigma2_01, theta_11, theta_02, W1) {
	sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
	mu_01_bar <- sigma2_01_bar*(mu_01/sigma2_01 + (theta_11 - theta_02)/W1)
	theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))
	return(theta_01)
}


# Used in: mh_montoril, pg_apf, sir_laplace, sir_collapsed
gibbs_sample_phi1 <- function(nu_01, eta_01, theta_01, theta1, theta_02, theta2, Tt) {
	nu_01_bar <- nu_01 + Tt/2
	dif1 <- theta1 - c(theta_01, theta1[-Tt])
	dif2 <- dif1 - c(theta_02, theta2[-Tt])
	eta_01_bar <- eta_01 + 0.5 * sum(dif2^2)
	phi1 <- rgamma(1, nu_01_bar, eta_01_bar)
	return(phi1)
}


# Used in: mh_montoril, pg_apf, sir_laplace, sir_collapsed
gibbs_sample_phi2 <- function(nu_02, eta_02, theta_02, theta2, Tt) {
	nu_02_bar <- nu_02 + Tt/2
	diffs2 <- theta2 - c(theta_02, theta2[-Tt])
	eta_02_bar <- eta_02 + 0.5 * sum(diffs2^2)
	phi2 <- rgamma(1, nu_02_bar, eta_02_bar)
	return(phi2)
}


# Used in: cwmh_sample_theta1 
# Full conditional log-posterior for theta_t1, t=1, ..., T-1
# theta_t1: theta_{t1}
# theta_tm11: theta_{t-1,1}
# theta_tp11: theta_{t+1,1}
# theta_t2: theta_{t2}
# theta_tm12: theta_{t-1,2}
logpost_theta_t1 <- function(theta_t1, theta_tm11, theta_tp11,
							 theta_t2, theta_tm12, yt, W1) {
	sigma2_star <- W1/2
	mu_star <- ((theta_tm11+theta_tm12)+(theta_tp11-theta_t2))/2
	
	p1 <- yt*theta_t1 - exp(theta_t1) # log-likelihood
	p2 <- -(theta_t1 - mu_star)^2/(2*sigma2_star)
	logp <- p1 + p2
	return(logp)
}


# Used in: cwmh_sample_theta1
# Full conditional log-posterior for theta_T1 (t=T)
logpost_theta_T1 <- function(theta_t1, theta_tm11, theta_tm12, yt, W1) {
	p1 <- yt*theta_t1 - exp(theta_t1) # log-likelihood
	p2 <- -(theta_t1 - theta_tm11 - theta_tm12)^2/(2*W1)
	logp <- p1+p2
	return(logp)
}


# Used in: cwmh_sample_theta1 
# Sample theta_t1 ~ logpost_theta_t1 (Metropolis step)
# final_t: boolean (0: t<T; 1: t=T)
sample_theta_t1_mh <- function(theta_t1_current, theta_tm11, theta_tp11,
							theta_t2, theta_tm12,
							yt, W1, varsigma2, final_t) {
	
	# proposed theta
	theta_t1_prop <- rnorm(1, mean=theta_t1_current, sd=sqrt(varsigma2))
	
	# acceptance/rejection step
	ac <- 0 # accepted flag
	logu <- log(runif(1))
	
	if (final_t) {
		logp1 <- logpost_theta_T1(theta_t1_prop, theta_tm11, theta_tm12, yt, W1)
		logp2 <- logpost_theta_T1(theta_t1_current, theta_tm11, theta_tm12, yt, W1)
		
	} else {
		logp1 <- logpost_theta_t1(theta_t1_prop, theta_tm11, theta_tp11,
								  theta_t2, theta_tm12, yt, W1)
		logp2 <- logpost_theta_t1(theta_t1_current, theta_tm11, theta_tp11,
								  theta_t2, theta_tm12, yt, W1)
	}
	logr <- logp1 - logp2
	
	# acceptance criteria
	if (logu < logr){
		theta_t1 <- theta_t1_prop
		ac <- 1
	} else {
		theta_t1 <- theta_t1_current
	}
	
	return(list(theta_t1=theta_t1, ac=ac))
}


# Used in: mh_montoril
# Sample theta1 using component-wise Metropolis (Random Walking)
cwmh_sample_theta1 <- function(y, theta_01, theta_02, 
							   theta1, theta2, W1, varsigma2, Tt) {
	
	ac <- numeric(Tt)
	for (t in 1:Tt) {
		
		if (t < Tt) {
			
			if (t==1) {
				# theta_t11
				res <- sample_theta_t1_mh(theta1[t], theta_01, theta1[t+1],
										  theta2[t], theta_02,
										  y[t], W1, varsigma2[t], final_t=FALSE)
				theta1[t] <- res$theta_t1
			} else {
				
				res <- sample_theta_t1_mh(theta1[t], theta1[t-1], theta1[t+1],
										  theta2[t], theta2[t-1],
										  y[t], W1, varsigma2[t], final_t=FALSE)
				theta1[t] <- res$theta_t1
			}
			
		} else {
			res <- sample_theta_t1_mh(theta1[t], theta1[t-1], NULL,
									  theta2[t], theta2[t-1],
									  y[t], W1, varsigma2[t], final_t=TRUE)
			theta1[t] <- res$theta_t1
		}
		
		ac[t] <- res$ac # flag: sample accepted(1) or not (0)
	}
	return(list(theta1 = theta1, ac = ac))
}	


# Used in: chan_build_static_objects, chan_build_static_objects_ext
# Shared low-level builder for the Chan sparse chain (RW1) precision structure.
# first_node_order: degree of the first node of the n-node chain
#   2: ANCHORED chain of n free nodes, with an externally fixed predecessor
#      (the original, T-dimensional blocks - theta_01/theta_02 held fixed)
#   1: FREE/unanchored chain of n nodes (pure random walk skeleton), used by
#      the extended, (T+1)-dimensional joint blocks (theta_01/theta_02 folded
#      in as node "0" of the chain itself)
# Note: c(2, rep(2,n-2), 1) == c(rep(2,n-1), 1), so first_node_order=2
# reproduces exactly the original (pre-extension) diagonal pattern.
chan_build_chain <- function(n, first_node_order) {
	
	# FIXED SPARSE STRUCTURES FOR CHAN METHOD ####
	# Base for the prior Precision Matrix K
	sub_diag_base <- rep(-1, n-1)
	main_diag_base <- c(first_node_order, rep(2, n-2), 1)
	
	# The unanchored/free chain (first_node_order=1) is singular on its own
	# (a pure RW skeleton has a free translation direction), so the initial
	# K0/Cholesky() call below needs a positive-definite placeholder
	# (main_diag_base + 1) just to extract the sparsity pattern; the real
	# (always PD, thanks to the prior + data terms added on top) values are
	# always filled in via update() before actual use, using main_diag_base
	# (without the placeholder) - see make_chan_theta1_smoother_ext/
	# make_chan_theta2_smoother_ext. The anchored chain (first_node_order=2)
	# is already PD as-is, so it keeps the exact main_diag_base: unlike the
	# extended blocks, its Ch0_factor is read directly (never update()'d) by
	# chan_log_det_K0, so K0 must equal the true precision matrix there.
	symbolic_diag <- if (first_node_order == 1) main_diag_base + 1 else main_diag_base
	K0 <- Matrix::bandSparse(n=n, k=c(0, -1),
							 diagonals=list(symbolic_diag, sub_diag_base),
							 symmetric = TRUE)
	# diagonal mask
	# @x: slot of the Sparce matrix (S4 object) that contains the non-zero values
	diag_pattern <- Matrix::bandSparse(n=n, k=c(0, -1),
									   diagonals=list(rep(TRUE, n), rep(FALSE, n-1)),
									   symmetric=TRUE)
	idx_diag <- which(diag_pattern@x) # index of subpattern@x which is non-zero
	
	# subdiagonal mask
	sub_pattern <- Matrix::bandSparse(n=n, k=c(0, -1),
									  diagonals=list(rep(FALSE, n), rep(TRUE, n-1)),
									  symmetric=TRUE)
	idx_sub <- which(sub_pattern@x)
	
	# Initial symbolic Cholesky factor
	Ch0_factor <- Matrix::Cholesky(K0, perm = FALSE, LDL = TRUE)
	
	return(list(
		K0 = K0, 
		Ch0_factor = Ch0_factor,
		main_diag_base = main_diag_base,
		sub_diag_base = sub_diag_base,
		idx_diag = idx_diag,
		idx_sub = idx_sub))
}


# Used in: sampler_sir_collapsed's make_chan_theta2_smoother (W2 marginal
# likelihood, theta_02 held fixed) - ANCHORED, T-dimensional block: T free
# nodes with an externally fixed predecessor (theta_01 or theta_02 held fixed)
chan_build_static_objects <- function(Tt) {
	chan_build_chain(Tt, first_node_order=2)
}


# Used in: make_chan_theta2_smoother_ext, make_chan_theta1_smoother_ext
# Extended, (T+1)-dimensional block: theta_01 (resp. theta_02) is folded in as
# node "0" of the chain, so the WHOLE chain (T+1 nodes) is FREE/unanchored
# (pure random walk skeleton) - node 0's own prior and cross-block coupling
# terms are added on top of this skeleton (see make_chan_theta1_smoother_ext/
# make_chan_theta2_smoother_ext)
chan_build_static_objects_ext <- function(Ttp1) {
	chan_build_chain(Ttp1, first_node_order=1)
}


# Used in: sir_collapsed
# log|K0| (constant, precomputed once - used in the exact W2 marginal likelihood)
chan_log_det_K0 <- function(Tt) {
	res <- chan_build_static_objects(Tt)
	log_det_K0 <- 2 * as.numeric(Matrix::determinant(res$Ch0_factor, logarithm = TRUE)$modulus)
	return(log_det_K0)
}


# Used in: sampler_amh_montoril, sampler_pg_apf, sampler_sir_laplace, sampler_sir_collapsed
# Sample (theta_02, theta2) JOINTLY via the extended, (T+1)-dimensional Chan
# block: theta_02 is folded in as node "0" of the chain. theta_02's own
# diagonal entry combines its Normal prior (1/sigma2_02) with the phi1
# contribution coming from theta1[1] = theta_01 + theta_02 + omega_1 (the
# channel into the OTHER extended block, theta1). Exact (no Poisson
# likelihood involved) - build the block here, then draw from it with
# chan_sample_from_build(build, Ttp1).
make_chan_theta2_smoother_ext <- function(Ttp1) {
	
	res <- chan_build_static_objects_ext(Ttp1)
	P2_matrix      <- res$K0
	Ch02_factor    <- res$Ch0_factor
	main_diag_rw   <- res$main_diag_base
	idx_diag_e     <- res$idx_diag
	idx_sub_e      <- res$idx_sub
	Tt <- Ttp1 - 1
	
	chan_smoothing_theta2 <- function(theta1, phi1, phi2, mu_02, sigma2_02, theta_01) {
		z <- diff(theta1)   # z_t = theta1[t+1] - theta1[t], t=1,...,T-1
		
		extra_diag <- c(1/sigma2_02 + phi1, rep(phi1, Tt-1), 0)
		P2_matrix@x[idx_diag_e] <<- (main_diag_rw*phi2) + extra_diag
		P2_matrix@x[idx_sub_e]  <<- -phi2
		
		Ch2_factor <- Matrix::update(Ch02_factor, P2_matrix)
		
		b <- numeric(Ttp1)
		b[1] <- mu_02/sigma2_02 + phi1*(theta1[1] - theta_01)
		b[2:Tt] <- z*phi1
		b[Ttp1] <- 0
		
		theta2_hat <- as.numeric(Matrix::solve(Ch2_factor, b, system="A"))
		return(list(theta_hat=theta2_hat, ch=Ch2_factor, z=z))
	}
	
	return(chan_smoothing_theta2)
}


# Used in: sampler_sir_laplace, sampler_sir_collapsed
# Laplace/IRLS approximation for (theta_01, theta1) JOINTLY via the extended,
# (T+1)-dimensional Chan block: theta_01 is folded in as node "0" of the
# chain. theta_01 needs no linearization itself (no likelihood) - only its
# exact Gaussian prior and the exact process link into theta1[1], both
# already captured by main_diag_rw*phi1 + extra_diag (no extra "+phi1" term
# at node 0, unlike theta_02, since theta_01 has no external channel into
# another block).
make_chan_theta1_smoother_ext <- function(Ttp1) {
	
	res <- chan_build_static_objects_ext(Ttp1)
	P1_matrix      <- res$K0
	Ch01_factor    <- res$Ch0_factor
	main_diag_rw   <- res$main_diag_base
	idx_diag_e     <- res$idx_diag
	idx_sub_e      <- res$idx_sub
	Tt <- Ttp1 - 1
	
	chan_smoothing_theta1 <- function(z_t, phi_V, phi1, mu_01, sigma2_01, theta_02, theta2) {
		extra_diag <- c(1/sigma2_01, phi_V)
		P1_matrix@x[idx_diag_e] <<- (main_diag_rw*phi1) + extra_diag
		P1_matrix@x[idx_sub_e]  <<- -phi1
		
		Ch1_factor <- Matrix::update(Ch01_factor, P1_matrix)
		
		RHS_ext <- c(theta_02, theta2[-Tt])
		Hb_ext <- numeric(Ttp1)
		Hb_ext[1] <- -RHS_ext[1]
		Hb_ext[2:Tt] <- RHS_ext[1:(Tt-1)] - RHS_ext[2:Tt]
		Hb_ext[Ttp1] <- RHS_ext[Tt]
		b <- c(mu_01/sigma2_01, phi_V*z_t) + phi1*Hb_ext
		
		theta1_hat <- as.numeric(Matrix::solve(Ch1_factor, b, system="A"))
		return(list(theta1_hat=theta1_hat, ch=Ch1_factor))
	}
	
	return(chan_smoothing_theta1)
}


# Used in: sampler_amh_montoril, sampler_pg_apf, sampler_sir_laplace, sampler_sir_collapsed
# It can be used with theta1 or theta2, ANCHORED (Tt) or extended (Ttp1) builds
chan_sample_from_build <- function(build, Tt) {
	
	ch <- build$ch
	theta_hat <- build[[1]] # theta1_hat or theta2_hat, always the first element
	
	d <- Matrix::diag(ch)
	u <- rnorm(Tt)
	w <- u / sqrt(d)
	x <- as.vector(Matrix::solve(ch, w, system="Lt"))
	
	return(theta_hat + x)
}
