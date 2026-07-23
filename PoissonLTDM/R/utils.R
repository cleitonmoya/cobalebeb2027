# Common functions for the Poisson LTDM samplers

# Used in: sampler_pg_apf, sampler_sir_laplace, sampler_sir_collapsed 
logsumexp <- function(x) {
	cc <- max(x)
	return(cc + log(sum(exp(x - cc))))
}


# Used in: sampler_pg_apf
# Log-likelihood
log_p_yt <- function(yt, theta_t1) {
	res <- yt * theta_t1 - exp(theta_t1)
	return(res)
}

# Used in: mh_cw, mh_montoril
gibbs_sample_theta01 <- function(mu_01, sigma2_01, theta_11, theta_02, W1) {
	sigma2_01_bar <- (1/sigma2_01 + 1/W1)^(-1)
	mu_01_bar <- sigma2_01_bar*(mu_01/sigma2_01 + (theta_11 - theta_02)/W1)
	theta_01 <- rnorm(1, mean=mu_01_bar, sd=sqrt(sigma2_01_bar))
	return(theta_01)
}


# Used in: mh_cw, mh_montoril
gibbs_sample_theta02 <- function(mu_02, sigma2_02, theta_01, theta_11,
								 theta_12, W1, W2) {
	
	sigma2_02_bar <- (1/sigma2_02 + 1/W1 + 1/W2)^(-1)
	mu_02_bar <- sigma2_02_bar*((theta_11 - theta_01)/W1 +
									theta_12/W2 + mu_02/sigma2_02)
	theta_02 <- rnorm(1, mean=mu_02_bar, sd=sqrt(sigma2_02_bar))
	return(theta_02)
}


# Used in: mh_cw, mh_montoril
gibbs_sample_phi1 <- function(nu_01, eta_01, theta_01, theta1, theta_02, theta2, Tt) {
	nu_01_bar <- nu_01 + Tt/2
	dif1 <- theta1 - c(theta_01, theta1[-Tt])
	dif2 <- dif1 - c(theta_02, theta2[-Tt])
	eta_01_bar <- eta_01 + 0.5 * sum(dif2^2)
	phi1 <- rgamma(1, nu_01_bar, eta_01_bar)
	return(phi1)
}


# Used in: mh_cw, mh_montoril
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


# Used in: mh_cw, mh_montoril
# Sample theta1 using component-wise Metropolis (Random Walking)
cwmh_sample_theta1 <- function(y, theta_01, theta_02, 
							   theta1, theta2, W1, varsigma2, Tt) {
	
	n_ac <- 0 # number of accepted samples
	
	for (t in 1:Tt) {
		
		if (t < Tt) {
			
			if (t==1) {
				# theta_t11
				res <- sample_theta_t1_mh(theta1[t], theta_01, theta1[t+1],
										  theta2[t], theta_02,
										  y[t], W1, varsigma2, final_t=FALSE)
				theta1[t] <- res$theta_t1
			} else {
				
				res <- sample_theta_t1_mh(theta1[t], theta1[t-1], theta1[t+1],
										  theta2[t], theta2[t-1],
										  y[t], W1, varsigma2, final_t=FALSE)
				theta1[t] <- res$theta_t1
			}
			
		} else {
			res <- sample_theta_t1_mh(theta1[t], theta1[t-1], NULL,
									  theta2[t], theta2[t-1],
									  y[t], W1, varsigma2, final_t=TRUE)
			theta1[t] <- res$theta_t1
		}
		
		ac <- res$ac # flag: sample accepted(1) or not (0)
		n_ac <- n_ac + ac
	}
	return(list(theta1 = theta1, n_ac = n_ac))
}	


# used in: sampler_mh_montoril
# Sample theta2 using the Chan Method (conjugated Normal)
make_chan_theta2_sampler <- function(Tt) {
	
	# FIXED SPARSE STRUCTURES FOR CHAN METHOD ####
	# Base for the prior Precision Matrix K
	sub_diag_base <- rep(-1, Tt-1)
	main_diag_base <- c(rep(2, Tt-1), 1)
	K0 <- bandSparse(n=Tt, k=c(0, -1),
					 diagonals=list(main_diag_base, sub_diag_base),
					 symmetric = TRUE)
	# diagonal mask
	# @x: slot of the Sparce matrix (S4 object) that contains the non-zero values
	diag_pattern <- bandSparse(n=Tt, k=c(0, -1),
							   diagonals=list(rep(TRUE, Tt), rep(FALSE, Tt-1)),
							   symmetric=TRUE)
	idx_diag <- which(diag_pattern@x) # index of subpattern@x which is non-zero
	
	# subdiagonal mask
	sub_pattern <- bandSparse(n=Tt, k=c(0, -1),
							  diagonals=list(rep(FALSE, Tt), rep(TRUE, Tt-1)),
							  symmetric=TRUE)
	idx_sub <- which(sub_pattern@x)
	
	# Initial symbolic Cholesky factor
	Ch02_factor <- Cholesky(K0, perm = FALSE, LDL = TRUE)
	
	# Work precision matrix (static)
	P2_matrix <- K0
	
	
	chan_sample_theta2 <- function(theta1, phi1, phi2, theta_02, Tt) {
		
		z <- diff(theta1)   # z_t = theta1[t+1] - theta1[t], t=1,...,T-1
		
		diag_obs <- c(rep(phi1, Tt-1), 0)
		P2_matrix@x[idx_diag] <<- (main_diag_base*phi2) + diag_obs
		P2_matrix@x[idx_sub]  <<- -phi2
		
		Ch2_factor <- update(Ch02_factor, P2_matrix)
		
		b <- numeric(Tt)
		b[1:(Tt-1)] <- z * phi1
		b[1] <- b[1] + theta_02 * phi2
		
		theta2_hat <- as.numeric(Matrix::solve(Ch2_factor, b, system="A"))
		d <- Matrix::diag(Ch2_factor)
		u <- rnorm(Tt)
		w <- u/sqrt(d)
		x <- as.vector(Matrix::solve(Ch2_factor, w, system="Lt"))
		
		return(theta2_hat + x)
	}
	
	return(chan_sample_theta2)
}
