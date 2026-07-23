# Poisson - 2nd Order Polynomial Dynamic Model
# Gibbs with Collapsed Samplers for W1 AND W2
# Strategy:
#   - Pre-run: standard Gibbs (R_prerun iterations) to collect phi1 and phi2 samples
#   - CE calibration: fit Gamma(c_hat, d_hat) proposals for both phi1 and phi2
#   - Main loop:
#       * theta_02, theta_01 : conjugated Normal
#       * W1  : collapsed MH with CE proposal, marginalizing theta1 via
#               Laplace approximation + importance sampling (Poisson likelihood
#               is non-conjugate, so the integrated likelihood is approximate)
#       * theta1* : SIR + Chan sampler given accepted phi1
#       * W2  : collapsed MH with CE proposal, marginalizing theta2 via
#               EXACT closed-form Gaussian marginal likelihood (no approximation
#               needed: theta2 never touches the Poisson likelihood, only the
#               linear-Gaussian evolution of theta1, so integrating it out is exact)
#       * theta2  : precision sampler (Chan), given accepted phi2
#
# Both collapsed samplers use the SAME mechanism (evaluate the integrated
# likelihood at the conditional mode via the Chan sparse Cholesky factor),
# but W1's version requires IS because the Poisson likelihood makes the
# posterior of theta1 non-Gaussian; W2's version is exact because the
# posterior of theta2 given theta1 is exactly Gaussian.
#
# Author: Cleiton Moya de Almeida

library(Matrix)
library(coda)

rm(list = ls())
options(error = function() traceback(2))
tp <- base::t
set.seed(42)

setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))

# Load the data
source_file <- "poisson_pol2_200"
data <- readRDS(paste("../data/", source_file, ".rds", sep = ""))
y <- data$y

Tt <- length(y)
if (Tt == 200) t_obs <- c(50, 100, 150, 175)
if (Tt == 2000) t_obs <- c(500, 1000, 1500, 1750)

theta1_present <- TRUE
theta2_present <- FALSE

if (theta1_present) {
    theta1_true <- data$theta
    lambda_true <- exp(theta1_true)
}
if (theta2_present)
    theta2_true <- data$theta2


#####
# Auxiliary functions

printf <- function(...) cat(paste(sprintf(...), "\n"))



#####
# Static sparse matrix for the Chan Method

start_time = proc.time()

# Base for the prior Precision Matrix K
sub_diag_base <- rep(-1, Tt-1)
main_diag_base <- c(rep(2, Tt-1), 1)
K0 <- bandSparse(n=Tt, k=c(0, -1),
                 diagonals=list(main_diag_base, sub_diag_base),
                 symmetric = TRUE)

# diagonal mask
diag_pattern <- bandSparse(n=Tt, k=c(0, -1),
                           diagonals=list(rep(TRUE, Tt), rep(FALSE, Tt-1)),
                           symmetric=TRUE)
idx_diag <- which(diag_pattern@x)

# subdiagonal mask
sub_pattern <- bandSparse(n=Tt, k=c(0, -1),
                          diagonals=list(rep(FALSE, Tt), rep(TRUE, Tt-1)),
                          symmetric=TRUE)
idx_sub <- which(sub_pattern@x)

# Initial symbolic Cholesky factors
Ch01_factor <- Cholesky(K0, perm = FALSE, LDL = TRUE)
Ch02_factor <- Cholesky(K0, perm = FALSE, LDL = TRUE)

# log|K0| (constant, precomputed once - used in the exact W2 marginal likelihood)
log_det_K0 <- 2 * as.numeric(determinant(Ch01_factor, logarithm = TRUE)$modulus)

# Work precision matrices (static)
P1_matrix <- K0
P2_matrix <- K0

building_time <- (proc.time()- start_time)[[1]]
printf("Sparse structures building: %.4f s", building_time)


#####
# Chan Method functions

chan_smoothing_theta1 <- function(y, phi_V, phi1, theta_01, theta_02, theta2) {
    Tt <- length(y)
    P1_matrix@x[idx_diag] <- (main_diag_base * phi1) + phi_V
    P1_matrix@x[idx_sub]  <- -phi1
    Ch1_factor <- update(Ch01_factor, P1_matrix)

    b <- y * phi_V
    Hb_theta2 <- numeric(Tt)
    Hb_theta2[1] <- -theta2[1]
    Hb_theta2[2:(Tt-1)] <- theta2[1:(Tt-2)] - theta2[2:(Tt-1)]
    Hb_theta2[Tt] <- theta2[Tt-1]
    b <- b + phi1*Hb_theta2
    b[1] <- b[1] + phi1*(theta_01 + theta_02)

    theta1_hat <- as.numeric(Matrix::solve(Ch1_factor, b, system="A"))
    list(theta1_hat=theta1_hat, ch=Ch1_factor)
}


chan_sample_theta1 <- function(build_res) {
    d <- Matrix::diag(build_res$ch)
    u <- rnorm(Tt)
    w <- u / sqrt(d)
    x <- as.vector(Matrix::solve(build_res$ch, w, system="Lt"))
    build_res$theta1_hat + x
}


# Builds the precision system for theta2 given phi1, phi2, theta_02.
# Returns theta2_hat (posterior mode) and the Cholesky factor (for sampling
# and for the exact log-determinant used in the collapsed W2 likelihood).
chan_smoothing_theta2 <- function(theta1, phi1, phi2, theta_02) {
    z <- diff(theta1)   # z_t = theta1[t+1] - theta1[t], t=1,...,T-1

    diag_obs <- c(rep(phi1, Tt-1), 0)
    P2_matrix@x[idx_diag] <- (main_diag_base*phi2) + diag_obs
    P2_matrix@x[idx_sub]  <- -phi2

    Ch2_factor <- update(Ch02_factor, P2_matrix)

    b <- numeric(Tt)
    b[1:(Tt-1)] <- z * phi1
    b[1] <- b[1] + theta_02 * phi2

    theta2_hat <- as.numeric(Matrix::solve(Ch2_factor, b, system="A"))
    list(theta2_hat = theta2_hat, ch = Ch2_factor, z = z)
}


chan_sample_theta2 <- function(build_res) {
    d <- Matrix::diag(build_res$ch)
    u <- rnorm(Tt)
    w <- u/sqrt(d)
    x <- as.vector(Matrix::solve(build_res$ch, w, system="Lt"))
    build_res$theta2_hat + x
}


# IRLS: build Laplace approximation around theta1_tilde
run_irls <- function(theta1_tilde, theta2, theta_01, theta_02, phi1,
                     y, tol, M_irls_max) {
    for (j in 1:M_irls_max) {
        f_t   <- exp(-theta1_tilde)
        phi_V <- 1 / f_t
        z_t   <- theta1_tilde + f_t * y - 1   # Poisson pseudo-observation

        res <- chan_smoothing_theta1(z_t, phi_V, phi1, theta_01, theta_02, theta2)
        theta1_tilde_new <- res$theta1_hat
        if (any(!is.finite(theta1_tilde_new))) break

        if (max(abs(theta1_tilde_new - theta1_tilde)) < tol) {
            theta1_tilde <- theta1_tilde_new
            break
        }
        theta1_tilde <- theta1_tilde_new
    }
    return(list(res = res, theta1_tilde = theta1_tilde, itr = j,
         f_t = f_t, phi_V = phi_V))
}


# ---------------------------------------------------------------------------
# IS estimator of log p(y | phi1, theta2)  [integrated over theta1]
# Approximate: the Poisson likelihood makes p(theta1|y,phi1,theta2) non-Gaussian,
# so this uses importance sampling around the Laplace/IRLS mode.
# ---------------------------------------------------------------------------
is_log_lik <- function(irls_res, y, phi1, theta2, theta_01, theta_02, M_is_lik) {
    res     <- irls_res$res
    eta_hat <- res$theta1_hat
    ch      <- res$ch
    W1 <- 1/phi1

    log_det_H <- 2 * as.numeric(determinant(ch, logarithm = TRUE)$modulus)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])
    log_norm_H <- -Tt / 2 * log(2 * pi) + 0.5 * log_det_H

    log_w <- numeric(M_is_lik)
    u_mat <- matrix(rnorm(M_is_lik * Tt), nrow = M_is_lik)
    d_ch  <- Matrix::diag(ch)

    for (i in 1:M_is_lik) {
        u <- u_mat[i, ]
        w <- u / sqrt(d_ch)
        x <- as.vector(solve(ch, w, system = "Lt"))
        th <- eta_hat + x

        log_py <- sum(y * th - exp(th))
        th_lag1 <- c(theta_01, th[-Tt])
        eps     <- th - th_lag1 - th_lag2_fixed
        log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
        log_q <- log_norm_H - 0.5 * sum(u^2)

        log_w[i] <- log_py + log_prior_th - log_q
    }

    log_lik <- logsumexp(log_w) - log(M_is_lik)
    w_norm  <- exp(log_w - logsumexp(log_w))
    ess_is  <- 1 / sum(w_norm^2)

    return(list(log_lik = log_lik, ess_is = ess_is))
}


# ---------------------------------------------------------------------------
# EXACT integrated likelihood log p(z | phi2, phi1, theta_02) [integrated over theta2]
# Exact: z|theta2 is Gaussian (linear evolution, no Poisson nonlinearity), so
# integrating theta2 out analytically has closed form. Evaluated via the
# Chan sparse Cholesky factor at the conditional mode theta2_hat:
#
#   log p(z|phi2) = log p(z|theta2_hat) + log p(theta2_hat|phi2) - log q(theta2_hat)
#
# where q(.) = N(theta2_hat, H^{-1}) is the EXACT conditional posterior
# (not an approximation), so this identity is exact for any phi2.
# ---------------------------------------------------------------------------
log_marginal_lik_w2 <- function(theta1, phi1, phi2, theta_02) {
    build <- chan_smoothing_theta2(theta1, phi1, phi2, theta_02)
    theta2_hat <- build$theta2_hat
    ch <- build$ch
    z  <- build$z

    log_pz <- -0.5*(Tt-1)*log(2*pi/phi1) - 0.5*phi1*sum((z - theta2_hat[1:(Tt-1)])^2)

    diffs2 <- theta2_hat - c(theta_02, theta2_hat[-Tt])
    log_p_theta2 <- -0.5*Tt*log(2*pi/phi2) + 0.5*log_det_K0 - 0.5*phi2*sum(diffs2^2)

    log_det_H <- 2 * as.numeric(determinant(ch, logarithm = TRUE)$modulus)
    log_q <- -0.5*Tt*log(2*pi) + 0.5*log_det_H

    list(log_lik = log_pz + log_p_theta2 - log_q, build = build)
}


# Cross-entropy calibration: fit Gamma(c, d) to samples via moment matching
calibrate_ce_gamma <- function(samples) {
    m1    <- mean(samples)
    m_log <- mean(log(samples))
    rhs   <- log(m1) - m_log

    c_hat <- if (rhs <= 0.5772) {
        (3 - rhs + sqrt((rhs - 3)^2 + 24 * rhs)) / (12 * rhs)
    } else {
        1 / rhs
    }

    for (k in 1:100) {
        f_val <- log(c_hat) - digamma(c_hat) - rhs
        fp    <- 1 / c_hat - trigamma(c_hat)
        c_new <- c_hat - f_val / fp
        if (!is.finite(c_new) || c_new <= 0) break
        if (abs(c_new - c_hat) < 1e-10) { c_hat <- c_new; break }
        c_hat <- c_new
    }
    d_hat <- c_hat / m1
    return(list(shape = c_hat, rate = d_hat))
}


# Collapsed MH step for W1 (independence chain with CE proposal, approximate IS likelihood)
mh_w1_collapsed <- function(phi1_cur, log_lik_cur, irls_cur,
                            ce_params, nu_01, eta_01,
                            theta2, theta_01, theta_02,
                            theta1_tilde, y, tol, M_irls_max, M_is_lik) {

    phi1_prop <- rgamma(1, shape = ce_params$shape, rate = ce_params$rate)

    irls_prop    <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                             phi1_prop, y, tol, M_irls_max)
    res_lik_prop <- is_log_lik(irls_prop, y, phi1_prop, theta2, theta_01, theta_02, M_is_lik)
    log_lik_prop <- res_lik_prop$log_lik

    log_prior <- function(phi) dgamma(phi, shape = nu_01, rate = eta_01, log = TRUE)
    log_prop  <- function(phi) dgamma(phi, shape = ce_params$shape, rate = ce_params$rate, log = TRUE)

    log_alpha <- (log_lik_prop + log_prior(phi1_prop) + log_prop(phi1_cur)) -
        (log_lik_cur  + log_prior(phi1_cur)  + log_prop(phi1_prop))

    if (log(runif(1)) < log_alpha) {
        return(list(phi1 = phi1_prop, log_lik = log_lik_prop,
                    irls_res = irls_prop, accepted = TRUE))
    } else {
        return(list(phi1 = phi1_cur, log_lik = log_lik_cur,
                    irls_res = irls_cur, accepted = FALSE))
    }
}


# Collapsed MH step for W2 (independence chain with CE proposal, EXACT likelihood)
mh_w2_collapsed <- function(phi2_cur, log_lik2_cur, theta2_build_cur,
                            ce2_params, nu_02, eta_02,
                            theta1_star, phi1, theta_02) {

    phi2_prop <- rgamma(1, shape = ce2_params$shape, rate = ce2_params$rate)

    res_prop <- log_marginal_lik_w2(theta1_star, phi1, phi2_prop, theta_02)
    log_lik2_prop <- res_prop$log_lik

    log_prior <- function(phi) dgamma(phi, shape = nu_02, rate = eta_02, log = TRUE)
    log_prop  <- function(phi) dgamma(phi, shape = ce2_params$shape, rate = ce2_params$rate, log = TRUE)

    log_alpha <- (log_lik2_prop + log_prior(phi2_prop) + log_prop(phi2_cur)) -
        (log_lik2_cur  + log_prior(phi2_cur)  + log_prop(phi2_prop))

    if (!is.finite(log_alpha)) {
        return(list(phi2 = phi2_cur, log_lik2 = log_lik2_cur,
                    build = theta2_build_cur, accepted = FALSE))
    }

    if (log(runif(1)) < log_alpha) {
        return(list(phi2 = phi2_prop, log_lik2 = log_lik2_prop,
                    build = res_prop$build, accepted = TRUE))
    } else {
        return(list(phi2 = phi2_cur, log_lik2 = log_lik2_cur,
                    build = theta2_build_cur, accepted = FALSE))
    }
}


# Sampling Importance Resampling for theta1
sir_theta1 <- function(irls_res, y, phi1, theta2, theta_01, theta_02, M_sir_theta1) {
    eta_hat <- irls_res$res$theta1_hat
    ch      <- irls_res$res$ch
    W1 <- 1/phi1

    log_det_H     <- 2 * as.numeric(determinant(ch, logarithm = TRUE)$modulus)
    th_lag2_fixed <- c(theta_02, theta2[-Tt])
    log_norm_H    <- -Tt / 2 * log(2 * pi) + 0.5 * log_det_H

    log_w  <- numeric(M_sir_theta1)
    draws  <- matrix(0, M_sir_theta1, Tt)
    d_ch   <- Matrix::diag(ch)

    for (i in 1:M_sir_theta1) {
        u  <- rnorm(Tt)
        w <- u / sqrt(d_ch)
        x  <- as.vector(Matrix::solve(ch, w, system = "Lt"))
        th <- eta_hat + x
        draws[i, ] <- th

        log_py       <- sum(y * th - exp(th))
        th_lag1      <- c(theta_01, th[-Tt])
        eps          <- th - th_lag1 - th_lag2_fixed
        log_prior_th <- -Tt / 2 * log(2 * pi * W1) - sum(eps^2) / (2 * W1)
        log_q        <- log_norm_H - 0.5 * sum(u^2)

        log_w[i] <- log_py + log_prior_th - log_q
    }

    w_norm  <- exp(log_w - logsumexp(log_w))
    ess_sir <- 1 / sum(w_norm^2)
    idx     <- sample.int(M_sir_theta1, size = 1, prob = w_norm)

    list(theta1 = draws[idx, ], ess = ess_sir)
}


# Prior hyperparameters
mu_01     <- 0
sigma2_01 <- 100

mu_02     <- 0
sigma2_02 <- 100

nu_01 <- 2
eta_01  <- 0.01

nu_02 <- 2
eta_02  <- 0.0001


# Simulation parameters
R_prerun   <- 3000   # pre-run iterations to calibrate CE proposals (phi1 and phi2)
N <- 10000            # Gibbs iterations
burnin <- 1000
M_is_lik <- 3         # Number of particles - IS for W1's integrated likelihood
M_sir_theta1 <- 3     # Number of particles - SIR of theta1
M_irls_max <- 20
tol <- 1e-4


# Initial values (shared for both pre-run and main run)
theta1_star <- numeric(Tt)
theta2 <- numeric(Tt)
theta1_tilde <- 0
theta_01 <- 0
theta_02 <- 0
W1 <- 0.01
phi1 <- 1/W1
W2 <- 0.01
phi2 <- 1/W2

# Auxiliary variables
W1_hist <- numeric(N)
W2_hist <- numeric(N)
theta_01_hist <- numeric(N)
theta_02_hist <- numeric(N)
theta1_hist <- matrix(0, N, Tt)
theta2_hist <- matrix(0, N, Tt)
accepted_hist  <- logical(N)   # W1 MH acceptance
accepted2_hist <- logical(N)   # W2 MH acceptance
itr_irls <- numeric(N)
ess_is_hist  <- numeric(N)     # IS of W1's integrated likelihood
ess_sir_hist <- numeric(N)     # SIR of theta1


#####
# PRE-RUN: standard Gibbs to collect phi1 AND phi2 samples for CE calibration

printf("Starting pre-run (%d iterations)...", R_prerun)
phi1_prerun <- numeric(R_prerun)
phi2_prerun <- numeric(R_prerun)
time_prerun <- proc.time()

for (n in 1:R_prerun) {

    # Sample theta_02 (conjugated Normal)
    sigma2_02_bar <- (1 / sigma2_02 + 1 / W1 + 1 / W2)^(-1)
    mu_02_bar <- sigma2_02_bar * ((theta1_star[1] - theta_01) / W1 +
                                      theta2[1] / W2 + mu_02 / sigma2_02)
    theta_02 <- rnorm(1, mean = mu_02_bar, sd = sqrt(sigma2_02_bar))

    # Sample theta_01 (conjugated Normal)
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01 / sigma2_01 +
                                      (theta1_star[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))

    # Sample phi1 (conjugated Gamma)
    dif1   <- theta1_star - c(theta_01, theta1_star[-Tt])
    diffs1 <- dif1 - c(theta_02, theta2[-Tt])
    nu_01_bar <- nu_01 + Tt / 2
    eta_01_bar  <- eta_01 + 0.5 * sum(diffs1^2)
    phi1 <- rgamma(1, shape = nu_01_bar, rate = eta_01_bar)
    W1 <- 1/phi1
    phi1_prerun[n] <- phi1

    # Sample phi2 (conjugated Gamma)
    diffs2 <- theta2 - c(theta_02, theta2[-Tt])
    nu_02_bar <- nu_02 + Tt / 2
    eta_02_bar  <- eta_02 + 0.5 * sum(diffs2^2)
    phi2 <- rgamma(1, shape = nu_02_bar, rate = eta_02_bar)
    W2 <- 1/phi2
    phi2_prerun[n] <- phi2

    # IRLS + Chan for theta1*
    irls_out <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                         phi1, y, tol, M_irls_max)
    theta1_tilde <- irls_out$theta1_tilde
    theta1_star  <- chan_sample_theta1(irls_out$res)

    # Chan sampler for theta2
    theta2_build <- chan_smoothing_theta2(theta1_star, phi1, phi2, theta_02)
    theta2 <- chan_sample_theta2(theta2_build)
}

elapsed_prerun <- (proc.time() - time_prerun)[[1]]
printf("Pre-run done in %.0f s", elapsed_prerun)

#####
# CE Calibration (both phi1 and phi2)
ce_params  <- calibrate_ce_gamma(phi1_prerun)
ce2_params <- calibrate_ce_gamma(phi2_prerun)

printf("CE Gamma proposal for phi1: shape = %.4f, rate = %.4f (mean phi1 = %.6f, mean W1 = %.6f)",
       ce_params$shape, ce_params$rate,
       ce_params$shape / ce_params$rate,
       ce_params$rate / ce_params$shape)

printf("CE Gamma proposal for phi2: shape = %.4f, rate = %.4f (mean phi2 = %.6f, mean W2 = %.6f)",
       ce2_params$shape, ce2_params$rate,
       ce2_params$shape / ce2_params$rate,
       ce2_params$rate / ce2_params$shape)


#####
# Gibbs Loop

# Compute initial IRLS and log integrated likelihood (W1)
irls_cur     <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                         phi1, y, tol, M_irls_max)
theta1_tilde <- irls_cur$theta1_tilde

res_lik <- is_log_lik(irls_cur, y, phi1, theta2, theta_01, theta_02, M_is_lik)
log_lik_cur <- res_lik$log_lik

# Compute initial exact marginal likelihood (W2)
res_lik2 <- log_marginal_lik_w2(theta1_star, phi1, phi2, theta_02)
log_lik2_cur <- res_lik2$log_lik
theta2_build_cur <- res_lik2$build

printf("Starting main MCMC (%d iterations)...", N)
start_time <- proc.time()

for (n in 1:N) {

    if (n %% 1000 == 0) {
        elapsed <- (proc.time() - start_time)[[1]]
        acc_rate  <- mean(accepted_hist[1:(n - 1)])
        acc_rate2 <- mean(accepted2_hist[1:(n - 1)])
        printf("Iter %d / %d | Elapsed: %.0f s | W1 accept: %.2f | W2 accept: %.2f",
               n, N, elapsed, acc_rate, acc_rate2)
    }


    # Sample theta_02 (conjugated Normal)
    sigma2_02_bar <- (1 / sigma2_02 + 1 / W1 + 1 / W2)^(-1)
    mu_02_bar <- sigma2_02_bar * ((theta1_star[1] - theta_01) / W1 +
                                      theta2[1] / W2 + mu_02 / sigma2_02)
    theta_02 <- rnorm(1, mean = mu_02_bar, sd = sqrt(sigma2_02_bar))


    # Sample theta_01 (conjugated Normal)
    sigma2_01_bar <- (1 / sigma2_01 + 1 / W1)^(-1)
    mu_01_bar <- sigma2_01_bar * (mu_01 / sigma2_01 +
                                      (theta1_star[1] - theta_02) / W1)
    theta_01 <- rnorm(1, mean = mu_01_bar, sd = sqrt(sigma2_01_bar))


    # Collapsed MH for W1 (approximate integrated likelihood, marginalizes theta1)
    mh_res <- mh_w1_collapsed(
        phi1_cur     = phi1,
        log_lik_cur  = log_lik_cur,
        irls_cur     = irls_cur,
        ce_params    = ce_params,
        nu_01     = nu_01,
        eta_01      = eta_01,
        theta2       = theta2,
        theta_01     = theta_01,
        theta_02     = theta_02,
        theta1_tilde = theta1_tilde,
        y            = y,
        tol          = tol,
        M_irls_max   = M_irls_max,
        M_is_lik     = M_is_lik
    )

    phi1 <- mh_res$phi1
    W1 <- 1 / phi1
    log_lik_cur <- mh_res$log_lik
    irls_cur <- mh_res$irls_res
    itr_irls[n] <- irls_cur$itr
    accepted_hist[n] <- mh_res$accepted


    # Collapsed MH for W2 (EXACT integrated likelihood, marginalizes theta2)
    # Recompute log_lik2_cur first: theta_02 and phi1 both may have changed
    # since the last time it was updated (theta_01/theta_02 sampled above,
    # phi1 updated by the W1 MH step just now).
    res_lik2_cur <- log_marginal_lik_w2(theta1_star, phi1, phi2, theta_02)
    log_lik2_cur <- res_lik2_cur$log_lik
    theta2_build_cur <- res_lik2_cur$build

    mh2_res <- mh_w2_collapsed(
        phi2_cur         = phi2,
        log_lik2_cur     = log_lik2_cur,
        theta2_build_cur = theta2_build_cur,
        ce2_params       = ce2_params,
        nu_02            = nu_02,
        eta_02           = eta_02,
        theta1_star      = theta1_star,
        phi1             = phi1,
        theta_02         = theta_02
    )

    phi2 <- mh2_res$phi2
    W2 <- 1 / phi2
    accepted2_hist[n] <- mh2_res$accepted


    # Sample theta2 via Chan sampler given accepted phi2
    # (must happen BEFORE SIR of theta1, so SIR sees the fresh theta2
    #  consistent with the just-accepted phi2 - matches the validated
    #  isolated test; doing SIR first would use a stale theta2 and
    #  weaken the phi2-theta1 coupling, degrading ESS)
    theta2_build <- chan_smoothing_theta2(theta1_star, phi1, phi2, theta_02)
    theta2 <- chan_sample_theta2(theta2_build)


    # Sample theta1* via SIR given accepted phi1 and the fresh theta2
    sir_res      <- sir_theta1(irls_cur, y, phi1, theta2, theta_01, theta_02, M_sir_theta1)
    theta1_star  <- sir_res$theta1
    theta1_tilde <- irls_cur$theta1_tilde
    ess_sir_hist[n] <- sir_res$ess


    # Update irls_cur and log_lik_cur for next iteration
    # (theta2, theta_01, theta_02 all changed in this iteration)
    irls_cur <- run_irls(theta1_tilde, theta2, theta_01, theta_02,
                             phi1, y, tol, M_irls_max)
    theta1_tilde <- irls_cur$theta1_tilde

    res_lik <- is_log_lik(irls_cur, y, phi1, theta2, theta_01, theta_02, M_is_lik)
    log_lik_cur <- res_lik$log_lik


    # Store the results
    theta_01_hist[n] <- theta_01
    theta_02_hist[n] <- theta_02
    W1_hist[n] <- W1
    W2_hist[n] <- W2
    theta1_hist[n, ] <- theta1_star
    theta2_hist[n, ] <- theta2
    ess_is_hist[n] <- res_lik$ess_is
}


#####
# Simulation summary

end_time <- proc.time()
sampling_time <- (end_time - time1)[[1]]
elapsed_time <- (end_time - start_time)[[1]]
printf("Sampling: %.2f s", sampling_time)
printf("Total CPU time: %.0f s", elapsed_time)

printf("W1 MH acceptance rate: %.3f", mean(accepted_hist))
printf("W2 MH acceptance rate: %.3f", mean(accepted2_hist))


# Posterior mean
theta1_mean <- colMeans(theta1_hist[-(1:burnin), ])
theta2_mean <- colMeans(theta2_hist[-(1:burnin), ])
lambda_mean <- exp(theta1_mean)

printf("W1 mean: %.5f", mean(W1_hist[-(1:burnin)]))
printf("W1 median: %.5f", median(W1_hist[-(1:burnin)]))
printf("W2 mean: %.5f", mean(W2_hist[-(1:burnin)]))
printf("W2 median: %.5f", median(W2_hist[-(1:burnin)]))


# Log-likelihood
loglik <- sum(dpois(y, lambda_mean, log=TRUE))
printf("Log-likelihood: %.2f", loglik)


# Effective sample size
printf("Effective Sample Size:")
ess_w1 <- effectiveSize(mcmc(W1_hist[-(1:burnin)]))
ess_w2 <- effectiveSize(mcmc(W2_hist[-(1:burnin)]))
printf("\tW1: %.0f", ess_w1)
printf("\tW2: %.0f", ess_w2)

ess_theta1 <- effectiveSize(mcmc(theta1_hist[-(1:burnin),]))
ess_theta2 <- effectiveSize(mcmc(theta2_hist[-(1:burnin),]))
printf("\ttheta1 (mean): %.2f", mean(ess_theta1))
printf("\ttheta2 (mean): %.2f", mean(ess_theta2))

# Effective sample size per second
printf("Effective Sample Size / second:")
printf("\tW1: %.2f", ess_w1/elapsed_time)
printf("\tW2: %.2f", ess_w2/elapsed_time)

ess_sec_theta1 <- ess_theta1/elapsed_time
ess_sec_theta2 <- ess_theta2/elapsed_time
printf("\ttheta1 (mean): %.2f", mean(ess_sec_theta1))
printf("\ttheta2 (mean): %.2f", mean(ess_sec_theta2))


# Geweke diagnostic: Z test for two mean difference
#   H0: segments same means -> chain has converged
printf("Geweke convergence diagnostic")
z_w1 <- unname(geweke.diag(W1_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
z_w2 <- unname(geweke.diag(W2_hist[-(1:burnin)], frac1=0.1, frac2=0.5)[[1]])
printf("\tz_w1: %.2f", z_w1)
printf("\tz_w2: %.2f", z_w2)

# Percent of instants in the H_0 rejection region:
z_theta1 <- unname(geweke.diag(theta1_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z_theta2 <- unname(geweke.diag(theta2_hist[-(1:burnin),], frac1=0.1, frac2=0.5)[[1]])
z1_out <- sum((z_theta1 < -1.96) | (z_theta1 > 1.96))/Tt
z2_out <- sum((z_theta2 < -1.96) | (z_theta2 > 1.96))/Tt
printf("\tPercent of theta1 out: %.3f", z1_out)
printf("\tPercent of theta2 out: %.3f", z2_out)


#####
# Plots
# y, lambda_true, lambda_estimated
x <- 1:Tt
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(x, y, type="l", xlab="t", ylab="", col="gray",
     main="Poisson Local Trend Polynomial Model")
points(x, y, pch = 20)
lines(x, lambda_mean, col="red", lwd=2)
if (theta1_present) {
    lines(x, lambda_true, col="blue", lwd=2)
    legend("topright",
           legend = expression(y[t], lambda[t], hat(lambda)[t]),
           col = c("black", "blue", "red"),
           lty = c(NA, 1, 1),
           lwd = c(NA, 2, 2),
           pch = c(20, NA, NA),
           bty = "n")
} else {
    legend("topright",
           legend = expression(y[t], hat(lambda)[t]),
           col = c("black", "red"),
           lty = c(NA, 1),
           lwd = c(NA, 2),
           pch = c(20, NA),
           bty = "n")
}


# theta1_true, theta1_mean
x <- 1:Tt
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
if (theta1_present) {
    ylim_range <- range(theta1_mean, theta1_true)
} else {
    ylim_range <- range(theta1_mean)
}
plot(x, theta1_mean, type="l", col="red", lwd=2, ylim=ylim_range,
     xlab="t", ylab="", main="theta_t1")
if (theta1_present) {
    lines(x, theta1_true, col="blue", lwd=2)
    legend("topright", legend=expression(hat(theta)[t1], theta[t1]),
           col=c("red","blue"), lwd=2, bty="n")
} else {
    legend("topright", legend=expression(hat(theta)[t1]),
           col="red", lwd=2, bty="n")
}


# theta2_true, theta2_mean
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
if (theta2_present) {
    ylim_range <- range(theta2_mean, theta2_true)
} else {
    ylim_range <- range(theta2_mean)
}
plot(x, theta2_mean, type="l", col="red", lwd=2, ylim=ylim_range,
     xlab="t", ylab="", main="theta_t2")
if (theta2_present) {
    lines(x, theta2_true, col="blue", lwd=2)
    legend("topright", legend=expression(hat(theta)[t2], theta[t2]),
           col=c("red","blue"), lwd=2, bty="n")
} else {
    legend("topright", legend=expression(hat(theta)[t2]),
           col="red", lwd=2, bty="n")
}


# Posterior distribution of theta_t1
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta1_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 1]),
         main = bquote("Posterior of " * theta[.(t) * "," * 1]))
    lines(density(theta1_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of theta_t2
par(mfrow = c(2, 2))
for (t in t_obs) {
    hist(theta2_hist[-(1:burnin), t], breaks = 50, freq = FALSE,
         xlab = bquote(theta[.(t) * "," * 2]),
         main = bquote("Posterior of " * theta[.(t) * "," * 2]))
    lines(density(theta2_hist[-(1:burnin), t]), col = "blue", lwd = 2)
}


# Posterior distribution of W1
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W1_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W1")
lines(density(W1_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Posterior distribution of W2
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex = 0.8)
hist(W2_hist[-(1:burnin)], breaks = 50, freq = FALSE, main ="Posterior of W2")
lines(density(W2_hist[-(1:burnin)]), col = "blue", lwd = 2)


# Traceplot for W1 and W2
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex = 0.8)
plot(W1_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W1")
abline(v=burnin, col="red")
plot(W2_hist[-(1:burnin)], type="l", xlab="n", ylab="W", main="Traceplot of W2")
abline(v=burnin, col="red")


# Traceplots for theta_t1
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta1_hist[, t], type="l", main=bquote(theta[.(t)*","*1]), xlab="", ylab="")
    abline(v=burnin, col="red")
}


# Traceplots for theta_t2
par(mfrow = c(2, 2))
for (t in t_obs) {
    plot(theta2_hist[, t], type="l", main=bquote(theta[.(t)*","*2]), xlab="", ylab="")
    abline(v=burnin, col="red")
}


# Effective sample size
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(ess_theta1, type="l", main=expression("Effective sample of " * theta[t1]), xlab="t")
plot(ess_theta2, type="l", main=expression("Effective sample of " * theta[t2]), xlab="t")

par(mfrow = c(1, 2), mar = c(4, 4, 2, 2), cex=0.8)
hist(ess_theta1)
hist(ess_theta2)


# Geweke diagnostic
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(z_theta1, type="l", main=expression("Geweke diagnostic for " * theta[t1]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")

plot(z_theta2, type="l", main=expression("Geweke diagnostic for " * theta[t2]),
     xlab="t", ylab="Z score")
abline(h=c(-1.96, 1.96), col="red")


# ACF for theta1 e theta2
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8)
for (t in t_obs) {
    acf(theta1_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*1]))
    acf(theta2_hist[-(1:burnin), t], main=bquote(theta[.(t)*","*2]))
}


# Prior vs posterior for phi1 and phi2 (comparing against the CE proposals)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
curve(dgamma(x, shape=nu_01, rate=eta_01), from=1e-6, to=max(1/W1_hist[-(1:burnin)]),
      main="phi1 prior vs. posterior", col="red", lwd=2)
lines(density(1/W1_hist[-(1:burnin)]), col="blue", lwd=2)
legend("topright", legend=c("Prior","Posterior"), col=c("red","blue"), lwd=2)

par(mfrow = c(1, 1), mar = c(4, 4, 2, 2), cex=0.8)
curve(dgamma(x, shape=nu_02, rate=eta_02), from=1e-6, to=max(1/W2_hist[-(1:burnin)]),
      main="phi2 prior vs. posterior", col="red", lwd=2)
lines(density(1/W2_hist[-(1:burnin)]), col="blue", lwd=2)
legend("topright", legend=c("Prior","Posterior"), col=c("red","blue"), lwd=2)


# Effective Sample Size - IS of W1's Integrated Likelihood
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_is_hist, type="l", main="ESS - IS Integrated Likelihood (W1)", xlab="n")


# Effective Sample Size - SIR of theta1
par(mfrow=c(1,1), mar=c(4,4,2,2), cex=0.8)
plot(ess_sir_hist, type="l", main="ESS - SIR theta1", xlab="n")


# W1 and W2 acceptance rate over time (rolling window)
roll_mean <- function(x, k) {
    n <- length(x)
    out <- rep(NA, n)
    for (i in k:n) out[i] <- mean(x[(i-k+1):i])
    out
}
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2), cex=0.8)
plot(roll_mean(as.numeric(accepted_hist), 200), type="l",
     main="W1 rolling acceptance rate (window=200)", xlab="n", ylab="rate")
plot(roll_mean(as.numeric(accepted2_hist), 200), type="l",
     main="W2 rolling acceptance rate (window=200)", xlab="n", ylab="rate")
