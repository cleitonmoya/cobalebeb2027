# tests/algorithm_config.R
#
# Per-algorithm hyperparameters, simulation parameters, and initial values,
# shared by test_prototype_R.R, test_cpp.R, and test_validation.R. Returns
# a named list ready to be combined with `y` and passed via do.call() to
# either the _r() or _cpp() version of the algorithm -- both have the same
# argument names, so the same config works for both.
#
# quick_N, if not NULL, overrides N (and R_prerun for sir_collapsed) with a
# smaller value -- useful for test_validation.R, where a full production
# run isn't needed to confirm R and C++ agree.

get_algorithm_config <- function(algorithm, y, burnin = 1000, quick_N = NULL) {

    Tt <- length(y)

    common_hyper <- list(
        mu_01 = 0, sigma2_01 = 100,
        mu_02 = 0, sigma2_02 = 100
    )

    if (algorithm == "pg_as") {
        N <- if (is.null(quick_N)) 3000 else quick_N
        theta1 <- log(y + 0.5)
        theta2 <- c(diff(theta1), 0)
        cfg <- c(common_hyper, list(
            nu_01 = 2, eta_01 = 0.1, nu_02 = 2, eta_02 = 0.01,
            K = 30, N = N,
            theta1 = theta1, theta2 = theta2,
            theta_01 = theta1[1], theta_02 = theta2[1],
            W1 = var(diff(theta1)), W2 = var(diff(theta2))
        ))

    } else if (algorithm == "amh_montoril") {
        N <- if (is.null(quick_N)) 10000 else quick_N
        cfg <- c(common_hyper, list(
            nu_01 = 2, eta_01 = 0.01, nu_02 = 2, eta_02 = 0.0001,
            N = N,
            theta1 = numeric(Tt), theta2 = numeric(Tt),
            theta_01 = 0, theta_02 = 0, W1 = 0.01, W2 = 0.01,
            ac_ref = 0.44, varsigma2 = rep(0.02, Tt)
        ))

    } else if (algorithm == "sir_laplace") {
        N <- if (is.null(quick_N)) 10000 else quick_N
        cfg <- c(common_hyper, list(
            nu_01 = 2, eta_01 = 0.01, nu_02 = 2, eta_02 = 0.0001,
            N = N,
            theta1 = numeric(Tt), theta2 = numeric(Tt),
            theta_01 = 0, theta_02 = 0, W1 = 0.01, W2 = 0.01,
            M_irls_max = 20, M_is = 3, tol = 1e-4
        ))

    } else if (algorithm == "sir_collapsed") {
        N <- if (is.null(quick_N)) 10000 else quick_N
        R_prerun <- if (is.null(quick_N)) 3000 else max(quick_N, 50)
        cfg <- c(common_hyper, list(
            nu_01 = 2, eta_01 = 0.01, nu_02 = 2, eta_02 = 0.0001,
            R_prerun = R_prerun, N = N,
            theta1 = numeric(Tt), theta2 = numeric(Tt), theta1_tilde = numeric(Tt),
            theta_01 = 0, theta_02 = 0, W1 = 0.01, W2 = 0.01,
            M_is_lik = 3, M_sir_theta1 = 3, M_irls_max = 20, tol = 1e-4
        ))

    } else {
        stop(sprintf("Unknown algorithm: '%s'. Use one of: pg_as, amh_montoril, sir_laplace, sir_collapsed.", algorithm))
    }

    attr(cfg, "burnin") <- burnin
    cfg
}
