// sir_collapsed.cpp
// Collapsed Gibbs sampler for the Poisson second-order polynomial DLM
// (local trend model), with BOTH W1 and W2 collapsed (marginalized) out of
// their respective Metropolis-Hastings steps:
//   - phi1 = 1/W1: collapsed MH with a Cross-Entropy-calibrated Gamma
//     independence proposal. The target log-likelihood log p(y | phi1,
//     theta2) [integrated over (theta_01, theta1) JOINTLY] is itself only
//     available via Importance Sampling on top of a Laplace/IRLS Gaussian
//     approximation (is_log_lik), since the Poisson likelihood makes the
//     exact integral intractable.
//   - phi2 = 1/W2: collapsed MH with a CE-calibrated Gamma independence
//     proposal. The target log-likelihood log p(z | phi1, phi2, theta_02)
//     [integrated over theta2, theta_02 held fixed] IS available in closed
//     form (log_marginal_lik_w2), since this half of the model is exactly
//     Gaussian -- no IS needed here.
//   - (theta_02, theta2): sampled jointly, exact, via the Chan-method
//     tridiagonal LDL^T block (same as in the other samplers).
//   - (theta_01, theta1): sampled jointly via SIR (Sampling Importance
//     Resampling) on top of the same Laplace/IRLS approximation used by
//     is_log_lik (sir_theta1).
//
// The Cross-Entropy Gamma proposals are calibrated once, from a short
// PRE-RUN of R_prerun iterations of a simple (non-collapsed) Gibbs sampler
// that draws phi1/phi2 from their ordinary conjugate full conditionals
// (gibbs_sample_phi1/phi2, from utils.h) -- see calibrate_ce_gamma.
//
// RNG-ORDER NOTE (important for R/C++ parity): is_log_lik draws its
// M_is_lik x Ttp1 standard normals ALL AT ONCE (mirroring R's `u_mat <-
// matrix(rnorm(M_is_lik * Ttp1), nrow = M_is_lik)`, filled COLUMN-MAJOR),
// so the draw order is: for each of the Ttp1 "columns", draw M_is_lik
// values in a row. sir_theta1, in contrast, draws a fresh block of Ttp1
// normals per importance sample (a plain sequential loop) -- this
// asymmetry exists in the original R code and is preserved here exactly,
// not "fixed" to be consistent, since the two functions must stay
// RNG-compatible with their R counterparts.
//
// Helpers shared with the other samplers (logsumexp, log_p_yt, log_dnorm,
// gibbs_sample_theta01/phi1/phi2, sample_index_from_logw, and the Chan-
// method tridiagonal LDL^T solver) live in utils.h and MUST stay identical
// across all converted algorithms -- do not fork/duplicate them here.
//
// Author of the original R code: Cleiton Moya de Almeida
// C++/Rcpp port for maximum performance.
//
// [[Rcpp::plugins(cpp11)]]
#include "utils.h"

// ---------------------------------------------------------------------
// IRLS: Laplace approximation for (theta_01, theta1) jointly, given a
// fixed phi1 and theta2 (T-node linearization, extended (T+1)-dim block).
// Mirrors chan_smoothing_theta1 + run_irls in sir_collapsed.R.
// ---------------------------------------------------------------------
struct IrlsResult {
    TriLDLT fac;                      // last factorization (size Ttp1)
    std::vector<double> theta1_hat;   // mode, length Ttp1 (node 0 = theta_01)
    std::vector<double> theta1_tilde; // updated linearization point, length Tt
    int itr;
};

static IrlsResult run_irls(std::vector<double> theta1_tilde, // by value: local working copy
                            const std::vector<double>& theta2, double theta_02, double phi1,
                            const NumericVector& y, double tol, int M_irls_max,
                            double mu_01, double sigma2_01, int Tt, int Ttp1) {
    std::vector<double> d1(Ttp1), e1(Ttp1 - 1), b1(Ttp1);
    std::vector<double> f_t(Tt), z_t(Tt);
    std::vector<double> RHS_ext(Tt), Hb_ext(Ttp1);
    std::vector<double> theta1_hat(Ttp1);
    TriLDLT fac1;

    int j = M_irls_max;
    for (int jj = 1; jj <= M_irls_max; ++jj) {
        for (int t = 0; t < Tt; ++t) {
            f_t[t] = std::exp(-theta1_tilde[t]);
            z_t[t] = theta1_tilde[t] + f_t[t] * y[t] - 1.0;
        }

        d1[0] = phi1 + 1.0 / sigma2_01;
        for (int i = 1; i < Tt; ++i) d1[i] = 2.0 * phi1 + 1.0 / f_t[i - 1];
        d1[Tt] = phi1 + 1.0 / f_t[Tt - 1];
        for (int i = 0; i < Ttp1 - 1; ++i) e1[i] = -phi1;

        RHS_ext[0] = theta_02;
        for (int i = 1; i < Tt; ++i) RHS_ext[i] = theta2[i - 1];
        Hb_ext[0] = -RHS_ext[0];
        for (int i = 1; i < Tt; ++i) Hb_ext[i] = RHS_ext[i - 1] - RHS_ext[i];
        Hb_ext[Tt] = RHS_ext[Tt - 1];

        b1[0] = mu_01 / sigma2_01 + phi1 * Hb_ext[0];
        for (int i = 1; i <= Tt; ++i) b1[i] = z_t[i - 1] / f_t[i - 1] + phi1 * Hb_ext[i];

        fac1 = tri_ldlt_factor(d1, e1);
        tri_solve_A(fac1, b1, theta1_hat);

        bool all_finite = true;
        for (int t = 0; t < Tt; ++t) {
            if (!std::isfinite(theta1_hat[t + 1])) { all_finite = false; break; }
        }
        if (!all_finite) { j = jj; break; } // theta1_tilde left UNCHANGED, matches R

        double maxdiff = 0.0;
        for (int t = 0; t < Tt; ++t) {
            double ad = std::fabs(theta1_hat[t + 1] - theta1_tilde[t]);
            if (ad > maxdiff) maxdiff = ad;
        }
        if (maxdiff < tol) {
            for (int t = 0; t < Tt; ++t) theta1_tilde[t] = theta1_hat[t + 1];
            j = jj;
            break;
        }
        for (int t = 0; t < Tt; ++t) theta1_tilde[t] = theta1_hat[t + 1];
    }

    IrlsResult out;
    out.fac = fac1;
    out.theta1_hat = theta1_hat;
    out.theta1_tilde = theta1_tilde;
    out.itr = j;
    return out;
}

// ---------------------------------------------------------------------
// IS estimator of log p(y | phi1, theta2) [integrated over (theta_01,
// theta1) JOINTLY]. Mirrors is_log_lik in sir_collapsed.R -- including its
// upfront, column-major-order draw of the M_is_lik x Ttp1 standard normals
// (see the RNG-ORDER NOTE above).
// ---------------------------------------------------------------------
struct LogLikResult { double log_lik; double ess_is; };

static LogLikResult is_log_lik(const IrlsResult& irls, const NumericVector& y, double phi1,
                                const std::vector<double>& theta2, double theta_02, int M_is_lik,
                                double mu_01, double sigma2_01, int Tt, int Ttp1) {
    const TriLDLT& fac = irls.fac;
    const std::vector<double>& eta_hat = irls.theta1_hat;
    double W1 = 1.0 / phi1;

    double log_det_H = 0.0;
    for (int i = 0; i < Ttp1; ++i) log_det_H += std::log(fac.D[i]);

    std::vector<double> th_lag2_fixed(Tt);
    th_lag2_fixed[0] = theta_02;
    for (int i = 1; i < Tt; ++i) th_lag2_fixed[i] = theta2[i - 1];

    double log_norm_H = -0.5 * Ttp1 * std::log(2.0 * M_PI) + 0.5 * log_det_H;

    // u_mat[i][j], column-major like R's matrix(rnorm(M_is_lik*Ttp1), nrow=M_is_lik):
    // draw order is "for each column j, draw M_is_lik values" (see file header note).
    std::vector<double> u_mat((size_t)M_is_lik * Ttp1);
    for (int j = 0; j < Ttp1; ++j)
        for (int i = 0; i < M_is_lik; ++i)
            u_mat[(size_t)j * M_is_lik + i] = R::norm_rand();

    std::vector<double> w(Ttp1), x(Ttp1), draw(Ttp1);
    std::vector<double> log_w(M_is_lik);

    for (int i = 0; i < M_is_lik; ++i) {
        double sumu2 = 0.0;
        for (int j = 0; j < Ttp1; ++j) {
            double uj = u_mat[(size_t)j * M_is_lik + i];
            sumu2 += uj * uj;
            w[j] = uj / std::sqrt(fac.D[j]);
        }
        tri_solve_Lt(fac, w, x);
        for (int k = 0; k < Ttp1; ++k) draw[k] = eta_hat[k] + x[k];

        double theta_01_i = draw[0];
        double log_py = 0.0;
        for (int t = 0; t < Tt; ++t) {
            double th_t = draw[t + 1];
            log_py += y[t] * th_t - std::exp(th_t);
        }
        double log_prior_01 = -0.5 * std::log(2.0 * M_PI * sigma2_01)
                               - (theta_01_i - mu_01) * (theta_01_i - mu_01) / (2.0 * sigma2_01);

        double sq_eps = 0.0;
        double th_lag1_prev = theta_01_i;
        for (int t = 0; t < Tt; ++t) {
            double th_t = draw[t + 1];
            double eps = th_t - th_lag1_prev - th_lag2_fixed[t];
            sq_eps += eps * eps;
            th_lag1_prev = th_t;
        }
        double log_prior_th = -0.5 * Tt * std::log(2.0 * M_PI * W1) - sq_eps / (2.0 * W1);
        double log_q = log_norm_H - 0.5 * sumu2;

        log_w[i] = log_py + log_prior_01 + log_prior_th - log_q;
    }

    double lse = logsumexp(log_w);
    double log_lik = lse - std::log((double)M_is_lik);

    double s2 = 0.0;
    for (int i = 0; i < M_is_lik; ++i) {
        double wn = std::exp(log_w[i] - lse);
        s2 += wn * wn;
    }
    LogLikResult out;
    out.log_lik = log_lik;
    out.ess_is = 1.0 / s2;
    return out;
}

// ---------------------------------------------------------------------
// SIR for (theta_01, theta1) jointly. Mirrors sir_theta1 in
// sir_collapsed.R -- a fresh block of Ttp1 normals is drawn PER importance
// sample (plain sequential loop, unlike is_log_lik's upfront matrix draw).
// ---------------------------------------------------------------------
struct SirResult { double theta_01; std::vector<double> theta1; double ess; };

static SirResult sir_theta1(const IrlsResult& irls, const NumericVector& y, double phi1,
                             const std::vector<double>& theta2, double theta_02, int M_sir,
                             double mu_01, double sigma2_01, int Tt, int Ttp1,
                             std::vector<double>& samp_w_buf) {
    const TriLDLT& fac = irls.fac;
    const std::vector<double>& eta_hat = irls.theta1_hat;
    double W1 = 1.0 / phi1;

    double log_det_H = 0.0;
    for (int i = 0; i < Ttp1; ++i) log_det_H += std::log(fac.D[i]);

    std::vector<double> th_lag2_fixed(Tt);
    th_lag2_fixed[0] = theta_02;
    for (int i = 1; i < Tt; ++i) th_lag2_fixed[i] = theta2[i - 1];

    double log_norm_H = -0.5 * Ttp1 * std::log(2.0 * M_PI) + 0.5 * log_det_H;

    std::vector<double> log_w(M_sir);
    std::vector<double> draws((size_t)M_sir * Ttp1);
    std::vector<double> u(Ttp1), w(Ttp1), x(Ttp1);

    for (int i = 0; i < M_sir; ++i) {
        double sumu2 = 0.0;
        for (int k = 0; k < Ttp1; ++k) {
            double uk = R::norm_rand();
            sumu2 += uk * uk;
            w[k] = uk / std::sqrt(fac.D[k]);
        }
        tri_solve_Lt(fac, w, x);
        for (int k = 0; k < Ttp1; ++k) draws[(size_t)i * Ttp1 + k] = eta_hat[k] + x[k];

        double theta_01_i = draws[(size_t)i * Ttp1 + 0];
        double log_py = 0.0;
        for (int t = 0; t < Tt; ++t) {
            double th_t = draws[(size_t)i * Ttp1 + t + 1];
            log_py += y[t] * th_t - std::exp(th_t);
        }
        double log_prior_01 = -0.5 * std::log(2.0 * M_PI * sigma2_01)
                               - (theta_01_i - mu_01) * (theta_01_i - mu_01) / (2.0 * sigma2_01);

        double sq_eps = 0.0;
        double th_lag1_prev = theta_01_i;
        for (int t = 0; t < Tt; ++t) {
            double th_t = draws[(size_t)i * Ttp1 + t + 1];
            double eps = th_t - th_lag1_prev - th_lag2_fixed[t];
            sq_eps += eps * eps;
            th_lag1_prev = th_t;
        }
        double log_prior_th = -0.5 * Tt * std::log(2.0 * M_PI * W1) - sq_eps / (2.0 * W1);
        double log_q = log_norm_H - 0.5 * sumu2;

        log_w[i] = log_py + log_prior_01 + log_prior_th - log_q;
    }

    double lse = logsumexp(log_w);
    double s2 = 0.0;
    for (int i = 0; i < M_sir; ++i) {
        double wn = std::exp(log_w[i] - lse);
        s2 += wn * wn;
    }

    int idx = sample_index_from_logw(log_w.data(), M_sir, samp_w_buf);

    SirResult out;
    out.theta_01 = draws[(size_t)idx * Ttp1 + 0];
    out.theta1.resize(Tt);
    for (int t = 0; t < Tt; ++t) out.theta1[t] = draws[(size_t)idx * Ttp1 + t + 1];
    out.ess = 1.0 / s2;
    return out;
}

// ---------------------------------------------------------------------
// EXACT integrated likelihood log p(z | phi1, phi2, theta_02) [integrated
// over theta2, theta_02 held FIXED]: a T-node (non-extended, "anchored")
// tridiagonal block -- mirrors log_marginal_lik_w2 in sir_collapsed.R.
// log_det_K0 is the precomputed (phi-independent) log|K0| constant for
// this T-node chain's unit-precision skeleton.
// ---------------------------------------------------------------------
static double log_marginal_lik_w2(const std::vector<double>& theta1, double phi1, double phi2,
                                   double theta_02, double log_det_K0, int Tt) {
    std::vector<double> z(Tt - 1);
    for (int i = 0; i < Tt - 1; ++i) z[i] = theta1[i + 1] - theta1[i];

    std::vector<double> d(Tt), e(Tt - 1), b(Tt, 0.0), theta2_hat(Tt);
    for (int i = 0; i < Tt - 1; ++i) d[i] = 2.0 * phi2 + phi1;
    d[Tt - 1] = phi2;
    for (int i = 0; i < Tt - 1; ++i) e[i] = -phi2;

    for (int i = 0; i < Tt - 1; ++i) b[i] = z[i] * phi1;
    b[0] += theta_02 * phi2;

    TriLDLT fac = tri_ldlt_factor(d, e);
    tri_solve_A(fac, b, theta2_hat);

    double sse_z = 0.0;
    for (int i = 0; i < Tt - 1; ++i) {
        double diff = z[i] - theta2_hat[i];
        sse_z += diff * diff;
    }
    double log_pz = -0.5 * (Tt - 1) * std::log(2.0 * M_PI / phi1) - 0.5 * phi1 * sse_z;

    double sse_theta2 = 0.0;
    double prev = theta_02;
    for (int i = 0; i < Tt; ++i) {
        double diffs2 = theta2_hat[i] - prev;
        sse_theta2 += diffs2 * diffs2;
        prev = theta2_hat[i];
    }
    double log_p_theta2 = -0.5 * Tt * std::log(2.0 * M_PI / phi2) + 0.5 * log_det_K0 - 0.5 * phi2 * sse_theta2;

    double log_det_H = 0.0;
    for (int i = 0; i < Tt; ++i) log_det_H += std::log(fac.D[i]);
    double log_q = -0.5 * Tt * std::log(2.0 * M_PI) + 0.5 * log_det_H;

    return log_pz + log_p_theta2 - log_q;
}

// ---------------------------------------------------------------------
// Cross-Entropy calibration of a Gamma(shape, rate) proposal from a
// vector of samples (moment/digamma matching + Newton-Raphson refinement).
// Mirrors calibrate_ce_gamma in sir_collapsed.R.
// ---------------------------------------------------------------------
struct GammaParams { double shape; double rate; };

static GammaParams calibrate_ce_gamma(const std::vector<double>& samples) {
    int n = (int) samples.size();
    double m1 = 0.0, m_log = 0.0;
    for (double s : samples) { m1 += s; m_log += std::log(s); }
    m1 /= n;
    m_log /= n;
    double rhs = std::log(m1) - m_log;

    double c_hat;
    if (rhs <= 0.5772) {
        c_hat = (3.0 - rhs + std::sqrt((rhs - 3.0) * (rhs - 3.0) + 24.0 * rhs)) / (12.0 * rhs);
    } else {
        c_hat = 1.0 / rhs;
    }

    for (int k = 0; k < 100; ++k) {
        double f_val = std::log(c_hat) - R::digamma(c_hat) - rhs;
        double fp = 1.0 / c_hat - R::trigamma(c_hat);
        double c_new = c_hat - f_val / fp;
        if (!std::isfinite(c_new) || c_new <= 0.0) break;
        if (std::fabs(c_new - c_hat) < 1e-10) { c_hat = c_new; break; }
        c_hat = c_new;
    }

    GammaParams out;
    out.shape = c_hat;
    out.rate = c_hat / m1;
    return out;
}

// ---------------------------------------------------------------------
// Main SIR-Collapsed sampler
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List sir_collapsed_cpp(NumericVector y,
                        int R_prerun, int N,
                        double mu_01, double sigma2_01,
                        double mu_02, double sigma2_02,
                        double nu_01, double eta_01,
                        double nu_02, double eta_02,
                        NumericVector theta1, NumericVector theta2,
                        NumericVector theta1_tilde,
                        double theta_01, double theta_02,
                        double W1, double W2,
                        int M_is_lik, int M_sir_theta1, int M_irls_max, double tol,
                        bool verbose = true, int print_every = 1000) {

    RNGScope rngScope; // sync with R's RNG / set.seed

    const int Tt = y.size();
    const int Ttp1 = Tt + 1;

    // ---- histories (row-major access via NumericMatrix, R-compatible) ----
    // Declared here (outer scope, not inside the nested block below) so
    // they remain accessible at the final return.
    NumericMatrix theta1_hist(N, Tt), theta2_hist(N, Tt);
    NumericVector theta_01_hist(N), theta_02_hist(N);
    NumericVector W1_hist(N), W2_hist(N);
    LogicalVector accepted_hist(N), accepted2_hist(N);
    IntegerVector itr_irls(N);
    NumericVector ess_is_hist(N), ess_sir_hist(N);
    double ce1_shape = 0.0, ce1_rate = 0.0, ce2_shape = 0.0, ce2_rate = 0.0;

    // ---- state ----
    // theta_01, theta_02, W1, W2 are already local (by-value) copies of the
    // R-supplied initial values, so they are mutated in place below -- no
    // separate "_init" parameter/copy pair needed.
    // theta1/theta2/theta1_tilde need a type conversion (NumericVector ->
    // std::vector<double>); see pg_as.cpp for why this needs a nested
    // block (a parameter's scope is the function's outer block, so a
    // same-named local there would be an illegal redeclaration, not mere
    // shadowing).
    {
    std::vector<double> theta1_tmp(theta1.begin(), theta1.end());
    std::vector<double> theta1_star = std::move(theta1_tmp); // R's "theta1_star"
    std::vector<double> theta2_tmp(theta2.begin(), theta2.end());
    std::vector<double> theta2_v = std::move(theta2_tmp);
    std::vector<double> theta1_tilde_tmp(theta1_tilde.begin(), theta1_tilde.end());
    std::vector<double> theta1_tilde_v = std::move(theta1_tilde_tmp);

    double phi1 = 1.0 / W1;
    double phi2 = 1.0 / W2;

    std::vector<double> samp_w_buf(std::max(M_sir_theta1, 2));

    // Precompute log_det_K0: T-node ANCHORED chain, unit precision (phi-
    // independent constant, used by log_marginal_lik_w2).
    std::vector<double> d_K0(Tt), e_K0(Tt - 1);
    for (int i = 0; i < Tt - 1; ++i) d_K0[i] = 2.0;
    d_K0[Tt - 1] = 1.0;
    for (int i = 0; i < Tt - 1; ++i) e_K0[i] = -1.0;
    TriLDLT fac_K0 = tri_ldlt_factor(d_K0, e_K0);
    double log_det_K0 = 0.0;
    for (int i = 0; i < Tt; ++i) log_det_K0 += std::log(fac_K0.D[i]);

    // Chan (theta_02, theta2) extended-block buffers (size Ttp1)
    std::vector<double> d2(Ttp1), e2(Ttp1 - 1), b2(Ttp1), theta2_hat(Ttp1);
    std::vector<double> u2(Ttp1), w2(Ttp1), xsol2(Ttp1), draw2(Ttp1);
    std::vector<double> zdiff(Tt - 1);

    // ---- PRE-RUN: simple (non-collapsed) Gibbs, to calibrate the CE
    // Gamma proposals for phi1 and phi2 ----
    if (verbose) Rcout << "Starting pre-run (" << R_prerun << " iterations)...\n";
    std::vector<double> phi1_prerun(R_prerun), phi2_prerun(R_prerun);

    for (int r = 0; r < R_prerun; ++r) {
        phi1 = gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1_star, theta_02, theta2_v, Tt);
        W1 = 1.0 / phi1;
        phi1_prerun[r] = phi1;

        phi2 = gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2_v, Tt);
        W2 = 1.0 / phi2;
        phi2_prerun[r] = phi2;

        IrlsResult irls_out = run_irls(theta1_tilde_v, theta2_v, theta_02, phi1, y, tol, M_irls_max,
                                        mu_01, sigma2_01, Tt, Ttp1);
        theta1_tilde_v = irls_out.theta1_tilde;
        SirResult res_sir = sir_theta1(irls_out, y, phi1, theta2_v, theta_02, M_sir_theta1,
                                        mu_01, sigma2_01, Tt, Ttp1, samp_w_buf);
        theta_01 = res_sir.theta_01;
        theta1_star = res_sir.theta1;

        // (theta_02, theta2) jointly via extended block
        for (int i = 0; i < Tt - 1; ++i) zdiff[i] = theta1_star[i + 1] - theta1_star[i];
        d2[0] = phi2 + (1.0 / sigma2_02 + phi1);
        for (int i = 1; i < Tt; ++i) d2[i] = 2.0 * phi2 + phi1;
        d2[Tt] = phi2;
        for (int i = 0; i < Ttp1 - 1; ++i) e2[i] = -phi2;
        b2[0] = mu_02 / sigma2_02 + phi1 * (theta1_star[0] - theta_01);
        for (int i = 1; i < Tt; ++i) b2[i] = zdiff[i - 1] * phi1;
        b2[Tt] = 0.0;
        TriLDLT fac2 = tri_ldlt_factor(d2, e2);
        tri_solve_A(fac2, b2, theta2_hat);
        for (int i = 0; i < Ttp1; ++i) u2[i] = R::norm_rand();
        for (int i = 0; i < Ttp1; ++i) w2[i] = u2[i] / std::sqrt(fac2.D[i]);
        tri_solve_Lt(fac2, w2, xsol2);
        for (int i = 0; i < Ttp1; ++i) draw2[i] = theta2_hat[i] + xsol2[i];
        theta_02 = draw2[0];
        for (int i = 0; i < Tt; ++i) theta2_v[i] = draw2[i + 1];
    }

    // ---- CE Calibration ----
    GammaParams ce_params = calibrate_ce_gamma(phi1_prerun);
    GammaParams ce2_params = calibrate_ce_gamma(phi2_prerun);
    ce1_shape = ce_params.shape; ce1_rate = ce_params.rate;
    ce2_shape = ce2_params.shape; ce2_rate = ce2_params.rate;
    if (verbose) {
        Rcout << "CE Gamma phi1: shape=" << ce_params.shape << " rate=" << ce_params.rate << "\n";
        Rcout << "CE Gamma phi2: shape=" << ce2_params.shape << " rate=" << ce2_params.rate << "\n";
    }

    // ---- Setup before main loop ----
    IrlsResult irls_cur = run_irls(theta1_tilde_v, theta2_v, theta_02, phi1, y, tol, M_irls_max,
                                    mu_01, sigma2_01, Tt, Ttp1);
    theta1_tilde_v = irls_cur.theta1_tilde;
    LogLikResult res_lik = is_log_lik(irls_cur, y, phi1, theta2_v, theta_02, M_is_lik,
                                       mu_01, sigma2_01, Tt, Ttp1);
    double log_lik_cur = res_lik.log_lik;

    if (verbose) Rcout << "Starting main MCMC (" << N << " iterations)...\n";

    // ---- MAIN LOOP ----
    for (int n = 0; n < N; ++n) {

        if (verbose && ((n + 1) % print_every == 0)) {
            Rcout << "Iter " << (n + 1) << " / " << N << "\n";
        }

        // ---- Collapsed MH for W1 (phi1) ----
        double phi1_prop = R::rgamma(ce_params.shape, 1.0 / ce_params.rate);
        IrlsResult irls_prop = run_irls(theta1_tilde_v, theta2_v, theta_02, phi1_prop, y, tol, M_irls_max,
                                         mu_01, sigma2_01, Tt, Ttp1);
        LogLikResult res_lik_prop = is_log_lik(irls_prop, y, phi1_prop, theta2_v, theta_02, M_is_lik,
                                                mu_01, sigma2_01, Tt, Ttp1);
        double log_lik_prop = res_lik_prop.log_lik;

        double log_alpha1 =
            (log_lik_prop + R::dgamma(phi1_prop, nu_01, 1.0 / eta_01, 1)
                          + R::dgamma(phi1, ce_params.shape, 1.0 / ce_params.rate, 1))
          - (log_lik_cur  + R::dgamma(phi1, nu_01, 1.0 / eta_01, 1)
                          + R::dgamma(phi1_prop, ce_params.shape, 1.0 / ce_params.rate, 1));

        bool accept1 = (std::log(unif_rand()) < log_alpha1);
        if (accept1) {
            phi1 = phi1_prop;
            log_lik_cur = log_lik_prop;
            irls_cur = irls_prop;
        }
        W1 = 1.0 / phi1;
        itr_irls[n] = irls_cur.itr;
        accepted_hist[n] = accept1;

        // ---- Collapsed MH for W2 (phi2) ----
        double log_lik2_cur = log_marginal_lik_w2(theta1_star, phi1, phi2, theta_02, log_det_K0, Tt);
        double phi2_prop = R::rgamma(ce2_params.shape, 1.0 / ce2_params.rate);
        double log_lik2_prop = log_marginal_lik_w2(theta1_star, phi1, phi2_prop, theta_02, log_det_K0, Tt);

        double log_alpha2 =
            (log_lik2_prop + R::dgamma(phi2_prop, nu_02, 1.0 / eta_02, 1)
                           + R::dgamma(phi2, ce2_params.shape, 1.0 / ce2_params.rate, 1))
          - (log_lik2_cur  + R::dgamma(phi2, nu_02, 1.0 / eta_02, 1)
                           + R::dgamma(phi2_prop, ce2_params.shape, 1.0 / ce2_params.rate, 1));

        bool accept2;
        if (!std::isfinite(log_alpha2)) {
            accept2 = false; // matches R: no unif_rand() draw consumed in this branch
        } else {
            accept2 = (std::log(unif_rand()) < log_alpha2);
        }
        if (accept2) phi2 = phi2_prop;
        W2 = 1.0 / phi2;
        accepted2_hist[n] = accept2;

        // ---- (theta_02, theta2) jointly via Chan method, using theta1_star
        // (still the value from the PREVIOUS iteration at this point) ----
        for (int i = 0; i < Tt - 1; ++i) zdiff[i] = theta1_star[i + 1] - theta1_star[i];
        d2[0] = phi2 + (1.0 / sigma2_02 + phi1);
        for (int i = 1; i < Tt; ++i) d2[i] = 2.0 * phi2 + phi1;
        d2[Tt] = phi2;
        for (int i = 0; i < Ttp1 - 1; ++i) e2[i] = -phi2;
        b2[0] = mu_02 / sigma2_02 + phi1 * (theta1_star[0] - theta_01);
        for (int i = 1; i < Tt; ++i) b2[i] = zdiff[i - 1] * phi1;
        b2[Tt] = 0.0;
        TriLDLT fac2 = tri_ldlt_factor(d2, e2);
        tri_solve_A(fac2, b2, theta2_hat);
        for (int i = 0; i < Ttp1; ++i) u2[i] = R::norm_rand();
        for (int i = 0; i < Ttp1; ++i) w2[i] = u2[i] / std::sqrt(fac2.D[i]);
        tri_solve_Lt(fac2, w2, xsol2);
        for (int i = 0; i < Ttp1; ++i) draw2[i] = theta2_hat[i] + xsol2[i];
        theta_02 = draw2[0];
        for (int i = 0; i < Tt; ++i) theta2_v[i] = draw2[i + 1];

        // ---- (theta_01, theta1) jointly via SIR, using irls_cur (built
        // with the OLD theta2, but reweighted with the FRESH theta2_v just
        // sampled above -- see the file header note on this ordering) ----
        SirResult res_sir = sir_theta1(irls_cur, y, phi1, theta2_v, theta_02, M_sir_theta1,
                                        mu_01, sigma2_01, Tt, Ttp1, samp_w_buf);
        theta_01 = res_sir.theta_01;
        theta1_star = res_sir.theta1;
        theta1_tilde_v = irls_cur.theta1_tilde;
        ess_sir_hist[n] = res_sir.ess;

        // ---- Refresh irls_cur/log_lik_cur for the NEXT iteration's phi1
        // MH step (uses the fresh theta1_tilde/theta2/theta_01/theta_02) ----
        irls_cur = run_irls(theta1_tilde_v, theta2_v, theta_02, phi1, y, tol, M_irls_max,
                             mu_01, sigma2_01, Tt, Ttp1);
        theta1_tilde_v = irls_cur.theta1_tilde;

        res_lik = is_log_lik(irls_cur, y, phi1, theta2_v, theta_02, M_is_lik,
                              mu_01, sigma2_01, Tt, Ttp1);
        log_lik_cur = res_lik.log_lik;
        ess_is_hist[n] = res_lik.ess_is;

        // ---- Store results ----
        theta_01_hist[n] = theta_01;
        theta_02_hist[n] = theta_02;
        W1_hist[n] = W1;
        W2_hist[n] = W2;
        for (int t = 0; t < Tt; ++t) {
            theta1_hist(n, t) = theta1_star[t];
            theta2_hist(n, t) = theta2_v[t];
        }
    }
    } // end of nested scope (theta1_star/theta2_v/theta1_tilde_v as std::vector<double>)

    return List::create(
        _["theta1_hist"] = theta1_hist,
        _["theta2_hist"] = theta2_hist,
        _["theta_01_hist"] = theta_01_hist,
        _["theta_02_hist"] = theta_02_hist,
        _["W1_hist"] = W1_hist,
        _["W2_hist"] = W2_hist,
        _["accepted_hist"] = accepted_hist,
        _["accepted2_hist"] = accepted2_hist,
        _["itr_irls"] = itr_irls,
        _["ess_is_hist"] = ess_is_hist,
        _["ess_sir_hist"] = ess_sir_hist,
        _["ce1_shape"] = ce1_shape,
        _["ce1_rate"] = ce1_rate,
        _["ce2_shape"] = ce2_shape,
        _["ce2_rate"] = ce2_rate
    );
}
