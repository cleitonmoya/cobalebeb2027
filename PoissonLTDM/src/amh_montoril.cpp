// amh_montoril.cpp
// Adaptive Metropolis-within-Gibbs (Roberts & Rosenthal, 2009) sampler for
// the Poisson second-order polynomial DLM (local trend model), following
// Montoril et al.'s component-wise scan for theta_t1:
//   - theta_t1, t = 1, ..., T: single-site random-walk Metropolis, with a
//     continuous Robbins-Monro update of log(varsigma_t) after each full
//     sweep (Andrieu & Thoms, 2008, Eq. 20/22), instead of the batch update
//     of Roberts & Rosenthal.
//   - (theta_02, theta2): sampled jointly, exact, via the Chan-method
//     tridiagonal LDL^T block (same as in the PG samplers).
//
// Unlike PG-AS/PG-APF, this is NOT a particle method: there is a single
// state path theta1 (no K particles, no resampling, no backward tracking).
//
// Helpers shared with the other samplers (logsumexp, log_p_yt, log_dnorm,
// gibbs_sample_theta01/phi1/phi2, logpost_theta_t1/logpost_theta_T1/
// sample_theta_t1_mh, and the Chan-method tridiagonal LDL^T solver) live in
// utils.h and MUST stay identical across all converted algorithms -- do not
// fork/duplicate them here.
//
// Author of the original R code: Cleiton Moya de Almeida
// C++/Rcpp port for maximum performance.
//
// [[Rcpp::plugins(cpp11)]]
#include "utils.h"

// ---------------------------------------------------------------------
// Main AMH-Montoril sampler
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List amh_montoril_cpp(NumericVector y,
                       int N,
                       double mu_01, double sigma2_01,
                       double mu_02, double sigma2_02,
                       double nu_01, double eta_01,
                       double nu_02, double eta_02,
                       NumericVector theta1, NumericVector theta2,
                       double theta_01, double theta_02,
                       double W1, double W2,
                       double ac_ref, NumericVector varsigma2,
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
    IntegerMatrix ac_hist(N, Tt);

    // ---- state ----
    // theta_01, theta_02, W1, W2 are already local (by-value) copies of the
    // R-supplied initial values, so they are mutated in place below -- no
    // separate "_init" parameter/copy pair needed.
    // theta1/theta2/varsigma2 need a type conversion (NumericVector ->
    // std::vector<double>); see pg_as.cpp for why this needs a nested
    // block (a parameter's scope is the function's outer block, so a
    // same-named local there would be an illegal redeclaration, not mere
    // shadowing).
    {
    std::vector<double> theta1_tmp(theta1.begin(), theta1.end());
    std::vector<double> theta1 = std::move(theta1_tmp);
    std::vector<double> theta2_tmp(theta2.begin(), theta2.end());
    std::vector<double> theta2 = std::move(theta2_tmp);
    std::vector<double> varsigma2_tmp(varsigma2.begin(), varsigma2.end());
    std::vector<double> varsigma2 = std::move(varsigma2_tmp);

    // Chan-method buffers (size Ttp1)
    std::vector<double> d(Ttp1), e(Ttp1 - 1), b(Ttp1), theta_hat(Ttp1);
    std::vector<double> u(Ttp1), w(Ttp1), xsol(Ttp1), draw2(Ttp1);
    std::vector<double> z(Tt - 1);

    for (int n = 0; n < N; ++n) {

        if (verbose && ((n + 1) % print_every == 0)) {
            Rcout << "Iteration " << (n + 1) << " / " << N << "\n";
        }

        // ---- Sample theta_01 (conjugate Normal, using the OLD W1) ----
        theta_01 = gibbs_sample_theta01(mu_01, sigma2_01, theta1[0], theta_02, W1);

        // ---- Sample phi1 -> W1 (using the NEW theta_01, OLD theta_02/theta2) ----
        double phi1 = gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1, theta_02, theta2, Tt);
        W1 = 1.0 / phi1;

        // ---- Sample phi2 -> W2 ----
        double phi2 = gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt);
        W2 = 1.0 / phi2;

        // ---- Sample theta_t1, t = 1, ..., Tt (component-wise random-walk MH) ----
        for (int t = 0; t < Tt; ++t) {
            MHStep res;
            if (t < Tt - 1) {
                if (t == 0) {
                    res = sample_theta_t1_mh(theta1[t], theta_01, theta1[t + 1],
                                              theta2[t], theta_02,
                                              y[t], W1, varsigma2[t], false);
                } else {
                    res = sample_theta_t1_mh(theta1[t], theta1[t - 1], theta1[t + 1],
                                              theta2[t], theta2[t - 1],
                                              y[t], W1, varsigma2[t], false);
                }
            } else {
                res = sample_theta_t1_mh(theta1[t], theta1[t - 1], 0.0,
                                          theta2[t], theta2[t - 1],
                                          y[t], W1, varsigma2[t], true);
            }
            theta1[t] = res.theta_t1;
            ac_hist(n, t) = res.ac;
        }

        // ---- Adaptive stage of varsigma2 (continuous Robbins-Monro update) ----
        double delta = std::min(0.01, 1.0 / std::sqrt((double)(n + 1)));
        for (int t = 0; t < Tt; ++t) {
            double ls = std::log(varsigma2[t]) / 2.0 + delta * ((double)ac_hist(n, t) - ac_ref);
            varsigma2[t] = std::exp(2.0 * ls);
        }

        // ---- (theta_02, theta2) jointly via Chan method (tridiagonal LDL^T) ----
        for (int i = 0; i < Tt - 1; ++i) z[i] = theta1[i + 1] - theta1[i];

        d[0] = phi2 + (1.0 / sigma2_02 + phi1);
        for (int i = 1; i < Tt; ++i) d[i] = 2.0 * phi2 + phi1;
        d[Tt] = phi2;
        for (int i = 0; i < Ttp1 - 1; ++i) e[i] = -phi2;

        b[0] = mu_02 / sigma2_02 + phi1 * (theta1[0] - theta_01);
        for (int i = 1; i < Tt; ++i) b[i] = z[i - 1] * phi1;
        b[Tt] = 0.0;

        TriLDLT fac = tri_ldlt_factor(d, e);
        tri_solve_A(fac, b, theta_hat);

        for (int i = 0; i < Ttp1; ++i) u[i] = R::norm_rand();
        for (int i = 0; i < Ttp1; ++i) w[i] = u[i] / std::sqrt(fac.D[i]);
        tri_solve_Lt(fac, w, xsol);
        for (int i = 0; i < Ttp1; ++i) draw2[i] = theta_hat[i] + xsol[i];

        theta_02 = draw2[0];
        for (int i = 0; i < Tt; ++i) theta2[i] = draw2[i + 1];

        // ---- Store results ----
        theta_01_hist[n] = theta_01;
        theta_02_hist[n] = theta_02;
        W1_hist[n] = W1;
        W2_hist[n] = W2;
        for (int t = 0; t < Tt; ++t) {
            theta1_hist(n, t) = theta1[t];
            theta2_hist(n, t) = theta2[t];
        }
    }
    } // end of nested scope (theta1/theta2/varsigma2 as std::vector<double>)

    return List::create(
        _["theta1_hist"] = theta1_hist,
        _["theta2_hist"] = theta2_hist,
        _["theta_01_hist"] = theta_01_hist,
        _["theta_02_hist"] = theta_02_hist,
        _["W1_hist"] = W1_hist,
        _["W2_hist"] = W2_hist,
        _["ac_hist"] = ac_hist
    );
}
