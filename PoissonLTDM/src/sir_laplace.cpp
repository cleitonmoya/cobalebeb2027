// sir_laplace.cpp
// Sampling Importance Resampling with a Laplace (Gaussian) proposal for the
// Poisson second-order polynomial DLM (local trend model):
//   - theta1: Importance Sampling with a Laplace/IRLS approximation as the
//     importance density. The mode/precision of the proposal are found by
//     iterating the Chan-method tridiagonal block solve on a pseudo-linear
//     Gaussian model (Schnatter/IRLS pseudo-observations z_t, heteroscedastic
//     pseudo-precision phi_V_t = 1/f_t), warm-started from the previous
//     Gibbs iteration's converged linearization point. M_is > 1 draws M_is
//     proposals and importance-resamples one; M_is == 1 skips reweighting
//     (plain draw from the Laplace approximation).
//   - theta2: sampled jointly with theta_02, exact, via the Chan-method
//     tridiagonal LDL^T block (same as in the other samplers).
//
// theta1's Chan block has a DIFFERENT diagonal structure than theta2's (a
// heteroscedastic, per-t precision phi_V_t from the IRLS linearization,
// instead of a constant phi1) -- this is not a shared helper in utils.h,
// just a different (d, e, b) assembly using the same generic TriLDLT solver.
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
// Main SIR-Laplace sampler
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List sir_laplace_cpp(NumericVector y,
                      int N,
                      double mu_01, double sigma2_01,
                      double mu_02, double sigma2_02,
                      double nu_01, double eta_01,
                      double nu_02, double eta_02,
                      NumericVector theta1, NumericVector theta2,
                      double theta_01, double theta_02,
                      double W1, double W2,
                      int M_irls_max, int M_is, double tol,
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
    NumericVector ess_is(N);
    IntegerVector itr_irls(N);

    // ---- state ----
    // theta_01, theta_02, W1, W2 are already local (by-value) copies of the
    // R-supplied initial values, so they are mutated in place below -- no
    // separate "_init" parameter/copy pair needed.
    // theta1/theta2 need a type conversion (NumericVector -> std::vector<
    // double>); see pg_as.cpp for why this needs a nested block (a
    // parameter's scope is the function's outer block, so a same-named
    // local there would be an illegal redeclaration, not mere shadowing).
    {
    std::vector<double> theta1_tmp(theta1.begin(), theta1.end());
    std::vector<double> theta1 = std::move(theta1_tmp);
    std::vector<double> theta2_tmp(theta2.begin(), theta2.end());
    std::vector<double> theta2 = std::move(theta2_tmp);

    // IRLS linearization point, persisted (warm-started) ACROSS Gibbs
    // iterations -- mirrors the R script, where theta1_tilde is declared
    // once before the main loop and never reset inside it (independent of
    // theta1's own initial value).
    std::vector<double> theta1_tilde(Tt, 0.0);
    std::vector<double> theta1_tilde_old(Tt, 0.0);

    // theta1's Chan block buffers (size Ttp1) -- kept across IRLS
    // iterations so the LAST factorization/mode is available for the
    // Importance Sampling step below.
    std::vector<double> d1(Ttp1), e1(Ttp1 - 1), b1(Ttp1), theta1_hat(Ttp1);
    std::vector<double> f_t(Tt), z_t(Tt);
    std::vector<double> RHS_ext(Tt), Hb_ext(Ttp1);

    // theta2's Chan block buffers (size Ttp1)
    std::vector<double> d2(Ttp1), e2(Ttp1 - 1), b2(Ttp1), theta2_hat(Ttp1);
    std::vector<double> zdiff(Tt - 1);

    // shared draw/solve scratch (size Ttp1)
    std::vector<double> u(Ttp1), w(Ttp1), xsol(Ttp1), draw(Ttp1);

    // Importance Sampling scratch (only used when M_is > 1)
    std::vector<double> log_w(M_is);
    std::vector<double> theta_01_traj(M_is);
    std::vector<double> theta1_traj(M_is * (size_t)Tt);
    std::vector<double> samp_w_buf(M_is);

    for (int n = 0; n < N; ++n) {

        if (verbose && ((n + 1) % print_every == 0)) {
            Rcout << "Iteration " << (n + 1) << " / " << N << "\n";
        }

        // ---- Sample phi1 -> W1 ----
        double phi1 = gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1, theta_02, theta2, Tt);
        W1 = 1.0 / phi1;

        // ---- Sample phi2 -> W2 ----
        double phi2 = gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt);
        W2 = 1.0 / phi2;

        // ---- IRLS: Laplace approximation for (theta_01, theta1) ----
        theta1_tilde_old = theta1_tilde;
        TriLDLT fac1;
        int j = M_irls_max;
        for (int jj = 1; jj <= M_irls_max; ++jj) {

            for (int t = 0; t < Tt; ++t) {
                f_t[t] = std::exp(-theta1_tilde[t]);
                z_t[t] = theta1_tilde[t] + f_t[t] * y[t] - 1.0;
            }

            // extended (T+1)-dimensional block: node 0 = theta_01, nodes
            // 1..Tt = theta_11, ..., theta_T1. main_diag_rw = (1, 2,...,2, 1);
            // extra_diag = (1/sigma2_01, phi_V_1, ..., phi_V_Tt), phi_V_t = 1/f_t[t].
            d1[0] = phi1 + 1.0 / sigma2_01;
            for (int i = 1; i < Tt; ++i) d1[i] = 2.0 * phi1 + 1.0 / f_t[i - 1];
            d1[Tt] = phi1 + 1.0 / f_t[Tt - 1];
            for (int i = 0; i < Ttp1 - 1; ++i) e1[i] = -phi1;

            // RHS_ext = (theta_02, theta2[1], ..., theta2[Tt-1]) (0-based:
            // theta2 indices 0..Tt-2), length Tt.
            RHS_ext[0] = theta_02;
            for (int i = 1; i < Tt; ++i) RHS_ext[i] = theta2[i - 1];
            Hb_ext[0] = -RHS_ext[0];
            for (int i = 1; i < Tt; ++i) Hb_ext[i] = RHS_ext[i - 1] - RHS_ext[i];
            Hb_ext[Tt] = RHS_ext[Tt - 1];

            b1[0] = mu_01 / sigma2_01 + phi1 * Hb_ext[0];
            for (int i = 1; i <= Tt; ++i) b1[i] = z_t[i - 1] / f_t[i - 1] + phi1 * Hb_ext[i];

            fac1 = tri_ldlt_factor(d1, e1);
            tri_solve_A(fac1, b1, theta1_hat);

            for (int t = 0; t < Tt; ++t) theta1_tilde[t] = theta1_hat[t + 1];

            double maxdiff = 0.0;
            for (int t = 0; t < Tt; ++t) {
                double ad = std::fabs(theta1_tilde[t] - theta1_tilde_old[t]);
                if (ad > maxdiff) maxdiff = ad;
            }
            if (maxdiff < tol) { j = jj; break; }
            theta1_tilde_old = theta1_tilde;
        }
        itr_irls[n] = j;

        // ---- Importance Sampling step for (theta_01, theta1) ----
        if (M_is > 1) {
            for (int i = 0; i < M_is; ++i) {
                for (int k = 0; k < Ttp1; ++k) u[k] = R::norm_rand();
                for (int k = 0; k < Ttp1; ++k) w[k] = u[k] / std::sqrt(fac1.D[k]);
                tri_solve_Lt(fac1, w, xsol);
                for (int k = 0; k < Ttp1; ++k) draw[k] = theta1_hat[k] + xsol[k];

                theta_01_traj[i] = draw[0];
                double log_p = 0.0, log_g = 0.0;
                for (int t = 0; t < Tt; ++t) {
                    double th = draw[t + 1];
                    theta1_traj[(size_t)i * Tt + t] = th;
                    log_p += y[t] * th - std::exp(th); // no isfinite guard, matches R exactly
                    log_g += log_dnorm(th, z_t[t], std::sqrt(f_t[t]));
                }
                log_w[i] = log_p - log_g;
            }

            double lse = logsumexp(log_w);
            double s2 = 0.0;
            for (int i = 0; i < M_is; ++i) {
                double wi = std::exp(log_w[i] - lse);
                s2 += wi * wi;
            }
            ess_is[n] = 1.0 / s2;

            int idx = sample_index_from_logw(log_w.data(), M_is, samp_w_buf);
            theta_01 = theta_01_traj[idx];
            for (int t = 0; t < Tt; ++t) theta1[t] = theta1_traj[(size_t)idx * Tt + t];
        } else {
            for (int k = 0; k < Ttp1; ++k) u[k] = R::norm_rand();
            for (int k = 0; k < Ttp1; ++k) w[k] = u[k] / std::sqrt(fac1.D[k]);
            tri_solve_Lt(fac1, w, xsol);
            for (int k = 0; k < Ttp1; ++k) draw[k] = theta1_hat[k] + xsol[k];
            theta_01 = draw[0];
            for (int t = 0; t < Tt; ++t) theta1[t] = draw[t + 1];
        }

        // ---- (theta_02, theta2) jointly via Chan method (tridiagonal LDL^T) ----
        for (int i = 0; i < Tt - 1; ++i) zdiff[i] = theta1[i + 1] - theta1[i];

        d2[0] = phi2 + (1.0 / sigma2_02 + phi1);
        for (int i = 1; i < Tt; ++i) d2[i] = 2.0 * phi2 + phi1;
        d2[Tt] = phi2;
        for (int i = 0; i < Ttp1 - 1; ++i) e2[i] = -phi2;

        b2[0] = mu_02 / sigma2_02 + phi1 * (theta1[0] - theta_01);
        for (int i = 1; i < Tt; ++i) b2[i] = zdiff[i - 1] * phi1;
        b2[Tt] = 0.0;

        TriLDLT fac2 = tri_ldlt_factor(d2, e2);
        tri_solve_A(fac2, b2, theta2_hat);

        for (int i = 0; i < Ttp1; ++i) u[i] = R::norm_rand();
        for (int i = 0; i < Ttp1; ++i) w[i] = u[i] / std::sqrt(fac2.D[i]);
        tri_solve_Lt(fac2, w, xsol);
        for (int i = 0; i < Ttp1; ++i) draw[i] = theta2_hat[i] + xsol[i];

        theta_02 = draw[0];
        for (int i = 0; i < Tt; ++i) theta2[i] = draw[i + 1];

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
    } // end of nested scope (theta1/theta2 as std::vector<double>)

    return List::create(
        _["theta1_hist"] = theta1_hist,
        _["theta2_hist"] = theta2_hist,
        _["theta_01_hist"] = theta_01_hist,
        _["theta_02_hist"] = theta_02_hist,
        _["W1_hist"] = W1_hist,
        _["W2_hist"] = W2_hist,
        _["ess_is"] = ess_is,
        _["itr_irls"] = itr_irls
    );
}
