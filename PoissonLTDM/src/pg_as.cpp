// pg_as.cpp
// Particle Gibbs with Ancestral Sampling (PG-AS): Bootstrap Particle Filter
// with Ancestral Sampling for theta_t1, reconstructed by simple Backward
// Tracking (genealogy lookup, no Backward Simulation), and exact block
// (Chan-method) sampling of (theta_02, theta2) via a tridiagonal LDL^T
// factorization (the extended precision matrix is tridiagonal, so no
// sparse-matrix machinery is needed).
//
// Helpers shared with the other samplers (logsumexp, log_p_yt, log_dnorm,
// gibbs_sample_theta01/phi1/phi2, sample_index_from_logw,
// sample_indices_from_logw, and the Chan-method tridiagonal LDL^T solver)
// live in utils.h and MUST stay identical across all converted algorithms --
// do not fork/duplicate them here.
//
// PERFORMANCE NOTES (specific to PG-AS):
//  - theta_1_k is stored as a flat row-major buffer (index t*K + k) instead
//    of std::vector<std::vector<double>>, removing one level of pointer
//    indirection and improving cache locality.
//  - Unlike PG-APF (which needs the full T x K weight history for Backward
//    Simulation), PG-AS only ever reads log_w_tilde at t-1 (to build the
//    resampling/ancestral-sampling probabilities) and at the final time T
//    (to seed the backward tracking). So the weights are kept as two K-sized
//    buffers (log_w_prev, log_w_curr) swapped every t, instead of a T x K
//    matrix -- an O(K) memory footprint instead of O(K*T), specific to this
//    algorithm's structure.
//  - The resampling helpers (from pg_common.h) take pre-allocated scratch
//    buffers instead of allocating a new std::vector internally on every
//    call (they were being called T*N times).
//  - sample_indices_from_logw performs SYSTEMATIC resampling (Kitagawa,
//    1996; Doucet, de Freitas & Gordon, 2001) instead of multinomial
//    resampling for the K-1 non-reference particles.
//
// Author of the original R code: Cleiton Moya de Almeida
// C++/Rcpp port for maximum performance.
//
// [[Rcpp::plugins(cpp11)]]
#include "utils.h"
#include <utility>

// ---------------------------------------------------------------------
// TEST-ONLY export: solve a tridiagonal SPD system A x = b given diag d
// and off-diag e, for validation against a dense reference in R.
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List test_tridiag_solve(NumericVector d, NumericVector e, NumericVector b,
                          NumericVector w_for_Lt) {
    int n = d.size();
    std::vector<double> dv(d.begin(), d.end()), ev(e.begin(), e.end()), bv(b.begin(), b.end());
    TriLDLT fac = tri_ldlt_factor(dv, ev);
    std::vector<double> x(n);
    tri_solve_A(fac, bv, x);

    std::vector<double> wv(w_for_Lt.begin(), w_for_Lt.end());
    std::vector<double> xlt(n);
    tri_solve_Lt(fac, wv, xlt);

    return List::create(
        _["x_A"] = NumericVector(x.begin(), x.end()),
        _["D"] = NumericVector(fac.D.begin(), fac.D.end()),
        _["L"] = NumericVector(fac.L.begin(), fac.L.end()),
        _["x_Lt"] = NumericVector(xlt.begin(), xlt.end())
    );
}

// ---------------------------------------------------------------------
// TEST-ONLY export: Chan-method matrix/vector assembly, for validation
// against Matrix::bandSparse in R.
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List test_chan_assembly(NumericVector theta1, double phi1, double phi2,
                          double mu_02, double sigma2_02, double theta_01) {
    int Tt = theta1.size();
    int Ttp1 = Tt + 1;
    std::vector<double> d(Ttp1), e(Ttp1 - 1), b(Ttp1);

    d[0] = phi2 + (1.0 / sigma2_02 + phi1);
    for (int i = 1; i < Tt; ++i) d[i] = 2.0 * phi2 + phi1;
    d[Tt] = phi2;
    for (int i = 0; i < Ttp1 - 1; ++i) e[i] = -phi2;

    b[0] = mu_02 / sigma2_02 + phi1 * (theta1[0] - theta_01);
    for (int i = 1; i < Tt; ++i) b[i] = (theta1[i] - theta1[i - 1]) * phi1;
    b[Tt] = 0.0;

    return List::create(_["d"] = NumericVector(d.begin(), d.end()),
                         _["e"] = NumericVector(e.begin(), e.end()),
                         _["b"] = NumericVector(b.begin(), b.end()));
}

// ---------------------------------------------------------------------
// Main PG-AS sampler
// ---------------------------------------------------------------------
// [[Rcpp::export]]
List pg_as_cpp(NumericVector y,
                  int K, int N,
                  double mu_01, double sigma2_01,
                  double mu_02, double sigma2_02,
                  double nu_01, double eta_01,
                  double nu_02, double eta_02,
                  NumericVector theta1, NumericVector theta2,
                  double theta_01, double theta_02,
                  double W1, double W2,
                  bool verbose = true, int print_every = 100) {

    RNGScope rngScope; // sync with R's RNG / set.seed

    const int Tt = y.size();
    const int Ttp1 = Tt + 1;

    // ---- histories (row-major access via NumericMatrix, R-compatible) ----
    // Declared here (outer scope, not inside the nested block below) so
    // they remain accessible at the final return.
    NumericMatrix theta1_hist(N, Tt), theta2_hist(N, Tt);
    NumericVector theta_01_hist(N), theta_02_hist(N);
    NumericVector W1_hist(N), W2_hist(N);
    NumericVector ess_smc(N);

    // ---- state ----
    // theta_01, theta_02, W1, W2 are already local (by-value) copies of the
    // R-supplied initial values, so they are mutated in place below -- no
    // separate "_init" parameter/copy pair needed.
    // theta1/theta2 need a type conversion (NumericVector -> std::vector<
    // double>). A parameter's scope is the function's outer block, so a
    // same-named local there would be an illegal redeclaration (not mere
    // shadowing) -- the nested block below opens a genuinely new scope,
    // where reusing the name theta1/theta2 for the std::vector<double>
    // state is legal C++ and lets the rest of the function read exactly
    // like before.
    {
    std::vector<double> theta1_tmp(theta1.begin(), theta1.end());
    std::vector<double> theta1 = std::move(theta1_tmp);
    std::vector<double> theta2_tmp(theta2.begin(), theta2.end());
    std::vector<double> theta2 = std::move(theta2_tmp);

    // particle-system buffers (reused across iterations)
    // theta_1_k is flat, row-major (index t*K + k) to avoid the
    // pointer-chasing/cache-miss cost of std::vector<std::vector<double> >.
    std::vector<double> theta_1_k(Tt * (size_t)K, 0.0);
    // log_w_tilde only needs the previous and current time slice (see the
    // performance note above) -- kept as two K-sized buffers, swapped every
    // t, instead of a T x K matrix.
    std::vector<double> log_w_prev(K), log_w_curr(K);
    // A_hist (ancestor indices, for Backward Tracking) is flat, row-major
    // (index t*K + k), for the same cache-locality reason as theta_1_k.
    std::vector<int> A_hist(Tt * (size_t)K, 0);
    std::vector<double> log_as(K);
    std::vector<int> A(K - 1); // ancestors of the K - 1 non-reference particles
    std::vector<double> samp_w_buf(K), samp_cum_buf(K); // scratch for the resampling helpers

    // Chan-method buffers (size Ttp1)
    std::vector<double> d(Ttp1), e(Ttp1 - 1), b(Ttp1), theta_hat(Ttp1);
    std::vector<double> u(Ttp1), w(Ttp1), xsol(Ttp1), draw2(Ttp1);
    std::vector<double> z(Tt - 1);

    for (int n = 0; n < N; ++n) {

        if (verbose && ((n + 1) % print_every == 0)) {
            Rcout << "Iteration " << (n + 1) << " / " << N << "\n";
        }

        // ---- Sample W1 ----
        double phi1 = gibbs_sample_phi1(nu_01, eta_01, theta_01, theta1, theta_02, theta2, Tt);
        W1 = 1.0 / phi1;
        double sd_W1 = std::sqrt(W1);

        // ---- Sample W2 ----
        double phi2 = gibbs_sample_phi2(nu_02, eta_02, theta_02, theta2, Tt);
        W2 = 1.0 / phi2;

        // ---- Sample theta_01 (conjugate Normal) ----
        theta_01 = gibbs_sample_theta01(mu_01, sigma2_01, theta1[0], theta_02, W1);

        // ---- Conditional SMC (Bootstrap PF + Ancestral Sampling) for theta_t1 ----
        // t = 0 (R: t=1)
        for (int k = 0; k < K; ++k) {
            theta_1_k[k] = R::norm_rand() * sd_W1 + (theta_01 + theta_02);
        }
        theta_1_k[K - 1] = theta1[0];
        for (int k = 0; k < K; ++k) log_w_curr[k] = log_p_yt(y[0], theta_1_k[k]);
        {
            double lse = logsumexp(log_w_curr);
            for (int k = 0; k < K; ++k) log_w_curr[k] -= lse;
        }
        std::swap(log_w_prev, log_w_curr); // log_w_prev now holds the t = 0 weights

        // t = 1 .. Tt-1 (R: t = 2..Tt)
        for (int t = 1; t < Tt; ++t) {
            const size_t base_prev = (size_t)(t - 1) * K;
            const size_t base_curr = (size_t)t * K;

            // resampling + propagation for the K - 1 non-reference particles
            // (bootstrap PF: resample directly from log_w_prev)
            sample_indices_from_logw(log_w_prev.data(), K, samp_w_buf, samp_cum_buf, A);
            for (int k = 0; k < K - 1; ++k) {
                double mean_k = theta_1_k[base_prev + A[k]] + theta2[t - 1];
                theta_1_k[base_curr + k] = R::norm_rand() * sd_W1 + mean_k;
            }

            // reference particle: fixed value, ancestor drawn via Ancestral Sampling
            theta_1_k[base_curr + (K - 1)] = theta1[t];
            for (int j = 0; j < K; ++j) {
                double mean_j = theta_1_k[base_prev + j] + theta2[t - 1];
                log_as[j] = log_w_prev[j] + log_dnorm(theta1[t], mean_j, sd_W1);
            }
            int a_ref = sample_index_from_logw(log_as.data(), K, samp_w_buf);

            for (int k = 0; k < K - 1; ++k) A_hist[base_curr + k] = A[k];
            A_hist[base_curr + (K - 1)] = a_ref;

            // updated weights (bootstrap: likelihood only)
            for (int k = 0; k < K; ++k)
                log_w_curr[k] = log_p_yt(y[t], theta_1_k[base_curr + k]);
            double lse = logsumexp(log_w_curr);
            for (int k = 0; k < K; ++k) log_w_curr[k] -= lse;

            std::swap(log_w_prev, log_w_curr); // log_w_prev now holds the weights at time t
        }

        // ESS at final time
        // (after the loop above, log_w_prev holds the weights at t = Tt - 1)
        {
            double s2 = 0.0;
            for (int k = 0; k < K; ++k) s2 += std::exp(2.0 * log_w_prev[k]);
            ess_smc[n] = 1.0 / s2;
        }

        // ---- Backward tracking for theta1 ----
        // PG-AS: the genealogy already encodes the backward path, so we
        // only need to look up the stored ancestor indices -- no Backward
        // Simulation (no density evaluations here).
        int b_idx = sample_index_from_logw(log_w_prev.data(), K, samp_w_buf);
        theta1[Tt - 1] = theta_1_k[(size_t)(Tt - 1) * K + b_idx];
        for (int t = Tt - 2; t >= 0; --t) {
            b_idx = A_hist[(size_t)(t + 1) * K + b_idx];
            theta1[t] = theta_1_k[(size_t)t * K + b_idx];
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
    } // end of nested scope (theta1/theta2 as std::vector<double>)

    return List::create(
        _["theta1_hist"] = theta1_hist,
        _["theta2_hist"] = theta2_hist,
        _["theta_01_hist"] = theta_01_hist,
        _["theta_02_hist"] = theta_02_hist,
        _["W1_hist"] = W1_hist,
        _["W2_hist"] = W2_hist,
        _["ess_smc"] = ess_smc
    );
}
