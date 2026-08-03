// utils.h
// Common helper functions and structures shared across ALL samplers for the
// Poisson second-order polynomial DLM (local trend model) -- not just the
// Particle Gibbs family (PG-AS, PG-APF, ...), but also the other algorithms
// in the project (Metropolis-within-Gibbs / CWMH, SIR-Laplace, SIR-collapsed,
// ...). Mirrors utils.R: anything defined there is meant to eventually live
// here too, so the R and C++/Rcpp implementations stay in lockstep.
//
// Centralizing this code means any correction or optimization here
// propagates automatically to every algorithm that includes this header --
// no forking, no per-file duplication, no drift between implementations.
// This also anticipates the future R-package modularization: the plan is
// for the package's C++ backend to reuse this file as-is.
//
// Naming convention: function names mirror utils.R 1:1 (e.g.
// gibbs_sample_phi1 here <-> gibbs_sample_phi1 in R) to keep the R and C++
// versions easy to cross-reference.
//
// Contents so far:
//   - logsumexp, log_p_yt, log_dnorm
//   - gibbs_sample_theta01, gibbs_sample_phi1, gibbs_sample_phi2
//   - logpost_theta_t1, logpost_theta_T1, sample_theta_t1_mh
//     (component-wise Metropolis-within-Gibbs step, Montoril et al.)
//   - sample_index_from_logw       (single-index draw from log-weights)
//   - sample_indices_from_logw     (systematic resampling, Kitagawa 1996)
//   - TriLDLT + tri_ldlt_factor/tri_solve_Lt/tri_solve_A
//     (Chan method: tridiagonal LDL^T factorization/solve)
//
#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ---------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------

// log-sum-exp of a std::vector<double>
static inline double logsumexp(const std::vector<double>& x) {
    double m = x[0];
    for (std::size_t i = 1; i < x.size(); ++i) if (x[i] > m) m = x[i];
    double s = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) s += std::exp(x[i] - m);
    return m + std::log(s);
}

// Poisson log-likelihood contribution: y*theta - exp(theta)
static inline double log_p_yt(double yt, double theta) {
    double res = yt * theta - std::exp(theta);
    if (!std::isfinite(res)) res = -std::numeric_limits<double>::infinity();
    return res;
}

// log-density of N(x; mean, sd)
static inline double log_dnorm(double x, double mean, double sd) {
    static const double LOG_SQRT_2PI = 0.9189385332046727;
    double z = (x - mean) / sd;
    return -LOG_SQRT_2PI - std::log(sd) - 0.5 * z * z;
}

// ---------------------------------------------------------------------
// Conjugate Gibbs steps for theta_01, phi1 = 1/W1, phi2 = 1/W2
// (mirrors gibbs_sample_theta01 / gibbs_sample_phi1 / gibbs_sample_phi2
// in utils.R -- same signatures, same return values: phi1/phi2 are
// PRECISIONS, not variances, exactly as in the R version).
// ---------------------------------------------------------------------

// theta_01 | ... ~ N(mu_01_bar, sigma2_01_bar)
static inline double gibbs_sample_theta01(double mu_01, double sigma2_01,
                                           double theta_11, double theta_02,
                                           double W1) {
    double sigma2_01_bar = 1.0 / (1.0 / sigma2_01 + 1.0 / W1);
    double mu_01_bar = sigma2_01_bar * (mu_01 / sigma2_01 + (theta_11 - theta_02) / W1);
    return R::norm_rand() * std::sqrt(sigma2_01_bar) + mu_01_bar;
}

// phi1 = 1/W1 | ... ~ Gamma(nu_01_bar, eta_01_bar)
// theta1, theta2: 0-based vectors of length Tt (R indices 1..Tt)
static inline double gibbs_sample_phi1(double nu_01, double eta_01, double theta_01,
                                        const std::vector<double>& theta1,
                                        double theta_02,
                                        const std::vector<double>& theta2,
                                        int Tt) {
    double ssq = 0.0;
    double prev1 = theta_01, prev2 = theta_02;
    for (int t = 0; t < Tt; ++t) {
        double dif1 = theta1[t] - prev1;
        double diffs = dif1 - prev2;
        ssq += diffs * diffs;
        prev1 = theta1[t];
        prev2 = theta2[t];
    }
    double nu_01_bar = nu_01 + Tt / 2.0;
    double eta_01_bar = eta_01 + 0.5 * ssq;
    return R::rgamma(nu_01_bar, 1.0 / eta_01_bar); // rate -> scale = 1/rate
}

// phi2 = 1/W2 | ... ~ Gamma(nu_02_bar, eta_02_bar)
// theta2: 0-based vector of length Tt (R indices 1..Tt)
static inline double gibbs_sample_phi2(double nu_02, double eta_02, double theta_02,
                                        const std::vector<double>& theta2,
                                        int Tt) {
    double ssq = 0.0;
    double prev = theta_02;
    for (int t = 0; t < Tt; ++t) {
        double diffs = theta2[t] - prev;
        ssq += diffs * diffs;
        prev = theta2[t];
    }
    double nu_02_bar = nu_02 + Tt / 2.0;
    double eta_02_bar = eta_02 + 0.5 * ssq;
    return R::rgamma(nu_02_bar, 1.0 / eta_02_bar); // rate -> scale = 1/rate
}

// ---------------------------------------------------------------------
// Component-wise Metropolis-within-Gibbs step for theta_t1 (Montoril et
// al.): random-walk proposal, full conditional log-posterior evaluated
// in closed form (mirrors logpost_theta_t1 / logpost_theta_T1 /
// sample_theta_t1_mh in utils.R).
// ---------------------------------------------------------------------

// Full conditional log-posterior for theta_t1, t = 1, ..., T-1 (0-based
// t = 0, ..., Tt-2): has both a predecessor and a successor time point.
static inline double logpost_theta_t1(double theta_t1, double theta_tm11, double theta_tp11,
                                       double theta_t2, double theta_tm12,
                                       double yt, double W1) {
    double sigma2_star = W1 / 2.0;
    double mu_star = ((theta_tm11 + theta_tm12) + (theta_tp11 - theta_t2)) / 2.0;
    double p1 = log_p_yt(yt, theta_t1);
    double diff = theta_t1 - mu_star;
    double p2 = -(diff * diff) / (2.0 * sigma2_star);
    return p1 + p2;
}

// Full conditional log-posterior for theta_T1 (t = T, the last time point:
// no successor term).
static inline double logpost_theta_T1(double theta_t1, double theta_tm11, double theta_tm12,
                                       double yt, double W1) {
    double p1 = log_p_yt(yt, theta_t1);
    double diff = theta_t1 - theta_tm11 - theta_tm12;
    double p2 = -(diff * diff) / (2.0 * W1);
    return p1 + p2;
}

// Result of a single Metropolis step for theta_t1: the (possibly
// unchanged) value and whether the proposal was accepted (mirrors R's
// sample_theta_t1_mh, which returns list(theta_t1, ac)).
struct MHStep {
    double theta_t1;
    int ac;
};

// Metropolis step for theta_t1 (random-walk proposal, variance varsigma2).
// final_t = true selects logpost_theta_T1 (t = T; theta_tp11 is then
// unused -- pass any value, e.g. 0.0).
static inline MHStep sample_theta_t1_mh(double theta_t1_current, double theta_tm11, double theta_tp11,
                                         double theta_t2, double theta_tm12,
                                         double yt, double W1, double varsigma2, bool final_t) {
    double theta_t1_prop = R::norm_rand() * std::sqrt(varsigma2) + theta_t1_current;

    double logu = std::log(unif_rand());
    double logp1, logp2;
    if (final_t) {
        logp1 = logpost_theta_T1(theta_t1_prop, theta_tm11, theta_tm12, yt, W1);
        logp2 = logpost_theta_T1(theta_t1_current, theta_tm11, theta_tm12, yt, W1);
    } else {
        logp1 = logpost_theta_t1(theta_t1_prop, theta_tm11, theta_tp11, theta_t2, theta_tm12, yt, W1);
        logp2 = logpost_theta_t1(theta_t1_current, theta_tm11, theta_tp11, theta_t2, theta_tm12, yt, W1);
    }
    double logr = logp1 - logp2;

    MHStep out;
    if (logu < logr) {
        out.theta_t1 = theta_t1_prop;
        out.ac = 1;
    } else {
        out.theta_t1 = theta_t1_current;
        out.ac = 0;
    }
    return out;
}

// ---------------------------------------------------------------------
// Draw one index in [0, K) from (unnormalized) log-weights, using a single
// uniform draw and the cumulative distribution (inverse-CDF sampling).
// logw is a raw pointer, so it can point directly into a scratch weight
// buffer without a copy. w_buf is a caller-provided scratch buffer of size
// >= K, reused across calls to avoid a heap allocation on every invocation.
// ---------------------------------------------------------------------
static inline int sample_index_from_logw(const double* logw, int K,
                                          std::vector<double>& w_buf) {
    double m = logw[0];
    for (int k = 1; k < K; ++k) if (logw[k] > m) m = logw[k];
    double tot = 0.0;
    for (int k = 0; k < K; ++k) { w_buf[k] = std::exp(logw[k] - m); tot += w_buf[k]; }
    double u = unif_rand() * tot;
    double cum = 0.0;
    for (int k = 0; k < K; ++k) {
        cum += w_buf[k];
        if (u <= cum) return k;
    }
    return K - 1;
}

// Draw out.size() indices (with replacement) in [0, K) via SYSTEMATIC
// resampling (Kitagawa, 1996): a single unif_rand() draw u0 ~ U(0, total/M)
// (M = out.size()) determines M equally-spaced points u_j = u0 + total*j/M,
// located by ONE forward sweep through the cumulative weights (u_j is
// increasing in j, so the search pointer never needs to go backward). This
// replaces multinomial resampling (M independent unif_rand() draws + M
// binary searches) with 1 draw + a single O(K) linear scan. Systematic
// resampling is unbiased and has lower variance than multinomial resampling
// (Doucet, de Freitas & Gordon, 2001), so this is a standard SMC choice,
// not merely a speed shortcut.
static inline void sample_indices_from_logw(const double* logw, int K,
                                             std::vector<double>& w_buf,
                                             std::vector<double>& cum_buf,
                                             std::vector<int>& out) {
    double m = logw[0];
    for (int k = 1; k < K; ++k) if (logw[k] > m) m = logw[k];
    double tot = 0.0;
    for (int k = 0; k < K; ++k) { w_buf[k] = std::exp(logw[k] - m); tot += w_buf[k]; }
    cum_buf[0] = w_buf[0];
    for (int k = 1; k < K; ++k) cum_buf[k] = cum_buf[k - 1] + w_buf[k];
    double total = cum_buf[K - 1];
    int ndraws = (int) out.size();
    double u0 = unif_rand() * (total / ndraws);
    int k = 0;
    for (int j = 0; j < ndraws; ++j) {
        double u = u0 + total * ((double) j / ndraws);
        while (cum_buf[k] < u) ++k;
        out[j] = k;
    }
}

// ---------------------------------------------------------------------
// Tridiagonal LDL^T factorization/solve (Chan method), 0-based, size n.
// A is symmetric tridiagonal: diagonal d[0..n-1], sub/super-diagonal e[0..n-2]
// (e[i] connects rows i and i+1). Produces L (unit lower-bidiagonal,
// L[i] = multiplier at (i, i-1), i = 1..n-1) and D (diagonal, D[0..n-1]).
// ---------------------------------------------------------------------
struct TriLDLT {
    std::vector<double> L; // size n, L[0] unused
    std::vector<double> D; // size n
};

static inline TriLDLT tri_ldlt_factor(const std::vector<double>& d,
                                       const std::vector<double>& e) {
    int n = (int) d.size();
    TriLDLT f;
    f.L.assign(n, 0.0);
    f.D.assign(n, 0.0);
    f.D[0] = d[0];
    for (int i = 1; i < n; ++i) {
        f.L[i] = e[i - 1] / f.D[i - 1];
        f.D[i] = d[i] - f.L[i] * f.L[i] * f.D[i - 1];
    }
    return f;
}

// Solve L^T x = w  (backward substitution; L upper-bidiagonal transpose)
static inline void tri_solve_Lt(const TriLDLT& f, const std::vector<double>& w,
                                 std::vector<double>& x) {
    int n = (int) f.D.size();
    x[n - 1] = w[n - 1];
    for (int i = n - 2; i >= 0; --i) x[i] = w[i] - f.L[i + 1] * x[i + 1];
}

// Solve A x = b  (L D L^T x = b): forward, diag, backward(=tri_solve_Lt)
static inline void tri_solve_A(const TriLDLT& f, const std::vector<double>& b,
                                std::vector<double>& x) {
    int n = (int) f.D.size();
    std::vector<double> z(n), w(n);
    z[0] = b[0];
    for (int i = 1; i < n; ++i) z[i] = b[i] - f.L[i] * z[i - 1];
    for (int i = 0; i < n; ++i) w[i] = z[i] / f.D[i];
    tri_solve_Lt(f, w, x);
}

#endif // UTILS_H
