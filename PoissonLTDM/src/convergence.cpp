// PoissonLTDM/src/convergence.cpp
//
// Reimplementação em C++ (Rcpp + RcppArmadillo + FFTW3) de
// posterior::rhat(), posterior::ess_bulk(), posterior::ess_tail(),
// fiel ao algoritmo de Vehtari et al. (2021), conforme extraído do
// código-fonte oficial (stan-dev/posterior, R/convergence.R e
// R/split_chains.R). Vetorizado sobre Tt "variáveis" (colunas na terceira
// dimensão de um draws_array) num loop C++ paralelizado com OpenMP,
// reutilizando o plano FFTW3 por thread (mesmo tamanho de FFT em todas as
// chamadas, já que n_post é fixo).
//
// Exportado como PoissonLTDM::rhat_ess_fast() após build/load_all() do
// pacote. Usado por PoissonLTDM::metrics_convergence() e
// PoissonLTDM::metrics_convergence_by_chain() (R/metrics.R), que montam o
// array (n_post, N_chains, Tt) a partir de result_list e chamam esta
// função -- não é chamado diretamente pelo código de diagnóstico.
//
// Input: array [[double]] com dimensões (n_post, N_chains, Tt) -- mesma
// convenção de posterior::as_draws_array().
// Output: lista com 3 vetores de tamanho Tt: rhat, ess_bulk, ess_tail.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <fftw3.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// ---- qnorm thread-safe (AS 241, Wichura 1988) ----
// R::qnorm() (Rf_qnorm5 via a API do R) NAO e thread-safe -- nao pode ser
// chamado de dentro de threads OpenMP. Esta e uma porta direta do
// qnorm5() do R (src/nmath/qnorm.c, Ihaka/R-core, baseado em Wichura AS 241),
// especializada para o unico caso usado aqui: mu=0, sigma=1, lower_tail=TRUE,
// log_p=FALSE. E uma funcao pura (sem estado global), portanto thread-safe.
// Precisao: ~1 parte em 10^16, identica ao qnorm() do R para este caso.
static inline double qnorm_std(double p) {
    double q = p - 0.5;
    double r, val;
    if (std::fabs(q) <= .425) {
        r = .180625 - q * q;
        val = q * (((((((r * 2509.0809287301226727 +
                   33430.575583588128105) * r + 67265.770927008700853) * r +
                 45921.953931549871457) * r + 13731.693765509461125) * r +
               1971.5909503065514427) * r + 133.14166789178437745) * r +
             3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                 28729.085735721942674) * r + 39307.89580009271061) * r +
               21213.794301586595867) * r + 5394.1960214247511077) * r +
             687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    } else {
        double lp = std::log((q > 0) ? (1.0 - p) : p);
        r = std::sqrt(-lp);
        if (r <= 5.) {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        } else if (r <= 27) {
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        } else {
            if (r >= 6.4e8) {
                val = r * M_SQRT2;
            } else {
                double s2 = -2.0 * lp;
                double x2 = s2 - std::log(2.0 * M_PI * s2);
                if (r < 36000.) {
                    x2 = s2 - std::log(2.0 * M_PI * x2) - 2. / (2. + x2);
                    if (r < 840.) {
                        x2 = s2 - std::log(2.0 * M_PI * x2) +
                             2 * std::log1p(-(1 - 1 / (4 + x2)) / (2. + x2));
                        if (r < 109.) {
                            x2 = s2 - std::log(2.0 * M_PI * x2) +
                                 2 * std::log1p(-(1 - (1 - 5 / (6 + x2)) / (4. + x2)) / (2. + x2));
                            if (r < 55.) {
                                x2 = s2 - std::log(2.0 * M_PI * x2) +
                                     2 * std::log1p(-(1 - (1 - (5 - 9 / (8. + x2)) / (6. + x2)) / (4. + x2)) / (2. + x2));
                            }
                        }
                    }
                }
                val = std::sqrt(x2);
            }
        }
        if (q < 0.0) val = -val;
    }
    return val;
}

// ---- rank() com ties.method="average", igual ao rank() do R ----
// Retorna ranks médios para empates, 1-indexado (como o R).
static arma::vec rank_average(const arma::vec& x) {
    int n = x.n_elem;
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return x[a] < x[b]; });

    arma::vec ranks(n);
    int i = 0;
    while (i < n) {
        int j = i;
        while (j + 1 < n && x[idx[j + 1]] == x[idx[i]]) ++j;
        // empate entre posições i..j (ordenadas), rank médio = média de (i+1)..(j+1)
        double avg_rank = (double)(i + 1 + j + 1) / 2.0;
        for (int k = i; k <= j; ++k) ranks[idx[k]] = avg_rank;
        i = j + 1;
    }
    return ranks;
}

// ---- backtransform_ranks: (r - c) / (S - 2c + 1), c = 3/8 ----
static arma::vec backtransform_ranks(const arma::vec& r, double c = 3.0 / 8.0) {
    double S = r.n_elem;
    return (r - c) / (S - 2.0 * c + 1.0);
}

// ---- z_scale: rank-normalização (rank -> backtransform -> qnorm) ----
// Opera sobre a matriz achatada inteira (equivalente a as.array(x) no R:
// rank() é calculado sobre TODOS os elementos da matriz junta, não por coluna)
static arma::mat z_scale(const arma::mat& x) {
    int n = x.n_rows, m = x.n_cols;
    arma::vec flat = arma::vectorise(x); // cópia segura (x é const)
    arma::vec r = rank_average(flat);
    arma::vec u = backtransform_ranks(r);
    arma::vec z(n * m);
    for (int i = 0; i < n * m; ++i) z[i] = qnorm_std(u[i]);
    arma::mat out(z.memptr(), n, m); // reconstrói como matriz (cópia segura)
    return out;
}

// ---- split_chains: divide cada cadeia (coluna) ao meio, dobra colunas ----
// Réplica fiel de posterior:::.split_chains: cbind(x[1:floor(half),],
// x[ceiling(half+1):niter,]), com half = niter/2 (double, não inteiro).
// Quando niter é ímpar, floor(half) e ceiling(half+1) produzem duas metades
// de MESMO tamanho, descartando automaticamente a linha do meio (sem erro,
// sem warning) -- confirmado empiricamente contra o R para niter=201.
static arma::mat split_chains(const arma::mat& x) {
    int niter = x.n_rows;
    if (niter == 1) return x;
    double half = niter / 2.0;
    int first_end = (int)std::floor(half);           // floor(half), 1-indexed length
    int second_start = (int)std::ceil(half + 1.0);    // ceiling(half+1), 1-indexed start
    // 0-indexed: primeira metade linhas [0, first_end), segunda [second_start-1, niter)
    arma::mat first = x.rows(0, first_end - 1);
    arma::mat second = x.rows(second_start - 1, niter - 1);
    return arma::join_rows(first, second);
}

// ---- fold_draws: |x - median(x)| ----
static arma::mat fold_draws(const arma::mat& x) {
    arma::vec flat = arma::vectorise(x);
    double med = arma::median(flat);
    return arma::abs(x - med);
}

// ---- .rhat: razão de variância entre/dentro de cadeias ----
static double rhat_core(const arma::mat& x) {
    int niterations = x.n_rows;
    int nchains = x.n_cols;
    arma::rowvec chain_mean = arma::mean(x, 0);
    arma::rowvec chain_var(nchains);
    for (int c = 0; c < nchains; ++c) {
        arma::vec col = x.col(c);
        chain_var[c] = arma::var(col);
    }
    double var_between = niterations * arma::var(chain_mean.t());
    double var_within = arma::mean(chain_var);
    return std::sqrt((var_between / var_within + niterations - 1) / (double)niterations);
}

// ---- Contexto FFTW3 reutilizável entre chamadas (mesmo N sempre) ----
struct FFTWContext {
    int N;      // tamanho original da série (niterations pós-split)
    int Mt2;    // tamanho com zero-padding (2 * nextn(N))
    fftw_plan plan_fwd;
    fftw_plan plan_inv;
    double* in;
    fftw_complex* freq;
    double* out;
    int freq_len;

    FFTWContext(int N_) : N(N_) {
        // nextn(N) no R usa fatores 2,3,5 por padrão; para simplicidade e
        // robustez usamos a próxima potência de 2 >= N, que FFTW lida melhor
        // (mesmo princípio de "zero padding acelera FFT", resultado
        // matematicamente idêntico já que ambos são zero-padding válido --
        // a normalização por ac[0] no autocovariance() torna o tamanho exato
        // do padding irrelevante para o resultado final).
        int M = 1;
        while (M < N) M <<= 1;
        Mt2 = 2 * M;
        freq_len = Mt2 / 2 + 1;

        in = (double*) fftw_malloc(sizeof(double) * Mt2);
        freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * freq_len);
        out = (double*) fftw_malloc(sizeof(double) * Mt2);

        plan_fwd = fftw_plan_dft_r2c_1d(Mt2, in, freq, FFTW_MEASURE);
        plan_inv = fftw_plan_dft_c2r_1d(Mt2, freq, out, FFTW_MEASURE);
    }

    ~FFTWContext() {
        fftw_destroy_plan(plan_fwd);
        fftw_destroy_plan(plan_inv);
        fftw_free(in);
        fftw_free(freq);
        fftw_free(out);
    }

    // autocovariance(x): réplica fiel de posterior:::autocovariance()
    // usando FFT real-to-complex (equivalente a Re(fft(abs(fft(yc))^2, inverse=TRUE)))
    arma::vec autocovariance(const arma::vec& x) {
        int n = x.n_elem; // == N
        double varx = arma::var(x);
        arma::vec ac(n);
        if (varx == 0.0) {
            ac.zeros();
            return ac;
        }
        double mean_x = arma::mean(x);
        std::fill(in, in + Mt2, 0.0);
        for (int i = 0; i < n; ++i) in[i] = x[i] - mean_x;

        fftw_execute(plan_fwd);
        // |fft(yc)|^2 no domínio de frequência
        for (int k = 0; k < freq_len; ++k) {
            double re = freq[k][0], im = freq[k][1];
            double mag2 = re * re + im * im;
            freq[k][0] = mag2;
            freq[k][1] = 0.0;
        }
        fftw_execute(plan_inv);
        // FFTW não normaliza a inversa (resultado vem multiplicado por Mt2)
        double scale = 1.0 / (double)Mt2;
        arma::vec ac_raw(n);
        for (int i = 0; i < n; ++i) ac_raw[i] = out[i] * scale;

        // normalização igual ao posterior: ac / ac[0] * varx * (N-1)/N
        double ac0 = ac_raw[0];
        for (int i = 0; i < n; ++i) {
            ac[i] = ac_raw[i] / ac0 * varx * (n - 1) / (double)n;
        }
        return ac;
    }
};

// ---- .ess: effective sample size via autocorrelação + soma de Geyer ----
static double ess_core(const arma::mat& x, FFTWContext& ctx) {
    int niterations = x.n_rows;
    int nchains = x.n_cols;
    if (niterations < 3) return NA_REAL;

    arma::mat acov(niterations, nchains);
    for (int c = 0; c < nchains; ++c) {
        arma::vec col = x.col(c);
        acov.col(c) = ctx.autocovariance(col);
    }
    arma::vec acov_means = arma::mean(acov, 1); // rowMeans

    double mean_var = acov_means[0] * niterations / (double)(niterations - 1);
    double var_plus = mean_var * (niterations - 1) / (double)niterations;
    if (nchains > 1) {
        arma::rowvec chain_means = arma::mean(x, 0);
        var_plus += arma::var(chain_means.t());
    }

    arma::vec rho_hat_t(niterations, arma::fill::zeros);
    int t = 0;
    double rho_hat_even = 1.0;
    rho_hat_t[t] = rho_hat_even; // rho_hat_t[t+1] no R (1-indexed) -> rho_hat_t[t] aqui (0-indexed)
    double rho_hat_odd = 1.0 - (mean_var - acov_means[t + 1]) / var_plus;
    rho_hat_t[t + 1] = rho_hat_odd;

    while (t < (int)acov.n_rows - 5 && !std::isnan(rho_hat_even + rho_hat_odd) &&
           (rho_hat_even + rho_hat_odd > 0)) {
        t += 2;
        rho_hat_even = 1.0 - (mean_var - acov_means[t]) / var_plus;
        rho_hat_odd = 1.0 - (mean_var - acov_means[t + 1]) / var_plus;
        if ((rho_hat_even + rho_hat_odd) >= 0) {
            rho_hat_t[t] = rho_hat_even;
            rho_hat_t[t + 1] = rho_hat_odd;
        }
    }
    int max_t = t;
    if (rho_hat_even > 0) rho_hat_t[max_t] = rho_hat_even; // rho_hat_t[max_t+1] no R 1-indexed

    t = 0;
    while (t <= max_t - 4) {
        t += 2;
        // R (1-indexed): rho_hat_t[t+1]+rho_hat_t[t+2] vs rho_hat_t[t-1]+rho_hat_t[t]
        // 0-indexed aqui: rho_hat_t[t]+rho_hat_t[t+1] vs rho_hat_t[t-2]+rho_hat_t[t-1]
        if (rho_hat_t[t] + rho_hat_t[t + 1] > rho_hat_t[t - 2] + rho_hat_t[t - 1]) {
            rho_hat_t[t] = (rho_hat_t[t - 2] + rho_hat_t[t - 1]) / 2.0;
            rho_hat_t[t + 1] = rho_hat_t[t];
        }
    }

    double ess = nchains * niterations;
    // sum(rho_hat_t[1:max_t]) no R (1-indexed, inclusive) -> índices 0..max_t-1 aqui
    double sum_rho = 0.0;
    for (int i = 0; i < max_t; ++i) sum_rho += rho_hat_t[i];
    double tau_hat = -1.0 + 2.0 * sum_rho + rho_hat_t[max_t]; // rho_hat_t[max_t+1] no R
    double tau_bound = 1.0 / std::log10(ess);
    if (tau_hat < tau_bound) tau_hat = tau_bound;
    ess = ess / tau_hat;
    return ess;
}

// ---- ess_quantile: ESS de indicadora (x <= quantile(x, prob)) ----
static double ess_quantile_core(const arma::mat& x, double prob, FFTWContext& ctx) {
    arma::vec flat = arma::vectorise(x);
    arma::vec sorted = arma::sort(flat);
    int n = sorted.n_elem;
    // quantile() tipo 7 do R (default): mesmo método usado por quantile(x, prob)
    double h = (n - 1) * prob;
    int lo = (int)std::floor(h);
    int hi = (int)std::ceil(h);
    double qval = sorted[lo] + (h - lo) * (sorted[hi] - sorted[lo]);

    arma::mat I(x.n_rows, x.n_cols);
    for (arma::uword i = 0; i < x.n_elem; ++i) I[i] = (x[i] <= qval) ? 1.0 : 0.0;
    arma::mat I_split = split_chains(I);
    return ess_core(I_split, ctx);
}

// [[Rcpp::export]]
List rhat_ess_fast(NumericVector arr_r, int n_threads = 0) {
    IntegerVector dims = arr_r.attr("dim");
    if (dims.size() != 3) stop("arr must be a 3D array (n_post x N_chains x Tt)");
    int n_post = dims[0], N_chains = dims[1], Tt = dims[2];

    arma::cube arr(arr_r.begin(), n_post, N_chains, Tt, false, true);

    NumericVector out_rhat(Tt), out_ess_bulk(Tt), out_ess_tail(Tt);

    int n_split = n_post / 2;

#ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    int nt = (n_threads > 0) ? std::min(n_threads, max_threads) : max_threads;
#else
    int nt = 1;
    if (n_threads > 1) {
        Rcpp::warning("Compilado sem suporte a OpenMP -- rodando single-thread. "
                       "Verifique se -fopenmp foi usado na compilacao.");
    }
#endif

    // Um FFTWContext por thread, criados SEQUENCIALMENTE (fora da regiao
    // paralela): fftw_plan_dft_*_1d nao e thread-safe durante o planning
    // (estado global interno da FFTW), mesmo que a EXECUCAO do plano via
    // fftw_execute() seja thread-safe quando cada thread usa seus proprios
    // buffers de entrada/saida. Criar os planos antes evita qualquer corrida
    // nessa etapa.
    std::vector<std::unique_ptr<FFTWContext>> contexts;
    contexts.reserve(nt);
    for (int i = 0; i < nt; ++i) {
        contexts.push_back(std::make_unique<FFTWContext>(n_split));
    }

#ifdef _OPENMP
    #pragma omp parallel for num_threads(nt) schedule(dynamic)
#endif
    for (int t = 0; t < Tt; ++t) {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        FFTWContext& ctx = *contexts[tid];

        arma::mat x = arr.slice(t); // n_post x N_chains

        // ---- rhat ----
        arma::mat x_split = split_chains(x);
        arma::mat x_split_z = z_scale(x_split);
        double rhat_bulk = rhat_core(x_split_z);

        arma::mat x_folded = fold_draws(x);
        arma::mat x_folded_split = split_chains(x_folded);
        arma::mat x_folded_split_z = z_scale(x_folded_split);
        double rhat_tail = rhat_core(x_folded_split_z);

        out_rhat[t] = std::max(rhat_bulk, rhat_tail);

        // ---- ess_bulk (reusa x_split_z já calculado acima) ----
        out_ess_bulk[t] = ess_core(x_split_z, ctx);

        // ---- ess_tail ----
        double ess_q05 = ess_quantile_core(x, 0.05, ctx);
        double ess_q95 = ess_quantile_core(x, 0.95, ctx);
        out_ess_tail[t] = std::min(ess_q05, ess_q95);
    }

    return List::create(
        Named("rhat") = out_rhat,
        Named("ess_bulk") = out_ess_bulk,
        Named("ess_tail") = out_ess_tail
    );
}
