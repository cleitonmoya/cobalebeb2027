# Poisson Local Trend Dynamic Model

Code, data, and instructions to reproduce the simulation study and real-data
application of **"Efficient Samplers for the Poisson Local Trend Dynamic
Model"**, submitted to the VIII Latin American Meeting on Bayesian Statistics
(VIII COBAL) / XVIII Brazilian Meeting of Bayesian Statistics (EBEB), 2027.

*Author information is withheld while this submission is under double-blind
review.*

## Overview

The **Poisson Local Trend Dynamic Model (PoissonLTDM)** is a second-order
polynomial Dynamic Generalized Linear Model with Poisson observations. The
state vector $(\theta_{t1}, \theta_{t2})$ has $\theta_{t1}$ as the level
(log-intensity, entering the observation equation via $\exp(\theta_{t1})$)
and $\theta_{t2}$ as the slope, following the Durbin & Koopman / Harvey
structural-time-series nomenclature.

The paper compares five Bayesian samplers for this model:

| Method            | Description                                                        |
|-------------------|---------------------------------------------------------------------|
| `amh_montoril`    | Adaptive Metropolis-Hastings (componentwise, Robbins–Monro tuning)  |
| `pg_apf` / `pg_as`| Particle Gibbs with Auxiliary Particle Filter / Ancestor Sampling   |
| `sir_laplace`     | Sampling-Importance-Resampling with Laplace/IRLS approximation      |
| `sir_collapsed`   | **Main methodological contribution**: collapsed Gibbs with a Cross-Entropy-calibrated Gamma MH proposal |
| `stan`            | HMC/NUTS via Stan — reference / gold-standard sampler                |

Two applications are reported:

1. A **simulation study** across a grid of series lengths and function types
   (piecewise constant, piecewise linear, piecewise quadratic, sinusoidal),
   comparing sampler efficiency and estimation accuracy.
2. A **real-data application** to the `campy` dataset (weekly campylobacteriosis
   counts, Quebec, 1990–2000; from the `tscount` package), benchmarked against
   observation-driven models fit with `tscount::tsglm()`.

## Repository structure

```
cobalebeb2027/
├── PoissonLTDM/                  # R package (Rcpp) — the samplers
│   ├── R/                        # metrics.R, sampler_stan.R, RcppExports.R
│   ├── src/                      # utils.h + one .cpp per sampler + RcppExports.cpp
│   └── inst/stan/                # poisson_ltdm.stan
│
├── R_prototypes/                 # Didactic, pure-R reference implementations
│                                  # (pg_as_r(), amh_montoril_r(), sir_laplace_r(),
│                                  #  sir_collapsed_r() — signatures mirror the *_cpp() ones)
│
├── tests/                        # Correctness / validation harness
│   ├── test_prototype_R.R
│   ├── test_cpp.R
│   ├── test_validation.R
│   └── plot_diagnostics.R
│
├── data/
│   ├── data_generation.R         # Generates the synthetic series used by both
│   │                              # calibration and the production simulation
│   └── simulated/                # <function>_<Tt>_<replica>.rds files
│
├── calibration/                  # Tuning phase (N / burnin / K per method)
│   ├── calibration_run.R
│   ├── calibration_aggregate.R
│   ├── calibration_pbs.tmpl
│   ├── check_calibration_progress.R
│   └── registry_calibration/     # batchtools registry (cluster runs)
│
├── simulation/                   # Production simulation study
│   ├── simulation_run.R
│   ├── simulation_grid_config.R  # run_mode / grid_subset — edit without touching simulation_run.R
│   ├── simulation_aggregate.R
│   ├── simulation_pbs.tmpl
│   ├── check_simulation_progress.R
│   └── registry_simulation/      # batchtools registry (cluster runs)
│
├── results/
│   ├── calibration/
│   └── partial/                  # per-task .rds results from the production run
│
├── real_data_run.R               # Fits all five samplers to `campy`
├── real_data_aggregate.R
├── figures_real_data.R           # article_real_data_fit.pdf, article_real_data_ess.pdf
├── predictive_loglik_pf.R        # Predictive log-likelihood via particle filter
├── tscount_literature_benchmark.R# INGARCH-family benchmarks via tscount::tsglm()
│
├── cache/                        # Compiled Stan model cache (not versioned)
└── latex/                        # Manuscript source
```

## The `PoissonLTDM` package

`PoissonLTDM` is the R/Rcpp package containing the production implementation
of all five samplers. Performance-critical inner loops (`amh_montoril`,
`pg_apf`/`pg_as`, `sir_laplace`, `sir_collapsed`) are implemented in C++
(`src/`), sharing common routines through `utils.h`. `sampler_stan.R` wraps
the Stan/`rstan` model in `inst/stan/poisson_ltdm.stan`.

Install locally with:

```r
# from the repository root
pkgload::load_all("PoissonLTDM")   # for development / interactive use
# or
devtools::install("PoissonLTDM")    # for a regular package install
```

<!-- TODO: list hard package dependencies (Rcpp, RcppArmadillo/Eigen if used,
     rstan, posterior, Matrix, ...) once DESCRIPTION is finalized -->

### `R_prototypes/`

Standalone, pure-R implementations of the same four non-Stan samplers
(`pg_as_r()`, `amh_montoril_r()`, `sir_laplace_r()`, `sir_collapsed_r()`),
kept outside the package. Their purpose is **didactic**: each function's
signature and return structure mirror the corresponding `PoissonLTDM::*_cpp()`
routine one-to-one, so a reader can follow the algorithm in plain R before
looking at the optimized C++ code. They also serve as the reference
implementation against which the C++ port is validated (see below).

## Tests and validation

The `tests/` directory is the correctness harness, not an automated
`testthat` suite bundled with the package — it is meant to be run manually.

- **`test_prototype_R.R`** / **`test_cpp.R`** — run a single sampler
  (selected via an `algorithm <- "..."` variable at the top of the file) on
  one dataset, with all hyperparameters and initial values exposed inline in
  the script (no shared config layer). Diagnostics and plots are produced by
  the shared `plot_diagnostics.R` helper (`print_and_plot_diagnostics()`),
  which adapts automatically to the fields returned by each algorithm
  (SMC-ESS, MH acceptance rate, IS-ESS, CE-calibration diagnostics, etc.).
- **`test_validation.R`** — runs the R prototype and the C++ port with the
  same seed and configuration and compares every output field with
  `testthat::expect_equal(tolerance = 1e-10)`. This is the check that the
  Rcpp port is numerically faithful to the didactic R implementation.

Note on exact reproducibility across implementations: bit-exact agreement
between the R and C++ RNG paths (`sample()`/`sample.int()` in R vs. the
`sample_one_from_logw()` routine in `utils.h`) is validated at moderate
$T \times N$; at very large $T \times K \times N$ combinations, floating-point
boundary sensitivity in systematic resampling causes small, expected
divergence — this is documented, not treated as a bug.

To run a validation check:

```r
setwd("tests")
source("test_validation.R")   # edit `algorithm` at the top to choose which sampler
```

## Data generation

`data/data_generation.R` generates the synthetic series shared by both the
calibration phase and the production simulation, saved under
`data/simulated/` following the naming convention
`<function>_<Tt>_<replica>.rds`, where `<function>` is one of the four
trend types, `<Tt>` the series length, and `<replica>` the replica index.

## Calibration

`calibration/` tunes the number of MCMC/SMC iterations (`N`), burn-in, and
(for `pg_apf`/`pg_as`) the number of particles `K` for each sampler, over a
grid of configurations, selecting the setting with the best efficiency/ESS
trade-off. The values used in the production study were:

| Method                                | N       | burnin | K   |
|----------------------------------------|---------|--------|-----|
| `amh_montoril`                         | 110,000 | 10,000 | —   |
| `pg_apf`/`pg_as`                       | 22,000  | 2,000  | 200 |
| `sir_laplace`, `sir_collapsed`, `stan` | 11,000  | 1,000  | —   |

These are stored per-method in `R_config` inside `simulation/simulation_run.R`.

## Production simulation

The full simulation study spans $T \in \{200, 400, 800, 2000\}$, four
function types, and $R = 200$ replicas for **all five methods**, including
`stan` (the `stan` replica count was raised mid-study from an initial
$R=50$ reference subsample to the full $R=200$, once compute budget allowed
it). It was run on an external HPC cluster (PBS Pro scheduler,
`batchtools` job arrays), with each job handling a chunk of tasks
(`simulation_run.R` + `simulation_pbs.tmpl`). `simulation_grid_config.R`
isolates `run_mode` and `grid_subset`, so the grid can be edited directly on
the cluster without resubmitting the whole script.

Aggregation (`simulation_aggregate.R`) reduces the per-task `.rds` files in
`results/partial/` into three CSVs at different granularities:

- `summary_replicas.csv` — one row per replica (scalar summaries only).
- `summary_by_t.csv` — one row per method × function × $T$ × time index $t$
  (bias, interval width, and dispersion band of the estimator over time).
- `summary_aggregated.csv` — one row per method × function × $T$, fully
  reduced.

Running this full grid requires an HPC-like environment and was **not**
designed to be reproduced end-to-end on a laptop; it took roughly a weekend
of largely unsupervised, cluster-parallel execution. A local reader who wants
to sanity-check the pipeline should reduce `grid_subset` in
`simulation_grid_config.R` to a handful of configurations and lower `R`.

## Real-data application

`real_data_run.R` fits all five samplers to the `campy` dataset;
`real_data_aggregate.R` collects the results; `figures_real_data.R` produces
the article figures (`article_real_data_fit.pdf`, a five-method overlay with
observed counts, HPD ribbon, and posterior mean lines; and
`article_real_data_ess.pdf`). `predictive_loglik_pf.R` computes out-of-sample
predictive log-likelihood via a particle filter (parallelized with
`parallel::mclapply`, FORK backend). `tscount_literature_benchmark.R` fits
four observation-driven INGARCH-family models via `tscount::tsglm()` for
comparison. All nine methods (5 PoissonLTDM samplers + 4 literature
benchmarks) are combined into a single results table,
`table_predictive_loglik_cputime.tex`.

This part of the pipeline is lightweight and runs on a standard laptop in
minutes — it is the recommended entry point for reproducing a concrete result
from the paper without cluster access.

## Reproducing the results

1. **Install the package** and its dependencies — see `PoissonLTDM/DESCRIPTION`
   for the authoritative list; the packages known to be required are `Rcpp`,
   `rstan`, `Matrix`, `posterior`, `coda`, `batchtools`, `tscount`,
   `this.path`, `testthat`, `parallel`, `foreach`, `doRNG`, and `ggplot2`
   (see [Requirements](#requirements) below).
2. **Sanity-check correctness**: run `tests/test_validation.R` to confirm the
   R prototype and the Rcpp port agree.
3. **Reproduce the real-data application** (fastest path to a paper figure):
   run `real_data_run.R`, then `real_data_aggregate.R` and
   `figures_real_data.R`.
4. **Reproduce a small slice of the simulation study locally**: generate a
   reduced dataset with `data/data_generation.R`, restrict
   `simulation_grid_config.R` to a few configurations, and run
   `simulation/simulation_run.R` directly (outside `batchtools`) for a single
   task.
5. **Reproduce the full simulation study**: requires a PBS Pro (or adaptable)
   HPC cluster; see `calibration/` and `simulation/` for the job submission
   templates and `batchtools` registries.

## Requirements

- **R** 4.5.1 (2025-06-13), `x86_64-conda-linux-gnu`, tested on Debian GNU/Linux
  13 (trixie); BLAS/LAPACK via OpenBLAS 0.3.30 / LAPACK 3.12.0 (conda
  environment). The production simulation additionally ran under R 4.4.3 on
  an external HPC cluster — both are compatible with the packages below.
- **R packages** (pinned versions from the development environment):

  | Package      | Version    |
  |--------------|------------|
  | `Rcpp`       | 1.1.1.1.1  |
  | `rstan`      | 2.32.7     |
  | `Matrix`     | 1.7.5      |
  | `posterior`  | 1.7.0      |
  | `coda`       | 0.19.4.1 <!-- status: possibly no longer imported directly — ESS/R-hat migrated to `posterior` on 2026-08-13; confirm with `grep -rn "coda::" --include="*.R" .` before pinning as a direct dependency --> |
  | `batchtools` | 0.9.18     |
  | `tscount`    | 1.4.3      |
  | `this.path`  | 2.8.0      |
  | `testthat`   | 3.3.2      |
  | `foreach`    | 1.5.2      |
  | `doRNG`      | 1.8.6.3    |
  | `ggplot2`    | 4.0.3      |

- **System**: a working C++ compiler toolchain for `Rcpp`/`rstan` compilation
  (`gcc`/`g++`); PBS Pro (or an adaptable scheduler) only if reproducing the
  full-grid production simulation on a cluster.

A pinned `renv.lock` (via `renv::init()` + `renv::snapshot()` at the
repository root) is the recommended way to freeze this environment for
long-term reproducibility, in place of the manual table above.

## Coming soon

- `CITATION.cff`
- Zenodo archival with a version-specific DOI (post-review)

## Citation

The paper is currently under review; no DOI or camera-ready citation exists
yet. A BibTeX entry will be added here once the paper is accepted, and a
Zenodo concept/version DOI will be added once the repository is archived
(see [Coming soon](#coming-soon)).

## License

MIT.
