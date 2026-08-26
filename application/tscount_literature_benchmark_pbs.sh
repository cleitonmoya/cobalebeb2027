#!/bin/bash
#PBS -N tscount_benchmark
#PBS -l select=1:ncpus=20
#PBS -l place=excl
#PBS -l walltime=08:00:00
#PBS -j oe
#PBS -o tscount_benchmark.log

# ncpus=20, NOT ncpus=1: routing to Parallel_longX (vs. the sequential
# Seq_longX "Fat Node" pool) depends on the requested ncpus matching a
# parallel-node resource profile, not just on walltime -- a ncpus=1
# request gets routed to Seq_* regardless of place=excl (confirmed: this
# job landed in Seq_longX with ncpus=1). 20 matches the physical core
# count of the parallel nodes used elsewhere in the project
# (simulation_run.R/calibration_run.R/real_data_run.R all reserve a
# Parallel_longX node via place=excl at their own ncpus). The R script
# itself still only uses 1 core (single-threaded, see below) -- the other
# 19 are requested purely to land on, and exclusively reserve, the
# correct node type for a fair comparison against real_data_run.R's
# timings, not because the script needs them.
#
# Plain PBS script, single-core R job -- but place=excl anyway (reserves
# the WHOLE node regardless of ncpus requested, same as real_data_pbs.sh),
# NOT the shared Seq_*/single-core pool that would otherwise fit a
# ncpus=1 resource profile. This is deliberate, not an oversight: the
# CPU-time comparison against the 5 PoissonLTDM methods (run in
# real_data_run.R, also under place=excl) is only fair if BOTH sides are
# measured free of contention from other users' co-scheduled jobs on the
# same physical node -- contention can inflate user.self itself (not just
# wall time), via cache/memory-bus pressure from neighboring processes,
# not merely scheduling wait. A shared Seq_* node would reintroduce
# exactly the noise source the repeated-block timing scheme
# (B_timing/R_timing in tscount_literature_benchmark.R) was designed to
# average out.
#
# walltime=8h: NOT a measured runtime (the script itself, R_timing=30 x
# B_timing=50 x 4 models = 6000 tsglm() calls, should take well under 10
# minutes) -- set purely to route to Parallel_longX rather than
# Parallel_shortX (walltime <= 6h routes to the short queue -- see
# project memory / real_data_pbs.sh's identical reasoning). The job exits
# as soon as the script finishes; walltime is only an upper bound, not
# time actually reserved/consumed.
#
# Single-threaded BLAS, explicit: tsglm()'s internal optimizer
# (constrOptim/BFGS, Hessian via numDeriv) can invoke multi-threaded BLAS
# if the cluster's R build links against one -- forcing 1 thread here
# keeps user.self measuring the algorithm's own serial CPU cost, matching
# the (single-threaded) Rcpp samplers on the other side of the
# comparison, and avoiding the "user.self can exceed wall time under
# multi-threaded BLAS" distortion already flagged elsewhere in this
# project.

source ~/setup_env.sh

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd ~/cobalebeb2027/application

Rscript -e 'source("tscount_literature_benchmark.R")'
