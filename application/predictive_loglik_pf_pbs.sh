#!/bin/bash
#PBS -N predictive_loglik_pf
#PBS -l select=1:ncpus=20
#PBS -l place=excl
#PBS -l walltime=08:00:00
#PBS -j oe
#PBS -o predictive_loglik_pf.log

# Plain PBS script, no batchtools -- same reasoning as real_data_pbs.sh:
# this is a single job (5 methods x 100 replicates = 500 independent
# tasks, all fitting in one ncpus=20 exclusive node), not a sweep.
# Submit with: qsub predictive_loglik_pf_pbs.sh
#
# walltime=8h (not a measured worst case -- the actual run should take
# ~5-10 minutes): walltime <= 6h routes to Parallel_shortX, > 6h routes to
# Parallel_longX (see project memory). The first run of this job, with
# walltime=1h, was routed to Parallel_shortX and failed with a PSOCK
# "invalid connection" error right after all 20 workers started -- real_
# data_pbs.sh (walltime=8h, Parallel_longX) has never had this problem
# with the same PSOCK pattern. Matching that queue is the likely actual
# fix, not the PSOCK->FORK switch made earlier (kept anyway as a more
# robust choice regardless of queue, but this may have been the real
# cause all along).
#
# No explicit mem= request: place=excl already reserves the whole node
# exclusively, so the job has access to all of the node's RAM regardless.
# Adding an explicit mem= value gives the scheduler an extra, unnecessary
# criterion to match against (potentially narrowing the set of candidate
# nodes or otherwise slowing queueing) without changing how much memory is
# actually available -- same convention as real_data_pbs.sh.

source ~/setup_env.sh

cd ~/cobalebeb2027/application

Rscript -e 'source("predictive_loglik_pf.R")'
