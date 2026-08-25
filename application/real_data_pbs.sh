#!/bin/bash
#PBS -N real_data_campy
#PBS -l select=1:ncpus=15
#PBS -l place=excl
#PBS -l walltime=08:00:00
#PBS -j oe
#PBS -o real_data.log

# Plain PBS script, no batchtools -- deliberately (see real_data_run.R's
# header comment). This is a single one-off job (5 methods x K_chains=3 =
# 15 chain-units, all fitting in one ncpus=15 exclusive node), not a sweep,
# so batchtools' registry/job-collection machinery (built for managing many
# jobs) adds nothing here. Submit with: qsub real_data_pbs.sh
#
# ncpus=15 above must match nrow(chain_grid) in real_data_run.R (5 methods
# x K_chains) -- if either changes, update both.
#
# walltime=8h: same "route to the long parallel queue" reasoning as
# simulation_pbs.tmpl/calibration_pbs.tmpl (any walltime <= 6h routes to
# the short queue, capped at ~2-4 concurrent place=excl jobs at last check
# with Euler support) -- not a measured worst-case runtime for this job.

source ~/setup_env.sh

cd ~/cobalebeb2027/application

Rscript -e 'source("real_data_run.R")'
