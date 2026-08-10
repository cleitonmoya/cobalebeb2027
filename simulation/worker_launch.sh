#!/bin/bash
# Wrapper used as the `rscript=` target for parallel::makePSOCKcluster()
# workers. Invoked as: env TASKSET_CPU=<cpu> worker_launch.sh <args...>
# Kept as a plain script (no nested quoting) because embedding a composed
# "bash -c '...'" string directly in rscript= breaks: parallel/ssh pass it
# through as a single literal token instead of letting a shell parse it,
# causing "No such file or directory" on the remote side.
source ~/setup_env.sh > /dev/null 2>&1
exec taskset -c "${TASKSET_CPU}" Rscript "$@"
