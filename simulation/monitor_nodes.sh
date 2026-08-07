#!/bin/bash
# monitor_nodes.sh
#
# Watches the given PBS job IDs and records which node each one ran on,
# as soon as it transitions to running (state R), by parsing exec_host
# from `qstat -f`. Appends to node_log.txt, one line per job, without
# duplicating a job already recorded. Stops automatically once none of
# the given jobs still appear in `qstat -u $USER` (all finished).
#
# Usage:
#   ./monitor_nodes.sh 4329173 4329174 4329175 ... 4329187
#   ./monitor_nodes.sh $(seq 4329173 4329187)
#
# Run this in the background (nohup ... &) or in a separate terminal/screen
# session while the batch is still in the queue -- it needs to catch each
# job while qstat -f still reports it (running or recently finished),
# since that information is not retained indefinitely.

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <job_id> [<job_id> ...]"
    echo "Example: $0 4329173 4329174 4329175"
    exit 1
fi

JOB_IDS=("$@")
LOG_FILE="node_log.txt"
POLL_INTERVAL=60  # seconds

echo "Monitoring ${#JOB_IDS[@]} job(s): ${JOB_IDS[*]}"
echo "Logging to: $LOG_FILE (polling every ${POLL_INTERVAL}s)"
echo "Press Ctrl+C to stop early (already-recorded jobs are kept)."
echo

touch "$LOG_FILE"

while true; do
    any_still_queued=false

    for jid in "${JOB_IDS[@]}"; do
        full_jid="${jid}.icex"

        # Skip jobs already recorded (avoid duplicate entries on repeated polls)
        if grep -q "^${full_jid} " "$LOG_FILE" 2>/dev/null; then
            continue
        fi

        qstat_out=$(qstat -f "$full_jid" 2>/dev/null)

        if [ -z "$qstat_out" ]; then
            # Job no longer known to the scheduler at all (finished and
            # already purged from qstat -f, or never existed) -- nothing
            # more to record for it.
            continue
        fi

        exec_host=$(echo "$qstat_out" | grep -oP '(?<=exec_host = ).*' | head -1)
        job_state=$(echo "$qstat_out" | grep -oP '(?<=job_state = ).*' | head -1)

        if [ -n "$exec_host" ]; then
            # Node assigned -- record it (works whether the job is
            # currently R or has already finished but qstat -f still
            # remembers it).
            node_name=$(echo "$exec_host" | cut -d'/' -f1)
            echo "${full_jid} node=${node_name} state=${job_state} exec_host=${exec_host} recorded_at=$(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
            echo "Recorded: ${full_jid} -> ${node_name} (state=${job_state})"
        elif [ "$job_state" = "Q" ]; then
            any_still_queued=true
        fi
    done

    # Stop condition: every job either recorded already, or no longer
    # known to qstat -f at all (finished and purged without us catching
    # exec_host -- rare if POLL_INTERVAL is short enough, but possible).
    still_pending=0
    for jid in "${JOB_IDS[@]}"; do
        full_jid="${jid}.icex"
        if ! grep -q "^${full_jid} " "$LOG_FILE" 2>/dev/null; then
            if qstat -f "$full_jid" >/dev/null 2>&1; then
                still_pending=$((still_pending + 1))
            fi
        fi
    done

    if [ "$still_pending" -eq 0 ]; then
        echo
        echo "All jobs recorded or no longer trackable. Stopping."
        break
    fi

    sleep "$POLL_INTERVAL"
done

echo
echo "=== Final node_log.txt ==="
cat "$LOG_FILE"
