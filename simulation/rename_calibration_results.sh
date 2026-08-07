#!/bin/bash
# rename_calibration_results.sh
#
# Renames existing calibration .rds (results/calibration/) and .pdf
# (results/calibration/plots/) files to the new naming convention that
# includes N/burnin/K, so they are not confused with future tasks that
# use a different configuration for the same (method, f, Tt).
#
# Old:  <method>_<f>_<Tt>.rds
# New:  <method>_<f>_<Tt>_N<N>_b<burnin>[_K<K>].rds
#
# The mapping below was generated directly from the current summary.csv
# (the source of truth for which N/burnin/K each existing .rds actually
# used) -- see the conversation this came from if you need to regenerate
# it for a different summary.csv.
#
# Run from the directory containing results/ (i.e. simulation/).
#
# SAFE BY DEFAULT: dry-run (only prints what would be renamed). Pass
# --apply to actually rename.
#
# Usage:
#   ./rename_calibration_results.sh            # dry-run, prints only
#   ./rename_calibration_results.sh --apply     # actually renames

set -e

APPLY=false
if [ "$1" == "--apply" ]; then
    APPLY=true
fi

RESULTS_DIR="results/calibration"
PLOTS_DIR="results/calibration/plots"

# old_name new_name (space-separated), one pair per line
MAPPING="
stan_constant_200 stan_constant_200_N11000_b1000
stan_quadratic_200 stan_quadratic_200_N11000_b1000
amh_montoril_constant_1600 amh_montoril_constant_1600_N55000_b5000
amh_montoril_constant_200 amh_montoril_constant_200_N55000_b5000
amh_montoril_quadratic_1600 amh_montoril_quadratic_1600_N55000_b5000
amh_montoril_quadratic_200 amh_montoril_quadratic_200_N55000_b5000
pg_as_constant_1600 pg_as_constant_1600_N11000_b1000_K100
pg_as_constant_200 pg_as_constant_200_N11000_b1000_K100
pg_as_quadratic_1600 pg_as_quadratic_1600_N11000_b1000_K100
pg_as_quadratic_200 pg_as_quadratic_200_N11000_b1000_K100
sir_collapsed_constant_1600 sir_collapsed_constant_1600_N11000_b1000
sir_collapsed_constant_200 sir_collapsed_constant_200_N11000_b1000
sir_collapsed_quadratic_1600 sir_collapsed_quadratic_1600_N11000_b1000
sir_collapsed_quadratic_200 sir_collapsed_quadratic_200_N11000_b1000
sir_laplace_constant_1600 sir_laplace_constant_1600_N11000_b1000
sir_laplace_constant_200 sir_laplace_constant_200_N11000_b1000
sir_laplace_quadratic_1600 sir_laplace_quadratic_1600_N11000_b1000
sir_laplace_quadratic_200 sir_laplace_quadratic_200_N11000_b1000
stan_constant_1600 stan_constant_1600_N11000_b1000
stan_quadratic_1600 stan_quadratic_1600_N11000_b1000
"

if [ "$APPLY" = false ]; then
    echo "=== DRY RUN (pass --apply to actually rename) ==="
fi
echo

n_renamed=0
n_missing=0

while read -r old_name new_name; do
    [ -z "$old_name" ] && continue

    old_rds="${RESULTS_DIR}/${old_name}.rds"
    new_rds="${RESULTS_DIR}/${new_name}.rds"
    old_pdf="${PLOTS_DIR}/${old_name}.pdf"
    new_pdf="${PLOTS_DIR}/${new_name}.pdf"

    if [ -f "$old_rds" ]; then
        if [ -f "$new_rds" ]; then
            echo "SKIP (target already exists): $new_rds"
        else
            echo "${old_rds} -> ${new_rds}"
            if [ "$APPLY" = true ]; then
                mv "$old_rds" "$new_rds"
            fi
            n_renamed=$((n_renamed + 1))
        fi
    else
        echo "MISSING (no .rds found, skipping): $old_rds"
        n_missing=$((n_missing + 1))
    fi

    if [ -f "$old_pdf" ]; then
        if [ -f "$new_pdf" ]; then
            echo "SKIP (target already exists): $new_pdf"
        else
            echo "${old_pdf} -> ${new_pdf}"
            if [ "$APPLY" = true ]; then
                mv "$old_pdf" "$new_pdf"
            fi
        fi
    fi

done <<< "$MAPPING"

echo
echo "=== Summary: ${n_renamed} .rds renamed, ${n_missing} not found ==="
if [ "$APPLY" = false ]; then
    echo "This was a dry run -- re-run with --apply to actually rename files."
fi
