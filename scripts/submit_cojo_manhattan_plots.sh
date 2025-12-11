#!/bin/bash
# Convenience wrapper to run the COJO Manhattan plots rule

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <target-analysis|all> [additional snakemake args...]" >&2
    exit 1
fi

TARGET_ARG="$1"
shift || true

if [[ "$TARGET_ARG" == "all" ]]; then
    TARGET_PATH="results/05_cojo/cojo_manhattan_plots_all.done"
else
    TARGET_PATH="results/05_cojo/${TARGET_ARG}/plots/manhattan_plots.done"
fi

snakemake --snakefile rules/cojo_manhattan_plots.smk "$TARGET_PATH" "$@"
