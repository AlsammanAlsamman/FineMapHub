#!/bin/bash
# Wrapper to run COJO Manhattan plotting helper

set -euo pipefail

TARGET_ANALYSIS=""
COJO_DIR=""
LOCI_DIR=""
OUTPUT_DIR=""
SIGNIFICANCE_THRESHOLD=""
R_MODULE=""
GENES_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --target-analysis)
            TARGET_ANALYSIS="$2"
            shift 2
            ;;
        --cojo-dir)
            COJO_DIR="$2"
            shift 2
            ;;
        --loci-dir)
            LOCI_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --significance-threshold)
            SIGNIFICANCE_THRESHOLD="$2"
            shift 2
            ;;
        --genes-file)
            GENES_FILE="$2"
            shift 2
            ;;
        --r-module)
            R_MODULE="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$TARGET_ANALYSIS" || -z "$COJO_DIR" || -z "$LOCI_DIR" || -z "$OUTPUT_DIR" || -z "$SIGNIFICANCE_THRESHOLD" || -z "$R_MODULE" || -z "$GENES_FILE" ]]; then
    echo "ERROR: Missing required arguments" >&2
    echo "Required: --target-analysis --cojo-dir --loci-dir --output-dir --significance-threshold --genes-file --r-module" >&2
    exit 1
fi

if [[ ! -d "$COJO_DIR" ]]; then
    echo "ERROR: COJO directory not found: $COJO_DIR" >&2
    exit 1
fi

if [[ ! -d "$LOCI_DIR" ]]; then
    echo "ERROR: Loci directory not found: $LOCI_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

module load "$R_MODULE"

Rscript "${SCRIPT_DIR}/cojo_manhattan_plots.R" \
    --target-analysis "$TARGET_ANALYSIS" \
    --cojo-dir "$COJO_DIR" \
    --loci-dir "$LOCI_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --significance-threshold "$SIGNIFICANCE_THRESHOLD" \
    --genes-file "$GENES_FILE"
