#!/bin/bash

# Helper script to run the COJO regional dependency experiment for
# hisp_euro_analysis_eur LOC_100.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

GWAS_FILE="${REPO_ROOT}/results/04_loci/hisp_euro_analysis_eur/LOC_100/LOC_100_gwas_ref_match.tsv"
BFILE_PREFIX="${REPO_ROOT}/results/04_loci/hisp_euro_analysis_eur/LOC_100/reference_EUR"
OUTPUT_DIR="${REPO_ROOT}/results/05_cojo/hisp_euro_analysis_eur/LOC_100/experiments/region_dependency"
SAMPLE_SIZE=25217

mkdir -p "${OUTPUT_DIR}"

python3 "${SCRIPT_DIR}/cojo_region_dependency.py" \
  --gwas-file "${GWAS_FILE}" \
  --bfile-prefix "${BFILE_PREFIX}" \
  --output-dir "${OUTPUT_DIR}" \
  --sample-size "${SAMPLE_SIZE}" \
  --region "BLK_REGION:8:11000000:12000000" \
  --region "PPP1R3B_REGION:8:8500000:9000000"

Rscript "${SCRIPT_DIR}/plot_region_dependency.R" \
  --summary "${OUTPUT_DIR}/region_dependency_summary.tsv" \
  --output-prefix "${OUTPUT_DIR}/region_dependency_scree"
