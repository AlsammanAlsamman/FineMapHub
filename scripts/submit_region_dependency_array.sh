#!/bin/bash
#SBATCH --job-name=cojo_region_dependency
#SBATCH --output=logs/region_dependency_array/cojo_region_dependency_%A_%a.out
#SBATCH --error=logs/region_dependency_array/cojo_region_dependency_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

# Usage:
#   sbatch --array=0-9 scripts/submit_region_dependency_array.sh
# The array range should cover ceil(total_selected_snps / CHUNK_SIZE) - 1.

set -euo pipefail

if command -v module &>/dev/null; then
  module load python || true
  module load gcta || true
fi

CHUNK_SIZE=${CHUNK_SIZE:-30}
PYTHON_BIN=${PYTHON_BIN:-python}

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  cd "${SLURM_SUBMIT_DIR}" || exit 1
fi

REPO_ROOT="$(pwd)"
SCRIPT_DIR="${REPO_ROOT}/scripts"

mkdir -p "${REPO_ROOT}/logs/region_dependency_array"

CHUNK_INDEX_VALUE="${SLURM_ARRAY_TASK_ID:-${CHUNK_INDEX:-}}"

if [[ -z "${CHUNK_INDEX_VALUE}" ]]; then
  cat >&2 <<'EOF'
This launcher expects a chunk index.
  • When submitting to Slurm:  sbatch --array=0-N scripts/submit_region_dependency_array.sh
  • For a manual dry run:      CHUNK_INDEX=0 bash scripts/submit_region_dependency_array.sh
Set either SLURM_ARRAY_TASK_ID (via --array) or the CHUNK_INDEX environment variable.
EOF
  exit 1
fi

${PYTHON_BIN} "${SCRIPT_DIR}/cojo_region_dependency.py" \
  --gwas-file "${REPO_ROOT}/results/04_loci/hisp_euro_analysis_eur/LOC_100/LOC_100_gwas_ref_match.tsv" \
  --bfile-prefix "${REPO_ROOT}/results/04_loci/hisp_euro_analysis_eur/LOC_100/reference_EUR" \
  --output-dir "${REPO_ROOT}/results/05_cojo/hisp_euro_analysis_eur/LOC_100/experiments/region_dependency" \
  --sample-size 25217 \
  --use-cojo-slct \
  --chunk-size "${CHUNK_SIZE}" \
  --chunk-index "${CHUNK_INDEX_VALUE}" \
  --region "BLK_REGION:8:11000000:12000000" \
  --region "PPP1R3B_REGION:8:8500000:9000000"
