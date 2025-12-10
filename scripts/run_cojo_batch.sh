#!/bin/bash

# COJO Batch Analysis Script for Snakemake Integration
# Processes all loci for a single target analysis
# Called by Snakemake rule: cojo_analysis

set -e

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --target-analysis)
            TARGET_ANALYSIS="$2"
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
        --population)
            POPULATION="$2"
            shift 2
            ;;
        --sample-size)
            SAMPLE_SIZE="$2"
            shift 2
            ;;
        --ref-panel)
            REF_PANEL="$2"
            shift 2
            ;;
        --max-iterations)
            MAX_ITERATIONS="$2"
            shift 2
            ;;
        --significance-threshold)
            SIGNIFICANCE_THRESHOLD="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            exit 1
            ;;
    esac
done

# Validate required parameters
if [ -z "$TARGET_ANALYSIS" ] || [ -z "$LOCI_DIR" ] || [ -z "$OUTPUT_DIR" ] || \
   [ -z "$POPULATION" ] || [ -z "$SAMPLE_SIZE" ] || [ -z "$REF_PANEL" ]; then
    echo "ERROR: Missing required parameters"
    echo "Usage: $0 --target-analysis <name> --loci-dir <dir> --output-dir <dir> --population <pop> --sample-size <n> --ref-panel <path> [--max-iterations <n>] [--significance-threshold <p>]"
    exit 1
fi

# Set defaults
MAX_ITERATIONS=${MAX_ITERATIONS:-20}
SIGNIFICANCE_THRESHOLD=${SIGNIFICANCE_THRESHOLD:-5e-8}

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================================"
echo "COJO Batch Analysis"
echo "============================================================"
echo "Target Analysis: ${TARGET_ANALYSIS}"
echo "Loci Directory: ${LOCI_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Population: ${POPULATION}"
echo "Sample Size: ${SAMPLE_SIZE}"
echo "Reference Panel: ${REF_PANEL}"
echo "Max Iterations: ${MAX_ITERATIONS}"
echo "Significance Threshold: ${SIGNIFICANCE_THRESHOLD}"
echo "============================================================"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Auto-detect GCTA and PLINK executables
# Try multiple possible names for GCTA
for gcta_cmd in gcta64 gcta gcta-1.94.1 gcta_1.94.1; do
    if command -v "$gcta_cmd" &> /dev/null; then
        GCTA_PATH="$gcta_cmd"
        break
    fi
done

# Try multiple possible names for PLINK
for plink_cmd in plink plink2 plink1.9; do
    if command -v "$plink_cmd" &> /dev/null; then
        PLINK_PATH="$plink_cmd"
        break
    fi
done

if [ -z "$GCTA_PATH" ]; then
    echo "ERROR: GCTA executable not found. Make sure gcta module is loaded."
    echo "Tried: gcta64, gcta, gcta-1.94.1, gcta_1.94.1"
    exit 1
fi

if [ -z "$PLINK_PATH" ]; then
    echo "ERROR: PLINK executable not found. Make sure plink module is loaded."
    echo "Tried: plink, plink2, plink1.9"
    exit 1
fi

echo "Using GCTA: ${GCTA_PATH}"
echo "Using PLINK: ${PLINK_PATH}"
echo "============================================================"

# Find all LOC_* directories
LOCI=()
for dir in "${LOCI_DIR}"/LOC_*; do
    [ -d "$dir" ] || continue
    locus_name=$(basename "$dir")
    
    # Skip merged loci (containing "-") and newLOC directories
    [[ "$locus_name" == *"-"* ]] && continue
    [[ "$locus_name" == newLOC* ]] && continue
    
    LOCI+=("$locus_name")
done

NUM_LOCI=${#LOCI[@]}
echo "Found ${NUM_LOCI} loci to process"

if [ $NUM_LOCI -eq 0 ]; then
    echo "WARNING: No loci found in ${LOCI_DIR}"
    exit 0
fi

# Initialize counters
SUCCESS_COUNT=0
SKIP_COUNT=0
FAIL_COUNT=0

# Process each locus
CURRENT_LOCUS=0
for locus_name in "${LOCI[@]}"; do
    CURRENT_LOCUS=$((CURRENT_LOCUS + 1))
    echo ""
    echo "------------------------------------------------------------"
    echo "Processing: ${locus_name} (${CURRENT_LOCUS}/${NUM_LOCI})"
    echo "------------------------------------------------------------"
    
    # Set paths
    LOCUS_INPUT_DIR="${LOCI_DIR}/${locus_name}"
    LOCUS_OUTPUT_DIR="${OUTPUT_DIR}/${locus_name}/cojo_results"
    
    GWAS_FILE="${LOCUS_INPUT_DIR}/${locus_name}_gwas_ref_match.tsv"
    BFILE_PREFIX="${LOCUS_INPUT_DIR}/reference_${POPULATION}"
    
    # Validate input files
    if [ ! -f "$GWAS_FILE" ]; then
        echo "ERROR: Matched GWAS file not found: $GWAS_FILE"
        echo "Make sure reference_match step has been completed for this locus"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi
    
    if [ ! -f "${BFILE_PREFIX}.bed" ]; then
        echo "ERROR: Reference panel not found: ${BFILE_PREFIX}.bed"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi
    
    # Check if already completed
    if [ -f "${LOCUS_OUTPUT_DIR}/final_cond_snplist.txt" ] && \
       [ -f "${LOCUS_OUTPUT_DIR}/final_cojo_results_merged.txt" ] && \
       [ -f "${LOCUS_OUTPUT_DIR}/iteration_summary.json" ]; then
        echo "COJO analysis already completed for ${locus_name}"
        echo "Skipping..."
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi
    
    # Create output directory
    mkdir -p "${LOCUS_OUTPUT_DIR}"
    
    # Check if locus has any significant SNPs
    echo "Checking for significant SNPs (p < ${SIGNIFICANCE_THRESHOLD})..."
    
    # Count significant SNPs - simpler and faster approach
    SIG_COUNT=$(awk -F'\t' -v threshold="$SIGNIFICANCE_THRESHOLD" '
        NR==1 { 
            for(i=1; i<=NF; i++) { 
                if($i=="P" || $i=="p" || $i=="PVAL" || $i=="pval" || $i=="P-value" || $i=="p-value") { 
                    p_col=i; 
                    break; 
                } 
            } 
            next; 
        } 
        p_col && $p_col != "NA" && $p_col != "" && $p_col+0 < threshold { count++ } 
        END { print count+0 }
    ' "$GWAS_FILE")
    
    # Check if SIG_COUNT is valid
    if ! [[ "$SIG_COUNT" =~ ^[0-9]+$ ]]; then
        echo "WARNING: Could not determine significant SNP count for ${locus_name}"
        echo "SIG_COUNT value: '${SIG_COUNT}'"
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi
    
    if [ "$SIG_COUNT" -eq 0 ]; then
        echo "SKIP: No significant SNPs found (p < ${SIGNIFICANCE_THRESHOLD}) for ${locus_name}"
        echo "Creating empty output files..."
        
        # Create empty output files
        mkdir -p "${LOCUS_OUTPUT_DIR}"
        touch "${LOCUS_OUTPUT_DIR}/final_cond_snplist.txt"
        echo -e "iteration\ttotal_iterations\tconverged\tfinal_snps\titerations" > "${LOCUS_OUTPUT_DIR}/iteration_summary.json"
        echo '{"iterations":[],"final_snps":[],"total_iterations":0,"converged":true}' > "${LOCUS_OUTPUT_DIR}/iteration_summary.json"
        touch "${LOCUS_OUTPUT_DIR}/final_cojo_results_merged.txt"
        echo "SKIPPED: No significant SNPs (threshold: ${SIGNIFICANCE_THRESHOLD})" > "${LOCUS_OUTPUT_DIR}/cojo_analysis.log"
        
        SKIP_COUNT=$((SKIP_COUNT + 1))
        continue
    fi
    
    echo "Found ${SIG_COUNT} significant SNPs, proceeding with COJO analysis..."
    
    # Run COJO analysis using the iterative Python script directly
    echo "Starting COJO analysis..."
    
    python3 "${SCRIPT_DIR}/run_cojo_iterative.py" \
        --gwas-file "$GWAS_FILE" \
        --bfile-prefix "$BFILE_PREFIX" \
        --output-dir "$LOCUS_OUTPUT_DIR" \
        --sample-size "$SAMPLE_SIZE" \
        --max-iterations "$MAX_ITERATIONS" \
        --significance-threshold "$SIGNIFICANCE_THRESHOLD" \
        --gcta-path "$GCTA_PATH" \
        --plink-path "$PLINK_PATH" \
        2>&1 | tee "${LOCUS_OUTPUT_DIR}/cojo_analysis.log"
    
    exit_code=${PIPESTATUS[0]}
    
    # Check if analysis completed successfully
    if [ $exit_code -eq 0 ] && [ -f "${LOCUS_OUTPUT_DIR}/final_cond_snplist.txt" ]; then
        NUM_SIGNALS=$(wc -l < "${LOCUS_OUTPUT_DIR}/final_cond_snplist.txt")
        echo "SUCCESS: Found ${NUM_SIGNALS} independent signals"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "WARNING: Analysis failed or output files missing for ${locus_name}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

# Summary
echo ""
echo "============================================================"
echo "COJO Batch Analysis Summary"
echo "============================================================"
echo "Target Analysis: ${TARGET_ANALYSIS}"
echo "Total Loci: ${NUM_LOCI}"
echo "Successful: ${SUCCESS_COUNT}"
echo "Skipped (already done): ${SKIP_COUNT}"
echo "Failed: ${FAIL_COUNT}"
echo "============================================================"

# Exit with error if all loci failed
if [ $SUCCESS_COUNT -eq 0 ] && [ $SKIP_COUNT -eq 0 ]; then
    echo "ERROR: All loci failed to process"
    exit 1
fi

exit 0
