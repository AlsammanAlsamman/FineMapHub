#!/bin/bash

# FINEMAP Batch Analysis Script for Snakemake Integration
# Processes all loci for a single target analysis
# Called by Snakemake rule: finemap_analysis

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
        --max-causal)
            MAX_CAUSAL="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
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
    echo "Usage: $0 --target-analysis <name> --loci-dir <dir> --output-dir <dir> --population <pop> --sample-size <n> --ref-panel <path> [--max-causal <n>] [--threads <n>] [--significance-threshold <p>]"
    exit 1
fi

# Set defaults
MAX_CAUSAL=${MAX_CAUSAL:-2}
THREADS=${THREADS:-4}
SIGNIFICANCE_THRESHOLD=${SIGNIFICANCE_THRESHOLD:-5e-8}

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "============================================================"
echo "FINEMAP Batch Analysis"
echo "============================================================"
echo "Target Analysis: ${TARGET_ANALYSIS}"
echo "Loci Directory: ${LOCI_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Population: ${POPULATION}"
echo "Sample Size: ${SAMPLE_SIZE}"
echo "Reference Panel: ${REF_PANEL}"
echo "Max Causal SNPs: ${MAX_CAUSAL}"
echo "Threads: ${THREADS}"
echo "Significance Threshold: ${SIGNIFICANCE_THRESHOLD}"
echo "============================================================"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Auto-detect executables
FINEMAP_PATH=$(command -v finemap || echo "")
PLINK_PATH=$(command -v plink2 || command -v plink || echo "")

if [ -z "$FINEMAP_PATH" ]; then
    echo "ERROR: FINEMAP executable not found. Make sure finemap module is loaded."
    exit 1
fi

if [ -z "$PLINK_PATH" ]; then
    echo "ERROR: PLINK executable not found. Make sure plink module is loaded."
    exit 1
fi

echo "Using FINEMAP: ${FINEMAP_PATH}"
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
for locus_name in "${LOCI[@]}"; do
    echo ""
    echo "------------------------------------------------------------"
    echo "Processing: ${locus_name}"
    echo "------------------------------------------------------------"
    
    # Set paths
    LOCUS_INPUT_DIR="${LOCI_DIR}/${locus_name}"
    LOCUS_OUTPUT_DIR="${OUTPUT_DIR}/${locus_name}/finemap_results"
    
    GWAS_FILE="${LOCUS_INPUT_DIR}/${locus_name}.loci.tsv"
    BFILE_PREFIX="${LOCUS_INPUT_DIR}/reference_${POPULATION}"
    
    # Validate input files
    if [ ! -f "$GWAS_FILE" ]; then
        echo "ERROR: GWAS file not found: $GWAS_FILE"
        ((FAIL_COUNT++))
        continue
    fi
    
    if [ ! -f "${BFILE_PREFIX}.bed" ]; then
        echo "ERROR: Reference panel not found: ${BFILE_PREFIX}.bed"
        ((FAIL_COUNT++))
        continue
    fi
    
    # Check if already completed
    if [ -f "${LOCUS_OUTPUT_DIR}/${locus_name}.config" ] && \
       [ -f "${LOCUS_OUTPUT_DIR}/${locus_name}.snp" ] && \
       [ -f "${LOCUS_OUTPUT_DIR}/${locus_name}_summary.txt" ]; then
        echo "FINEMAP analysis already completed for ${locus_name}"
        echo "Skipping..."
        ((SKIP_COUNT++))
        continue
    fi
    
    # Create output directory
    mkdir -p "${LOCUS_OUTPUT_DIR}"
    
    # Check if locus has any significant SNPs
    echo "Checking for significant SNPs (p < ${SIGNIFICANCE_THRESHOLD})..."
    
    # Count significant SNPs using awk
    SIG_COUNT=$(awk -v threshold="$SIGNIFICANCE_THRESHOLD" '
        NR==1 { 
            for(i=1; i<=NF; i++) { 
                if($i=="P" || $i=="p" || $i=="PVAL" || $i=="pval") { 
                    p_col=i; 
                    break; 
                } 
            } 
            if(p_col=="") { 
                print "ERROR: P-value column not found" > "/dev/stderr"; 
                exit 1; 
            } 
            next; 
        } 
        { 
            if($p_col != "NA" && $p_col != "" && $p_col+0 < threshold) { 
                count++; 
            } 
        } 
        END { 
            print count+0; 
        }
    ' "$GWAS_FILE")
    
    if [ "$SIG_COUNT" -eq 0 ]; then
        echo "SKIP: No significant SNPs found (p < ${SIGNIFICANCE_THRESHOLD}) for ${locus_name}"
        echo "Creating empty output files..."
        
        # Create empty marker files
        mkdir -p "${LOCUS_OUTPUT_DIR}"
        touch "${LOCUS_OUTPUT_DIR}/${locus_name}.config"
        touch "${LOCUS_OUTPUT_DIR}/${locus_name}.snp"
        echo "SKIPPED: No significant SNPs (threshold: ${SIGNIFICANCE_THRESHOLD})" > "${LOCUS_OUTPUT_DIR}/${locus_name}_summary.txt"
        echo "SKIPPED: No significant SNPs (threshold: ${SIGNIFICANCE_THRESHOLD})" > "${LOCUS_OUTPUT_DIR}/finemap_analysis.log"
        
        ((SKIP_COUNT++))
        continue
    fi
    
    echo "Found ${SIG_COUNT} significant SNPs, proceeding with FINEMAP analysis..."
    
    # Run FINEMAP analysis using the pipeline script
    echo "Starting FINEMAP analysis..."
    
    bash "${SCRIPT_DIR}/finemap_analysis_pipeline.sh" \
        "${BFILE_PREFIX}" \
        "$GWAS_FILE" \
        "${locus_name}" \
        "$SAMPLE_SIZE" \
        --output-dir "$LOCUS_OUTPUT_DIR" \
        --max-causal "$MAX_CAUSAL" \
        --threads "$THREADS" \
        --no-modules \
        2>&1 | tee "${LOCUS_OUTPUT_DIR}/finemap_analysis.log"
    
    exit_code=${PIPESTATUS[0]}
    
    # Check if analysis completed successfully
    if [ $exit_code -eq 0 ] && [ -f "${LOCUS_OUTPUT_DIR}/${locus_name}.config" ]; then
        echo "SUCCESS: FINEMAP analysis completed"
        ((SUCCESS_COUNT++))
    else
        echo "WARNING: Analysis failed or output files missing for ${locus_name}"
        ((FAIL_COUNT++))
    fi
done

# Summary
echo ""
echo "============================================================"
echo "FINEMAP Batch Analysis Summary"
echo "============================================================"
echo "Target Analysis: ${TARGET_ANALYSIS}"
echo "Total Loci: ${NUM_LOCI}"
echo "Successful: ${SUCCESS_COUNT}"
echo "Skipped (no significant SNPs or already done): ${SKIP_COUNT}"
echo "Failed: ${FAIL_COUNT}"
echo "============================================================"

# Exit with error if all loci failed
if [ $SUCCESS_COUNT -eq 0 ] && [ $SKIP_COUNT -eq 0 ]; then
    echo "ERROR: All loci failed to process"
    exit 1
fi

exit 0
