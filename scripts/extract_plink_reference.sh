#!/bin/bash

# Extract PLINK reference data for all loci in a target analysis
# Usage: extract_plink_reference.sh <loci_dir> <ref_panel> <population> <analysis_name>

set -euo pipefail

LOCI_DIR="$1"
REF_PANEL="$2"
POPULATION="$3"
ANALYSIS_NAME="$4"

echo "=== PLINK REFERENCE DATA EXTRACTION ==="
echo "Analysis: $ANALYSIS_NAME"
echo "Population: $POPULATION"
echo "Reference panel: $REF_PANEL"
echo "Loci directory: $LOCI_DIR"
echo "======================================="

# Validate reference panel files exist
if [ ! -f "${REF_PANEL}.bed" ] || [ ! -f "${REF_PANEL}.bim" ] || [ ! -f "${REF_PANEL}.fam" ]; then
    echo "ERROR: Reference panel files not found: ${REF_PANEL}.{bed,bim,fam}"
    exit 1
fi
echo "✓ Reference panel files validated"

# Check PLINK is available
if ! command -v plink >/dev/null 2>&1; then
    echo "ERROR: PLINK not found in PATH"
    exit 1
fi
echo "✓ PLINK available: $(which plink)"

# Find all locus directories
LOCUS_DIRS=($(find "$LOCI_DIR" -name "*.loci.tsv" -exec dirname {} \; | sort -u))

if [ ${#LOCUS_DIRS[@]} -eq 0 ]; then
    echo "ERROR: No locus directories found in $LOCI_DIR"
    exit 1
fi

echo "Found ${#LOCUS_DIRS[@]} loci to process"

# Process each locus
SUCCESS_COUNT=0
FAILED_COUNT=0

for locus_dir in "${LOCUS_DIRS[@]}"; do
    locus_name=$(basename "$locus_dir")
    echo "Processing locus: $locus_name"
    
    # Find the loci TSV file
    LOCI_FILE=$(find "$locus_dir" -name "*.loci.tsv" | head -1)
    
    if [ ! -f "$LOCI_FILE" ]; then
        echo "  ❌ No .loci.tsv file found in $locus_dir"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        continue
    fi
    
    echo "  📁 Loci file: $(basename "$LOCI_FILE")"
    
    # Create SNP list from loci file
    SNP_LIST="$locus_dir/snps_${POPULATION}.txt"
    
    # Extract SNPID column (try different possible column names)
    if head -1 "$LOCI_FILE" | grep -qi "SNPID"; then
        SNP_COL="SNPID"
    elif head -1 "$LOCI_FILE" | grep -qi "rsid"; then
        SNP_COL="rsid"
    elif head -1 "$LOCI_FILE" | grep -qi "SNP"; then
        SNP_COL="SNP"
    elif head -1 "$LOCI_FILE" | grep -qi "ID"; then
        SNP_COL="ID"
    else
        echo "  ❌ No SNP ID column found in $LOCI_FILE"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        continue
    fi
    
    # Extract SNPs using awk
    awk -F'\t' -v col="$SNP_COL" '
    NR==1 {
        for(i=1; i<=NF; i++) {
            if($i==col) snp_col=i
        }
        if(snp_col==0) {
            print "ERROR: Column " col " not found" > "/dev/stderr"
            exit 1
        }
    }
    NR>1 && $snp_col!="" && $snp_col!="NA" {
        print $snp_col
    }' "$LOCI_FILE" | sort -u > "$SNP_LIST"
    
    if [ $? -ne 0 ] || [ ! -s "$SNP_LIST" ]; then
        echo "  ❌ Failed to create SNP list for $locus_name"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        continue
    fi
    
    SNP_COUNT=$(wc -l < "$SNP_LIST")
    echo "  📊 SNPs to extract: $SNP_COUNT"
    
    # Extract SNPs using PLINK
    OUTPUT_PREFIX="$locus_dir/reference_${POPULATION}"
    
    echo "  🧬 Extracting SNPs with PLINK..."
    plink \
        --bfile "$REF_PANEL" \
        --extract "$SNP_LIST" \
        --make-bed \
        --out "$OUTPUT_PREFIX" \
        --allow-extra-chr \
        --memory 8000 \
        --threads 1 \
        --silent
    
    # Check if extraction was successful
    if [ -f "${OUTPUT_PREFIX}.bed" ] && [ -f "${OUTPUT_PREFIX}.bim" ] && [ -f "${OUTPUT_PREFIX}.fam" ]; then
        EXTRACTED_COUNT=$(wc -l < "${OUTPUT_PREFIX}.bim")
        echo "  ✅ Successfully extracted $EXTRACTED_COUNT SNPs for $locus_name"
        
        # Create extraction summary
        cat > "$locus_dir/extraction_summary_${POPULATION}.txt" << EOF
Reference Data Extraction Summary
================================
Locus: $locus_name
Analysis: $ANALYSIS_NAME
Population: $POPULATION
Reference Panel: $REF_PANEL
Completed: $(date)

SNP Counts:
- SNPs in loci file: $SNP_COUNT
- Successfully extracted: $EXTRACTED_COUNT

Output Files:
- PLINK bed file: reference_${POPULATION}.bed
- PLINK bim file: reference_${POPULATION}.bim
- PLINK fam file: reference_${POPULATION}.fam
- SNP list: snps_${POPULATION}.txt
EOF
        
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "  ❌ PLINK extraction failed for $locus_name"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi
    
    echo ""
done

echo "=== EXTRACTION SUMMARY ==="
echo "Total loci: ${#LOCUS_DIRS[@]}"
echo "Successfully processed: $SUCCESS_COUNT"
echo "Failed: $FAILED_COUNT"
echo "=========================="

if [ $SUCCESS_COUNT -eq 0 ]; then
    echo "ERROR: No loci were successfully processed"
    exit 1
fi

echo "✅ Reference data extraction completed for $ANALYSIS_NAME"