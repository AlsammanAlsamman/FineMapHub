#!/bin/bash

# Simple test to verify find logic
LOCI_DIR="results/04_loci/hisp_euro_chi_analysis_eur"
POPULATION="EUR"

echo "Testing loci directory discovery:"
echo "LOCI_DIR: $LOCI_DIR"

# Test the find command
echo "=== Testing find command ==="
LOCUS_DIRS=$(find "$LOCI_DIR" -maxdepth 1 -type d \( -name "LOC_*" -o -name "newLOC_*" \) | grep -v "/results/")
echo "Found directories:"
echo "$LOCUS_DIRS"
echo ""
echo "Count: $(echo "$LOCUS_DIRS" | wc -l)"
echo ""

# Test file detection for first few loci
echo "=== Testing file detection ==="
for LOCUS_DIR_PATH in $LOCUS_DIRS; do
    if [ -z "$LOCUS_DIR_PATH" ]; then
        continue
    fi
    
    LOCUS_NAME=$(basename "$LOCUS_DIR_PATH")
    MATCHED_LOCI_FILE="$LOCUS_DIR_PATH/${LOCUS_NAME}_gwas_ref_match.tsv"
    ORIGINAL_LOCI_FILE="$LOCUS_DIR_PATH/${LOCUS_NAME}.loci.tsv"
    
    echo "Checking $LOCUS_NAME:"
    echo "  Matched file exists: $(test -f "$MATCHED_LOCI_FILE" && echo "YES" || echo "NO")"
    echo "  Original file exists: $(test -f "$ORIGINAL_LOCI_FILE" && echo "YES" || echo "NO")"
    echo ""
    
    # Just test first 3 loci
    count=$((count + 1))
    if [ "$count" -ge 3 ]; then
        break
    fi
done