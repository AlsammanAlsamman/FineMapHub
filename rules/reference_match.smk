#!/usr/bin/env python3

import sys
sys.path.append("utils")
from bioconfigme import get_software_module, get_results_dir, get_analysis_value

# Get results directory at module level
RESULTS_DIR = get_results_dir()

rule reference_match:
    """
    Match GWAS summary statistics to PLINK BIM files for each locus
    Creates filtered and ordered GWAS files for downstream fine-mapping
    """
    input:
        extract_done = "results/04_loci/{target_analysis}/extract_reference_data.done"
    output:
        done = "results/04_loci/{target_analysis}/reference_match.done"
    params:
        r_module = get_software_module("r"),
        population = lambda wildcards: get_analysis_value(["target_analysis", wildcards.target_analysis, "population"]),
        ref_panel_prefix = lambda wildcards: get_analysis_value(["ref_panel", get_analysis_value(["target_analysis", wildcards.target_analysis, "population"])])
    log:
        "results/log/reference_match/{target_analysis}.log"
    resources:
        mem_mb = 32000,
        time = "00:30:00",
        cores = 2
    shell:
        """
        # Redirect all output to log file for debugging
        exec > {log} 2>&1
        
        echo "=== REFERENCE MATCHING DEBUG LOG ===" 
        echo "Timestamp: $(date)"
        echo "Target analysis: {wildcards.target_analysis}"
        echo "Population: {params.population}"
        echo "Reference panel prefix: {params.ref_panel_prefix}"
        echo "Working directory: $(pwd)"
        echo "Shell: $BASH_VERSION"
        echo "Environment modules available: $(module avail 2>&1 | head -5 || echo 'Module command not available')"
        
        # Test basic commands
        echo "Testing basic commands..."
        echo "ls command: $(which ls || echo 'ls not found')"
        echo "basename command: $(which basename || echo 'basename not found')"
        echo "Rscript command: $(which Rscript || echo 'Rscript not found')"
        
        # Load required modules
        echo "Loading R module: {params.r_module}"
        if module load {params.r_module} 2>&1; then
            echo "✓ R module loaded successfully"
        else
            echo "✗ Failed to load R module: {params.r_module}"
            exit 1
        fi
        
        # Check R availability after module load
        echo "R version after module load: $(Rscript --version 2>&1 || echo 'Rscript still not available')"
        
        # Process all loci directories - only direct children, not nested
        LOCI_DIR="results/04_loci/{wildcards.target_analysis}"
        echo "Loci directory: $LOCI_DIR"
        echo "Loci directory exists: $(test -d "$LOCI_DIR" && echo 'YES' || echo 'NO')"
        
        if [ ! -d "$LOCI_DIR" ]; then
            echo "ERROR: Loci directory does not exist: $LOCI_DIR"
            exit 1
        fi
        
        echo "Contents of loci directory:"
        ls -la "$LOCI_DIR" || echo "Failed to list directory contents"
        
        TOTAL_LOCI=0
        SUCCESSFUL_LOCI=0
        
        echo "Starting reference matching for {wildcards.target_analysis}"
        
        # Get list of locus directories - bash compatible approach
        echo "Looking for LOC_* and newLOC_* directories..."
        
        # Test shell capabilities
        echo "Testing shopt command: $(shopt -s nullglob 2>&1 || echo 'shopt not available')"
        
        # Use shell globbing instead of find
        shopt -s nullglob 2>/dev/null || echo "Warning: nullglob not supported"
        LOCUS_DIRS=("$LOCI_DIR"/LOC_* "$LOCI_DIR"/newLOC_*)
        shopt -u nullglob 2>/dev/null || echo "Warning: nullglob disable not supported"
        
        echo "Found locus directories: ${{#LOCUS_DIRS[@]}} items"
        for dir in "${{LOCUS_DIRS[@]}}"; do
            echo "  - $dir"
        done
        
        for LOCUS_DIR_PATH in "${{LOCUS_DIRS[@]}}"; do
            echo "Processing directory: $LOCUS_DIR_PATH"
            
            # Skip if no matches (empty array)
            if [ ! -d "$LOCUS_DIR_PATH" ]; then
                echo "  Skipping - not a directory: $LOCUS_DIR_PATH"
                continue
            fi
            
            # Skip nested results directories
            if [[ "$LOCUS_DIR_PATH" == *"/results/"* ]]; then
                echo "  Skipping - nested results directory: $LOCUS_DIR_PATH"
                continue
            fi
            
            LOCUS_NAME=$(basename "$LOCUS_DIR_PATH")
            echo "  Locus name: $LOCUS_NAME"
            
            BIM_FILE="$LOCUS_DIR_PATH/reference_{params.population}.bim"
            GWAS_FILE="$LOCUS_DIR_PATH/$LOCUS_NAME.loci.tsv"
            OUTPUT_FILE="$LOCUS_DIR_PATH/${{LOCUS_NAME}}_gwas_ref_match.tsv"
            
            echo "  Checking files:"
            echo "    BIM file: $BIM_FILE (exists: $(test -f "$BIM_FILE" && echo 'YES' || echo 'NO'))"
            echo "    GWAS file: $GWAS_FILE (exists: $(test -f "$GWAS_FILE" && echo 'YES' || echo 'NO'))"
            echo "    Output file: $OUTPUT_FILE"
            
            if [[ -f "$BIM_FILE" && -f "$GWAS_FILE" ]]; then
                echo "  Processing $LOCUS_NAME..."
                TOTAL_LOCI=$((TOTAL_LOCI + 1))
                
                echo "  Running R script: Rscript scripts/match_gwas_to_bim.R \"$BIM_FILE\" \"$GWAS_FILE\" \"$OUTPUT_FILE\""
                if Rscript scripts/match_gwas_to_bim.R "$BIM_FILE" "$GWAS_FILE" "$OUTPUT_FILE" 2>&1; then
                    echo "  ✓ Success: $LOCUS_NAME"
                    SUCCESSFUL_LOCI=$((SUCCESSFUL_LOCI + 1))
                else
                    echo "  ✗ Failed: $LOCUS_NAME (R script returned error code $?)"
                fi
            else
                echo "  Skipping $LOCUS_NAME - missing required files"
                if [ ! -f "$BIM_FILE" ]; then
                    echo "    Missing: $BIM_FILE"
                fi
                if [ ! -f "$GWAS_FILE" ]; then
                    echo "    Missing: $GWAS_FILE"
                fi
            fi
            echo ""
        done
        
        echo "=== PROCESSING SUMMARY ==="
        echo "Reference matching completed: $SUCCESSFUL_LOCI/$TOTAL_LOCI loci processed successfully"
        echo "Analysis: {wildcards.target_analysis}"
        echo "End timestamp: $(date)"
        
        # Create done marker
        echo "Reference matching completed for {wildcards.target_analysis}" > {output.done}
        echo "Successful loci: $SUCCESSFUL_LOCI" >> {output.done}
        echo "Total loci: $TOTAL_LOCI" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        
        echo "Done marker created at: {output.done}"
        echo "Script completed successfully"
        """

rule reference_match_all:
    """
    Aggregate rule to run reference matching for all target analyses in parallel
    Creates individual .done files for each analysis to enable selective rerun
    """
    input:
        lambda wildcards: [
            f"{RESULTS_DIR}/04_loci/{analysis_name}/reference_match.done"
            for analysis_name in get_analysis_value(["target_analysis"]).keys()
        ]
    output:
        done="results/04_loci/reference_match_all.done"
    log:
        "results/log/reference_match/reference_match_all.log"
    params:
        max_jobs=get_analysis_value(["max_parallel_jobs"])
    resources:
        mem_mb=8000,
        time="00:10:00",
        cores=1
    shell:
        """
        echo "All reference matching jobs completed successfully" > {output.done}
        echo "Processed analyses:" >> {output.done}
        echo "$(ls {RESULTS_DIR}/04_loci/*/reference_match.done | wc -l) analyses completed" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """