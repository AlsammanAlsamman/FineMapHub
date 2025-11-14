"""
FINEMAP rule
Step 6: Fine-mapping analysis using FINEMAP for target analyses configured with finemap type
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module, get_software_params, get_gwas_sample_size

# Get results directory at module level
RESULTS_DIR = get_results_dir()

rule finemap:
    """
    Run FINEMAP analysis for all loci in a target analysis
    Only processes target analyses that have 'finemap' in their type list
    """
    input:
        ref_done=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}/extract_reference_data.done",
        match_done=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}/reference_match.done"
    output:
        done=f"{RESULTS_DIR}/05_finemapping/{{target_analysis}}/finemap/finemap.done"
    log:
        f"{RESULTS_DIR}/log/finemap/{{target_analysis}}.log"
    params:
        analysis_name=lambda wildcards: wildcards.target_analysis,
        loci_dir=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}",
        output_dir=lambda wildcards: f"{RESULTS_DIR}/05_finemapping/{wildcards.target_analysis}/finemap",
        population=lambda wildcards: get_analysis_value(["target_analysis", wildcards.target_analysis, "population"]),
        sample_size=lambda wildcards: get_gwas_sample_size(wildcards.target_analysis),
        finemap_module=get_software_module("finemap"),
        plink_module="plink2/1.90b3w",
        r_module=get_software_module("R"),
        max_causal_snps=lambda wildcards: get_software_params("finemap").get("max_causal_snps", "2"),
        min_pval=lambda wildcards: get_software_params("finemap").get("min_pval", "0.05")
    resources:
        mem_mb=32000,
        time="02:00:00",
        cores=lambda wildcards: int(get_software_params("finemap").get("threads", "4"))
    shell:
        """
        # Check if this target analysis should run FINEMAP
        analysis_types=$(python3 -c "
import sys
sys.path.append('utils')
from bioconfigme import get_analysis_value
types = get_analysis_value(['target_analysis', '{wildcards.target_analysis}', 'type'])
if isinstance(types, list):
    print(' '.join(types))
else:
    print(types)
")
        
        if [[ ! "$analysis_types" =~ "finemap" ]]; then
            echo "Skipping {wildcards.target_analysis} - finemap not in analysis types: $analysis_types"
            echo "FINEMAP analysis skipped - not configured for this target analysis" > {output.done}
            echo "Analysis types: $analysis_types" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
            exit 0
        fi
        
        mkdir -p {RESULTS_DIR}/log/finemap
        mkdir -p {params.output_dir}
        
        # Load required modules - modules are loaded inside finemap_analysis.sh
        # module load {params.finemap_module}
        # module load {params.plink_module}
        # module load {params.r_module}
        
        # Find the reference PLINK files and loci file for this analysis
        LOCI_DIR="{params.loci_dir}"
        POPULATION="{params.population}"
        
        # Look for reference PLINK files (created by extract_reference_data rule)
        BFILE=$(find "$LOCI_DIR" -name "reference_${{POPULATION}}.bed" | head -1 | sed 's/\.bed$//')
        
        if [ -z "$BFILE" ]; then
            echo "ERROR: No reference PLINK files found for population ${{POPULATION}} in ${{LOCI_DIR}}" >> {log}
            exit 1
        fi
        
        # Find the loci TSV file
        LOCIFILE=$(find "$LOCI_DIR" -name "*.loci.tsv" | head -1)
        
        if [ -z "$LOCIFILE" ]; then
            echo "ERROR: No loci TSV file found in ${{LOCI_DIR}}" >> {log}
            exit 1
        fi
        
        echo "Found BFILE: $BFILE" >> {log}
        echo "Found LOCIFILE: $LOCIFILE" >> {log}
        
        # Run FINEMAP analysis with the new script
        bash scripts/finemap_analysis.sh \
            "$BFILE" \
            "$LOCIFILE" \
            "{wildcards.target_analysis}" \
            "{params.sample_size}" \
            "{params.output_dir}" \
            "{params.max_causal_snps}" \
            "{params.min_pval}" \
            > {log} 2>&1
        
        # Create done marker
        if [ $? -eq 0 ]; then
            echo "FINEMAP analysis completed successfully for {wildcards.target_analysis}" > {output.done}
            echo "Population: {params.population}" >> {output.done}
            echo "Sample size: {params.sample_size}" >> {output.done}
            echo "Max causal SNPs: {params.max_causal_snps}" >> {output.done}
            echo "Min p-value: {params.min_pval}" >> {output.done}
            echo "Threads: {resources.cores}" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
        else
            echo "ERROR: FINEMAP analysis failed" >> {log}
            exit 1
        fi
        """

rule finemap_all:
    """
    Aggregate rule to run FINEMAP analysis for all target analyses in parallel
    Creates individual .done files for each analysis to enable selective rerun
    """
    input:
        lambda wildcards: [
            f"{RESULTS_DIR}/05_finemapping/{analysis_name}/finemap/finemap.done"
            for analysis_name in get_analysis_value(["target_analysis"]).keys()
            if "finemap" in get_analysis_value(["target_analysis", analysis_name, "type"])
        ]
    output:
        done="results/05_finemapping/finemap_all.done"
    log:
        "results/log/finemap/finemap_all.log"
    params:
        max_jobs=get_analysis_value(["max_parallel_jobs"])
    resources:
        mem_mb=8000,
        time="00:10:00",
        cores=1
    shell:
        """
        echo "All FINEMAP analyses completed successfully" > {output.done}
        echo "Processed analyses:" >> {output.done}
        echo "$(ls {RESULTS_DIR}/05_finemapping/*/finemap/finemap.done | wc -l) analyses completed" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """