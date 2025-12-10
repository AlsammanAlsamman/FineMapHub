"""
FINEMAP analysis rule
Step 5: Perform fine-mapping analysis on extracted loci
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module, get_software_value

# Get results directory at module level
RESULTS_DIR = get_results_dir()

rule finemap_analysis:
    """
    FINEMAP analysis for each target analysis
    Processes all loci for a target analysis in a single job
    Only runs for target analyses where type includes "finemap"
    """
    input:
        reference_match_done=f"{RESULTS_DIR}/04_loci/{{target_analysis}}/reference_match.done"
    output:
        done=f"{RESULTS_DIR}/05_finemap/{{target_analysis}}/finemap_analysis.done"
    log:
        f"{RESULTS_DIR}/log/finemap/{{target_analysis}}.log"
    params:
        target_analysis=lambda wildcards: wildcards.target_analysis,
        loci_dir=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}",
        output_dir=lambda wildcards: f"{RESULTS_DIR}/05_finemap/{wildcards.target_analysis}",
        population=lambda wildcards: get_analysis_value(["target_analysis", wildcards.target_analysis, "population"]),
        gwas_table=lambda wildcards: get_analysis_value(["target_analysis", wildcards.target_analysis, "gwas_table"]),
        sample_size=lambda wildcards: next(
            (g["sample_size"] for g in get_analysis_value(["gwastables"]) 
             if g["name"] == get_analysis_value(["target_analysis", wildcards.target_analysis, "gwas_table"])),
            None
        ),
        ref_panel=lambda wildcards: get_analysis_value(["ref_panel", get_analysis_value(["target_analysis", wildcards.target_analysis, "population"])]),
        finemap_module=get_software_module("finemap"),
        plink_module=get_software_module("plink1"),
        r_module=get_software_module("r"),
        max_causal_snps=get_software_value(["finemap", "params", "max_causal_snps"]),
        credible_threshold=get_software_value(["finemap", "params", "credible_threshold"]),
        threads=get_software_value(["finemap", "params", "threads"]),
        significance_threshold=get_analysis_value(["significance_threshold"])
    resources:
        mem_mb=32000,
        time="04:00:00",
        cores=2
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {RESULTS_DIR}/log/finemap
        
        # Validate input directory
        echo "Validating input directory..." > {log}
        if [ ! -d "{params.loci_dir}" ]; then
            echo "ERROR: Loci directory not found: {params.loci_dir}" >> {log}
            exit 1
        fi
        echo "✓ Loci directory validated" >> {log}
        
        # Load required modules
        echo "Loading modules..." >> {log}
        module load {params.finemap_module} >> {log} 2>&1
        module load {params.plink_module} >> {log} 2>&1
        module load {params.r_module} >> {log} 2>&1
        
        # Run FINEMAP analysis for all loci
        echo "Starting FINEMAP analysis for {params.target_analysis}..." >> {log}
        bash scripts/run_finemap_batch.sh \
            --target-analysis {params.target_analysis} \
            --loci-dir {params.loci_dir} \
            --output-dir {params.output_dir} \
            --population {params.population} \
            --sample-size {params.sample_size} \
            --ref-panel {params.ref_panel} \
            --max-causal {params.max_causal_snps} \
            --threads {params.threads} \
            --significance-threshold {params.significance_threshold} \
            >> {log} 2>&1
        
        exit_code=$?
        echo "FINEMAP batch script exit code: $exit_code" >> {log}
        
        # Create done marker
        if [ $exit_code -eq 0 ]; then
            echo "FINEMAP analysis completed successfully for {wildcards.target_analysis}" > {output.done}
            echo "Output directory: {params.output_dir}" >> {output.done}
            echo "Number of loci processed: $(find {params.output_dir} -name "*.config" 2>/dev/null | wc -l)" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
        else
            echo "ERROR: FINEMAP analysis failed - exit code: $exit_code" >> {log}
            exit 1
        fi
        """

rule finemap_analysis_all:
    """
    Aggregate rule to run FINEMAP analysis for all target analyses where type includes "finemap"
    Creates individual .done files for each target analysis to enable selective rerun
    """
    input:
        lambda wildcards: [
            f"{RESULTS_DIR}/05_finemap/{analysis_name}/finemap_analysis.done"
            for analysis_name, analysis_config in get_analysis_value(["target_analysis"]).items()
            if "finemap" in analysis_config.get("type", [])
        ]
    output:
        done=f"{RESULTS_DIR}/05_finemap/finemap_analysis_all.done"
    log:
        f"{RESULTS_DIR}/log/finemap/finemap_analysis_all.log"
    params:
        max_jobs=get_analysis_value(["max_parallel_jobs"])
    resources:
        mem_mb=8000,
        time="00:10:00",
        cores=1
    shell:
        """
        echo "All FINEMAP analysis jobs completed successfully" > {output.done}
        echo "Processed target analyses:" >> {output.done}
        echo "$(ls {RESULTS_DIR}/05_finemap/*/finemap_analysis.done 2>/dev/null | wc -l) analyses completed" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """
