"""
COJO conditional analysis rule
Step 5: Perform iterative GCTA-COJO conditional analysis on extracted loci
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module, get_software_value

# Get results directory at module level
RESULTS_DIR = get_results_dir()

rule cojo_analysis:
    """
    GCTA-COJO conditional analysis for each target analysis
    Processes all loci for a target analysis in a single job
    Only runs for target analyses where type includes "cojo"
    """
    input:
        reference_match_done=f"{RESULTS_DIR}/04_loci/{{target_analysis}}/reference_match.done"
    output:
        done=f"{RESULTS_DIR}/05_cojo/{{target_analysis}}/cojo_analysis.done"
    log:
        f"{RESULTS_DIR}/log/cojo/{{target_analysis}}.log"
    params:
        target_analysis=lambda wildcards: wildcards.target_analysis,
        loci_dir=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}",
        output_dir=lambda wildcards: f"{RESULTS_DIR}/05_cojo/{wildcards.target_analysis}",
        population=lambda wildcards: get_analysis_value(["target_analysis", wildcards.target_analysis, "population"]),
        gwas_table=lambda wildcards: get_analysis_value(["target_analysis", wildcards.target_analysis, "gwas_table"]),
        sample_size=lambda wildcards: next(
            (g["sample_size"] for g in get_analysis_value(["gwastables"]) 
             if g["name"] == get_analysis_value(["target_analysis", wildcards.target_analysis, "gwas_table"])),
            None
        ),
        ref_panel=lambda wildcards: get_analysis_value(["ref_panel", get_analysis_value(["target_analysis", wildcards.target_analysis, "population"])]),
        gcta_module=get_software_module("gcta"),
        plink_module=get_software_module("plink1"),
        python_module=get_software_module("python"),
        max_iterations=get_software_value(["gcta", "cojo_max_iterations"]),
        significance_threshold=get_analysis_value(["significance_threshold"])
    resources:
        mem_mb=32000,
        time="04:00:00",
        cores=2
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {RESULTS_DIR}/log/cojo
        
        # Validate input directory
        echo "Validating input directory..." > {log}
        if [ ! -d "{params.loci_dir}" ]; then
            echo "ERROR: Loci directory not found: {params.loci_dir}" >> {log}
            exit 1
        fi
        echo "✓ Loci directory validated" >> {log}
        
        # Load required modules
        echo "Loading modules..." >> {log}
        module load {params.gcta_module} >> {log} 2>&1
        module load {params.plink_module} >> {log} 2>&1
        module load {params.python_module} >> {log} 2>&1
        
        # Run COJO analysis for all loci
        echo "Starting COJO analysis for {params.target_analysis}..." >> {log}
        bash scripts/run_cojo_batch.sh \
            --target-analysis {params.target_analysis} \
            --loci-dir {params.loci_dir} \
            --output-dir {params.output_dir} \
            --population {params.population} \
            --sample-size {params.sample_size} \
            --ref-panel {params.ref_panel} \
            --max-iterations {params.max_iterations} \
            --significance-threshold {params.significance_threshold} \
            >> {log} 2>&1
        
        exit_code=$?
        echo "COJO batch script exit code: $exit_code" >> {log}
        
        # Create done marker
        if [ $exit_code -eq 0 ]; then
            echo "COJO analysis completed successfully for {wildcards.target_analysis}" > {output.done}
            echo "Output directory: {params.output_dir}" >> {output.done}
            echo "Number of loci processed: $(find {params.output_dir} -name "final_cond_snplist.txt" 2>/dev/null | wc -l)" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
        else
            echo "ERROR: COJO analysis failed - exit code: $exit_code" >> {log}
            exit 1
        fi
        """

rule cojo_analysis_all:
    """
    Aggregate rule to run COJO analysis for all target analyses where type includes "cojo"
    Creates individual .done files for each target analysis to enable selective rerun
    """
    input:
        lambda wildcards: [
            f"{RESULTS_DIR}/05_cojo/{analysis_name}/cojo_analysis.done"
            for analysis_name, analysis_config in get_analysis_value(["target_analysis"]).items()
            if "cojo" in analysis_config.get("type", [])
        ]
    output:
        done=f"{RESULTS_DIR}/05_cojo/cojo_analysis_all.done"
    log:
        f"{RESULTS_DIR}/log/cojo/cojo_analysis_all.log"
    params:
        max_jobs=get_analysis_value(["max_parallel_jobs"])
    resources:
        mem_mb=8000,
        time="00:10:00",
        cores=1
    shell:
        """
        echo "All COJO analysis jobs completed successfully" > {output.done}
        echo "Processed target analyses:" >> {output.done}
        echo "$(ls {RESULTS_DIR}/05_cojo/*/cojo_analysis.done 2>/dev/null | wc -l) analyses completed" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """
