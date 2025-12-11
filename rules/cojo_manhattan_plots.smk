"""
COJO Manhattan Plot Generation
Creates before/after conditional Manhattan plots for each COJO iteration
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module

RESULTS_DIR = get_results_dir()

def _cojo_targets():
    return [
        name for name, config in get_analysis_value(["target_analysis"]).items()
        if "cojo" in config.get("type", [])
    ]

rule cojo_manhattan_plots:
    """
    Generate locus-level Manhattan plots before and after COJO conditioning
    """
    input:
        cojo_done=f"{RESULTS_DIR}/05_cojo/{{target_analysis}}/cojo_analysis.done"
    output:
        done=f"{RESULTS_DIR}/05_cojo/{{target_analysis}}/plots/manhattan_plots.done"
    log:
        f"{RESULTS_DIR}/log/cojo_manhattan/{{target_analysis}}.log"
    params:
        target_analysis=lambda wildcards: wildcards.target_analysis,
        cojo_dir=lambda wildcards: f"{RESULTS_DIR}/05_cojo/{wildcards.target_analysis}",
        loci_dir=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}",
        output_dir=lambda wildcards: f"{RESULTS_DIR}/05_cojo/{wildcards.target_analysis}/plots",
        significance_threshold=get_analysis_value(["significance_threshold"]),
        r_module=get_software_module("R"),
        genes_file=get_analysis_value(["genes_loc"])
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2
    shell:
        """
        mkdir -p {RESULTS_DIR}/log/cojo_manhattan
        mkdir -p {params.output_dir}

        bash scripts/run_cojo_manhattan_plots.sh \
            --target-analysis {params.target_analysis} \
            --cojo-dir {params.cojo_dir} \
            --loci-dir {params.loci_dir} \
            --output-dir {params.output_dir} \
            --significance-threshold {params.significance_threshold} \
            --genes-file {params.genes_file} \
            --r-module {params.r_module} \
            > {log} 2>&1

        status=$?
        if [ $status -eq 0 ]; then
            echo "COJO Manhattan plots completed for {params.target_analysis}" > {output.done}
            echo "Target analysis: {params.target_analysis}" >> {output.done}
            echo "Plots directory: {params.output_dir}" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
        else
            exit $status
        fi
        """

rule cojo_manhattan_plots_all:
    """
    Aggregate completion marker for all COJO Manhattan plot jobs
    """
    input:
        lambda wildcards: [
            f"{RESULTS_DIR}/05_cojo/{target}/plots/manhattan_plots.done"
            for target in _cojo_targets()
        ]
    output:
        done=f"{RESULTS_DIR}/05_cojo/cojo_manhattan_plots_all.done"
    log:
        f"{RESULTS_DIR}/log/cojo_manhattan/cojo_manhattan_plots_all.log"
    resources:
        mem_mb=8000,
        time="00:10:00",
        cores=1
    shell:
        """
        echo "COJO Manhattan plots completed for all target analyses" > {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """
