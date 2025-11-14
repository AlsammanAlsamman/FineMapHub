"""
Population-based harmonization rule
Step 3: Harmonize GWAS data with population-specific reference panels
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module

# Get results directory at module level
RESULTS_DIR = get_results_dir()

rule population_harmonize_gwas:
    """
    Population-based harmonization using reference .bim files
    Takes harmonized GWAS from step 02 and performs population-specific harmonization
    """
    input:
        harmonized_file=lambda wildcards: f"{RESULTS_DIR}/02_harmonized/{wildcards.gwas_table}.harmonized.tsv"
    output:
        popharmonized=f"{RESULTS_DIR}/03_popbased_harmonized/{{gwas_table}}_{{population}}.popharmonized.tsv",
        report=f"{RESULTS_DIR}/03_popbased_harmonized/{{gwas_table}}_{{population}}.popharmonization_report.txt",
        done=f"{RESULTS_DIR}/03_popbased_harmonized/{{gwas_table}}_{{population}}.popharmonized.done"
    log:
        f"{RESULTS_DIR}/log/popbased_harmonization/{{gwas_table}}_{{population}}.log"
    params:
        reference_bim=lambda wildcards: get_analysis_value(["ref_panel", wildcards.population]) + ".bim",
        r_module=get_software_module("r")
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2
    shell:
        """
        mkdir -p {RESULTS_DIR}/03_popbased_harmonized
        mkdir -p {RESULTS_DIR}/log/popbased_harmonization
        
        # Load R module
        module load {params.r_module}
        
        # Run population-based harmonization
        Rscript scripts/population_harmonize.R \
            {input.harmonized_file} \
            {params.reference_bim} \
            {RESULTS_DIR}/03_popbased_harmonized/{wildcards.gwas_table}_{wildcards.population} \
            > {log} 2>&1
        
        # Create done marker
        if [ -f "{output.popharmonized}" ] && [ -f "{output.report}" ]; then
            echo "Population harmonization completed successfully for {wildcards.gwas_table} with {wildcards.population} population" > {output.done}
            echo "Files created: {output.popharmonized}, {output.report}" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
        else
            echo "ERROR: Population harmonization failed - output files not created" >> {log}
            exit 1
        fi
        """

rule population_harmonize_all:
    """
    Aggregate rule to harmonize all target analysis combinations in parallel
    Creates individual .done files for each combination to enable selective rerun
    """
    input:
        lambda wildcards: [
            f"{RESULTS_DIR}/03_popbased_harmonized/{analysis_config['gwas_table']}_{analysis_config['population']}.popharmonized.done"
            for analysis_name, analysis_config in get_analysis_value(["target_analysis"]).items()
        ]
    output:
        done=f"{RESULTS_DIR}/03_popbased_harmonized/population_harmonize_all.done"
    log:
        f"{RESULTS_DIR}/log/popbased_harmonization/population_harmonize_all.log"
    params:
        max_jobs=get_analysis_value(["max_parallel_jobs"])
    resources:
        mem_mb=8000,
        time="00:10:00",
        cores=1
    shell:
        """
        echo "All population harmonization jobs completed successfully" > {output.done}
        echo "Processed combinations:" >> {output.done}
        echo "$(ls {RESULTS_DIR}/03_popbased_harmonized/*.popharmonized.done | wc -l) files created" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """