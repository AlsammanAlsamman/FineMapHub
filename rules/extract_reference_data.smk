"""
Reference data extraction rule
Step 5: Extract SNPs from population reference panels for each locus
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module

# Get results directory at module level
RESULTS_DIR = get_results_dir()

rule extract_reference_data:
    """
    Extract SNPs from population reference panels for all loci in a target analysis
    Creates PLINK format reference data files in each locus folder
    """
    input:
        loci_done=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}/extract_loci.done"
    output:
        done=f"{RESULTS_DIR}/04_loci/{{target_analysis}}/extract_reference_data.done"
    log:
        f"{RESULTS_DIR}/log/reference_extraction/{{target_analysis}}.log"
    params:
        analysis_name=lambda wildcards: wildcards.target_analysis,
        loci_dir=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}",
        population=lambda wildcards: get_analysis_value(["target_analysis", wildcards.target_analysis, "population"]),
        ref_panel=lambda wildcards: get_analysis_value(["ref_panel", get_analysis_value(["target_analysis", wildcards.target_analysis, "population"])]),
        plink_module=get_software_module("plink2")
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2
    shell:
        """
        mkdir -p {RESULTS_DIR}/log/reference_extraction
        
        # Load PLINK module
        module load {params.plink_module}
        
        # Run reference data extraction
        bash scripts/extract_plink_reference.sh \
            {params.loci_dir} \
            {params.ref_panel} \
            {params.population} \
            {params.analysis_name} \
            > {log} 2>&1
        
        # Create done marker
        if [ $? -eq 0 ]; then
            echo "Reference data extraction completed successfully for {wildcards.target_analysis}" > {output.done}
            echo "Population: {params.population}" >> {output.done}
            echo "Reference panel: {params.ref_panel}" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
        else
            echo "ERROR: Reference data extraction failed" >> {log}
            exit 1
        fi
        """

rule extract_reference_data_all:
    """
    Aggregate rule to run reference data extraction for all target analyses in parallel
    Creates a single done marker after all target analyses are processed
    """
    input:
        expand(f"{RESULTS_DIR}/04_loci/{{target_analysis}}/extract_reference_data.done",
               target_analysis=[analysis for analysis in get_analysis_value(["target_analysis"]).keys()])
    output:
        done=f"{RESULTS_DIR}/04_loci/extract_reference_data_all.done"
    log:
        f"{RESULTS_DIR}/log/reference_extraction/extract_reference_data_all.log"
    resources:
        mem_mb=8000,
        time="00:05:00",
        cores=1
    shell:
        """
        mkdir -p {RESULTS_DIR}/log/reference_extraction
        
        echo "All reference data extraction jobs completed successfully" > {output.done}
        echo "Total target analyses processed: $(echo {input} | wc -w)" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        
        # Log completion details
        echo "Reference data extraction aggregate completed at $(date)" > {log}
        echo "Input files processed:" >> {log}
        for file in {input}; do
            echo "  - $file" >> {log}
        done
        echo "Output marker: {output.done}" >> {log}
        """