"""
Loci extraction rule
Step 4: Extract genomic loci from population-harmonized GWAS data
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module

# Get results directory at module level
RESULTS_DIR = get_results_dir()

rule extract_loci:
    """
    Extract genomic loci from population-harmonized GWAS data
    Creates individual locus folders with extracted GWAS data for each target analysis
    """
    input:
        popharmonized_file=lambda wildcards: f"{RESULTS_DIR}/03_popbased_harmonized/{get_analysis_value(['target_analysis', wildcards.target_analysis, 'gwas_table'])}_{get_analysis_value(['target_analysis', wildcards.target_analysis, 'population'])}.popharmonized.tsv",
        loci_file=lambda wildcards: get_analysis_value(["target_analysis", wildcards.target_analysis, "snplist"])
    output:
        done=f"{RESULTS_DIR}/04_loci/{{target_analysis}}/extract_loci.done"
    log:
        f"{RESULTS_DIR}/log/loci_extraction/{{target_analysis}}.log"
    params:
        analysis_name=lambda wildcards: wildcards.target_analysis,
        output_dir=lambda wildcards: f"{RESULTS_DIR}/04_loci/{wildcards.target_analysis}",
        python_module=get_software_module("python"),
        numpy_module=get_software_module("numpy")
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p {RESULTS_DIR}/log/loci_extraction
        
        # Validate input files exist
        echo "Validating input files..." > {log}
        if [ ! -f "{input.popharmonized_file}" ]; then
            echo "ERROR: Population harmonized file not found: {input.popharmonized_file}" >> {log}
            exit 1
        fi
        if [ ! -f "{input.loci_file}" ]; then
            echo "ERROR: Loci file not found: {input.loci_file}" >> {log}
            exit 1
        fi
        echo "✓ Input files validated" >> {log}
        
        # Load Python module
        echo "Loading Python module: {params.python_module}" >> {log}
        module load {params.python_module}
        module load {params.numpy_module}
        
        # Check Python and pandas availability
        echo "Checking Python environment..." >> {log}
        python --version >> {log} 2>&1
        python -c "import pandas; print(f'pandas version: {{pandas.__version__}}')" >> {log} 2>&1
        
        # Run loci extraction with detailed logging
        echo "Starting loci extraction..." >> {log}
        python scripts/extract_loci_regions.py \
            --gwas-file {input.popharmonized_file} \
            --loci {input.loci_file} \
            --output-dir {params.output_dir} \
            --analysis-name {params.analysis_name} \
            >> {log} 2>&1
        
        exit_code=$?
        echo "Python script exit code: $exit_code" >> {log}
        
        # Create done marker
        if [ $exit_code -eq 0 ] && [ -d "{params.output_dir}" ] && [ "$(find {params.output_dir} -name "*.loci.tsv" | wc -l)" -gt 0 ]; then
            echo "Loci extraction completed successfully for {wildcards.target_analysis}" > {output.done}
            echo "Output directory: {params.output_dir}" >> {output.done}
            echo "Number of loci extracted: $(find {params.output_dir} -name "*.loci.tsv" | wc -l)" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
        else
            echo "ERROR: Loci extraction failed - exit code: $exit_code" >> {log}
            exit 1
        fi
        """

rule extract_loci_all:
    """
    Aggregate rule to extract loci for all target analysis combinations in parallel
    Creates individual .done files for each combination to enable selective rerun
    """
    input:
        lambda wildcards: [
            f"{RESULTS_DIR}/04_loci/{analysis_name}/extract_loci.done"
            for analysis_name in get_analysis_value(["target_analysis"]).keys()
        ]
    output:
        done=f"{RESULTS_DIR}/04_loci/extract_loci_all.done"
    log:
        f"{RESULTS_DIR}/log/loci_extraction/extract_loci_all.log"
    params:
        max_jobs=get_analysis_value(["max_parallel_jobs"])
    resources:
        mem_mb=8000,
        time="00:10:00",
        cores=1
    shell:
        """
        echo "All loci extraction jobs completed successfully" > {output.done}
        echo "Processed target analyses:" >> {output.done}
        echo "$(ls {RESULTS_DIR}/04_loci/*/extract_loci.done | wc -l) analyses completed" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """