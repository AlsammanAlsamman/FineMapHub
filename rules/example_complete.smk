"""
Example Snakemake rule demonstrating bioconfigme integration
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_default_resources

# Get configuration values
RESULTS_DIR = get_results_dir()
DEFAULT_RESOURCES = get_default_resources()

rule example_complete:
    """
    Example rule that processes sample data and creates a completion marker.
    Demonstrates proper bioconfigme integration and CLI parameter passing.
    """
    input:
        sample_file="inputs/{sample}.txt"
    output:
        done_marker=f"{RESULTS_DIR}/example/{{sample}}/example_complete.done"
    params:
        sample_name="{sample}",
        python_module=lambda wildcards: get_software_module("python"),
        output_dir=f"{RESULTS_DIR}/example/{{sample}}"
    log:
        f"{RESULTS_DIR}/log/example/{{sample}}/example_complete.log"
    resources:
        mem_mb=DEFAULT_RESOURCES.get("mem_mb", 32000),
        time=DEFAULT_RESOURCES.get("time", "00:30:00"),
        cores=DEFAULT_RESOURCES.get("cores", 2)
    shell:
        """
        # Load Python module
        module load {params.python_module}
        
        # Create output directory
        mkdir -p {params.output_dir}
        mkdir -p $(dirname {log})
        
        # Run helper script with CLI arguments
        python scripts/example_helper.py \\
            --input {input.sample_file} \\
            --sample {params.sample_name} \\
            --output-dir {params.output_dir} \\
            --done-marker {output.done_marker} \\
            > {log} 2>&1
        """