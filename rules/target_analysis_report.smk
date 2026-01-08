"""
Target Analysis Excel Report Generation
Generates an Excel report for each target analysis containing COJO results
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_software_module, get_analysis_value

# Get configuration values
RESULTS_DIR = get_results_dir()

rule target_analysis_report:
    """
    Generate Excel report for a single target analysis with COJO results
    Creates a comprehensive report with multiple sheets for loci and SNP details
    """
    input:
        cojo_done=f"{RESULTS_DIR}/05_cojo/{{target_analysis}}/cojo_analysis.done"
    output:
        report=f"{RESULTS_DIR}/05_cojo/{{target_analysis}}/{{target_analysis}}_report.xlsx",
        done=f"{RESULTS_DIR}/05_cojo/{{target_analysis}}/{{target_analysis}}_report.done"
    params:
        target_analysis="{target_analysis}",
        cojo_dir=f"{RESULTS_DIR}/05_cojo/{{target_analysis}}",
        loci_dir=f"{RESULTS_DIR}/04_loci/{{target_analysis}}",
        python_module=get_software_module("python"),
        loci_info_file=lambda wildcards: get_analysis_value(["loci_info"], required=False, default="inputs/LOCI_info.txt")
    log:
        f"{RESULTS_DIR}/log/reports/{{target_analysis}}_report.log"
    resources:
        mem_mb=32000,
        time="00:30:00",
        cores=2
    shell:
        """
        # Load Python module
        module load {params.python_module}
        
        # Create output directory
        mkdir -p $(dirname {output.report})
        mkdir -p $(dirname {log})
        
        # Run report generation script
        python scripts/generate_target_analysis_report.py \\
            --target-analysis {params.target_analysis} \\
            --cojo-dir {params.cojo_dir} \\
            --loci-dir {params.loci_dir} \\
            --loci-info-file {params.loci_info_file} \\
            --output {output.report} \\
            > {log} 2>&1
        
        # Create completion marker
        echo "Target analysis report generation completed successfully" > {output.done}
        echo "Target analysis: {params.target_analysis}" >> {output.done}
        echo "Output file: {output.report}" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """

rule target_analysis_report_all:
    """
    Aggregate rule to generate Excel reports for all target analyses with COJO in type
    """
    input:
        lambda wildcards: [
            f"{RESULTS_DIR}/05_cojo/{analysis_name}/{analysis_name}_report.done"
            for analysis_name, analysis_config in get_analysis_value(["target_analysis"]).items()
            if "cojo" in analysis_config.get("type", [])
        ]
    output:
        done=f"{RESULTS_DIR}/05_cojo/target_analysis_reports_all.done"
    log:
        f"{RESULTS_DIR}/log/reports/target_analysis_reports_all.log"
    resources:
        mem_mb=1000,
        time="00:05:00",
        cores=1
    shell:
        """
        echo "All target analysis reports generated successfully" > {output.done}
        echo "Number of reports: $(echo {input} | wc -w)" >> {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        
        # Log completion
        echo "Target analysis report generation aggregate completed at $(date)" > {log}
        for file in {input}; do
            echo "  - $file" >> {log}
        done
        """
