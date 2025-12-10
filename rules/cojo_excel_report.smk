"""
COJO Excel Report Generation
Generates a comprehensive Excel report summarizing COJO conditional analysis results
"""

import sys
sys.path.append("utils")
from bioconfigme import get_results_dir, get_analysis_value, get_software_module

# Get results directory at module level
RESULTS_DIR = get_results_dir()

rule cojo_excel_report:
    """
    Generate Excel report summarizing COJO results across all target analyses
    Creates multi-sheet report with locus summary, SNP details, and distribution statistics
    """
    input:
        cojo_done=f"{RESULTS_DIR}/05_cojo/cojo_analysis_all.done"
    output:
        excel=f"{RESULTS_DIR}/05_cojo/cojo_summary_report.xlsx",
        done=f"{RESULTS_DIR}/05_cojo/cojo_excel_report.done"
    log:
        f"{RESULTS_DIR}/log/cojo/cojo_excel_report.log"
    params:
        cojo_dir=f"{RESULTS_DIR}/05_cojo",
        loci_info_file=get_analysis_value(["loci_info_file"]),
        target_analyses=lambda wildcards: [
            name for name, config in get_analysis_value(["target_analysis"]).items()
            if "cojo" in config.get("type", [])
        ],
        populations=lambda wildcards: str({
            name: config.get("population")
            for name, config in get_analysis_value(["target_analysis"]).items()
            if "cojo" in config.get("type", [])
        }).replace("'", '"'),
        python_module=get_software_module("python")
    resources:
        mem_mb=16000,
        time="01:00:00",
        cores=1
    shell:
        """
        mkdir -p {RESULTS_DIR}/log/cojo
        
        # Load Python module
        echo "Loading Python module..." > {log}
        module load {params.python_module} >> {log} 2>&1
        
        # Generate Excel report
        echo "Generating COJO Excel report..." >> {log}
        
        # Create populations JSON string
        populations_json=$(cat <<'EOF'
{params.populations}
EOF
)
        
        python3 scripts/generate_cojo_excel_report.py \
            --cojo-dir {params.cojo_dir} \
            --loci-info-file {params.loci_info_file} \
            --target-analyses {params.target_analyses} \
            --populations "$populations_json" \
            --output {output.excel} \
            >> {log} 2>&1
        
        exit_code=$?
        echo "Excel report generation exit code: $exit_code" >> {log}
        
        # Create done marker
        if [ $exit_code -eq 0 ] && [ -f {output.excel} ]; then
            echo "COJO Excel report generated successfully" > {output.done}
            echo "Report location: {output.excel}" >> {output.done}
            echo "Timestamp: $(date)" >> {output.done}
        else
            echo "ERROR: Excel report generation failed - exit code: $exit_code" >> {log}
            exit 1
        fi
        """

rule cojo_excel_report_all:
    """
    Alias rule to generate COJO Excel report
    """
    input:
        f"{RESULTS_DIR}/05_cojo/cojo_excel_report.done"
    output:
        done=f"{RESULTS_DIR}/05_cojo/cojo_reports_all.done"
    shell:
        """
        echo "COJO Excel report completed" > {output.done}
        echo "Timestamp: $(date)" >> {output.done}
        """
