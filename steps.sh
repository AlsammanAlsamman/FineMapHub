#!/bin/bash
"""
FineMapHub Pipeline Step-by-Step Execution Guide

This script contains the sequential commands for running the FineMapHub pipeline.
Each step should be executed manually to ensure proper completion before proceeding.
"""

set -euo pipefail

# =============================================================================
# EXAMPLE STEP - Template for future pipeline steps
# =============================================================================

# Example processing step
# ./submit.sh --snakefile rules/example_complete.smk results/example/test_sample/example_complete.done

# Example with specific sample
# ./submit.sh --snakefile rules/example_complete.smk --cores 2 results/example/my_sample/example_complete.done

# =============================================================================
# FUTURE PIPELINE STEPS (to be added as rules are created)
# =============================================================================

# Step 1: Standardize GWAS data
# ./submit.sh --snakefile rules/standardize_gwas.smk results/01_standardized/analysis_standardized.done

# Step 2: Soft harmonization against 1000G reference
# ./submit.sh --snakefile rules/soft_harmonization.smk results/02_harmonized/analysis_harmonized.done

# Step 3: Population-based harmonization with reference panels
./submit.sh --snakefile rules/population_harmonize.smk --cores 1 results/03_popbased_harmonized/hisp_AMR.popharmonized.done

# Step 3 (All): Population-based harmonization for all target analyses
# 6 parallel jobs as per Target analysis
./submit.sh --snakefile rules/population_harmonize.smk --cores 4 --jobs 6 results/03_popbased_harmonized/population_harmonize_all.done

# Step 4: Loci extraction for target analysis
./submit.sh --snakefile rules/extract_loci.smk --cores 1 results/04_loci/hisp_euro_chi_analysis_eur/extract_loci.done

# Step 4 (All): Loci extraction for all target analyses
./submit.sh --snakefile rules/extract_loci.smk --cores 4 results/04_loci/extract_loci_all.done

# Step 5: Reference data extraction for target analysis
./submit.sh --snakefile rules/extract_reference_data.smk --cores 1 results/04_loci/hisp_euro_chi_analysis_eur/extract_reference_data.done

# Step 5 (All): Reference data extraction for all target analyses
./submit.sh --snakefile rules/extract_reference_data.smk --cores 4 --jobs 6 results/04_loci/extract_reference_data_all.done

# Step 5.5: Reference matching - match GWAS to PLINK BIM files
./submit.sh --snakefile rules/reference_match.smk --cores 2 results/04_loci/hisp_euro_chi_analysis_eur/reference_match.done

# Step 5.5 (All): Reference matching for all target analyses
./submit.sh --snakefile rules/reference_match.smk --cores 4 --jobs 6 results/04_loci/reference_match_all.done




