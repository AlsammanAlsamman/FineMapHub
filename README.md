# FineMapHub - Fine-mapping Pipeline

A simplified, modular GWAS fine-mapping pipeline with step-by-step execution and centralized configuration management.

## Project Goals

### Primary Objectives
- **Simplified Configuration**: Replace extremely complex configuration files with just `analysis.yml` and `software.yml`
- **Independent Rules**: Each Snakemake rule should be self-contained and runnable independently
- **Step-by-Step Execution**: Clear sequential workflow with individual `submit.sh` calls
- **Maintainability**: Easier to maintain and adjust compared to previous versions
- **Manual Control**: Run each step manually initially to ensure correctness before automation

### Key Improvements Over Previous Versions
1. **Reduced Complexity**: Streamlined configuration management
2. **Modular Design**: Independent rules that can be run step-by-step
3. **Clear Workflow**: Sequential execution with explicit dependencies
4. **Better Documentation**: Each step clearly documented and explained

## Workflow Architecture

### Sequential Pipeline Steps

The pipeline follows a clear 7-step process:

#### 1. Standardize GWAS Data
- **Purpose**: Convert different GWAS file formats to consistent column naming scheme
- **Input**: Raw GWAS summary statistics files
- **Output**: Standardized GWAS files with consistent column names
- **Rule**: `rules/01_standardize_gwas.smk`

#### 2. Soft Harmonization
- **Purpose**: Harmonize GWAS data against 1000 Genomes reference
- **Input**: Standardized GWAS files
- **Reference**: 1000G VCF and FASTA files
- **Output**: Soft harmonized GWAS files
- **Rule**: `rules/02_soft_harmonization.smk`

#### 3. Population Harmonization
- **Purpose**: Harmonize GWAS against specific population PLINK files
- **Input**: Soft harmonized GWAS files
- **Reference**: Population-specific PLINK files (EUR, AMR, EAS, AFR, etc.)
- **Output**: Population harmonized GWAS files
- **Rule**: `rules/03_population_harmonization.smk`

#### 4. Extract Loci
- **Purpose**: Extract significant regions/loci around target SNPs
- **Input**: Population harmonized GWAS + SNP lists
- **Output**: Loci-specific GWAS files
- **Rule**: `rules/04_extract_loci.smk`

#### 5. Calculate LD Matrices
- **Purpose**: Compute LD matrices for each locus using population reference panels
- **Input**: Loci GWAS files + population reference panels
- **Output**: LD matrices for fine-mapping
- **Rule**: `rules/05_calculate_ld.smk`

#### 6. Fine-mapping Analysis
- **Purpose**: Run multiple fine-mapping methods
- **Methods**: SuSiE-R, SuSiE-inf, FINEMAP, GCTA-COJO
- **Input**: Loci GWAS + LD matrices
- **Output**: Fine-mapping results for each method
- **Rules**: 
  - `rules/06a_susie_finemapping.smk`
  - `rules/06b_finemap_analysis.smk`
  - `rules/06c_cojo_analysis.smk`

#### 7. Reporting and Comparison
- **Purpose**: Generate comparison reports and visualization
- **Input**: All fine-mapping results
- **Output**: Comparison tables, plots, and summary reports
- **Rule**: `rules/07_reporting.smk`

## Configuration Structure

### Simplified Configuration Files

#### `analysis.yml`
```yaml
# Reference files for harmonization
harmonization:
  ref_fasta: "/path/to/human_g1k_v37.fasta"
  ref_vcf: "/path/to/1000g.vcf.gz"

# Population reference panels  
ref_panel:
  EUR: "/path/to/EUR_plink"
  AMR: "/path/to/AMR_plink"
  EAS: "/path/to/EAS_plink"
  AFR: "/path/to/AFR_plink"

# GWAS datasets
gwastables:
  - name: hisp
    samples: 8519
    file: "/path/to/Hisp.tsv"

# Analysis targets
target_analysis:
  hisp_snps_AMR:
    gwas_table: hisp
    snplist: "input/hisp.txt"
    population: AMR
    region: 1000000
    type: ["cojo-cond", "susieR", "susieinf"]
```

#### `software.yml`
```yaml
# Software modules and paths
modules:
  r: "R/4.3.0"
  python: "python/3.9"
  plink: "plink/1.9"
  
# Tool-specific parameters
tools:
  plink:
    memory: 32000
    threads: 4
  susieR:
    max_causal: 10
    coverage: 0.95
  susieinf:
    max_causal: 5
    coverage: 0.95
```

## Execution Model

### Step-by-Step Execution
Each step can be run independently using the main `submit.sh` script:

```bash
#!/bin/bash

# Step 1: Standardize GWAS data
./submit.sh --snakefile rules/01_standardize_gwas.smk results/01_standardized/analysis_standardized.done

# Step 2: Soft harmonization
./submit.sh --snakefile rules/02_soft_harmonization.smk results/02_soft_harmonized/analysis_soft_harmonized.done

# Step 3: Population harmonization  
./submit.sh --snakefile rules/03_population_harmonization.smk results/03_pop_harmonized/analysis_pop_harmonized.done

# Step 4: Extract loci
./submit.sh --snakefile rules/04_extract_loci.smk results/04_loci/analysis_loci_extracted.done

# Step 5: Calculate LD matrices
./submit.sh --snakefile rules/05_calculate_ld.smk results/05_ld_matrices/analysis_ld_calculated.done

# Step 6: Fine-mapping analysis
./submit.sh --snakefile rules/06a_susie_finemapping.smk results/06_finemapping/analysis_susie_done.done
./submit.sh --snakefile rules/06b_finemap_analysis.smk results/06_finemapping/analysis_finemap_done.done

# Step 7: Reporting
./submit.sh --snakefile rules/07_reporting.smk results/07_reports/analysis_reports.done
```

### Design Principles

1. **Independent Rules**: Each rule is self-contained and can be tested independently
2. **Marker Files**: Each step creates a `.done` marker file to track completion
3. **CLI-Only Scripts**: Helper scripts accept all parameters via command-line arguments
4. **No YAML in Scripts**: Scripts never read configuration files directly
5. **Clear Dependencies**: Each step depends on outputs from previous steps

## Project Structure

```
FineMapHub/
├── configs/
│   ├── analysis.yml          # Main analysis configuration
│   └── software.yml          # Software and tool configuration
├── rules/
│   ├── 01_standardize_gwas.smk
│   ├── 02_soft_harmonization.smk
│   ├── 03_population_harmonization.smk
│   ├── 04_extract_loci.smk
│   ├── 05_calculate_ld.smk
│   ├── 06a_susie_finemapping.smk
│   ├── 06b_finemap_analysis.smk
│   ├── 06c_cojo_analysis.smk
│   └── 07_reporting.smk
├── scripts/
│   └── [helper scripts for each rule]
├── utils/
│   └── bioconfigme.py        # Configuration management utility
├── inputs/
│   └── [SNP lists and sample data]
├── results/
│   └── [step-wise outputs]
├── submit.sh                 # Main submission script
├── steps.sh                  # Sequential step execution guide
└── README.md                 # This file
```

## Development Status

### Current Phase
- **Planning and Design**: Defining architecture and requirements
- **Requirements Gathering**: Identifying specific needs and constraints

### Next Steps
1. Answer design questions for formal specification
2. Create project skeleton with simplified configuration
3. Implement standardization rule as first step
4. Progressive implementation of remaining steps
5. Testing and validation of each step
6. Integration and end-to-end testing

## Key Design Questions (Pending)

The following questions need to be answered to proceed with implementation:

1. Primary programming language for helper scripts (Python recommended)
2. CLI-only argument passing to scripts (recommended: yes)
3. Include requirements.txt for dependencies
4. Include example invocations and test data
5. Include unit/smoke tests
6. Independent script execution capability
7. Input validation and error handling level

## Legacy Code Reference

This project builds upon extensive previous work including:
- Multiple Snakemake workflows for different fine-mapping methods
- SuSiE-R, SuSiE-inf, FINEMAP, and GCTA-COJO implementations
- LD calculation methods (PLINK and bigSNPr approaches)
- GWAS harmonization and standardization tools
- Complex configuration management systems

The goal is to simplify and improve upon these existing implementations while preserving their core functionality.

## Contact and Maintenance

**Author**: Alsamman M. Alsamman  
**Project**: FineMapHub - Improved Fine-mapping Pipeline  
**Status**: In Development  
**Last Updated**: November 5, 2025

---

*This README will be updated as the project progresses and requirements are refined.*