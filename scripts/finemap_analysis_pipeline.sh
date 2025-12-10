#!/bin/bash
#==============================================================================
# FINEMAP Analysis Pipeline - Integrated for Snakemake
# Purpose: Run FINEMAP on already-matched GWAS and reference data
# Uses: Pre-matched data from reference_match step
#==============================================================================

set -e

# Default parameters
DEFAULT_OUTPUT_DIR="finemap_results"
DEFAULT_MAX_CAUSAL="2"
DEFAULT_THREADS="4"
DEFAULT_QC_MIN_MAF="0.001"
DEFAULT_QC_MAX_MAF="0.5"

LOAD_MODULES=false
KEEP_TEMP=false
VERBOSE=false

# Usage
usage() {
    cat << EOF
Usage: $0 <BFILE> <GWAS_FILE> <PREFIX> <NSAMPLES> [OPTIONS]

REQUIRED:
  BFILE         PLINK file prefix (already matched with GWAS)
  GWAS_FILE     GWAS summary statistics (already matched with reference)
  PREFIX        Output prefix
  NSAMPLES      Sample size

OPTIONS:
  --output-dir DIR      Output directory (default: $DEFAULT_OUTPUT_DIR)
  --max-causal INT      Max causal SNPs (default: $DEFAULT_MAX_CAUSAL)
  --threads INT         Threads (default: $DEFAULT_THREADS)
  --keep-temp           Keep temporary files
  --verbose             Verbose output
  --no-modules          Skip module loading (tools in PATH)
EOF
    exit 0
}

# Parse arguments
if [ $# -lt 4 ]; then
    usage
fi

BFILE="$1"
GWAS_FILE="$2"
PREFIX="$3"
NSAMPLES="$4"
shift 4

OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
MAX_CAUSAL="$DEFAULT_MAX_CAUSAL"
THREADS="$DEFAULT_THREADS"
QC_MIN_MAF="$DEFAULT_QC_MIN_MAF"
QC_MAX_MAF="$DEFAULT_QC_MAX_MAF"

while [[ $# -gt 0 ]]; do
    case $1 in
        --output-dir) OUTPUT_DIR="$2"; shift 2;;
        --max-causal) MAX_CAUSAL="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --keep-temp) KEEP_TEMP=true; shift;;
        --verbose) VERBOSE=true; shift;;
        --no-modules) LOAD_MODULES=false; shift;;
        --help|-h) usage;;
        *) echo "Unknown option: $1"; usage;;
    esac
done

# Validate files
if [ ! -f "${BFILE}.bed" ] || [ ! -f "${BFILE}.bim" ] || [ ! -f "${BFILE}.fam" ]; then
    echo "ERROR: PLINK files not found: ${BFILE}.{bed,bim,fam}"
    exit 1
fi

if [ ! -f "$GWAS_FILE" ]; then
    echo "ERROR: GWAS file not found: $GWAS_FILE"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
SCRIPT_DIR="$(dirname "${BASH_SOURCE[0]}")"

echo "============================================================"
echo "FINEMAP Pipeline - Using Pre-Matched Data"
echo "============================================================"
echo "BFILE: $BFILE"
echo "GWAS: $GWAS_FILE"
echo "PREFIX: $PREFIX"
echo "NSAMPLES: $NSAMPLES"
echo "OUTPUT: $OUTPUT_DIR"
echo "MAX_CAUSAL: $MAX_CAUSAL"
echo "THREADS: $THREADS"
echo "============================================================"

cd "$OUTPUT_DIR"

#==============================================================================
# Step 1: Calculate MAF from PLINK
#==============================================================================
echo "[Step 1] Calculating MAF from PLINK data..."
plink --bfile "$BFILE" --freq --out "${PREFIX}_maf"

#==============================================================================
# Step 2: Compute LD matrix
#==============================================================================
echo "[Step 2] Computing LD matrix..."
plink --bfile "$BFILE" --r square --out "${PREFIX}_ldmatrix"

#==============================================================================
# Step 3: Format LD matrix
#==============================================================================
echo "[Step 3] Formatting LD matrix..."
cat > format_ld.R << 'EOFLD'
args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]
ld_file <- args[2]
prefix <- args[3]

bim <- read.table(bim_file, header=FALSE, stringsAsFactors=FALSE)
n_snps <- nrow(bim)
ld_data <- read.table(ld_file, header=FALSE)

ld_matrix <- as.matrix(ld_data)
if(nrow(ld_matrix) != n_snps || ncol(ld_matrix) != n_snps) {
    ld_matrix <- ld_matrix[1:n_snps, 1:n_snps]
}

# Fix diagonal and NA values
ld_matrix[is.na(ld_matrix)] <- 0
diag(ld_matrix) <- 1

write.table(ld_matrix, paste0(prefix, ".ld"),
            row.names=FALSE, col.names=FALSE, quote=FALSE)
cat("LD matrix formatted:", nrow(ld_matrix), "x", ncol(ld_matrix), "\n")
EOFLD

Rscript format_ld.R "${BFILE}.bim" "${PREFIX}_ldmatrix.ld" "$PREFIX"

#==============================================================================
# Step 4: Prepare Z file and SNP file
#==============================================================================
echo "[Step 4] Preparing Z and SNP files..."
cat > prepare_finemap_input.R << 'EOF'
args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]
gwas_file <- args[2]
maf_file <- args[3]
prefix <- args[4]
nsamples <- as.numeric(args[5])
qc_min_maf <- as.numeric(args[6])
qc_max_maf <- as.numeric(args[7])

# Read files
bim <- read.table(bim_file, header=FALSE, stringsAsFactors=FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")
gwas <- read.table(gwas_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
maf_data <- read.table(paste0(maf_file, ".frq"), header=TRUE, stringsAsFactors=FALSE)

# Merge data
merged <- merge(bim, gwas, by.x="SNP", by.y="SNPID", all.x=TRUE)
merged <- merge(merged, maf_data[,c("SNP","MAF")], by="SNP", all.x=TRUE)

# Order by BIM
merged <- merged[match(bim$SNP, merged$SNP),]

# Clean MAF values
merged$MAF[is.na(merged$MAF) | merged$MAF <= 0] <- qc_min_maf
merged$MAF[merged$MAF > qc_max_maf] <- qc_max_maf

# Create Z file (beta/se)
z_data <- data.frame(
    rsid = merged$SNP,
    chromosome = merged$CHR,
    position = merged$POS,
    allele1 = merged$EA,
    allele2 = merged$NEA,
    maf = merged$MAF,
    beta = merged$BETA,
    se = merged$SE
)

# Calculate z-score
z_data$z <- z_data$beta / z_data$se

# Remove rows with NA
z_data <- z_data[complete.cases(z_data),]

# Write Z file
write.table(z_data[,c("rsid","chromosome","position","allele1","allele2","maf","beta","se")],
            paste0(prefix, ".z"), row.names=FALSE, quote=FALSE)

# Write SNP file
snp_data <- z_data[,c("rsid","chromosome","position","allele1","allele2","maf")]
write.table(snp_data, paste0(prefix, ".snp"), row.names=FALSE, quote=FALSE)

cat("Prepared", nrow(z_data), "SNPs for FINEMAP\n")
cat("MAF range:", min(z_data$maf), "-", max(z_data$maf), "\n")
EOF

Rscript prepare_finemap_input.R "${BFILE}.bim" "$GWAS_FILE" "${PREFIX}_maf" "$PREFIX" "$NSAMPLES" "$QC_MIN_MAF" "$QC_MAX_MAF"

#==============================================================================
# Step 5: Create master file
#==============================================================================
echo "[Step 5] Creating FINEMAP master file..."
echo "z;ld;snp;config;cred;n_samples;log" > "${PREFIX}.master"
echo "${PREFIX}.z;${PREFIX}.ld;${PREFIX}.snp;${PREFIX}.config;${PREFIX}.cred;${NSAMPLES};${PREFIX}.log" >> "${PREFIX}.master"

#==============================================================================
# Step 6: Run FINEMAP
#==============================================================================
echo "[Step 6] Running FINEMAP..."
finemap --sss --in-files "${PREFIX}.master" --n-threads "$THREADS" --n-causal-snps "$MAX_CAUSAL"

#==============================================================================
# Step 7: Run post-processing analysis (both quick and best)
#==============================================================================
echo "[Step 7] Running post-processing analysis..."

if [ -f "${PREFIX}.config" ] && [ -f "${PREFIX}.snp" ]; then
    # Run quick analysis (all k combined)
    echo "Running quick analysis..."
    Rscript "${SCRIPT_DIR}/analyze_finemap_results.R" quick "$PREFIX" "${GWAS_FILE}" "${PREFIX}.ld" \
        --ld-high 0.8 --ld-low 0.2 --top-snps 5 || echo "Warning: Quick analysis failed"
    
    # Run best k analysis
    echo "Running best k analysis..."
    Rscript "${SCRIPT_DIR}/analyze_finemap_results.R" best "$PREFIX" "${GWAS_FILE}" "${PREFIX}.ld" \
        --max-k "$MAX_CAUSAL" --ld-high 0.8 --ld-low 0.2 --top-snps 5 || echo "Warning: Best k analysis failed"
else
    echo "Warning: FINEMAP output files not found, skipping post-processing"
fi

#==============================================================================
# Cleanup
#==============================================================================
if [ "$KEEP_TEMP" = false ]; then
    rm -f format_ld.R prepare_finemap_input.R
    rm -f "${PREFIX}_ldmatrix.ld" "${PREFIX}_maf.frq"
fi

cd - > /dev/null

echo "============================================================"
echo "FINEMAP Analysis Complete"
echo "============================================================"
echo "Results in: $OUTPUT_DIR"
echo "Config: ${OUTPUT_DIR}/${PREFIX}.config"
echo "SNP results: ${OUTPUT_DIR}/${PREFIX}.snp"
echo "Quick analysis: ${OUTPUT_DIR}/${PREFIX}_summary.txt"
echo "Best k analysis: ${OUTPUT_DIR}/${PREFIX}_best_credible_analysis.txt"
echo "============================================================"

exit 0
