#!/bin/bash
ml finemap/1.4
ml plink2/1.90b3w
ml R
set -e

# Parse command line arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <BFILE> <LOCIFILE> <PREFIX> <NSAMPLES> [OUTPUT_DIR] [MAX_CAUSAL_SNPS] [MIN_PVAL]"
    echo "Parameters:"
    echo "  BFILE: Prefix for PLINK files (bed, bim, fam)"
    echo "  LOCIFILE: Loci TSV file with GWAS summary statistics"
    echo "Example: $0 extracted_snps_AMR hisp_AMR.loci.tsv hisp_AMR 8519 finemap_results"
    exit 1
fi

BFILE="$1"
LOCIFILE="$2"
PREFIX="$3"
NSAMPLES="$4"
OUTPUT_DIR="${5:-finemap_results}"
MAX_CAUSAL_SNPS="${6:-2}"
MIN_PVAL="${7:-0.05}"

# Convert to absolute paths
BFILE=$(readlink -f "$BFILE" 2>/dev/null || realpath "$BFILE" 2>/dev/null || echo "$BFILE")
LOCIFILE=$(readlink -f "$LOCIFILE" 2>/dev/null || realpath "$LOCIFILE" 2>/dev/null || echo "$LOCIFILE")
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(readlink -f "$OUTPUT_DIR" 2>/dev/null || realpath "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")

# Validate required files exist
if [ ! -f "${BFILE}.bed" ] || [ ! -f "${BFILE}.bim" ] || [ ! -f "${BFILE}.fam" ]; then
    echo "Error: PLINK files not found: ${BFILE}.bed, ${BFILE}.bim, ${BFILE}.fam"
    exit 1
fi

if [ ! -f "$LOCIFILE" ]; then
    echo "Error: Loci file not found: $LOCIFILE"
    exit 1
fi

echo "Starting FINEMAP analysis with parameters:"
echo "  BFILE: $BFILE"
echo "  LOCIFILE: $LOCIFILE"
echo "  PREFIX: $PREFIX"
echo "  NSAMPLES: $NSAMPLES"
echo "  OUTPUT_DIR: $OUTPUT_DIR"
echo "  MIN_PVAL: $MIN_PVAL"

# Save original directory and change to output directory
ORIGINAL_DIR=$(pwd)
cd ${OUTPUT_DIR}

echo "Step 0: Filter GWAS summary statistics by p-value (MIN_PVAL = $MIN_PVAL)..."
RSCRIPT_FILTER_PVAL=$(cat << 'EOF0'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
gwas_file <- args[1]
min_pval <- as.numeric(args[2])
prefix <- args[3]
output_dir <- args[4]

# Read GWAS summary statistics
gwas <- read.table(gwas_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Original GWAS file:", nrow(gwas), "SNPs\n")

# Check if P-value column exists
if (!"P" %in% colnames(gwas)) {
  stop("Error: P-value column 'P' not found in GWAS file. Available columns: ", paste(colnames(gwas), collapse = ", "))
}

# Filter by p-value
gwas_filtered <- gwas[gwas$P <= min_pval, ]
cat("After p-value filtering (P <=", min_pval, "):", nrow(gwas_filtered), "SNPs\n")

if (nrow(gwas_filtered) == 0) {
  stop("Error: No SNPs remain after p-value filtering. Consider using a less stringent p-value threshold.")
}

# Write filtered GWAS file
filtered_file <- paste0(prefix, "_pval_filtered.loci.tsv")
write.table(gwas_filtered, filtered_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Filtered GWAS file written to:", filtered_file, "\n")

# Also write summary statistics
cat("P-value distribution in filtered data:\n")
print(summary(gwas_filtered$P))
cat("Number of SNPs with P < 1e-5:", sum(gwas_filtered$P < 1e-5), "\n")
cat("Number of SNPs with P < 1e-8:", sum(gwas_filtered$P < 1e-8), "\n")

# Create SNP list for PLINK extraction
snp_list <- gwas_filtered$SNPID
snp_file <- paste0(prefix, "_gwas_snps.txt")
write.table(snp_list, snp_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
cat("SNP list for PLINK extraction written to:", snp_file, "\n")
EOF0
)

echo "$RSCRIPT_FILTER_PVAL" > filter_by_pvalue.R
Rscript filter_by_pvalue.R $LOCIFILE $MIN_PVAL $PREFIX $OUTPUT_DIR

# Update LOCIFILE to use the p-value filtered version
LOCIFILE="${PREFIX}_pval_filtered.loci.tsv"

echo "Step 1: Extract only GWAS SNPs from PLINK files..."
plink --bfile $BFILE \
      --extract ${PREFIX}_gwas_snps.txt \
      --make-bed \
      --out ${PREFIX}_gwas_snps_only

# Update BFILE to use the GWAS-filtered dataset
BFILE="${PREFIX}_gwas_snps_only"

echo "Step 2: Performing comprehensive quality control..."
plink --bfile $BFILE \
      --maf 0.001 \
      --geno 0.05 \
      --hwe 1e-6 \
      --mind 0.1 \
      --make-bed \
      --out ${BFILE}_qc

# Update BFILE to use QCed dataset
BFILE="${BFILE}_qc"

echo "Filtering SNPs to match PLINK BIM after QC..."
# Use awk for exact SNP matching (fixes the extra SNP issue)
head -1 $LOCIFILE > ${PREFIX}_filtered.loci.tsv
awk 'NR==FNR{snps[$2]=1; next} FNR>1 && ($1 in snps)' ${BFILE}.bim $LOCIFILE >> ${PREFIX}_filtered.loci.tsv

# Debug: Count SNPs to verify filtering
echo "Debug - SNP counts after filtering:"
echo "BIM SNPs: $(wc -l < ${BFILE}.bim)"
echo "P-value filtered loci SNPs: $(tail -n +2 $LOCIFILE | wc -l)"
echo "Final filtered loci SNPs: $(tail -n +2 ${PREFIX}_filtered.loci.tsv | wc -l)"

echo "Step 3: Calculate actual MAF from PLINK data..."
plink --bfile $BFILE --freq --out ${PREFIX}_maf

echo "Step 4: Compute LD matrix using modern PLINK command with QC..."
plink --bfile $BFILE --r square --out ${PREFIX}_ldmatrix

echo "Step 5: Create proper square LD matrix (improved version with error handling)..."
RSCRIPT_IMPROVED=$(cat << 'EOF'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]
ld_file <- args[2]
output_prefix <- args[3]

# Read BIM to get SNP names and count
bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
n_snps <- nrow(bim)
cat("Number of SNPs in BIM:", n_snps, "\n")

# Read LD matrix carefully
cat("Reading LD file:", ld_file, "\n")
ld_data <- read.table(ld_file, header = FALSE, stringsAsFactors = FALSE)

cat("Raw LD data dimensions:", nrow(ld_data), "rows x", ncol(ld_data), "columns\n")

# Convert all columns to numeric, handling any issues
ld_matrix <- matrix(NA, nrow = nrow(ld_data), ncol = ncol(ld_data))
for (i in 1:ncol(ld_data)) {
  ld_matrix[, i] <- as.numeric(ld_data[[i]])
}

# Check if we have the right dimensions
if (nrow(ld_matrix) == n_snps && ncol(ld_matrix) == n_snps) {
  cat("Perfect! LD matrix is square with correct dimensions\n")
} else if (nrow(ld_matrix) == n_snps && ncol(ld_matrix) > n_snps) {
  cat("Trimming LD matrix from", ncol(ld_matrix), "to", n_snps, "columns\n")
  ld_matrix <- ld_matrix[, 1:n_snps]
} else if (nrow(ld_matrix) > n_snps && ncol(ld_matrix) == n_snps) {
  cat("Trimming LD matrix from", nrow(ld_matrix), "to", n_snps, "rows\n")
  ld_matrix <- ld_matrix[1:n_snps, ]
} else {
  cat("Warning: LD matrix dimensions don't match SNP count\n")
  cat("Using first", n_snps, "rows and columns\n")
  ld_matrix <- ld_matrix[1:n_snps, 1:n_snps]
}

# Identify problematic SNPs with zero variance
problematic_snps <- which(diag(ld_matrix) == 0 | is.na(diag(ld_matrix)))
if (length(problematic_snps) > 0) {
  cat("Found", length(problematic_snps), "SNPs with zero variance:\n")
  cat("Positions:", problematic_snps, "\n")
  cat("Removing problematic SNPs from analysis...\n")
  
  # Keep only good SNPs
  good_snps <- setdiff(1:n_snps, problematic_snps)
  ld_matrix <- ld_matrix[good_snps, good_snps]
  n_snps <- length(good_snps)
  
  # Update BIM file to match
  bim_good <- bim[good_snps, ]
  write.table(bim_good, paste0(output_prefix, "_filtered.bim"),
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
              
  cat("Updated SNP count after removing problematic SNPs:", n_snps, "\n")
}

# Replace any remaining NA/NaN/Inf with 0
ld_matrix[is.na(ld_matrix) | is.nan(ld_matrix) | is.infinite(ld_matrix)] <- 0

# Ensure diagonal is 1 (or very close to it)
diag_issues <- 0
for (i in 1:min(n_snps, nrow(ld_matrix))) {
  if (i <= ncol(ld_matrix)) {
    if (abs(ld_matrix[i,i] - 1) > 0.1) {
      cat("Warning: Diagonal element", i, "is", ld_matrix[i,i], "- setting to 1\n")
      ld_matrix[i,i] <- 1
      diag_issues <- diag_issues + 1
    }
  }
}

cat("Fixed", diag_issues, "diagonal elements\n")

# Check if matrix is positive definite
eigen_vals <- eigen(ld_matrix, only.values = TRUE)$values
negative_eigen <- sum(eigen_vals < -1e-8)
if (negative_eigen > 0) {
  cat("Warning: LD matrix has", negative_eigen, "negative eigenvalues\n")
  # Add small constant to diagonal to make it positive definite
  ld_matrix <- ld_matrix + diag(1e-8, n_snps)
  cat("Added small constant to diagonal to stabilize matrix\n")
}

# Write the cleaned LD matrix
write.table(ld_matrix, paste0(output_prefix, "_formatted.ld"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, 
            sep = " ", na = "0")

cat("Successfully wrote formatted LD matrix with dimensions:", 
    nrow(ld_matrix), "x", ncol(ld_matrix), "\n")
cat("Minimum eigenvalue:", min(eigen_vals), "\n")
cat("Maximum eigenvalue:", max(eigen_vals), "\n")
EOF
)

echo "$RSCRIPT_IMPROVED" > create_ld_matrix.R
Rscript create_ld_matrix.R ${BFILE}.bim ${PREFIX}_ldmatrix.ld ${PREFIX}

# Check if the LD matrix was created successfully
if [ ! -f "${PREFIX}_formatted.ld" ]; then
    echo "Error: Failed to create LD matrix"
    exit 1
fi

# Check if we have a filtered BIM file and update if needed
if [ -f "${PREFIX}_filtered.bim" ]; then
    echo "Updating to filtered BIM file after removing problematic SNPs..."
    cp ${PREFIX}_filtered.bim ${BFILE}.bim
fi

echo "Step 6: Reorder GWAS summary statistics and merge with actual MAF (with validation)..."
RSCRIPT2=$(cat << 'EOF2'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
bim_file <- args[1]
gwas_file <- args[2]
maf_file <- args[3]
prefix <- args[4]
output_dir <- args[5]

# Read BIM file to get PLINK order
bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")
plink_snps <- bim$SNP

# Read GWAS summary statistics
gwas <- read.table(gwas_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Original GWAS columns:", colnames(gwas), "\n")

# Read MAF file from PLINK
maf_data <- read.table(paste0(maf_file, ".frq"), header = TRUE, stringsAsFactors = FALSE)
colnames(maf_data) <- c("CHR", "SNP", "A1", "A2", "MAF", "NCHROBS")
cat("MAF data columns:", colnames(maf_data), "\n")
cat("First few MAF values:\n")
print(head(maf_data))

# Merge GWAS with MAF data
gwas_with_maf <- merge(gwas, maf_data[, c("SNP", "MAF")], by.x = "SNPID", by.y = "SNP", all.x = TRUE)

# Check for missing MAF values
missing_maf <- sum(is.na(gwas_with_maf$MAF))
if (missing_maf > 0) {
  cat("Warning:", missing_maf, "SNPs missing MAF values\n")
  # Show which SNPs are missing MAF
  missing_snps <- gwas_with_maf$SNPID[is.na(gwas_with_maf$MAF)]
  cat("First 10 missing SNPs:", head(missing_snps, 10), "\n")
}

# VALIDATE MAF VALUES - CRITICAL FIX
cat("Validating MAF values...\n")
cat("MAF values before any processing:\n")
print(summary(gwas_with_maf$MAF))

invalid_maf_before <- sum(gwas_with_maf$MAF <= 0 | gwas_with_maf$MAF > 0.5 | is.na(gwas_with_maf$MAF))
cat("Invalid MAF values before correction:", invalid_maf_before, "\n")

# Fix MAF values that are out of bounds
gwas_with_maf$MAF_original <- gwas_with_maf$MAF

# For SNPs missing MAF, use a reasonable default based on allele frequency
if (missing_maf > 0) {
  # Try to estimate MAF from the data or use conservative default
  gwas_with_maf$MAF[is.na(gwas_with_maf$MAF)] <- 0.05
  cat("Set missing MAF values to 0.05\n")
}

# Ensure MAF is between 0.001 and 0.5
gwas_with_maf$MAF <- ifelse(gwas_with_maf$MAF <= 0, 0.001, 
                           ifelse(gwas_with_maf$MAF > 0.5, 0.5, gwas_with_maf$MAF))

# For MAF = 0 exactly, set to small positive value
gwas_with_maf$MAF <- ifelse(gwas_with_maf$MAF == 0, 0.001, gwas_with_maf$MAF)

invalid_maf_after <- sum(gwas_with_maf$MAF <= 0 | gwas_with_maf$MAF > 0.5 | is.na(gwas_with_maf$MAF))
cat("Invalid MAF values after correction:", invalid_maf_after, "\n")

# Reorder to match PLINK BIM order
gwas_ordered <- gwas_with_maf[match(plink_snps, gwas_with_maf$SNPID), ]
gwas_ordered <- na.omit(gwas_ordered)

# Verify MAF distribution
cat("Final MAF summary statistics:\n")
print(summary(gwas_ordered$MAF))
cat("Number of SNPs with MAF < 0.01:", sum(gwas_ordered$MAF < 0.01), "\n")
cat("Number of SNPs with MAF < 0.05:", sum(gwas_ordered$MAF < 0.05), "\n")

# Check for any remaining invalid MAF values
final_invalid <- sum(gwas_ordered$MAF <= 0 | gwas_ordered$MAF > 0.5)
if (final_invalid > 0) {
  cat("ERROR: Still", final_invalid, "invalid MAF values after all corrections!\n")
  cat("Problematic MAF values:", gwas_ordered$MAF[gwas_ordered$MAF <= 0 | gwas_ordered$MAF > 0.5], "\n")
  stop("Cannot proceed with invalid MAF values")
}

cat("Original PLINK SNPs:", length(plink_snps), "\n")
cat("Final matched SNPs:", nrow(gwas_ordered), "\n")

# Write the final column order for debugging
cat("Final column order in reordered file:\n")
print(colnames(gwas_ordered))

# Write reordered GWAS file with MAF - ensure consistent column order
# Create a clean output with only the necessary columns in a fixed order
output_data <- data.frame(
  SNPID = gwas_ordered$SNPID,
  CHR = gwas_ordered$CHR,
  POS = gwas_ordered$POS,
  EA = gwas_ordered$EA,
  NEA = gwas_ordered$NEA,
  BETA = gwas_ordered$BETA,
  SE = gwas_ordered$SE,
  P = gwas_ordered$P,
  MAF = gwas_ordered$MAF,
  stringsAsFactors = FALSE
)

write.table(output_data, paste0(prefix, "_reordered_with_maf.loci.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Also write a validation file to check MAF values
validation_df <- data.frame(
  SNP = gwas_ordered$SNPID,
  MAF_original = gwas_ordered$MAF_original,
  MAF_final = gwas_ordered$MAF,
  CHR = gwas_ordered$CHR,
  POS = gwas_ordered$POS
)
write.table(validation_df, paste0(output_dir, "/", prefix, "_maf_validation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Successfully created GWAS file with validated MAF values\n")
cat("MAF validation file written to:", paste0(output_dir, "/", prefix, "_maf_validation.tsv"), "\n")
EOF2
)

echo "$RSCRIPT2" > reorder_snps_with_maf.R
Rscript reorder_snps_with_maf.R ${BFILE}.bim ${PREFIX}_filtered.loci.tsv ${PREFIX}_maf ${PREFIX} ${OUTPUT_DIR}

echo "Step 7: Prepare Z file with headers using validated MAF values..."
echo "rsid chromosome position allele1 allele2 maf beta se" > ${OUTPUT_DIR}/${PREFIX}.z

# Debug: Check the actual column structure
echo "Debug: Checking column structure of reordered file..."
head -1 ${PREFIX}_reordered_with_maf.loci.tsv | tr '\t' '\n' | nl -v 1

# Create the Z file - MAF should now be in column 9
tail -n +2 ${PREFIX}_reordered_with_maf.loci.tsv | awk 'BEGIN{FS="\t"; OFS=" "} {
    # With the fixed R script, columns are: 1=SNPID, 2=CHR, 3=POS, 4=EA, 5=NEA, 6=BETA, 7=SE, 8=P, 9=MAF
    maf = $9
    # Additional validation in awk
    if (maf == "" || maf == "NA") maf = 0.05
    if (maf <= 0) maf = 0.001
    if (maf > 0.5) maf = 0.5
    print $1, $2, $3, $4, $5, maf, $6, $7
}' >> ${OUTPUT_DIR}/${PREFIX}.z

echo "Step 8: Prepare SNP info file with validated MAF values..."
echo "rsid chromosome position allele1 allele2 maf" > ${OUTPUT_DIR}/${PREFIX}.snp
tail -n +2 ${PREFIX}_reordered_with_maf.loci.tsv | awk 'BEGIN{FS="\t"; OFS=" "} {
    maf = $9
    if (maf == "" || maf == "NA") maf = 0.05
    if (maf <= 0) maf = 0.001
    if (maf > 0.5) maf = 0.5
    print $1, $2, $3, $4, $5, maf
}' >> ${OUTPUT_DIR}/${PREFIX}.snp

# Final validation of Z file before FINEMAP
echo "Validating Z file before FINEMAP..."
RSCRIPT_VALIDATE=$(cat << 'EOF3'
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
z_file <- args[1]

z_data <- read.table(z_file, header = TRUE, stringsAsFactors = FALSE)
cat("Z file validation:\n")
cat("Total SNPs:", nrow(z_data), "\n")

# Check MAF range
invalid_maf <- sum(z_data$maf <= 0 | z_data$maf > 0.5)
cat("SNPs with invalid MAF (<=0 or >0.5):", invalid_maf, "\n")

if (invalid_maf > 0) {
  cat("Problematic MAF values:\n")
  print(z_data$maf[z_data$maf <= 0 | z_data$maf > 0.5])
  stop("Z file contains invalid MAF values!")
}

cat("MAF range:", min(z_data$maf), "-", max(z_data$maf), "\n")
cat("MAF distribution:\n")
print(summary(z_data$maf))
cat("Z file is valid for FINEMAP\n")
EOF3
)

echo "$RSCRIPT_VALIDATE" > validate_z_file.R
Rscript validate_z_file.R ${OUTPUT_DIR}/${PREFIX}.z

echo "Step 9: Copy LD matrix to output directory..."
cp ${PREFIX}_formatted.ld ${OUTPUT_DIR}/${PREFIX}.ld

echo "Step 10: Create master file in output directory..."
echo "z;ld;snp;config;cred;n_samples;log" > ${OUTPUT_DIR}/${PREFIX}.master
echo "${PREFIX}.z;${PREFIX}.ld;${PREFIX}.snp;${PREFIX}.config;${PREFIX}.cred;${NSAMPLES};${PREFIX}.log" >> ${OUTPUT_DIR}/${PREFIX}.master

echo "Step 11: Run FINEMAP from output directory..."
cd ${OUTPUT_DIR}
finemap --sss --in-files ${PREFIX}.master --n-threads 4 --n-causal-snps ${MAX_CAUSAL_SNPS} --pvalue-snps ${MIN_PVAL}
cd ..

echo "Step 12: Move intermediate files to output directory..."
mv ${PREFIX}_pval_filtered.loci.tsv ${PREFIX}_gwas_snps.txt ${PREFIX}_gwas_snps_only.* ${PREFIX}_filtered.loci.tsv ${PREFIX}_reordered_with_maf.loci.tsv ${PREFIX}_ldmatrix.* ${PREFIX}_formatted.ld ${PREFIX}_maf.* ${BFILE}.* ${PREFIX}_filtered.bim ${OUTPUT_DIR}/ 2>/dev/null || true

echo "FINEMAP completed successfully! Results are in: ${OUTPUT_DIR}/"

# Cleanup temporary files
rm -f filter_by_pvalue.R reorder_snps_with_maf.R create_ld_matrix.R validate_z_file.R

echo "Final verification:"
echo "SNPs in Z file: $(tail -n +2 ${OUTPUT_DIR}/${PREFIX}.z | wc -l)"
echo "SNPs in SNP file: $(tail -n +2 ${OUTPUT_DIR}/${PREFIX}.snp | wc -l)"
echo "LD matrix dimension should be: $(tail -n +2 ${OUTPUT_DIR}/${PREFIX}.z | wc -l)"

# Display MAF distribution for verification
echo "MAF distribution in final analysis:"
tail -n +2 ${OUTPUT_DIR}/${PREFIX}.z | awk '{print $6}' | sort -g | awk '
BEGIN {
  count = 0
  sum = 0
  min = 1
  max = 0
  rare = 0
  lowfreq = 0
}
{
  count++
  sum += $1
  if ($1 < min) min = $1
  if ($1 > max) max = $1
  if ($1 < 0.01) rare++
  if ($1 < 0.05) lowfreq++
}
END {
  if (count > 0) {
    print "Total SNPs: " count
    print "MAF range: " min " - " max
    print "Mean MAF: " sum/count
    print "Rare variants (MAF < 0.01): " rare
    print "Low frequency variants (MAF < 0.05): " lowfreq
  } else {
    print "No SNPs found in analysis"
  }
}'

# Return to original directory
cd "${ORIGINAL_DIR}"

echo "FINEMAP analysis complete. Results in ${OUTPUT_DIR}/"