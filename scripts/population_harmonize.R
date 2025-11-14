#!/usr/bin/env Rscript

# Title: population_harmonize.R
# Description: Population-based GWAS harmonization with reference .bim file using data.table

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  cat("Usage: Rscript population_harmonize.R <gwas_file> <reference_bim> <output_prefix>\n")
  cat("Example: Rscript population_harmonize.R 02_harmonized/hisp.harmonized.tsv EUR_subset.bim 03_popbased_harmonized/hisp_EUR\n")
  quit(status = 1)
}

gwas_file <- args[1]
reference_bim <- args[2]
output_prefix <- args[3]

cat("=== POPULATION-BASED HARMONIZATION ===\n")
cat("GWAS file:", gwas_file, "\n")
cat("Reference .bim:", reference_bim, "\n")
cat("Output prefix:", output_prefix, "\n")
cat("=======================================\n")

# Read GWAS data
cat("Reading GWAS data...\n")
if (!file.exists(gwas_file)) {
  cat("ERROR: GWAS file not found:", gwas_file, "\n")
  quit(status = 1)
}

# Check file size and provide progress estimate
file_info <- file.info(gwas_file)
file_size_mb <- round(file_info$size / 1024^2, 2)
cat("File size:", file_size_mb, "MB\n")
if (file_size_mb > 100) {
  cat("Large file detected - this may take several minutes...\n")
}

# Read with optimized fread parameters and progress monitoring
cat("Loading data with optimized settings...\n")
gwas_dt <- fread(gwas_file, verbose = TRUE, showProgress = TRUE, nThread = 2)
cat("GWAS variants loaded:", nrow(gwas_dt), "\n")
cat("Memory usage: ", round(object.size(gwas_dt) / 1024^2, 2), "MB\n")

# Read reference .bim file
cat("Reading reference .bim file...\n")
if (!file.exists(reference_bim)) {
  cat("ERROR: Reference .bim file not found:", reference_bim, "\n")
  quit(status = 1)
}
bim_dt <- fread(reference_bim, 
                col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
                verbose = FALSE)
cat("Reference variants:", nrow(bim_dt), "\n")

# Convert alleles to uppercase for consistency
gwas_dt[, `:=`(EA = toupper(EA), NEA = toupper(NEA))]
bim_dt[, `:=`(A1 = toupper(A1), A2 = toupper(A2))]

# Merge GWAS with reference
cat("Merging GWAS with reference...\n")
merged_dt <- merge(gwas_dt, bim_dt[, .(SNP, A1, A2)], by.x = "SNPID", by.y = "SNP", all.x = TRUE)

# Initialize harmonization status
merged_dt[, harmonization_status := "not_in_reference"]

# Filter to SNPs found in reference
ref_snps <- merged_dt[!is.na(A1)]
cat("SNPs found in reference:", nrow(ref_snps), "/", nrow(merged_dt), "\n")

if (nrow(ref_snps) == 0) {
  cat("ERROR: No SNPs found in reference!\n")
  quit(status = 1)
}

# Harmonization logic using data.table operations
cat("Performing harmonization...\n")

# Case 1: Already aligned (EA=A1, NEA=A2)
aligned_idx <- ref_snps[, EA == A1 & NEA == A2]
ref_snps[aligned_idx, harmonization_status := "aligned"]

# Case 2: Need swapping (EA=A2, NEA=A1)
swap_idx <- ref_snps[, EA == A2 & NEA == A1]
ref_snps[swap_idx, `:=`(
  EA = A1,
  NEA = A2,
  BETA = -BETA,
  harmonization_status = "swapped"
)]

# Case 3: Ambiguous SNPs (A/T or C/G pairs)
# Create helper function to check ambiguous pairs
is_ambiguous <- function(a1, a2) {
  paste0(a1, "_", a2) %in% c("A_T", "T_A", "C_G", "G_C")
}

ambiguous_idx <- ref_snps[, is_ambiguous(A1, A2)]
ref_snps[ambiguous_idx & harmonization_status == "not_in_reference", harmonization_status := "ambiguous"]

# Case 4: Filter to compatible alleles
compatible_idx <- ref_snps[, harmonization_status %in% c("aligned", "swapped", "ambiguous")]
final_dt <- ref_snps[compatible_idx]

# Remove temporary reference columns
final_dt[, `:=`(A1 = NULL, A2 = NULL)]

# Report harmonization results
cat("\n=== HARMONIZATION RESULTS ===\n")
status_counts <- final_dt[, .N, by = harmonization_status][order(-N)]
print(status_counts)

total_input <- nrow(merged_dt)
total_output <- nrow(final_dt)
cat("Total input variants:", total_input, "\n")
cat("Total output variants:", total_output, "\n")
cat("Removal rate:", round(100 * (1 - total_output/total_input), 2), "%\n")

# Save harmonized data
output_file <- paste0(output_prefix, ".popharmonized.tsv")
cat("Writing population harmonized data to:", output_file, "\n")
fwrite(final_dt, output_file, sep = "\t", verbose = FALSE)

# Create harmonization report
report_file <- paste0(output_prefix, ".popharmonization_report.txt")
cat("Writing harmonization report to:", report_file, "\n")

sink(report_file)
cat("Population-Based GWAS Harmonization Report\n")
cat("=========================================\n\n")
cat("Input GWAS file:", gwas_file, "\n")
cat("Reference .bim file:", reference_bim, "\n")
cat("Processing date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Input Statistics:\n")
cat("- Total GWAS variants:", total_input, "\n")
cat("- Reference variants:", nrow(bim_dt), "\n")
cat("- SNPs found in reference:", nrow(ref_snps), "\n\n")

cat("Harmonization Results:\n")
print(status_counts)
cat("\n")

cat("Output Statistics:\n")
cat("- Final harmonized variants:", total_output, "\n")
cat("- Variants removed:", total_input - total_output, "\n")
cat("- Removal rate:", round(100 * (1 - total_output/total_input), 2), "%\n\n")

cat("Harmonization Details:\n")
cat("- Aligned: Effect allele matches reference A1\n")
cat("- Swapped: Effect allele matches reference A2 (BETA sign flipped)\n")
cat("- Ambiguous: A/T or C/G pairs (kept with original orientation)\n\n")

cat("Files Created:\n")
cat("- Population harmonized GWAS:", output_file, "\n")
cat("- This report:", report_file, "\n")
sink()

cat("✅ SUCCESS: Population-based harmonization completed!\n")
cat("Harmonized variants:", nrow(final_dt), "\n")
cat("Output file:", output_file, "\n")
cat("Report file:", report_file, "\n")