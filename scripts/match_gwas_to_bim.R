#!/usr/bin/env Rscript

# Simple GWAS to BIM matching script
# Usage: Rscript match_gwas_to_bim.R <bim_file> <gwas_file> <output_file>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    cat("Usage: Rscript match_gwas_to_bim.R <bim_file> <gwas_file> <output_file>\n")
    cat("Arguments:\n")
    cat("  bim_file   - PLINK BIM file path\n")
    cat("  gwas_file  - GWAS summary statistics file (TSV)\n")
    cat("  output_file - Output matched GWAS file path\n")
    quit(status = 1)
}

bim_file <- args[1]
gwas_file <- args[2]
output_file <- args[3]

# Validate input files exist
if (!file.exists(bim_file)) {
    cat("ERROR: BIM file not found:", bim_file, "\n")
    quit(status = 1)
}

if (!file.exists(gwas_file)) {
    cat("ERROR: GWAS file not found:", gwas_file, "\n")
    quit(status = 1)
}

# Read BIM file to get PLINK SNP order
bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")
plink_snps <- bim$SNP

# Read GWAS summary statistics
gwas <- read.table(gwas_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Check if SNPID column exists
if (!"SNPID" %in% colnames(gwas)) {
    cat("ERROR: SNPID column not found in GWAS file\n")
    quit(status = 1)
}

# Match and reorder GWAS to match PLINK BIM order
gwas_matched <- gwas[match(plink_snps, gwas$SNPID), ]

# Remove rows with NA (SNPs not found in GWAS)
gwas_matched <- gwas_matched[!is.na(gwas_matched$SNPID), ]

# Write matched GWAS file
write.table(gwas_matched, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Matched", nrow(gwas_matched), "SNPs from", length(plink_snps), "BIM SNPs\n")