#!/usr/bin/env Rscript
#==============================================================================
# FINEMAP Results Analyzer - Improved & Integrated Version
# Purpose: Analyze FINEMAP results with flexible modes
# Modes: 
#   - quick: Combine all k values (fast overview)
#   - best: Select and analyze best k based on config probabilities
# Date: 2025-11-14
#==============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default parameters
DEFAULT_LD_HIGH <- 0.8
DEFAULT_LD_LOW <- 0.2
DEFAULT_TOP_SNPS <- 5
DEFAULT_MODE <- "best"

# Help message
print_usage <- function() {
  cat("Usage: Rscript analyze_finemap_results.R <mode> <prefix> <gwas_file> <ld_file> [OPTIONS]\n\n")
  cat("REQUIRED PARAMETERS:\n")
  cat("  mode        Analysis mode: 'quick' (all k combined) or 'best' (best k selection)\n")
  cat("  prefix      Output prefix used in FINEMAP\n")
  cat("  gwas_file   Filtered GWAS summary statistics file\n")
  cat("  ld_file     LD matrix file\n\n")
  cat("OPTIONAL PARAMETERS:\n")
  cat("  --max-k INT           Maximum k value (required for 'best' mode)\n")
  cat("  --ld-high FLOAT       High LD threshold (default: 0.8)\n")
  cat("  --ld-low FLOAT        Low LD threshold (default: 0.2)\n")
  cat("  --top-snps INT        Number of top SNPs to report (default: 5)\n")
  cat("  --help                Show this help message\n\n")
  cat("EXAMPLES:\n")
  cat("  # Quick analysis (all k combined)\n")
  cat("  Rscript analyze_finemap_results.R quick PREFIX gwas.tsv ld.txt\n\n")
  cat("  # Best k selection\n")
  cat("  Rscript analyze_finemap_results.R best PREFIX gwas.tsv ld.txt --max-k 5\n\n")
  quit(status = 0)
}

# Check minimum arguments
if (length(args) < 4 || "--help" %in% args) {
  print_usage()
}

# Required parameters
mode <- tolower(args[1])
prefix <- args[2]
gwas_file <- args[3]
ld_file <- args[4]

# Validate mode
if (!mode %in% c("quick", "best")) {
  cat("Error: mode must be 'quick' or 'best', got:", mode, "\n")
  quit(status = 1)
}

# Parse optional arguments
max_k <- NULL
ld_high <- DEFAULT_LD_HIGH
ld_low <- DEFAULT_LD_LOW
top_snps <- DEFAULT_TOP_SNPS

i <- 5
while (i <= length(args)) {
  if (args[i] == "--max-k") {
    max_k <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--ld-high") {
    ld_high <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--ld-low") {
    ld_low <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--top-snps") {
    top_snps <- as.numeric(args[i + 1])
    i <- i + 2
  } else {
    cat("Warning: Unknown option:", args[i], "\n")
    i <- i + 1
  }
}

# Validate parameters
if (mode == "best" && is.null(max_k)) {
  cat("Error: --max-k is required for 'best' mode\n")
  quit(status = 1)
}

cat("=== FINEMAP RESULTS ANALYZER ===\n\n")
cat("Mode:", mode, "\n")
cat("Prefix:", prefix, "\n")
cat("GWAS file:", gwas_file, "\n")
cat("LD file:", ld_file, "\n")
cat("LD thresholds: high =", ld_high, ", low =", ld_low, "\n")
cat("Top SNPs to report:", top_snps, "\n")
if (!is.null(max_k)) cat("Max k:", max_k, "\n")
cat("\n")

# Validate files exist
if (!file.exists(gwas_file)) {
  stop("GWAS file not found: ", gwas_file)
}
if (!file.exists(ld_file)) {
  stop("LD file not found: ", ld_file)
}

#==============================================================================
# Utility Functions
#==============================================================================

# Function to read credible set and extract all SNPs
read_credible_snps <- function(cred_file, k) {
  if (!file.exists(cred_file)) {
    return(data.frame())
  }
  
  tryCatch({
    cred_data <- read.table(cred_file, header = TRUE, stringsAsFactors = FALSE)
    snp_list <- list()
    
    # Extract SNPs from all causal configurations
    for (i in 1:10) {
      cred_col <- paste0("cred", i)
      prob_col <- paste0("prob", i)
      
      if (cred_col %in% colnames(cred_data) && prob_col %in% colnames(cred_data)) {
        snps <- na.omit(cred_data[[cred_col]])
        probs <- na.omit(cred_data[[prob_col]])
        
        if (length(snps) > 0) {
          snp_list[[i]] <- data.frame(
            SNP = snps,
            PIP = probs,
            Causal_Config = i,
            Credible_Set = paste0("cs", k),
            stringsAsFactors = FALSE
          )
        }
      } else {
        break
      }
    }
    
    if (length(snp_list) == 0) {
      if ("rsid" %in% colnames(cred_data) && "prob" %in% colnames(cred_data)) {
        return(data.frame(
          SNP = cred_data$rsid,
          PIP = cred_data$prob,
          Causal_Config = 1,
          Credible_Set = paste0("cs", k),
          stringsAsFactors = FALSE
        ))
      } else {
        return(data.frame())
      }
    }
    
    do.call(rbind, snp_list)
  }, error = function(e) {
    cat("Error reading credible set file:", cred_file, "\n")
    cat("Error message:", e$message, "\n")
    return(data.frame())
  })
}

# Function to calculate LD with leading SNP
calculate_ld_with_leader <- function(snps, ld_matrix, snp_list, leader_idx) {
  ld_results <- numeric(length(snps))
  names(ld_results) <- snps
  
  for (i in 1:length(snps)) {
    snp <- snps[i]
    if (snp %in% snp_list) {
      snp_idx <- which(snp_list == snp)
      if (length(snp_idx) > 0 && snp_idx <= nrow(ld_matrix) && leader_idx <= nrow(ld_matrix)) {
        ld_r <- ld_matrix[leader_idx, snp_idx]
        ld_results[i] <- round(ld_r^2, 4)
      } else {
        ld_results[i] <- NA
      }
    } else {
      ld_results[i] <- NA
    }
  }
  
  ld_results
}

# Function to read config probabilities
read_config_probabilities <- function(prefix, max_k) {
  config_probs <- data.frame(K = integer(), Probability = numeric())
  
  config_file_to_read <- paste0(prefix, ".config")
  
  if (!file.exists(config_file_to_read)) {
    cat("Error: Config file not found:", config_file_to_read, "\n")
    return(config_probs)
  }
  
  cat("Reading configuration from:", config_file_to_read, "\n")
  
  # Read the config file
  config_lines <- readLines(config_file_to_read)
  config_lines <- config_lines[!grepl("^#", config_lines)]
  
  if (length(config_lines) <= 1) {
    cat("Warning: Config file is empty or has no data\n")
    return(config_probs)
  }
  
  # Parse config data
  config_data <- read.table(text = paste(config_lines, collapse = "\n"), 
                           header = TRUE, stringsAsFactors = FALSE)
  
  # Verify required columns exist
  if (!("k" %in% colnames(config_data)) || !("prob" %in% colnames(config_data))) {
    cat("Error: Config file missing 'k' or 'prob' columns\n")
    cat("Available columns:", paste(colnames(config_data), collapse=", "), "\n")
    return(config_probs)
  }
  
  cat("\nConfig file contains", nrow(config_data), "configurations\n")
  
  # For each k value from 1 to max_k, sum all probabilities for configurations with that k
  for (k in 1:max_k) {
    matching_rows <- which(config_data$k == k)
    
    if (length(matching_rows) > 0) {
      # Sum probabilities for all configurations with this k value
      total_prob <- sum(config_data$prob[matching_rows])
      config_probs <- rbind(config_probs, data.frame(
        K = k,
        Probability = total_prob
      ))
      cat("  k=", k, ": ", length(matching_rows), " configuration(s), total probability = ", 
          round(total_prob, 6), "\n", sep="")
    } else {
      # No configurations for this k value
      config_probs <- rbind(config_probs, data.frame(
        K = k,
        Probability = 0
      ))
      cat("  k=", k, ": No configurations (probability = 0)\n", sep="")
    }
  }
  
  cat("\n")
  return(config_probs)
}

#==============================================================================
# Main Analysis
#==============================================================================

tryCatch({
  # Read GWAS data
  gwas <- read.table(gwas_file, header = TRUE, stringsAsFactors = FALSE)
  ld_matrix <- as.matrix(read.table(ld_file, header = FALSE))
  snp_list <- gwas$SNPID
  
  # Find leading SNP
  pval_col <- c("P", "PVAL", "P_VALUE")[c("P", "PVAL", "P_VALUE") %in% colnames(gwas)][1]
  if (is.na(pval_col)) {
    stop("No p-value column found in GWAS file. Available columns: ", paste(colnames(gwas), collapse = ", "))
  }
  
  leader_idx <- which.min(gwas[[pval_col]])
  leading_snp <- gwas$SNPID[leader_idx]
  leading_pval <- gwas[[pval_col]][leader_idx]
  
  cat("Leading GWAS SNP:", leading_snp, "(p =", format(leading_pval, scientific = TRUE, digits = 3), ")\n\n")
  
  #============================================================================
  # MODE: QUICK ANALYSIS (All k combined)
  #============================================================================
  if (mode == "quick") {
    cat("=== QUICK ANALYSIS MODE (All k Combined) ===\n\n")
    
    # Read config file to choose best k
    config_file <- paste0(prefix, ".config")
    if (!file.exists(config_file)) {
      stop("Config file not found: ", config_file)
    }
    
    # Read config file with proper parsing
    config_lines <- readLines(config_file)
    config_lines <- config_lines[!grepl("^#", config_lines)]
    
    if (length(config_lines) <= 1) {
      stop("Config file is empty or has no data: ", config_file)
    }
    
    config <- read.table(text = paste(config_lines, collapse = "\n"), 
                        header = TRUE, stringsAsFactors = FALSE)
    
    if (nrow(config) == 0) {
      stop("No configurations found in: ", config_file)
    }
    
    best_k <- config$k[which.max(config$prob)]
    best_prob <- max(config$prob)
    
    cat("Configuration probabilities:\n")
    for (i in 1:nrow(config)) {
      cat("  k=", config$k[i], ": ", round(config$prob[i], 5), "\n", sep="")
    }
    cat("\nBest k (from config):", best_k, "(Probability:", round(best_prob, 3), ")\n\n")
    
    # Read all credible sets and combine
    all_credible_snps <- data.frame()
    for (k in 1:5) {
      cred_file <- paste0(prefix, ".cred", k)
      if (file.exists(cred_file)) {
        cred_snps <- read_credible_snps(cred_file, k)
        if (nrow(cred_snps) > 0) {
          cat("Found", nrow(cred_snps), "SNPs in credible set k =", k, "\n")
          all_credible_snps <- rbind(all_credible_snps, cred_snps)
        }
      }
    }
    
    if (nrow(all_credible_snps) == 0) {
      stop("No credible set SNPs found")
    }
    
    # Remove duplicates, keep highest PIP
    all_credible_snps <- all_credible_snps[order(all_credible_snps$PIP, decreasing = TRUE), ]
    all_credible_snps <- all_credible_snps[!duplicated(all_credible_snps$SNP), ]
    
    cat("\nTotal unique SNPs across all credible sets:", nrow(all_credible_snps), "\n")
    
    # Calculate LD with leading SNP
    ld_with_leader <- calculate_ld_with_leader(all_credible_snps$SNP, ld_matrix, snp_list, leader_idx)
    all_credible_snps$LD_with_Leader <- ld_with_leader
    
    # Merge with GWAS data
    merged_data <- merge(all_credible_snps, gwas, by.x = "SNP", by.y = "SNPID", all.x = TRUE)
    merged_data$is_leading <- merged_data$SNP == leading_snp
    merged_data$Relationship <- ifelse(
      merged_data$is_leading, "top",
      paste0("LD=", round(merged_data$LD_with_Leader, 3))
    )
    
    merged_data <- merged_data[order(merged_data$PIP, decreasing = TRUE), ]
    
    # Write results
    output_file <- paste0(prefix, "_final_results.tsv")
    write.table(merged_data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Results written to:", output_file, "\n")
    
    # Create summary
    top_n <- head(merged_data, min(3, nrow(merged_data)))
    summary_parts <- sapply(1:nrow(top_n), function(i) {
      paste0(top_n$SNP[i], "(PIP=", round(top_n$PIP[i], 2), ",", 
            top_n$Credible_Set[i], ",", top_n$Relationship[i], ")")
    })
    summary_string <- paste(summary_parts, collapse = ", ")
    
    # Write summary file
    summary_file <- paste0(prefix, "_summary.txt")
    sink(summary_file)
    cat("=== FINEMAP QUICK ANALYSIS SUMMARY ===\n\n")
    cat("Locus:", prefix, "\n")
    cat("Best k:", best_k, "(Probability:", round(best_prob, 3), ")\n")
    cat("Leading SNP:", leading_snp, "(p =", format(leading_pval, scientific = TRUE, digits = 3), ")\n\n")
    cat("TOP", min(top_snps, nrow(merged_data)), "SNPS:\n")
    for (i in 1:min(top_snps, nrow(merged_data))) {
      cat(i, ". ", merged_data$SNP[i], " (PIP=", round(merged_data$PIP[i], 3),
          ", ", merged_data$Credible_Set[i], ", ", merged_data$Relationship[i],
          ", p=", format(merged_data[[pval_col]][i], scientific = TRUE, digits = 2), ")\n", sep = "")
    }
    cat("\nSummary string:", summary_string, "\n")
    sink()
    
    cat("\nSummary written to:", summary_file, "\n")
    cat("Quick analysis complete!\n")
  }
  
  #============================================================================
  # MODE: BEST K ANALYSIS
  #============================================================================
  if (mode == "best") {
    cat("=== BEST K SELECTION MODE ===\n\n")
    
    # Read config probabilities
    config_probs <- read_config_probabilities(prefix, max_k)
    
    if (nrow(config_probs) == 0) {
      stop("No config files found")
    }
    
    cat("Configuration probabilities:\n")
    for (i in 1:nrow(config_probs)) {
      cat("  k=", config_probs$K[i], ": ", round(config_probs$Probability[i], 5), "\n", sep="")
    }
    cat("\n")
    
    best_k <- config_probs$K[which.max(config_probs$Probability)]
    best_prob <- max(config_probs$Probability)
    
    cat("\nBest k:", best_k, "(Probability:", round(best_prob, 3), ")\n\n")
    
    # Read ONLY the best k credible set
    best_cred_file <- paste0(prefix, ".cred", best_k)
    cat("Reading best credible set from:", best_cred_file, "\n")
    
    # Read credible sets for best k
    best_credible_snps <- data.frame()
    if (file.exists(best_cred_file)) {
      cred_lines <- readLines(best_cred_file)
      data_lines <- cred_lines[!grepl("^#", cred_lines) | grepl("^index", cred_lines)]
      
      if (length(data_lines) > 1) {
        cred_data <- read.table(text = paste(data_lines, collapse = "\n"), 
                               header = TRUE, stringsAsFactors = FALSE)
        
        for (config in 1:best_k) {
          cred_col <- paste0("cred", config)
          prob_col <- paste0("prob", config)
          
          if (cred_col %in% colnames(cred_data) && prob_col %in% colnames(cred_data)) {
            valid_indices <- !is.na(cred_data[[cred_col]]) & !is.na(cred_data[[prob_col]])
            valid_snps <- cred_data[[cred_col]][valid_indices]
            valid_probs <- cred_data[[prob_col]][valid_indices]
            
            if (length(valid_snps) > 0) {
              best_credible_snps <- rbind(best_credible_snps, data.frame(
                SNP = valid_snps,
                PIP = valid_probs,
                Causal_Config = config,
                Credible_Set = paste0("cs", config),
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
    
    if (nrow(best_credible_snps) == 0) {
      stop("No SNPs found in best credible set")
    }
    
    cat("SNPs in best k =", best_k, ":", nrow(best_credible_snps), "\n\n")
    
    # Calculate LD with leading SNP using the main LD matrix
    if (exists("ld_matrix") && exists("snp_list") && exists("leader_idx")) {
      best_credible_snps$LD_with_Leader <- calculate_ld_with_leader(
        best_credible_snps$SNP, ld_matrix, snp_list, leader_idx
      )
    } else {
      cat("Warning: Main LD matrix or SNP list not available for LD calculation.\n")
      best_credible_snps$LD_with_Leader <- NA
    }
    
    # Add relationship information
    best_credible_snps$is_leading <- best_credible_snps$SNP == leading_snp
    best_credible_snps$Relationship <- ifelse(
      best_credible_snps$is_leading, "top",
      paste0("LD=", round(best_credible_snps$LD_with_Leader, 3))
    )
    
    # Merge with GWAS data
    final_results <- merge(best_credible_snps, gwas, by.x = "SNP", by.y = "SNPID", all.x = TRUE)
    final_results <- final_results[order(final_results$PIP, decreasing = TRUE), ]
    
    # Write results
    output_file <- paste0(prefix, "_best_credible_results.tsv")
    write.table(final_results, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Results written to:", output_file, "\n")
    
    # Get top SNP from each credible set
    top_snps_by_cs <- data.frame()
    for (cs in 1:best_k) {
      cs_snps <- best_credible_snps[best_credible_snps$Causal_Config == cs, ]
      if (nrow(cs_snps) > 0) {
        top_cs_snp <- cs_snps[which.max(cs_snps$PIP), ]
        top_snps_by_cs <- rbind(top_snps_by_cs, top_cs_snp)
      }
    }
    
    top_snps_by_cs <- merge(top_snps_by_cs, final_results[, c("SNP", pval_col)], 
                           by = "SNP", all.x = TRUE)
    top_snps_by_cs <- top_snps_by_cs[order(top_snps_by_cs$Causal_Config), ]
    
    # Create summary string
    summary_parts <- sapply(1:nrow(top_snps_by_cs), function(i) {
      paste0(top_snps_by_cs$SNP[i], "(PIP=", round(top_snps_by_cs$PIP[i], 2), ",", 
            top_snps_by_cs$Credible_Set[i], ",", top_snps_by_cs$Relationship[i], ")")
    })
    summary_string <- paste(summary_parts, collapse = ", ")
    
    # Write summary file
    summary_file <- paste0(prefix, "_best_credible_analysis.txt")
    sink(summary_file)
    cat("=== FINEMAP BEST K ANALYSIS SUMMARY ===\n\n")
    cat("Locus:", prefix, "\n")
    cat("Best k:", best_k, "(Probability:", round(best_prob, 3), ")\n")
    cat("Leading SNP:", leading_snp, "(p =", format(leading_pval, scientific = TRUE, digits = 3), ")\n\n")
    cat("TOP SNPS FROM EACH CREDIBLE SET (k =", best_k, "):\n")
    for (i in 1:nrow(top_snps_by_cs)) {
      cat(i, ". ", top_snps_by_cs$SNP[i], " (PIP=", round(top_snps_by_cs$PIP[i], 3),
          ", ", top_snps_by_cs$Credible_Set[i], ", ", top_snps_by_cs$Relationship[i],
          ", p=", format(top_snps_by_cs[[pval_col]][i], scientific = TRUE, digits = 2), ")\n", sep = "")
    }
    cat("\nSummary string:", summary_string, "\n")
    sink()
    
    # Write best_top.txt
    best_top_file <- paste0(prefix, "_best_top.txt")
    writeLines(summary_string, best_top_file)
    
    cat("\nSummary written to:", summary_file, "\n")
    cat("Best top written to:", best_top_file, "\n")
    cat("Best k analysis complete!\n")
  }
  
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  quit(status = 1)
})
