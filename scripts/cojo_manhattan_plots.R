#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for argument: %s", key))
    }
    value <- args[[i + 1]]
    parsed[[substring(key, 3)]] <- value
    i <- i + 2
  }
  parsed
}

coerce_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

save_plot <- function(plot, file_path) {
  ggsave(filename = file_path, plot = plot, width = 10, height = 5, dpi = 300)
}

save_plot_pdf <- function(plot, file_path) {
  ggsave(filename = file_path, plot = plot, width = 10, height = 5, device = "pdf")
}

create_placeholder_plot <- function(message_text, x_range = NULL, line_values = NULL, line_colors = NULL) {
  max_y <- if (is.null(line_values)) 1 else max(c(line_values, 1), na.rm = TRUE) + 1

  plot_obj <- ggplot() +
    annotate("text", x = mean(if (is.null(x_range)) c(0, 1) else x_range), y = max_y * 0.6,
             label = message_text, size = 4) +
    ylim(0, max_y) +
    labs(x = "Genomic position", y = "-log10(p)") +
    theme_minimal()

  if (!is.null(x_range) && length(x_range) == 2 && all(is.finite(x_range))) {
    plot_obj <- plot_obj +
      scale_x_continuous(limits = x_range, expand = c(0, 0))
  } else {
    plot_obj <- plot_obj + xlim(0, 1)
  }

  if (!is.null(line_values)) {
    for (i in seq_along(line_values)) {
      color_i <- if (!is.null(line_colors) && length(line_colors) >= i) line_colors[i] else "#444444"
      plot_obj <- plot_obj +
        geom_hline(yintercept = line_values[i], color = color_i, linetype = "dashed")
    }
  }

  plot_obj
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_keys <- c("target-analysis", "cojo-dir", "loci-dir", "output-dir", "significance-threshold", "genes-file")
missing_keys <- setdiff(required_keys, names(args))
if (length(missing_keys) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing_keys, collapse = ", ")))
}

target_analysis <- args[["target-analysis"]]
cojo_dir <- args[["cojo-dir"]]
loci_dir <- args[["loci-dir"]]
output_dir <- args[["output-dir"]]
significance_threshold <- coerce_numeric(args[["significance-threshold"]])
genes_file <- args[["genes-file"]]

if (is.null(genes_file) || genes_file == "") {
  stop("genes-file argument is required")
}

if (is.na(significance_threshold) || significance_threshold <= 0 || significance_threshold >= 1) {
  stop("significance-threshold must be a numeric value between 0 and 1")
}

if (!dir.exists(cojo_dir)) {
  stop(sprintf("COJO directory not found: %s", cojo_dir))
}

if (!dir.exists(loci_dir)) {
  stop(sprintf("Loci directory not found: %s", loci_dir))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(genes_file)) {
  fallback_candidates <- c(
    file.path(dirname(genes_file), "NCBI37.3.gene.loc"),
    "resources/NCBI37.3.gene.loc",
    "resources/NCBI37.3.gene.txt"
  )
  fallback_existing <- fallback_candidates[file.exists(fallback_candidates)]
  if (length(fallback_existing) == 0) {
    stop(sprintf("Genes file not found: %s (checked fallbacks: %s)",
                 genes_file,
                 paste(fallback_candidates, collapse = ", ")))
  }
  message(sprintf("Using fallback genes file: %s", fallback_existing[1]))
  genes_file <- fallback_existing[1]
}

genes_dt <- fread(genes_file, col.names = c("gene_id", "chr", "start", "end", "strand", "symbol"))
genes_dt[, chr_numeric := coerce_numeric(chr)]

build_gene_plot <- function(chr_value, region_start, region_end) {
  if (is.na(chr_value) || is.na(region_start) || is.na(region_end)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No region info", size = 3) +
             xlim(0, 1) + ylim(0, 1) +
             labs(x = "Genomic position", y = NULL) +
             theme_void())
  }

  region_genes <- genes_dt[chr_numeric == chr_value & end >= region_start & start <= region_end]

  base_plot <- ggplot() +
    scale_x_continuous(limits = c(region_start, region_end), expand = c(0, 0)) +
    labs(x = "Genomic position", y = NULL) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 2, r = 10, b = 2, l = 10)
    )

  if (nrow(region_genes) == 0) {
    return(base_plot +
             annotate("text", x = (region_start + region_end) / 2, y = 0.5,
                      label = "No annotated genes in region", size = 3))
  }

  setorder(region_genes, start, end)
  region_genes[, level := 0L]
  track_ends <- numeric()
  for (i in seq_len(nrow(region_genes))) {
    placed <- FALSE
    if (length(track_ends) == 0) {
      track_ends <- c(region_genes$end[i])
      region_genes$level[i] <- 1L
      next
    }
    for (track_index in seq_along(track_ends)) {
      if (region_genes$start[i] > (track_ends[track_index] + 1e3)) {
        region_genes$level[i] <- track_index
        track_ends[track_index] <- region_genes$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      track_ends <- c(track_ends, region_genes$end[i])
      region_genes$level[i] <- length(track_ends)
    }
  }

  region_genes[, mid := (start + end) / 2]

  max_level <- max(region_genes$level)

  base_plot +
    geom_segment(
      data = region_genes,
      aes(x = start, xend = end, y = level, yend = level,
          color = strand),
      lineend = "round",
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    geom_point(
      data = region_genes,
      aes(x = start, y = level),
      size = 1,
      color = "#636363"
    ) +
    geom_text(
      data = region_genes,
      aes(x = mid, y = level + 0.25, label = symbol),
      size = 2,
      fontface = "bold",
      color = "#1a1a1a"
    ) +
    scale_color_manual(values = c("+" = "#2C7BB6", "-" = "#D7191C"), na.translate = TRUE) +
    scale_y_continuous(limits = c(0.5, max_level + 0.75), expand = expansion(mult = c(0.05, 0.15)))
}

locus_dirs <- list.dirs(cojo_dir, full.names = TRUE, recursive = FALSE)
locus_dirs <- locus_dirs[startsWith(basename(locus_dirs), "LOC_")]

if (length(locus_dirs) == 0) {
  warning(sprintf("No locus directories found in %s", cojo_dir))
  quit(status = 0)
}

manifest_entries <- list()
genome_line <- -log10(5e-8)
suggestive_line <- -log10(5e-5)
line_values <- c(genome_line, suggestive_line)
line_colors <- c("#2C7BB6", "#D7191C")

for (locus_path in sort(locus_dirs)) {
  locus_name <- basename(locus_path)
  message(sprintf("Processing locus %s", locus_name))

  cojo_results_dir <- file.path(locus_path, "cojo_results")
  if (!dir.exists(cojo_results_dir)) {
    warning(sprintf("Skipping %s: cojo_results directory not found", locus_name))
    next
  }

  gwas_file <- file.path(loci_dir, locus_name, sprintf("%s_gwas_ref_match.tsv", locus_name))
  if (!file.exists(gwas_file)) {
    gwas_file <- file.path(loci_dir, locus_name, sprintf("%s.loci.tsv", locus_name))
  }

  if (!file.exists(gwas_file)) {
    warning(sprintf("Skipping %s: GWAS data file not found", locus_name))
    next
  }

  gwas_dt <- tryCatch(
    fread(gwas_file, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(gwas_dt)) {
    warning(sprintf("Skipping %s: failed to read GWAS data", locus_name))
    next
  }

  if (!all(c("POS", "P") %in% names(gwas_dt))) {
    warning(sprintf("Skipping %s: GWAS data missing POS/P columns", locus_name))
    next
  }

  gwas_dt[, POS := coerce_numeric(POS)]
  gwas_dt[, P := coerce_numeric(P)]
  before_dt <- gwas_dt[is.finite(POS) & is.finite(P) & P > 0 & P <= 1, .(pos = POS, logp = -log10(P))]

  if (nrow(before_dt) == 0) {
    warning(sprintf("Skipping %s: GWAS dataset has no valid SNPs", locus_name))
    next
  }

  cma_files <- list.files(cojo_results_dir, pattern = "^iteration_[0-9]+\\.cma\\.cojo$", full.names = TRUE)
  if (length(cma_files) == 0) {
    warning(sprintf("Skipping %s: no iteration files found", locus_name))
    next
  }

  for (cma_file in sort(cma_files)) {
    iter_id <- sub("^iteration_([0-9]+)\\.cma\\.cojo$", "\\1", basename(cma_file))
    given_file <- file.path(cojo_results_dir, sprintf("iteration_%s.given.cojo", iter_id))

    if (!file.exists(given_file)) {
      warning(sprintf("Skipping %s iteration %s: missing given file", locus_name, iter_id))
      next
    }

    cma_dt <- tryCatch(fread(cma_file, showProgress = FALSE), error = function(e) NULL)
    given_dt <- tryCatch(fread(given_file, showProgress = FALSE), error = function(e) NULL)

    if (is.null(cma_dt) || is.null(given_dt)) {
      warning(sprintf("Skipping %s iteration %s: failed to read iteration files", locus_name, iter_id))
      next
    }

    if (!all(c("bp", "p", "pC") %in% names(cma_dt))) {
      warning(sprintf("Skipping %s iteration %s: expected columns missing", locus_name, iter_id))
      next
    }

    combined_dt <- rbindlist(list(
      cma_dt[, source := "cma"],
      given_dt[, source := "given"]
    ), use.names = TRUE, fill = TRUE)

    combined_dt[, bp := coerce_numeric(bp)]
    combined_dt[, p := coerce_numeric(p)]
    combined_dt[, pC := coerce_numeric(pC)]

    after_dt <- combined_dt[is.finite(bp) & is.finite(pC) & pC > 0 & pC <= 1, .(pos = bp, logp = -log10(pC))]

    chr_candidates <- unique(coerce_numeric(gwas_dt$CHR))
    chr_candidates <- chr_candidates[is.finite(chr_candidates)]
    chr_value <- if (length(chr_candidates) > 0) chr_candidates[1] else NA_real_

    region_start <- suppressWarnings(min(before_dt$pos, na.rm = TRUE))
    region_end <- suppressWarnings(max(before_dt$pos, na.rm = TRUE))
    if (nrow(after_dt) > 0) {
      region_start <- min(region_start, suppressWarnings(min(after_dt$pos, na.rm = TRUE)))
      region_end <- max(region_end, suppressWarnings(max(after_dt$pos, na.rm = TRUE)))
    }

    if (!is.finite(region_start) || !is.finite(region_end) || region_start >= region_end) {
      region_start <- suppressWarnings(min(coerce_numeric(gwas_dt$POS), na.rm = TRUE))
      region_end <- suppressWarnings(max(coerce_numeric(gwas_dt$POS), na.rm = TRUE))
    }

    if (!is.finite(region_start) || !is.finite(region_end) || region_start >= region_end) {
      region_start <- 0
      region_end <- 1
    }

    before_gene_plot <- build_gene_plot(chr_value, region_start, region_end)
    after_gene_plot <- build_gene_plot(chr_value, region_start, region_end)

    output_iteration_dir <- file.path(output_dir, locus_name, sprintf("iteration_%s", iter_id))
    dir.create(output_iteration_dir, recursive = TRUE, showWarnings = FALSE)

    before_plot <- ggplot(before_dt, aes(x = pos, y = logp)) +
      geom_point(color = "#2C7BB6", size = 0.9) +
      geom_hline(yintercept = genome_line, color = "#2C7BB6", linetype = "dashed") +
      geom_hline(yintercept = suggestive_line, color = "#D7191C", linetype = "dashed") +
      labs(title = sprintf("%s: Before COJO (iteration %s)", locus_name, iter_id),
           x = "Genomic position", y = "-log10(P)") +
      scale_x_continuous(limits = c(region_start, region_end), expand = c(0.01, 0.01)) +
      theme_bw()

    if (nrow(after_dt) == 0) {
      warning(sprintf("%s iteration %s: no conditional p-values available; writing placeholder", locus_name, iter_id))
      after_plot <- create_placeholder_plot("Conditional results not available",
                                           c(region_start, region_end),
                                           line_values, line_colors)
    } else {
      after_plot <- ggplot(after_dt, aes(x = pos, y = logp)) +
        geom_point(color = "#1A9641", size = 0.9) +
        geom_hline(yintercept = genome_line, color = "#2C7BB6", linetype = "dashed") +
        geom_hline(yintercept = suggestive_line, color = "#D7191C", linetype = "dashed") +
        labs(title = sprintf("%s: After COJO (iteration %s)", locus_name, iter_id),
             x = "Genomic position", y = "-log10(pC)") +
        scale_x_continuous(limits = c(region_start, region_end), expand = c(0.01, 0.01)) +
        theme_bw()
    }

    conditioned_snps <- combined_dt[source == "given" & !is.na(SNP)]
    conditioned_snps <- unique(conditioned_snps$SNP)
    conditioned_label <- if (length(conditioned_snps) > 0) {
      paste("Conditioned SNPs:", paste(conditioned_snps, collapse = ", "))
    } else {
      "Conditioned SNPs: none"
    }

    combined_plot <- (before_plot / before_gene_plot) | (after_plot / after_gene_plot)
    combined_plot <- combined_plot +
      plot_layout(widths = c(1, 1), heights = c(3, 1), guides = "collect") +
      plot_annotation(
        title = sprintf("%s Iteration %s", locus_name, iter_id),
        subtitle = conditioned_label
      )

    combined_png <- file.path(output_iteration_dir, "before_after_with_genes.png")
    combined_pdf <- file.path(output_iteration_dir, "before_after_with_genes.pdf")
    save_plot(combined_plot, combined_png)
    save_plot_pdf(combined_plot, combined_pdf)

    manifest_entries[[length(manifest_entries) + 1]] <- data.table(
      target_analysis = target_analysis,
      locus = locus_name,
      iteration = iter_id,
      combined_png = combined_png
    )
  }
}

if (length(manifest_entries) > 0) {
  manifest <- rbindlist(manifest_entries, use.names = TRUE, fill = TRUE)
  fwrite(manifest, file.path(output_dir, "cojo_manhattan_manifest.tsv"), sep = "\t")
} else {
  warning("No plots were generated; manifest not written")
}
