#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
  library(scales)
})

parser <- ArgumentParser(description = "Plot COJO regional dependency scree lines")
parser$add_argument("--summary", required = TRUE, help = "TSV produced by cojo_region_dependency.py")
parser$add_argument("--output-prefix", required = TRUE, help = "Prefix for output files (PNG/PDF)")
parser$add_argument("--title", default = "COJO regional dependency scree plot", help = "Plot title")
parser$add_argument("--significance-threshold", type = "double", default = 5e-8, help = "Genome-wide threshold used to judge independence (default: 5e-8)")
args <- parser$parse_args()

data <- fread(args$summary, sep = "\t", na.strings = c("NA", ""))
if (!"mode" %in% names(data)) {
  stop("Summary file missing 'mode' column")
}

numeric_cols <- c("step", "total_conditioned", "top_other_p", "top_other_neg_logp")
for (col in numeric_cols) {
  if (col %in% names(data)) {
    data[[col]] <- as.numeric(data[[col]])
  }
}

if (!"conditioning_region" %in% names(data) || !"other_region" %in% names(data)) {
  stop("Summary file missing conditioning_region or other_region columns")
}

data[, panel := paste0(conditioning_region, " conditioned -> ", other_region)]

seq_data <- data[mode %in% c("baseline", "sequential")]
if (nrow(seq_data) == 0) {
  stop("No sequential entries found in summary file")
}

seq_data <- seq_data[order(panel, step)]
seq_data[, total_conditioned := fcoalesce(total_conditioned, step)]
seq_data[, color_key := conditioning_region]

palette_vec <- c("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b")
region_levels <- unique(seq_data$conditioning_region)
palette_map <- setNames(rep_len(palette_vec, length(region_levels)), region_levels)

threshold_value <- args$significance_threshold
genome_line <- -log10(threshold_value)
if (!is.finite(genome_line)) {
  genome_line <- NA_real_
}
plot_obj <- ggplot(seq_data, aes(x = total_conditioned, y = top_other_neg_logp, color = color_key, group = panel)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = genome_line, linetype = "dashed", color = "#7f7f7f", linewidth = 0.3) +
  facet_wrap(~panel, scales = "free_y") +
  scale_color_manual(values = palette_map, guide = "none") +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    title = args$title,
    x = "Number of conditioned SNPs",
    y = expression(-log[10](italic(p))~"of top SNP in other region")
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

all_data <- data[mode == "all_at_once"]
if (nrow(all_data) > 0) {
  all_data[, panel := paste0(conditioning_region, " conditioned -> ", other_region)]
  plot_obj <- plot_obj +
    geom_point(
      data = all_data,
      aes(x = total_conditioned, y = top_other_neg_logp),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      color = "#111111",
      size = 3
    )
}

png_file <- paste0(args$output_prefix, ".png")
pdf_file <- paste0(args$output_prefix, ".pdf")

ggsave(png_file, plot = plot_obj, width = 10, height = 5, dpi = 300)
ggsave(pdf_file, plot = plot_obj, width = 10, height = 5)

direction_summary <- data[, {
  baseline_rows <- .SD[mode == "baseline"]
  baseline_p <- if (nrow(baseline_rows) > 0) baseline_rows$top_other_p[1] else NA_real_
  baseline_neg <- if (nrow(baseline_rows) > 0) baseline_rows$top_other_neg_logp[1] else NA_real_

  final_rows <- .SD[mode == "all_at_once"][order(-total_conditioned)]
  final_mode <- "all_at_once"
  if (nrow(final_rows) == 0) {
    final_rows <- .SD[mode == "sequential"][order(-total_conditioned)]
    final_mode <- "sequential"
  }

  final_p <- if (nrow(final_rows) > 0) final_rows$top_other_p[1] else NA_real_
  final_neg <- if (nrow(final_rows) > 0) final_rows$top_other_neg_logp[1] else NA_real_
  final_count <- if (nrow(final_rows) > 0) final_rows$total_conditioned[1] else NA_real_

  p_to_check <- final_p
  if (is.na(p_to_check) || !is.finite(p_to_check)) {
    p_to_check <- baseline_p
  }

  status <- "independent"
  if (!is.na(p_to_check) && p_to_check < threshold_value) {
    status <- "dependent"
  }

  list(
    baseline_p = baseline_p,
    baseline_neg_logp = baseline_neg,
    final_p = final_p,
    final_neg_logp = final_neg,
    conditioned_snp_count = final_count,
    final_mode = if (nrow(final_rows) > 0) final_mode else NA_character_,
    independence_status = status
  )
}, by = .(conditioning_region, other_region)]

summary_file <- paste0(args$output_prefix, "_independence_summary.tsv")
fwrite(direction_summary, file = summary_file, sep = "\t", quote = FALSE, na = "NA")

cat(sprintf("Saved plots to %s and %s\n", png_file, pdf_file))
cat(sprintf("Wrote independence summary to %s\n", summary_file))
if (nrow(direction_summary) > 0) {
  for (idx in seq_len(nrow(direction_summary))) {
    row <- direction_summary[idx]
    final_label <- if (!is.na(row$final_mode)) row$final_mode else "none"
    final_p_text <- if (!is.na(row$final_p)) sprintf("%.3e", row$final_p) else "NA"
    cat(
      sprintf(
        "%s conditioned -> %s: %s (final mode = %s, final p = %s)\n",
        row$conditioning_region,
        row$other_region,
        row$independence_status,
        final_label,
        final_p_text
      )
    )
  }
}
