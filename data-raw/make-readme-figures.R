# Generate every figure used by README.md.
#
# Run from the package root after `devtools::load_all()`:
#   Rscript data-raw/make-readme-figures.R
#
# Figures are written to man/figures/. Nothing here is part of the installed
# package, so `patchwork` is only needed to rebuild the documentation.

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

ext <- function(...) system.file("extdata", ..., package = "metaplot")
fig_dir <- "man/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(plot, file, width = 7, height = 4.2, dpi = 150) {
  ggsave(file.path(fig_dir, file), plot, width = width, height = height,
         dpi = dpi, bg = "white")
  message("  -> ", file)
}

# ==============================================================================
# annotation (cached: prepare_features() on the bundled GTF takes a while)
# ==============================================================================
cache <- "data-raw/cache/features_anno.rds"
dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)

gtf_file <- ext("hg38.ncbiRefSeq.gtf.gz")

if (file.exists(cache)) {
  features_anno <- readRDS(cache)
} else {
  features_anno <- prepare_features(gtf_file, features_name = c("5UTR", "CDS", "3UTR"))
  saveRDS(features_anno, cache)
}

m6a_bed    <- ext("m6a_peaks.bed")
m1a_bed    <- ext("m1a.sorted2.bed")
helako_bed <- ext("helaA30ko-1.m1a.sorted2.bed")
m6a_fa     <- ext("m6a_peaks.fa")

# ==============================================================================
# 1. prepare_features(): what the annotation looks like
# ==============================================================================
feat_df <- as.data.frame(features_anno)

feat_summary <- feat_df %>%
  group_by(type) %>%
  summarise(regions = n(),
            transcripts = n_distinct(transcript_id),
            median_len = median(f_len),
            .groups = "drop") %>%
  mutate(type = factor(type, levels = c("UTR5", "CDS", "UTR3")))

print(feat_summary)
cat("total regions:", nrow(feat_df), "\n")
cat("total transcripts:", n_distinct(feat_df$transcript_id), "\n")
cat("mcols:", paste(names(S4Vectors::mcols(features_anno)), collapse = ", "), "\n")

p_counts <- ggplot(feat_summary, aes(type, regions, fill = type)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(regions)), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(UTR5 = "#4C72B0", CDS = "#DD8452", UTR3 = "#55A868")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = scales::comma) +
  labs(x = NULL, y = "Annotated regions", title = "Region count") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

p_len <- ggplot(feat_summary, aes(type, median_len / 1000, fill = type)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = round(median_len / 1000, 2)), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c(UTR5 = "#4C72B0", CDS = "#DD8452", UTR3 = "#55A868")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Median length (kb)", title = "Median region length") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

save_fig(p_counts + p_len, "features-overview.png", width = 7.5, height = 3.6)

# ==============================================================================
# 2. calculate_relative_position(): the data behind the plot
# ==============================================================================
rel_df <- calculate_relative_position(m6a_bed, features_anno, scale_region = TRUE)

cat("\ncalculate_relative_position() output:\n")
print(utils::head(rel_df, 4))
cat("rows:", nrow(rel_df), " rel_pos range:",
    paste(round(range(rel_df$rel_pos), 3), collapse = " ~ "), "\n")

rel_plot_df <- rel_df %>%
  mutate(type = factor(type, levels = c("UTR5", "CDS", "UTR3")))

p_rel <- ggplot(rel_plot_df, aes(rel_pos, fill = type, colour = type)) +
  geom_histogram(aes(y = after_stat(density)), bins = 120,
                 position = "identity", alpha = 0.45, linewidth = 0.2) +
  geom_vline(xintercept = c(1, 2), linetype = "dashed", colour = "grey30") +
  annotate("text", x = c(0.5, 1.5, 2.5), y = Inf, label = c("5'UTR", "CDS", "3'UTR"),
           vjust = 1.6, fontface = "bold", size = 3.4) +
  scale_fill_manual(values = c(UTR5 = "#4C72B0", CDS = "#DD8452", UTR3 = "#55A868")) +
  scale_colour_manual(values = c(UTR5 = "#4C72B0", CDS = "#DD8452", UTR3 = "#55A868")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Relative position (scaled)", y = "Density",
       title = "Peak midpoints mapped onto the 0-1 / 1-2 / 2-3 meta-gene") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

save_fig(p_rel, "relative-position.png", width = 7, height = 3.8)

# ==============================================================================
# 3. metaPlot(): every plotting option
# ==============================================================================
# 3.1 single sample
res_single <- metaPlot(bed_file = list(m6a_bed), group_names = "m6A",
                       features_anno = features_anno, scale_region = TRUE)
save_fig(res_single$plot, "metaplot-single-sample.png", width = 6, height = 4)

# 3.2 scale_region on/off
res_unscaled <- metaPlot(bed_file = list(m6a_bed), group_names = "m6A",
                         features_anno = features_anno, scale_region = FALSE)

p_scale <- (res_unscaled$plot + labs(subtitle = "scale_region = FALSE")) +
  (res_single$plot + labs(subtitle = "scale_region = TRUE")) +
  plot_annotation(tag_levels = "A") &
  theme(plot.subtitle = element_text(size = 9, colour = "grey30"))
save_fig(p_scale, "metaplot-scale-region.png", width = 9.5, height = 3.8)

# 3.3 cut_ratio
res_auto <- metaPlot(bed_file = list(m6a_bed), group_names = "m6A",
                     features_anno = features_anno, scale_region = TRUE)
res_cut  <- metaPlot(bed_file = list(m6a_bed), group_names = "m6A",
                     features_anno = features_anno, scale_region = TRUE,
                     cut_ratio = c(0.1, 0.8, 0.1))

p_cut <- (res_auto$plot + labs(subtitle = "estimated from the annotation")) +
  (res_cut$plot + labs(subtitle = "cut_ratio = c(0.1, 0.8, 0.1)")) +
  plot_annotation(tag_levels = "A") &
  theme(plot.subtitle = element_text(size = 9, colour = "grey30"))
save_fig(p_cut, "metaplot-cut-ratio.png", width = 9.5, height = 3.8)

# 3.4 several samples overlaid
res_multi <- metaPlot(bed_file = list(m6a_bed, m1a_bed, helako_bed),
                      group_names = c("m6A", "m1A", "m1A_helaA30KO"),
                      features_anno = features_anno, scale_region = TRUE)
save_fig(res_multi$plot, "metaplot-multi-sample.png", width = 7, height = 4)

# 3.5 faceted
res_facet <- metaPlot(bed_file = list(m6a_bed, m1a_bed, helako_bed),
                      group_names = c("m6A", "m1A", "m1A_helaA30KO"),
                      features_anno = features_anno, scale_region = TRUE,
                      facet = TRUE)
save_fig(res_facet$plot, "metaplot-facet.png", width = 10, height = 3.4)

# 3.6 styling arguments
p_default <- metaPlot(list(m6a_bed, m1a_bed), c("m6A", "m1A"), features_anno,
                      scale_region = TRUE)$plot +
  labs(subtitle = "default (col_alpha = 0.5, line_width = 1)")

p_alpha <- metaPlot(list(m6a_bed, m1a_bed), c("m6A", "m1A"), features_anno,
                    scale_region = TRUE, col_alpha = 0.15,
                    line_width = 1.6)$plot +
  labs(subtitle = "col_alpha = 0.15, line_width = 1.6")

p_bw_small <- metaPlot(list(m6a_bed, m1a_bed), c("m6A", "m1A"), features_anno,
                       scale_region = TRUE, geom_density_bw = 0.01)$plot +
  labs(subtitle = "geom_density_bw = 0.01 (undersmoothed)")

p_bw_big <- metaPlot(list(m6a_bed, m1a_bed), c("m6A", "m1A"), features_anno,
                     scale_region = TRUE, geom_density_bw = 0.15)$plot +
  labs(subtitle = "geom_density_bw = 0.15 (oversmoothed)")

p_style <- (p_default | p_alpha) / (p_bw_small | p_bw_big) +
  plot_annotation(tag_levels = "A") &
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        legend.position = "bottom")
save_fig(p_style, "metaplot-style.png", width = 9.5, height = 6.6)

# 3.7 single-sample colour
p_line <- metaPlot(list(m6a_bed), "m6A", features_anno, scale_region = TRUE,
                   line_color = "#2C7FB8", col_alpha = 0.3)$plot +
  labs(subtitle = 'line_color = "#2C7FB8", col_alpha = 0.3')
save_fig(p_line, "metaplot-line-color.png", width = 6, height = 4)

# ==============================================================================
# 4. motif helpers
# ==============================================================================
motif_hits <- findMotif(peaksfa = m6a_fa)
cat("\nfindMotif(): ", length(motif_hits), " hits, ",
    length(unique(motif_hits)), " distinct sequences\n", sep = "")
print(utils::head(sort(table(toupper(motif_hits)), decreasing = TRUE), 8))

motif_tab <- as.data.frame(table(toupper(motif_hits)), stringsAsFactors = FALSE) %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 12) %>%
  mutate(Var1 = factor(Var1, levels = rev(Var1)))

p_hits <- ggplot(motif_tab, aes(Freq, Var1)) +
  geom_col(fill = "#CC0033", width = 0.7) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = scales::comma) +
  labs(x = "Occurrences", y = NULL, title = "Top 12 matched motif sequences") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
save_fig(p_hits, "motif-hits.png", width = 5.2, height = 3.8)

rel_dist <- CalculatePeaksRelativeDist(peaksfa = m6a_fa)
cat("\nCalculatePeaksRelativeDist(): ", length(rel_dist), " positions, range ",
    paste(round(range(rel_dist), 1), collapse = " ~ "), "\n", sep = "")

p_motif <- plotRelMotif(relpos = list(rel_dist), exp = "m6A") +
  theme(legend.position = "none")
save_fig(p_motif, "motif-relative-position.png", width = 6, height = 4)

message("\nall figures written to ", fig_dir)
