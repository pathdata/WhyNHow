# ============================================================================
# COMPREHENSIVE QC ANALYSIS FROM 10X (CHROMIUM) H5 FILE
# Including MT genes detection and distribution plots
# ============================================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
# ============================================================================
# LOAD DATA FROM H5 FILE
# ============================================================================

cat("LOADING 10X GENOMICS H5 FILE\n")

# Read 10X H5 file
h5_file <- "data\\GEX\\CAR-T10191180\\filtered_feature_bc_matrix.h5"
data_matrix <- Read10X_h5(h5_file)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = data_matrix, 
                                 project = "CAR-T_QC",
                                 min.cells = 0,
                                 min.features = 0)

# Get dimensions
n_cells <- ncol(seurat_obj)
n_genes <- nrow(seurat_obj)

cat("✅ Data loaded successfully!\n")
cat("n_obs (cells):  ", n_cells, "\n")
cat("n_vars (genes): ", n_genes, "\n\n")

# ============================================================================
# FIND MITOCHONDRIAL GENES
# ============================================================================

cat(strrep("=", 70), "\n")
cat("DETECTING MITOCHONDRIAL GENES\n")
cat(strrep("=", 70), "\n\n")

# Find MT genes (human: MT-, mouse: mt- or Mt-)
mt_genes_human <- rownames(seurat_obj)[grepl("^MT-", rownames(seurat_obj), ignore.case = FALSE)]
mt_genes_mouse <- rownames(seurat_obj)[grepl("^mt-", rownames(seurat_obj), ignore.case = FALSE)]
mt_genes <- unique(c(mt_genes_human, mt_genes_mouse))

cat("Mitochondrial genes found: ", length(mt_genes), "\n")
cat("Percentage of total genes: ", round(length(mt_genes)/n_genes*100, 2), "%\n\n")

if(length(mt_genes) > 0) {
  cat("MT genes detected:\n")
  for(i in 1:length(mt_genes)) {
    cat(sprintf("  %2d. %s\n", i, mt_genes[i]))
  }
  
  # Save MT genes list
  write.table(mt_genes, "mitochondrial_genes_list.txt", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat("\n✅ Saved MT genes list to: mitochondrial_genes_list.txt\n")
} else {
  cat("⚠️  No mitochondrial genes detected!\n")
  cat("    First 20 gene names:\n")
  print(head(rownames(seurat_obj), 20))
}

# ============================================================================
# CALCULATE QC METRICS
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("CALCULATING QC METRICS\n")
cat(strrep("=", 70), "\n\n")

# Calculate mitochondrial percentage
if(length(mt_genes) > 0) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
} else {
  # If no MT genes, set to 0
  seurat_obj[["percent.mt"]] <- 0
}

# Print summary statistics
cat("UMI Counts per cell:\n")
cat("  Mean:  ", round(mean(seurat_obj$nCount_RNA)), "\n")
cat("  Median:", round(median(seurat_obj$nCount_RNA)), "\n")
cat("  Range: ", min(seurat_obj$nCount_RNA), "-", max(seurat_obj$nCount_RNA), "\n\n")

cat("Genes detected per cell:\n")
cat("  Mean:  ", round(mean(seurat_obj$nFeature_RNA)), "\n")
cat("  Median:", round(median(seurat_obj$nFeature_RNA)), "\n")
cat("  Range: ", min(seurat_obj$nFeature_RNA), "-", max(seurat_obj$nFeature_RNA), "\n\n")

if(length(mt_genes) > 0) {
  cat("Mitochondrial percentage:\n")
  cat("  Mean:  ", round(mean(seurat_obj$percent.mt), 2), "%\n")
  cat("  Median:", round(median(seurat_obj$percent.mt), 2), "%\n")
  cat("  Range: ", round(min(seurat_obj$percent.mt), 2), "% -", 
      round(max(seurat_obj$percent.mt), 2), "%\n")
}

# ============================================================================
# CREATE QC PLOTS
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("GENERATING QC PLOTS\n")
cat(strrep("=", 70), "\n\n")

# Plot 1: Violin plot - Genes detected
p1 <- ggplot(seurat_obj@meta.data, aes(x = "", y = nFeature_RNA)) +
  geom_violin(fill = "#66d9ef", color = "#2d3748", alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = 16, outlier.size = 0.5) +
  labs(title = "Distribution of Genes per Cell",
       x = "", y = "Genes Detected (nFeature_RNA)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  ) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 4, 
               fill = "red", color = "darkred") +
  annotate("text", x = 1, y = Inf, 
           label = paste0("n = ", nrow(seurat_obj@meta.data), " cells\n",
                          "Median = ", round(median(seurat_obj@meta.data$nFeature_RNA))),
           hjust = 0.5, vjust = 2, size = 4, fontface = "bold")

ggsave("qc_plots\\80\\violin_genes_detected.png", p1, width = 6, height = 8, dpi = 300, bg = "white")
cat("✅ Saved: violin_genes_detected.png\n")

# Plot 2: Cell ranking by UMI counts
rank_df <- data.frame(
  rank = 1:n_cells,
  nCount_RNA = sort(seurat_obj@meta.data$nCount_RNA, decreasing = TRUE)
)

p2 <- ggplot(rank_df, aes(x = rank, y = nCount_RNA)) +
  geom_line(color = "#667eea", linewidth = 1) +
  geom_point(aes(color = nCount_RNA), size = 0.8, alpha = 0.6) +
  scale_color_gradient(low = "#feca57", high = "#ee5a6f", name = "UMI\nCounts") +
  scale_y_log10(labels = scales::comma) +
  geom_hline(yintercept = median(seurat_obj@meta.data$nCount_RNA), 
             linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Cell Ranking by Transcript Count",
       x = "Cell Index (Ranked by UMI)", 
       y = "UMI Counts (log scale)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  ) +
  annotate("text", x = n_cells * 0.8, y = median(seurat_obj@meta.data$nCount_RNA) * 1.5,
           label = paste0("Median = ", round(median(seurat_obj@meta.data$nCount_RNA))),
           color = "red", size = 4, fontface = "bold")

ggsave("qc_plots\\80\\cell_ranking_by_umi.png", p2, width = 10, height = 6, dpi = 300, bg = "white")
cat("✅ Saved: cell_ranking_by_umi.png\n")

# Plot 3: UMI count histogram
p3 <- ggplot(seurat_obj@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = median(seurat_obj@meta.data$nCount_RNA), 
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Median = ", round(median(seurat_obj@meta.data$nCount_RNA))),
           hjust = 1.1, vjust = 2, color = "red", size = 4, fontface = "bold") +
  labs(title = "UMI Count Distribution",
       x = "UMI Counts (nCount_RNA)", y = "Number of Cells") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave("qc_plots\\80\\umi_distribution.png", p3, width = 8, height = 6, dpi = 300, bg = "white")
cat("✅ Saved: umi_distribution.png\n")

# Plot 4: Genes detected histogram
p4 <- ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100, fill = "forestgreen", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = median(seurat_obj@meta.data$nFeature_RNA), 
             color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Median = ", round(median(seurat_obj@meta.data$nFeature_RNA))),
           hjust = 1.1, vjust = 2, color = "red", size = 4, fontface = "bold") +
  labs(title = "Genes Detected Distribution",
       x = "Genes Detected (nFeature_RNA)", y = "Number of Cells") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave("qc_plots\\80\\genes_distribution.png", p4, width = 8, height = 6, dpi = 300, bg = "white")
cat("✅ Saved: genes_distribution.png\n")

# Plot 5: Mitochondrial % distribution
if(length(mt_genes) > 0) {
  p5 <- ggplot(seurat_obj@meta.data, aes(x = percent.mt)) +
    geom_histogram(bins = 100, fill = "coral", color = "black", linewidth = 0.3) +
    geom_vline(xintercept = median(seurat_obj@meta.data$percent.mt), 
               color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Median = ", round(median(seurat_obj@meta.data$percent.mt), 2), "%"),
             hjust = 1.1, vjust = 2, color = "red", size = 4, fontface = "bold") +
    labs(title = "Mitochondrial Gene % Distribution",
         x = "Mitochondrial %", y = "Number of Cells") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
    )
  
  ggsave("qc_plots\\80\\mt_percent_distribution.png", p5, width = 8, height = 6, dpi = 300, bg = "white")
  cat("✅ Saved: mt_percent_distribution.png\n")
  
  # Plot 6: Violin plot - MT%
  p6 <- ggplot(seurat_obj@meta.data, aes(x = "", y = percent.mt)) +
    geom_violin(fill = "#fa709a", color = "#2d3748", alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.2, fill = "white", outlier.shape = 16, outlier.size = 0.5) +
    labs(title = "Distribution of Mitochondrial %",
         x = "", y = "Mitochondrial % (percent.mt)") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    ) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 4, 
                 fill = "red", color = "darkred") +
    annotate("text", x = 1, y = Inf, 
             label = paste0("Median = ", round(median(seurat_obj@meta.data$percent.mt), 2), "%"),
             hjust = 0.5, vjust = 2, size = 4, fontface = "bold")
  
  ggsave("qc_plots\\80\\violin_mt_percent.png", p6, width = 6, height = 8, dpi = 300, bg = "white")
  cat("✅ Saved: violin_mt_percent.png\n")
}

# Plot 7: UMI vs Genes scatter
p7 <- ggplot(seurat_obj@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = nFeature_RNA)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_viridis_c(option = "viridis", name = "Genes\nDetected") +
  labs(title = "UMI Counts vs Genes Detected",
       x = "UMI Counts (nCount_RNA)", y = "Genes Detected (nFeature_RNA)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave("qc_plots\\80\\umi_vs_genes_scatter.png", p7, width = 8, height = 6, dpi = 300, bg = "white")
cat("✅ Saved: umi_vs_genes_scatter.png\n")

# Combined plot
if(length(mt_genes) > 0) {
  combined <- (p3 | p4 | p5) / (p7 | p1 | p6)
  combined <- combined + plot_annotation(
    title = "Comprehensive QC Metrics from H5 File",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
} else {
  combined <- (p3 | p4 | p7) / (p1 | p2 | plot_spacer())
  combined <- combined + plot_annotation(
    title = "QC Metrics from H5 File (No MT genes detected)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )
}

ggsave("qc_plots\\80\\combined_qc_plots.png", combined, width = 18, height = 10, dpi = 300, bg = "white")
cat("✅ Saved: combined_qc_plots.png\n")

# ============================================================================
# SAVE SUMMARY REPORT
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("SAVING SUMMARY REPORT\n")
cat(strrep("=", 70), "\n\n")

sink("qc_plots\\80\\qc_summary_report.txt")
cat(strrep("=", 70), "\n")
cat("SINGLE-CELL QC SUMMARY FROM H5 FILE\n")
cat(strrep("=", 70), "\n\n")

cat("Input file: ", h5_file, "\n\n")

cat("Dataset dimensions:\n")
cat("  n_obs (cells):   ", n_cells, "\n")
cat("  n_vars (genes):  ", n_genes, "\n")
cat("  Total datapoints:", n_cells * n_genes, "\n\n")

cat("Mitochondrial genes:\n")
cat("  Count:      ", length(mt_genes), "\n")
cat("  Percentage: ", round(length(mt_genes)/n_genes*100, 2), "% of total genes\n")
if(length(mt_genes) > 0) {
  cat("  Genes:      ", paste(mt_genes, collapse = ", "), "\n")
} else {
  cat("  Genes:       None detected\n")
}
cat("\n")

cat("UMI Counts per cell:\n")
cat("  Mean:   ", round(mean(seurat_obj$nCount_RNA)), "\n")
cat("  Median: ", round(median(seurat_obj$nCount_RNA)), "\n")
cat("  SD:     ", round(sd(seurat_obj$nCount_RNA)), "\n")
cat("  Range:  ", min(seurat_obj$nCount_RNA), "-", max(seurat_obj$nCount_RNA), "\n\n")

cat("Genes detected per cell:\n")
cat("  Mean:   ", round(mean(seurat_obj$nFeature_RNA)), "\n")
cat("  Median: ", round(median(seurat_obj$nFeature_RNA)), "\n")
cat("  SD:     ", round(sd(seurat_obj$nFeature_RNA)), "\n")
cat("  Range:  ", min(seurat_obj$nFeature_RNA), "-", max(seurat_obj$nFeature_RNA), "\n\n")

if(length(mt_genes) > 0) {
  cat("Mitochondrial percentage:\n")
  cat("  Mean:   ", round(mean(seurat_obj$percent.mt), 2), "%\n")
  cat("  Median: ", round(median(seurat_obj$percent.mt), 2), "%\n")
  cat("  SD:     ", round(sd(seurat_obj$percent.mt), 2), "%\n")
  cat("  Range:  ", round(min(seurat_obj$percent.mt), 2), "% -", 
      round(max(seurat_obj$percent.mt), 2), "%\n\n")
}

cat("Generated files:\n")
cat("  - violin_genes_detected.png\n")
cat("  - cell_ranking_by_umi.png\n")
cat("  - umi_distribution.png\n")
cat("  - genes_distribution.png\n")
if(length(mt_genes) > 0) {
  cat("  - mt_percent_distribution.png\n")
  cat("  - violin_mt_percent.png\n")
}
cat("  - umi_vs_genes_scatter.png\n")
cat("  - combined_qc_plots.png\n")
if(length(mt_genes) > 0) {
  cat("  - qc_plots\\80\\mitochondrial_genes_list.txt\n")
}
cat("  - qc_plots\\80\\qc_summary_report.txt\n")

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 70), "\n")
sink()

cat("✅ Saved: qc_plots\\80\\qc_summary_report.txt\n")

cat("\n", strrep("=", 70), "\n")
cat("ALL DONE! 🎉\n")
cat(strrep("=", 70), "\n")

