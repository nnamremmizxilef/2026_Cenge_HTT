# ==============================================================================
# 04_05_te_annotation_classification.R
# Transposable element annotation visualization for Cenococcum geophilum
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggpubr)

# --- Paths --------------------------------------------------------------------
data_dir    <- "data/04_05_te_annotation_classification"
results_dir <- "results/04_05_te_annotation_classification"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ----------------------------------------------------------------
te_summary   <- read.delim(file.path(data_dir, "te_annotation_summary.tsv"), 
                           stringsAsFactors = FALSE)
te_class     <- read.delim(file.path(data_dir, "te_class_summary.tsv"), 
                           stringsAsFactors = FALSE)
te_filtering <- read.delim(file.path(data_dir, "te_filtering_summary.tsv"), 
                           stringsAsFactors = FALSE)

# Remove TOTAL row from filtering if present
te_filtering <- te_filtering %>% filter(genome != "TOTAL")

cat("Loaded TE annotation data for", nrow(te_summary), "genomes\n")
cat("Total TEs after filtering:", sum(te_summary$te_count), "\n")
cat("Total TE content:", round(sum(te_summary$total_bp) / 1e9, 2), "Gb\n")

# --- Summary statistics -------------------------------------------------------
cat("\n=== Summary Statistics ===\n")
cat("TEs per genome: mean =", round(mean(te_summary$te_count)), 
    ", range =", min(te_summary$te_count), "-", max(te_summary$te_count), "\n")
cat("TE content per genome: mean =", round(mean(te_summary$total_bp) / 1e6, 1), "Mb",
    ", range =", round(min(te_summary$total_bp) / 1e6, 1), "-", 
    round(max(te_summary$total_bp) / 1e6, 1), "Mb\n")
cat("Filtering removed:", round(mean(te_filtering$pct_removed), 1), "% on average\n")

# --- Plot 1: TE count per genome (stacked by class) ---------------------------
te_class_long <- te_class %>%
  pivot_longer(cols = -genome, names_to = "class", values_to = "count")

# Order genomes by total TE count
genome_order <- te_summary %>%
  arrange(te_count) %>%
  pull(genome)

te_class_long$genome <- factor(te_class_long$genome, levels = genome_order)

# Define class colors
class_colors <- c(
  "LTR" = "#264653",
  "DNA" = "#2A9D8F",
  "LINE" = "#E9C46A",
  "Unknown" = "#CCCCCC",
  "RC" = "#F4A261",
  "PLE" = "#E76F51",
  "SINE" = "#A23B72"
)

# Order classes for stacking
te_class_long$class <- factor(te_class_long$class, 
                              levels = c("LTR", "DNA", "LINE", "RC", "PLE", "SINE", "Unknown"))

p1 <- ggplot(te_class_long, aes(x = genome, y = count, fill = class)) +
  geom_col(alpha = 0.9) +
  geom_hline(yintercept = mean(te_summary$te_count), linetype = "dashed", 
             color = "#E63946", linewidth = 0.6) +
  scale_fill_manual(values = class_colors, name = "TE Class") +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Genome",
    y = "Number of TEs"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "te_count_per_genome.pdf"), p1, width = 12, height = 6)
ggsave(file.path(results_dir, "te_count_per_genome.png"), p1, width = 12, height = 6, dpi = 300)

# --- Plot 2: TE content (bp) per genome ---------------------------------------
te_summary_bp <- te_summary %>%
  arrange(total_bp) %>%
  mutate(
    genome = factor(genome, levels = genome),
    total_mb = total_bp / 1e6
  )

p2 <- ggplot(te_summary_bp, aes(x = genome, y = total_mb)) +
  geom_col(fill = "#2A9D8F", alpha = 0.8) +
  geom_hline(yintercept = mean(te_summary$total_bp) / 1e6, linetype = "dashed", 
             color = "#E63946", linewidth = 0.6) +
  labs(
    x = "Genome",
    y = "Total TE content (Mb)"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "te_content_per_genome.pdf"), p2, width = 10, height = 6)
ggsave(file.path(results_dir, "te_content_per_genome.png"), p2, width = 10, height = 6, dpi = 300)

# --- Plot 3: TE class composition (stacked bar) -------------------------------
te_class_long <- te_class %>%
  pivot_longer(cols = -genome, names_to = "class", values_to = "count") %>%
  group_by(genome) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Order genomes by LTR proportion
genome_order <- te_class_long %>%
  filter(class == "LTR") %>%
  arrange(proportion) %>%
  pull(genome)

te_class_long$genome <- factor(te_class_long$genome, levels = genome_order)

# Define class colors
class_colors <- c(
  "LTR" = "#264653",
  "DNA" = "#2A9D8F",
  "LINE" = "#E9C46A",
  "Unknown" = "#CCCCCC",
  "RC" = "#F4A261",
  "PLE" = "#E76F51",
  "SINE" = "#A23B72"
)

# Order classes for stacking
te_class_long$class <- factor(te_class_long$class, 
                              levels = c("LTR", "DNA", "LINE", "RC", "PLE", "SINE", "Unknown"))

p3 <- ggplot(te_class_long, aes(x = genome, y = proportion, fill = class)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = class_colors, name = "TE Class") +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Genome",
    y = "Proportion of TEs"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "te_class_composition.pdf"), p3, width = 12, height = 6)
ggsave(file.path(results_dir, "te_class_composition.png"), p3, width = 12, height = 6, dpi = 300)

# --- Plot 4: TE class totals across all genomes (pie/bar) ---------------------
te_class_totals <- te_class %>%
  select(-genome) %>%
  colSums() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("class") %>%
  setNames(c("class", "count")) %>%
  mutate(
    proportion = count / sum(count),
    label = paste0(class, "\n", round(proportion * 100, 1), "%")
  ) %>%
  arrange(desc(count))

te_class_totals$class <- factor(te_class_totals$class, levels = te_class_totals$class)

p4 <- ggplot(te_class_totals, aes(x = class, y = count, fill = class)) +
  geom_col(alpha = 0.9) +
  scale_fill_manual(values = class_colors) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "TE Class",
    y = "Total count across all genomes"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "te_class_totals.pdf"), p4, width = 8, height = 6)
ggsave(file.path(results_dir, "te_class_totals.png"), p4, width = 8, height = 6, dpi = 300)

# --- Plot 5: Filtering effect (raw vs filtered) -------------------------------
te_filtering_long <- te_filtering %>%
  select(genome, raw_count, filtered_count) %>%
  pivot_longer(cols = c(raw_count, filtered_count), 
               names_to = "type", values_to = "count") %>%
  mutate(type = factor(type, levels = c("raw_count", "filtered_count"),
                       labels = c("Raw", "Filtered (≥300bp, no satellites)")))

# Order by filtered count
genome_order_filt <- te_filtering %>%
  arrange(filtered_count) %>%
  pull(genome)

te_filtering_long$genome <- factor(te_filtering_long$genome, levels = genome_order_filt)

p5 <- ggplot(te_filtering_long, aes(x = genome, y = count, fill = type)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Raw" = "#CCCCCC", "Filtered (≥300bp, no satellites)" = "#2E86AB"),
                    name = "") +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Genome",
    y = "Number of TEs"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "te_filtering_effect.pdf"), p5, width = 12, height = 6)
ggsave(file.path(results_dir, "te_filtering_effect.png"), p5, width = 12, height = 6, dpi = 300)

# --- Plot 6: Percent removed by filtering -------------------------------------
te_filtering_ordered <- te_filtering %>%
  arrange(pct_removed) %>%
  mutate(genome = factor(genome, levels = genome))

p6 <- ggplot(te_filtering_ordered, aes(x = genome, y = pct_removed)) +
  geom_col(fill = "#E76F51", alpha = 0.8) +
  geom_hline(yintercept = mean(te_filtering$pct_removed), linetype = "dashed", 
             color = "#264653", linewidth = 0.6) +
  labs(
    x = "Genome",
    y = "Percent removed by filtering"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "te_percent_removed.pdf"), p6, width = 10, height = 6)
ggsave(file.path(results_dir, "te_percent_removed.png"), p6, width = 10, height = 6, dpi = 300)

# --- Plot 7: TE count vs TE content correlation -------------------------------
p7 <- ggplot(te_summary, aes(x = te_count, y = total_bp / 1e6)) +
  geom_point(size = 3, color = "#2E86AB", alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#E63946", linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = genome), size = 2.5, max.overlaps = 20) +
  labs(
    x = "Number of TEs",
    y = "Total TE content (Mb)"
  ) +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

# Calculate correlation
cor_test <- cor.test(te_summary$te_count, te_summary$total_bp)
cat("\nCorrelation TE count vs content: r =", round(cor_test$estimate, 3), 
    ", p =", format(cor_test$p.value, digits = 3), "\n")

ggsave(file.path(results_dir, "te_count_vs_content.pdf"), p7, width = 8, height = 6)
ggsave(file.path(results_dir, "te_count_vs_content.png"), p7, width = 8, height = 6, dpi = 300)

# --- Summary table ------------------------------------------------------------
summary_table <- te_summary %>%
  left_join(te_filtering %>% select(genome, raw_count, pct_removed), by = "genome") %>%
  mutate(
    total_mb = round(total_bp / 1e6, 1),
    mean_te_length = round(total_bp / te_count)
  ) %>%
  select(genome, raw_count, te_count, pct_removed, total_mb, mean_te_length)

write.table(summary_table, file.path(results_dir, "te_summary_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Class summary table ------------------------------------------------------
class_summary <- te_class_totals %>%
  select(class, count, proportion) %>%
  mutate(proportion = round(proportion * 100, 1))

write.table(class_summary, file.path(results_dir, "te_class_summary_table.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Output -------------------------------------------------------------------
cat("\n")
cat("Plots saved to:", results_dir, "\n")
cat("  - te_count_per_genome.pdf/png\n")
cat("  - te_content_per_genome.pdf/png\n")
cat("  - te_class_composition.pdf/png\n")
cat("  - te_class_totals.pdf/png\n")
cat("  - te_filtering_effect.pdf/png\n")
cat("  - te_percent_removed.pdf/png\n")
cat("  - te_count_vs_content.pdf/png\n")
cat("\n")
cat("Tables saved to:", results_dir, "\n")
cat("  - te_summary_table.tsv\n")
cat("  - te_class_summary_table.tsv\n")

# --- Session info -------------------------------------------------------------
sessionInfo()