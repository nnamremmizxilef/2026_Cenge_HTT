#!/usr/bin/env Rscript
# ==============================================================================
# Flanking region analysis - HTT validation
# ==============================================================================
# Compares TE similarity and flanking synteny between HTT and non-HTT pairs

library(tidyverse)
library(ggpubr)
library(scales)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
data_dir <- "data/07_flankid"
out_dir  <- "results/07_flankid"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
cat("Loading data...\n")

df    <- read_tsv(file.path(data_dir, "flankid_results.tsv"), show_col_types = FALSE)
stats <- read_tsv(file.path(data_dir, "summary_stats.tsv"), show_col_types = FALSE)

cat("  HTT pairs:", sum(df$label == "HTT"), "\n")
cat("  non-HTT pairs:", sum(df$label == "nonHTT"), "\n")

cat("Loading filtered HTT set...\n")

df_filt <- read_tsv(
  file.path(data_dir, "htt_candidates_flank_0pct_teid_gt75.tsv"),
  show_col_types = FALSE
)

# Ensure filtered set is HTT only (defensive)
if ("label" %in% colnames(df_filt)) {
  df_filt <- df_filt %>% filter(label == "HTT")
}

# df_filt may be either:
#  - a strict flankid results subset (same columns as df), OR
#  - a candidates list; in that case we join by pair_id to pull metrics from df.
need_metrics <- c("mean_pid", "flank_synt_cov_mean")
if (!all(need_metrics %in% colnames(df_filt))) {
  if (!("pair_id" %in% colnames(df_filt))) {
    stop("Filtered HTT file must contain either (mean_pid & flank_synt_cov_mean) OR pair_id.")
  }
  df_filt <- df %>%
    semi_join(df_filt %>% select(pair_id) %>% distinct(), by = "pair_id")
}

# Build 3-group dataset:
# - "HTT candidates (unfiltered)" from df where label == HTT
# - "HTT candidates (filtered)" from df_filt
# - "Non-HTT controls" from df where label == nonHTT
df3 <- bind_rows(
  df %>%
    filter(label == "HTT") %>%
    mutate(label3 = "HTT candidates (unfiltered)"),
  df_filt %>%
    mutate(label3 = "HTT candidates (filtered)"),
  df %>%
    filter(label == "nonHTT") %>%
    mutate(label3 = "Non-HTT controls")
) %>%
  mutate(
    label3 = factor(
      label3,
      levels = c("HTT candidates (unfiltered)", "HTT candidates (filtered)", "Non-HTT controls")
    )
  )

cols3 <- c(
  "HTT candidates (unfiltered)" = "#E64B35",
  "HTT candidates (filtered)"   = "#3C5488",
  "Non-HTT controls"            = "#4DBBD5"
)

cat("  Unfiltered HTT pairs:", sum(df3$label3 == "HTT candidates (unfiltered)"), "\n")
cat("  Filtered HTT pairs:  ", sum(df3$label3 == "HTT candidates (filtered)"), "\n")
cat("  Non-HTT pairs:       ", sum(df3$label3 == "Non-HTT controls"), "\n")

# Clean up labels for plotting (2-group plots)
df <- df %>%
  mutate(
    label = factor(label, levels = c("HTT", "nonHTT"),
                   labels = c("HTT candidates", "Non-HTT controls"))
  )

# Color palette (2-group)
cols <- c("HTT candidates" = "#E64B35", "Non-HTT controls" = "#4DBBD5")

# ------------------------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------------------------
cat("\n=== Summary Statistics ===\n\n")

summary_df <- df %>%
  group_by(label) %>%
  summarise(
    n = n(),
    te_cov_median = median(te_cov_mean),
    te_cov_mean   = mean(te_cov_mean),
    te_cov_sd     = sd(te_cov_mean),
    flank_median  = median(flank_synt_cov_mean),
    flank_mean    = mean(flank_synt_cov_mean),
    flank_sd      = sd(flank_synt_cov_mean),
    flank_zero_pct = 100 * sum(flank_synt_cov_mean == 0) / n(),
    htt_score_median = median(htt_score_te_minus_flank),
    htt_score_mean   = mean(htt_score_te_minus_flank),
    .groups = "drop"
  )

print(summary_df)
write_tsv(summary_df, file.path(out_dir, "summary_statistics.tsv"))

# Statistical tests (2-group, as before)
cat("\n=== Statistical Tests ===\n\n")

te_test <- wilcox.test(te_cov_mean ~ label, data = df)
cat("TE coverage (Wilcoxon):\n")
cat("  W =", te_test$statistic, ", p =", format.pval(te_test$p.value), "\n\n")

flank_test <- wilcox.test(flank_synt_cov_mean ~ label, data = df)
cat("Flank synteny (Wilcoxon):\n")
cat("  W =", flank_test$statistic, ", p =", format.pval(flank_test$p.value), "\n\n")

# ------------------------------------------------------------------------------
# Plot 1: TE coverage distribution (violin + boxplot) [unchanged]
# ------------------------------------------------------------------------------
p1 <- ggplot(df, aes(x = label, y = te_cov_mean, fill = label)) +
  geom_violin(alpha = 0.6, trim = FALSE, width = 0.8) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, alpha = 0.9) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.05)) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x.npc = "center", vjust = -0.5) +
  labs(
    x = NULL,
    y = "TE sequence coverage",
    title = "TE body similarity between paired elements"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(out_dir, "01_te_coverage_violin.png"), p1,
       width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "01_te_coverage_violin.pdf"), p1,
       width = 6, height = 5)

# ------------------------------------------------------------------------------
# Plot 2: Flank synteny distribution (violin + boxplot) [unchanged]
# ------------------------------------------------------------------------------
p2 <- ggplot(df, aes(x = label, y = flank_synt_cov_mean, fill = label)) +
  geom_violin(alpha = 0.6, trim = FALSE, width = 0.8) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, alpha = 0.9) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.05)) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x.npc = "center", vjust = -0.5) +
  labs(
    x = NULL,
    y = "Flanking region synteny",
    title = "Synteny-aware flanking similarity (±10 kb)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(out_dir, "02_flank_synteny_violin.png"), p2,
       width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "02_flank_synteny_violin.pdf"), p2,
       width = 6, height = 5)

# ------------------------------------------------------------------------------
# Plot 3: Combined panel (TE coverage + Flank) - PUBLICATION FIGURE [unchanged]
# ------------------------------------------------------------------------------
df_long <- df %>%
  select(pair_id, label, te_cov_mean, flank_synt_cov_mean) %>%
  pivot_longer(
    cols = c(te_cov_mean, flank_synt_cov_mean),
    names_to = "metric",
    values_to = "coverage"
  ) %>%
  mutate(
    metric = factor(metric,
                    levels = c("te_cov_mean", "flank_synt_cov_mean"),
                    labels = c("TE body", "Flanking regions (±10 kb)"))
  )

p3 <- ggplot(df_long, aes(x = label, y = coverage, fill = label)) +
  geom_violin(alpha = 0.6, trim = FALSE, width = 0.8) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, alpha = 0.9) +
  facet_wrap(~metric, scales = "fixed") +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.05)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     label.x.npc = "center", vjust = 1.5, size = 4) +
  labs(
    x = NULL,
    y = "Sequence similarity (coverage)",
    title = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10)
  )

ggsave(file.path(out_dir, "03_combined_panel_PUBLICATION.png"), p3,
       width = 8, height = 5, dpi = 300)
ggsave(file.path(out_dir, "03_combined_panel_PUBLICATION.pdf"), p3,
       width = 8, height = 5)

# ------------------------------------------------------------------------------
# Plot 3b: Combined panel (TE IDENTITY + Flank) with FILTERED HTTs added
#         - This plot shows EXACTLY what was filtered on:
#           mean_pid > 75% and flank_synt_cov_mean == 0
# ------------------------------------------------------------------------------
df3_long <- df3 %>%
  select(pair_id, label3, mean_pid, flank_synt_cov_mean) %>%
  pivot_longer(
    cols = c(mean_pid, flank_synt_cov_mean),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric,
                    levels = c("mean_pid", "flank_synt_cov_mean"),
                    labels = c("TE identity", "Flanking regions (±10 kb)"))
  )

p3b <- ggplot(df3_long, aes(x = label3, y = value, fill = label3)) +
  geom_violin(alpha = 0.6, trim = FALSE, width = 0.85) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, alpha = 0.9) +
  facet_wrap(~metric, scales = "fixed") +
  geom_hline(
    data = tibble(metric = "TE identity"),
    aes(yintercept = 75),
    linetype = "dashed",
    color = "gray40"
  ) +
  scale_fill_manual(values = cols3) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = NULL,
    y = "Sequence similarity",
    title = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 20, hjust = 1)
  )

ggsave(file.path(out_dir, "03b_identity_panel_WITH_FILTERED_PUBLICATION.png"), p3b,
       width = 9, height = 5, dpi = 300)
ggsave(file.path(out_dir, "03b_identity_panel_WITH_FILTERED_PUBLICATION.pdf"), p3b,
       width = 9, height = 5)

# ------------------------------------------------------------------------------
# Plot 4: Scatter plot TE coverage vs Flank [unchanged]
# ------------------------------------------------------------------------------
p4 <- ggplot(df, aes(x = te_cov_mean, y = flank_synt_cov_mean, color = label)) +
  geom_point(alpha = 0.3, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = cols) +
  scale_x_continuous(labels = percent_format(), limits = c(0, 1.05)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.05)) +
  labs(
    x = "TE body similarity",
    y = "Flanking region synteny",
    color = NULL,
    title = "TE similarity vs flanking synteny"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.75, 0.85),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(out_dir, "04_scatter_te_vs_flank.png"), p4,
       width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "04_scatter_te_vs_flank.pdf"), p4,
       width = 7, height = 6)

# ------------------------------------------------------------------------------
# Plot 5: HTT score distribution [unchanged]
# ------------------------------------------------------------------------------
p5 <- ggplot(df, aes(x = htt_score_te_minus_flank, fill = label)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(labels = percent_format()) +
  labs(
    x = "HTT score (TE coverage − flank synteny)",
    y = "Density",
    fill = NULL,
    title = "Distribution of HTT validation scores"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.2, 0.85),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(out_dir, "05_htt_score_density.png"), p5,
       width = 8, height = 5, dpi = 300)
ggsave(file.path(out_dir, "05_htt_score_density.pdf"), p5,
       width = 8, height = 5)

# ------------------------------------------------------------------------------
# Plot 6: Proportion with zero flanking synteny [unchanged]
# ------------------------------------------------------------------------------
zero_flank <- df %>%
  mutate(has_flank_synteny = flank_synt_cov_mean > 0) %>%
  group_by(label) %>%
  summarise(
    n = n(),
    n_zero = sum(!has_flank_synteny),
    n_nonzero = sum(has_flank_synteny),
    pct_zero = 100 * n_zero / n,
    .groups = "drop"
  )

cat("\n=== Flanking synteny breakdown ===\n")
print(zero_flank)

zero_flank_long <- zero_flank %>%
  select(label, `No synteny` = n_zero, `Some synteny` = n_nonzero) %>%
  pivot_longer(-label, names_to = "category", values_to = "count") %>%
  mutate(category = factor(category, levels = c("Some synteny", "No synteny")))

p6 <- ggplot(zero_flank_long, aes(x = label, y = count, fill = category)) +
  geom_col(position = "fill", width = 0.6, alpha = 0.9) +
  scale_fill_manual(values = c("No synteny" = "#2c3e50", "Some synteny" = "#95a5a6")) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = NULL,
    y = "Proportion of pairs",
    fill = NULL,
    title = "Proportion of pairs with detectable flanking synteny"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(out_dir, "06_zero_flank_proportion.png"), p6,
       width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "06_zero_flank_proportion.pdf"), p6,
       width = 6, height = 5)

# ------------------------------------------------------------------------------
# Plot 7: By TE superfamily [unchanged]
# ------------------------------------------------------------------------------
top_families <- df %>%
  count(family, sort = TRUE) %>%
  head(8) %>%
  pull(family)

df_top <- df %>%
  filter(family %in% top_families) %>%
  mutate(family = factor(family, levels = top_families))

p7 <- ggplot(df_top, aes(x = family, y = flank_synt_cov_mean, fill = label)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, position = position_dodge(0.8)) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = NULL,
    y = "Flanking region synteny",
    fill = NULL,
    title = "Flanking synteny by TE superfamily"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(out_dir, "07_flank_by_family.png"), p7,
       width = 10, height = 5, dpi = 300)
ggsave(file.path(out_dir, "07_flank_by_family.pdf"), p7,
       width = 10, height = 5)

# ------------------------------------------------------------------------------
# Output summary
# ------------------------------------------------------------------------------
cat("\n=== Output files ===\n")
cat("  ", file.path(out_dir, "summary_statistics.tsv"), "\n")
cat("  ", file.path(out_dir, "01_te_coverage_violin.png/pdf"), "\n")
cat("  ", file.path(out_dir, "02_flank_synteny_violin.png/pdf"), "\n")
cat("  ", file.path(out_dir, "03_combined_panel_PUBLICATION.png/pdf"), "\n")
cat("  ", file.path(out_dir, "03b_identity_panel_WITH_FILTERED_PUBLICATION.png/pdf"), " <- USE THIS (plots filtering criteria)\n")
cat("  ", file.path(out_dir, "04_scatter_te_vs_flank.png/pdf"), "\n")
cat("  ", file.path(out_dir, "05_htt_score_density.png/pdf"), "\n")
cat("  ", file.path(out_dir, "06_zero_flank_proportion.png/pdf"), "\n")
cat("  ", file.path(out_dir, "07_flank_by_family.png/pdf"), "\n")

cat("\nDone.\n")
