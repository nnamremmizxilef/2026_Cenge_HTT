# ==============================================================================
# 06_htt_candidates.R
# Horizontal transposon transfer (HTT) detection visualization for C. geophilum
# Including quantile sensitivity analysis
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(data.table)
library(ape)
library(ggtree)
library(stringr)
library(circlize)
library(RColorBrewer)
library(scales)

# --- Paths --------------------------------------------------------------------
data_dir    <- "data/06_htt_candidates"
results_dir <- "results/06_htt_candidates"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ----------------------------------------------------------------
cat("Loading data...\n")

# Use data.table for large files
te_ks <- fread(file.path(data_dir, "te_ks_values.tsv"), stringsAsFactors = FALSE) %>% as.data.frame()
busco_ks <- read.delim(file.path(data_dir, "node_ks_summary.tsv"), stringsAsFactors = FALSE)
busco_raw <- fread(file.path(data_dir, "ks_values_by_node.tsv"), stringsAsFactors = FALSE) %>% as.data.frame()
clade_assign <- read.delim(file.path(data_dir, "clade_assignment.tsv"), stringsAsFactors = FALSE)

# Load HTT candidates (0.5% quantile)
htt_file <- file.path(data_dir, "htt_candidates.tsv")
if (file.exists(htt_file) && file.size(htt_file) > 0) {
  htt_candidates <- fread(htt_file, stringsAsFactors = FALSE) %>% as.data.frame()
} else {
  htt_candidates <- data.frame()
}

# Load non-HTT pairs (contains MRCA info for all other pairs)
non_htt_file <- file.path(data_dir, "non_htt_pairs.tsv")
if (file.exists(non_htt_file) && file.size(non_htt_file) > 0) {
  cat("Loading non-HTT pairs (this may take a moment)...\n")
  non_htt_pairs <- fread(non_htt_file, stringsAsFactors = FALSE) %>% as.data.frame()
} else {
  non_htt_pairs <- data.frame()
}

# Load clade pair summary
clade_pair_file <- file.path(data_dir, "htt_by_clade_pair.tsv")
if (file.exists(clade_pair_file) && file.size(clade_pair_file) > 0) {
  htt_by_clade <- read.delim(clade_pair_file, stringsAsFactors = FALSE)
} else {
  htt_by_clade <- data.frame()
}

# Combine HTT and non-HTT pairs to get all pairs with MRCA info
cat("Combining HTT and non-HTT pairs...\n")
if (nrow(htt_candidates) > 0 && nrow(non_htt_pairs) > 0) {
  all_pairs <- bind_rows(
    htt_candidates %>% select(pair_id, mrca_node, te_ks, ka, qclade, sclade, qgenome, sgenome),
    non_htt_pairs %>% select(pair_id, mrca_node, te_ks, ka, qclade, sclade, qgenome, sgenome)
  )
} else if (nrow(htt_candidates) > 0) {
  all_pairs <- htt_candidates %>% select(pair_id, mrca_node, te_ks, ka, qclade, sclade, qgenome, sgenome)
} else if (nrow(non_htt_pairs) > 0) {
  all_pairs <- non_htt_pairs %>% select(pair_id, mrca_node, te_ks, ka, qclade, sclade, qgenome, sgenome)
} else {
  stop("No HTT or non-HTT pair data found!")
}

# Rename te_ks to ks for consistency
if ("te_ks" %in% colnames(all_pairs)) {
  all_pairs$ks <- all_pairs$te_ks
}

# Filter valid Ks values
all_pairs_valid <- all_pairs %>% filter(is.finite(ks) & ks > 0 & ks < 10)
busco_raw_valid <- busco_raw %>% filter(is.finite(ks) & ks > 0)

# Mark original HTT status (from htt_candidates.tsv) WITHOUT overwriting later
all_pairs_valid$is_htt_0.5_original <- all_pairs_valid$pair_id %in% htt_candidates$pair_id

cat("Total pairs with MRCA info:", nrow(all_pairs_valid), "\n")
cat("HTT candidates (0.5% file):", sum(all_pairs_valid$is_htt_0.5_original), "\n")

# Free memory
rm(non_htt_pairs)
gc()

# --- Quantile Sensitivity Analysis --------------------------------------------

cat("\n=== Quantile Sensitivity Analysis ===\n")

# --- Calculate all quantile thresholds from raw BUSCO Ks ----------------------
quantiles_to_test <- c(0.005, 0.01, 0.025, 0.05, 0.075, 0.10)
quantile_labels <- c("0.5%", "1%", "2.5%", "5%", "7.5%", "10%")

quantile_thresholds <- busco_raw_valid %>%
  group_by(node) %>%
  summarise(
    n_genes = n(),
    q_0.5 = quantile(ks, 0.005, na.rm = TRUE),
    q_1 = quantile(ks, 0.01, na.rm = TRUE),
    q_2.5 = quantile(ks, 0.025, na.rm = TRUE),
    q_5 = quantile(ks, 0.05, na.rm = TRUE),
    q_7.5 = quantile(ks, 0.075, na.rm = TRUE),
    q_10 = quantile(ks, 0.10, na.rm = TRUE),
    median = median(ks, na.rm = TRUE),
    .groups = "drop"
  )

cat("Calculated quantile thresholds for", nrow(quantile_thresholds), "nodes\n")

# Save quantile thresholds
write.table(quantile_thresholds, file.path(results_dir, "busco_quantile_thresholds.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Classify HTT for each quantile threshold ---------------------------------
cat("Classifying HTT for each quantile...\n")

# Create threshold lookup for each quantile
threshold_cols <- c("q_0.5", "q_1", "q_2.5", "q_5", "q_7.5", "q_10")

# Add thresholds for each quantile to all_pairs_valid
all_pairs_classified <- all_pairs_valid %>%
  left_join(quantile_thresholds %>% select(node, all_of(threshold_cols)),
            by = c("mrca_node" = "node"))

# Report pairs without thresholds
na_thresh <- sum(is.na(all_pairs_classified$q_0.5))
if (na_thresh > 0) {
  cat("Pairs with no BUSCO threshold (NA q_0.5):", na_thresh, "\n")
}

# Classify HTT for each quantile (keep your existing column names)
all_pairs_classified <- all_pairs_classified %>%
  mutate(
    is_htt_0.5 = !is.na(q_0.5) & ks < q_0.5,
    is_htt_1 = !is.na(q_1) & ks < q_1,
    is_htt_2.5 = !is.na(q_2.5) & ks < q_2.5,
    is_htt_5 = !is.na(q_5) & ks < q_5,
    is_htt_7.5 = !is.na(q_7.5) & ks < q_7.5,
    is_htt_10 = !is.na(q_10) & ks < q_10
  )

# Compare original file vs recomputed 0.5%
cmp_tbl <- data.frame(
  metric = c("original_file_0.5%", "recomputed_0.5%", "overlap", "file_only", "recomputed_only"),
  value = c(
    sum(all_pairs_classified$is_htt_0.5_original, na.rm = TRUE),
    sum(all_pairs_classified$is_htt_0.5, na.rm = TRUE),
    sum(all_pairs_classified$is_htt_0.5_original & all_pairs_classified$is_htt_0.5, na.rm = TRUE),
    sum(all_pairs_classified$is_htt_0.5_original & !all_pairs_classified$is_htt_0.5, na.rm = TRUE),
    sum(!all_pairs_classified$is_htt_0.5_original & all_pairs_classified$is_htt_0.5, na.rm = TRUE)
  )
)
cat("\n0.5% HTT comparison (file vs recomputed):\n")
print(cmp_tbl)
write.table(cmp_tbl, file.path(results_dir, "htt_0.5_file_vs_recomputed.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Summarize results for each quantile
quantile_results <- data.frame(
  quantile = quantile_labels,
  quantile_col = threshold_cols,
  htt_count = c(
    sum(all_pairs_classified$is_htt_0.5, na.rm = TRUE),
    sum(all_pairs_classified$is_htt_1, na.rm = TRUE),
    sum(all_pairs_classified$is_htt_2.5, na.rm = TRUE),
    sum(all_pairs_classified$is_htt_5, na.rm = TRUE),
    sum(all_pairs_classified$is_htt_7.5, na.rm = TRUE),
    sum(all_pairs_classified$is_htt_10, na.rm = TRUE)
  ),
  total_pairs = nrow(all_pairs_classified %>% filter(!is.na(q_0.5))),
  stringsAsFactors = FALSE
) %>%
  mutate(
    htt_percentage = 100 * htt_count / total_pairs,
    threshold_median = c(
      median(quantile_thresholds$q_0.5, na.rm = TRUE),
      median(quantile_thresholds$q_1, na.rm = TRUE),
      median(quantile_thresholds$q_2.5, na.rm = TRUE),
      median(quantile_thresholds$q_5, na.rm = TRUE),
      median(quantile_thresholds$q_7.5, na.rm = TRUE),
      median(quantile_thresholds$q_10, na.rm = TRUE)
    )
  )

# Add Ks statistics for HTT candidates at each threshold
for (i in 1:nrow(quantile_results)) {
  htt_col <- paste0("is_htt_", gsub("%", "", quantile_results$quantile[i]))
  htt_ks <- all_pairs_classified$ks[all_pairs_classified[[htt_col]] == TRUE]
  
  quantile_results$htt_ks_min[i] <- ifelse(length(htt_ks) > 0, min(htt_ks, na.rm = TRUE), NA)
  quantile_results$htt_ks_max[i] <- ifelse(length(htt_ks) > 0, max(htt_ks, na.rm = TRUE), NA)
  quantile_results$htt_ks_median[i] <- ifelse(length(htt_ks) > 0, median(htt_ks, na.rm = TRUE), NA)
}

quantile_results$quantile <- factor(quantile_results$quantile, levels = quantile_labels)

cat("\nHTT detection by quantile threshold:\n")
print(quantile_results %>% select(quantile, htt_count, htt_percentage, threshold_median))

# Save quantile results
write.table(quantile_results %>% select(-quantile_col),
            file.path(results_dir, "htt_quantile_sensitivity.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Plots --------------------------------------------------------------------

# Plot 1: HTT count by quantile (bar plot)
p_q1 <- ggplot(quantile_results, aes(x = quantile, y = htt_count)) +
  geom_col(fill = "#2A9D8F", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = comma(htt_count)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "BUSCO Ks quantile threshold",
    y = "Number of HTT candidates"
  ) +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

ggsave(file.path(results_dir, "htt_quantile_count.pdf"), p_q1, width = 8, height = 6)
ggsave(file.path(results_dir, "htt_quantile_count.png"), p_q1, width = 8, height = 6, dpi = 300)

# Plot 2: HTT percentage by quantile
p_q2 <- ggplot(quantile_results, aes(x = quantile, y = htt_percentage)) +
  geom_col(fill = "#E9C46A", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.3f%%", htt_percentage)), vjust = -0.5, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "BUSCO Ks quantile threshold",
    y = "HTT percentage"
  ) +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

ggsave(file.path(results_dir, "htt_quantile_percentage.pdf"), p_q2, width = 8, height = 6)
ggsave(file.path(results_dir, "htt_quantile_percentage.png"), p_q2, width = 8, height = 6, dpi = 300)

# Plot 3: Sensitivity curve
quantile_results_plot <- quantile_results %>%
  mutate(quantile_num = as.numeric(gsub("%", "", as.character(quantile))))

p_q3 <- ggplot(quantile_results_plot, aes(x = quantile_num, y = htt_count)) +
  geom_line(color = "#264653", linewidth = 1.2) +
  geom_point(color = "#264653", size = 3) +
  geom_text(aes(label = comma(htt_count)), vjust = -1, size = 3) +
  scale_x_continuous(
    breaks = quantile_results_plot$quantile_num,
    labels = paste0(quantile_results_plot$quantile_num, "%")
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "BUSCO Ks quantile threshold",
    y = "Number of HTT candidates"
  ) +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

ggsave(file.path(results_dir, "htt_sensitivity_curve.pdf"), p_q3, width = 8, height = 6)
ggsave(file.path(results_dir, "htt_sensitivity_curve.png"), p_q3, width = 8, height = 6, dpi = 300)

# Plot 4: Ks distribution with all quantile thresholds
threshold_lines <- data.frame(
  quantile = factor(quantile_labels, levels = quantile_labels),
  threshold = quantile_results$threshold_median
)

p_q4 <- ggplot(all_pairs_classified, aes(x = ks)) +
  geom_histogram(bins = 100, fill = "#AAAAAA", alpha = 0.7) +
  geom_vline(data = threshold_lines, aes(xintercept = threshold, color = quantile),
             linetype = "dashed", linewidth = 0.8) +
  scale_x_log10(labels = scales::scientific) +
  scale_color_viridis_d(option = "D", name = "Quantile\n(median)") +
  labs(
    x = "TE Ks (log scale)",
    y = "Number of TE pairs"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 9)
  )

ggsave(file.path(results_dir, "htt_ks_all_thresholds.pdf"), p_q4, width = 10, height = 6)
ggsave(file.path(results_dir, "htt_ks_all_thresholds.png"), p_q4, width = 10, height = 6, dpi = 300)

# Plot 5: Fold change from 0.5% baseline
baseline_htt <- quantile_results$htt_count[quantile_results$quantile == "0.5%"]

quantile_results_fc <- quantile_results %>%
  mutate(
    fold_change = htt_count / baseline_htt,
    quantile_num = as.numeric(gsub("%", "", as.character(quantile)))
  )

p_q5 <- ggplot(quantile_results_fc, aes(x = quantile_num, y = fold_change)) +
  geom_line(color = "#E76F51", linewidth = 1.2) +
  geom_point(color = "#E76F51", size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#666666") +
  geom_text(aes(label = sprintf("%.1fx", fold_change)), vjust = -1, size = 3) +
  scale_x_continuous(
    breaks = quantile_results_fc$quantile_num,
    labels = paste0(quantile_results_fc$quantile_num, "%")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "BUSCO Ks quantile threshold",
    y = "Fold change relative to 0.5%"
  ) +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

ggsave(file.path(results_dir, "htt_fold_change.pdf"), p_q5, width = 8, height = 6)
ggsave(file.path(results_dir, "htt_fold_change.png"), p_q5, width = 8, height = 6, dpi = 300)

# Plot 6: HTT by node for multiple quantiles
htt_by_node_quantile <- all_pairs_classified %>%
  group_by(mrca_node) %>%
  summarise(
    `0.5%` = sum(is_htt_0.5, na.rm = TRUE),
    `1%` = sum(is_htt_1, na.rm = TRUE),
    `2.5%` = sum(is_htt_2.5, na.rm = TRUE),
    `5%` = sum(is_htt_5, na.rm = TRUE),
    `7.5%` = sum(is_htt_7.5, na.rm = TRUE),
    `10%` = sum(is_htt_10, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = -mrca_node, names_to = "quantile", values_to = "htt_count") %>%
  left_join(quantile_thresholds %>% select(node, q_0.5), by = c("mrca_node" = "node"))

htt_by_node_quantile$quantile <- factor(htt_by_node_quantile$quantile, levels = quantile_labels)

node_order <- quantile_thresholds %>% arrange(q_0.5) %>% pull(node)
htt_by_node_quantile$mrca_node <- factor(htt_by_node_quantile$mrca_node, levels = node_order)

p_q6 <- ggplot(htt_by_node_quantile %>% filter(!is.na(mrca_node)),
               aes(x = mrca_node, y = htt_count, fill = quantile)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_viridis_d(option = "D", name = "Quantile") +
  scale_y_continuous(labels = comma) +
  labs(
    x = "MRCA node (ordered by divergence)",
    y = "Number of HTT candidates"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "bottom"
  )

ggsave(file.path(results_dir, "htt_by_node_quantiles.pdf"), p_q6, width = 14, height = 6)
ggsave(file.path(results_dir, "htt_by_node_quantiles.png"), p_q6, width = 14, height = 6, dpi = 300)

write.table(htt_by_node_quantile %>% select(-q_0.5) %>% pivot_wider(names_from = quantile, values_from = htt_count),
            file.path(results_dir, "htt_by_node_quantiles.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Standard HTT Visualizations (0.5% quantile) ------------------------------
cat("\n=== Standard HTT Analysis (0.5% quantile) ===\n")

# Use is_htt_0.5 as the standard HTT classification (recomputed)
te_ks_valid <- all_pairs_classified %>%
  mutate(is_htt = is_htt_0.5)

# --- HTT aggregated by clade pair (filter within-clade) -----------------------
clade_map <- clade_assign %>% select(sample, clade)

te_ks_valid <- te_ks_valid %>%
  left_join(clade_map, by = c("qgenome" = "sample")) %>%
  rename(qclade_map = clade) %>%
  left_join(clade_map, by = c("sgenome" = "sample")) %>%
  rename(sclade_map = clade)

# fall back to original columns if mapping missing
te_ks_valid$qclade_final <- ifelse(is.na(te_ks_valid$qclade_map), te_ks_valid$qclade, te_ks_valid$qclade_map)
te_ks_valid$sclade_final <- ifelse(is.na(te_ks_valid$sclade_map), te_ks_valid$sclade, te_ks_valid$sclade_map)

# keep only between-clade HTT (collapsed analysis unit)
htt_by_clade <- te_ks_valid %>%
  filter(is_htt == TRUE) %>%
  filter(qclade_final != sclade_final) %>%
  count(qclade = qclade_final, sclade = sclade_final, name = "htt_count")

# make symmetric (undirected) summary so A-B and B-A are combined
htt_by_clade_sym <- bind_rows(
  htt_by_clade,
  htt_by_clade %>% rename(qclade = sclade, sclade = qclade)
) %>%
  group_by(qclade, sclade) %>%
  summarise(htt_count = sum(htt_count), .groups = "drop")

write.table(htt_by_clade, file.path(results_dir, "htt_by_clade_pair_from_strain_calls.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(htt_by_clade_sym, file.path(results_dir, "htt_by_clade_pair_from_strain_calls_symmetric.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("HTT by clade pair (between-clade only):", nrow(htt_by_clade), "pairs\n")

# --- Summary statistics -------------------------------------------------------
cat("Total TE pairs with Ks:", nrow(te_ks_valid), "\n")
cat("HTT candidates:", sum(te_ks_valid$is_htt, na.rm = TRUE),
    "(", round(100 * sum(te_ks_valid$is_htt, na.rm = TRUE) / nrow(te_ks_valid), 4), "%)\n")

if (sum(te_ks_valid$is_htt, na.rm = TRUE) > 0) {
  htt_ks <- te_ks_valid$ks[te_ks_valid$is_htt == TRUE]
  cat("HTT Ks range:", round(min(htt_ks), 4), "-", round(max(htt_ks), 4), "\n")
  cat("HTT Ks median:", round(median(htt_ks), 4), "\n")
}

non_htt_ks <- te_ks_valid$ks[te_ks_valid$is_htt == FALSE]
cat("Non-HTT Ks range:", round(min(non_htt_ks), 4), "-", round(max(non_htt_ks), 4), "\n")
cat("Non-HTT Ks median:", round(median(non_htt_ks), 4), "\n")

busco_valid <- busco_ks %>% filter(pass_0.01 == TRUE)
cat("\nBUSCO Ks 0.5% quantile range:",
    round(min(busco_valid$ks_0.5_quantile, na.rm = TRUE), 4), "-",
    round(max(busco_valid$ks_0.5_quantile, na.rm = TRUE), 4), "\n")

# Plot 1: TE Ks distribution (HTT vs non-HTT)
min_threshold <- min(busco_valid$ks_0.5_quantile, na.rm = TRUE)

p1 <- ggplot(te_ks_valid, aes(x = ks, fill = is_htt)) +
  geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = min_threshold, linetype = "dashed",
             color = "#E63946", linewidth = 0.8) +
  scale_x_log10(labels = scales::scientific) +
  scale_fill_manual(
    values = c("FALSE" = "#AAAAAA", "TRUE" = "#2A9D8F"),
    labels = c("Non-HTT", "HTT candidate"),
    name = ""
  ) +
  annotate("text", x = min_threshold * 1.5, y = Inf, vjust = 2, hjust = 0,
           label = paste0("Min BUSCO 0.5% = ", round(min_threshold, 4)),
           size = 3, color = "#E63946") +
  labs(
    x = "Ks (log scale)",
    y = "Number of TE pairs"
  ) +
  theme_pubr() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "htt_ks_distribution.pdf"), p1, width = 10, height = 6)
ggsave(file.path(results_dir, "htt_ks_distribution.png"), p1, width = 10, height = 6, dpi = 300)

# Plot 2: Ks density comparison
p2 <- ggplot(te_ks_valid, aes(x = ks, fill = is_htt, color = is_htt)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  scale_x_log10(labels = scales::scientific) +
  scale_fill_manual(
    values = c("FALSE" = "#AAAAAA", "TRUE" = "#2A9D8F"),
    labels = c("Non-HTT", "HTT candidate"),
    name = ""
  ) +
  scale_color_manual(
    values = c("FALSE" = "#666666", "TRUE" = "#1D7A6E"),
    labels = c("Non-HTT", "HTT candidate"),
    name = ""
  ) +
  labs(
    x = "Ks (log scale)",
    y = "Density"
  ) +
  theme_pubr() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "htt_ks_density.pdf"), p2, width = 10, height = 6)
ggsave(file.path(results_dir, "htt_ks_density.png"), p2, width = 10, height = 6, dpi = 300)

# Plot 2b: Ks density with HTT colored by TE class
extract_class <- function(te_id) {
  classification <- str_extract(te_id, "(?<=\\|)[^:]+(?=::)")
  str_extract(classification, "^[^/]+")
}

# Join TE IDs from htt_candidates (te_ks only has HTT pairs)
te_ks_valid <- te_ks_valid %>%
  left_join(htt_candidates %>% select(pair_id, qseqid, sseqid), by = "pair_id")

te_ks_valid <- te_ks_valid %>%
  mutate(
    q_class = extract_class(qseqid),
    s_class = extract_class(sseqid)
  )

# Prepare data for plot
htt_with_class <- te_ks_valid %>% filter(is_htt == TRUE & !is.na(q_class))
non_htt <- te_ks_valid %>% filter(is_htt == FALSE)

class_counts <- htt_with_class %>% count(q_class, sort = TRUE)
cat("\nHTT candidates by TE class:\n")
print(class_counts)

p2b <- ggplot() +
  # Non-HTT as grey background
  geom_density(data = non_htt, aes(x = ks), 
               fill = "#AAAAAA", color = "#666666", alpha = 0.4, linewidth = 0.6) +
  # HTT colored by TE class
  geom_density(data = htt_with_class, aes(x = ks, fill = q_class, color = q_class), 
               alpha = 0.5, linewidth = 0.6) +
  scale_x_log10(labels = scales::scientific) +
  scale_fill_brewer(palette = "Set1", name = "HTT TE Class") +
  scale_color_brewer(palette = "Set1", name = "HTT TE Class") +
  labs(
    x = "Ks (log scale)",
    y = "Density"
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "htt_ks_density_by_class.pdf"), 
       p2b, width = 10, height = 6)
ggsave(file.path(results_dir, "htt_ks_density_by_class.png"), 
       p2b, width = 10, height = 6, dpi = 300)

write.table(class_counts, file.path(results_dir, "htt_by_te_class.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot 3: Ka vs Ks scatter
te_ks_scatter <- te_ks_valid %>% filter(is.finite(ka) & ka > 0 & ka < 10)

p3 <- ggplot(te_ks_scatter, aes(x = ks, y = ka, color = is_htt)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#666666") +
  scale_x_log10(labels = scales::scientific) +
  scale_y_log10(labels = scales::scientific) +
  scale_color_manual(
    values = c("FALSE" = "#AAAAAA", "TRUE" = "#2A9D8F"),
    labels = c("Non-HTT", "HTT candidate"),
    name = ""
  ) +
  labs(
    x = "Ks (log scale)",
    y = "Ka (log scale)"
  ) +
  theme_pubr() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "htt_ka_vs_ks.pdf"), p3, width = 8, height = 7)
ggsave(file.path(results_dir, "htt_ka_vs_ks.png"), p3, width = 8, height = 7, dpi = 300)

# Plot 4: HTT by clade pair heatmap
if (nrow(htt_by_clade_sym) > 0) {
  all_clades <- unique(c(htt_by_clade_sym$qclade, htt_by_clade_sym$sclade))
  
  complete_grid <- expand.grid(qclade = all_clades, sclade = all_clades,
                               stringsAsFactors = FALSE)
  
  htt_heatmap_data <- complete_grid %>%
    left_join(htt_by_clade_sym, by = c("qclade", "sclade")) %>%
    mutate(htt_count = replace_na(htt_count, 0))
  
  htt_heatmap_sym <- htt_heatmap_data %>%
    rename(c1 = qclade, c2 = sclade) %>%
    bind_rows(
      htt_heatmap_data %>% rename(c1 = sclade, c2 = qclade)
    ) %>%
    group_by(c1, c2) %>%
    summarize(htt_count = sum(htt_count), .groups = "drop")
  
  p4 <- ggplot(htt_heatmap_sym, aes(x = c1, y = c2, fill = htt_count)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(
      low = "white", high = "#264653",
      name = "HTT count",
      trans = "log1p",
      breaks = c(0, 10, 100, 1000),
      labels = c("0", "10", "100", "1k")
    ) +
    labs(x = "Clade", y = "Clade") +
    theme_pubr() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      legend.position = "right",
      panel.grid = element_blank()
    ) +
    coord_fixed()
  
  ggsave(file.path(results_dir, "htt_clade_heatmap.pdf"), p4, width = 10, height = 9)
  ggsave(file.path(results_dir, "htt_clade_heatmap.png"), p4, width = 10, height = 9, dpi = 300)
}

# Plot 5: HTT Ks vs BUSCO threshold scatter
htt_pairs <- te_ks_valid %>% filter(is_htt == TRUE)

if (nrow(htt_pairs) > 0) {
  p5 <- ggplot(htt_pairs, aes(x = q_0.5, y = ks)) +
    geom_point(alpha = 0.3, size = 1, color = "#2A9D8F") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#E63946") +
    scale_x_log10(labels = scales::scientific) +
    scale_y_log10(labels = scales::scientific) +
    labs(
      x = "BUSCO Ks threshold (0.5% quantile)",
      y = "TE Ks",
      caption = "Points below diagonal = HTT candidates (TE Ks < BUSCO threshold)"
    ) +
    theme_pubr() +
    theme(
      axis.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(results_dir, "htt_te_vs_busco_ks.pdf"), p5, width = 8, height = 7)
  ggsave(file.path(results_dir, "htt_te_vs_busco_ks.png"), p5, width = 8, height = 7, dpi = 300)
}

# Plot 6: Boxplot of Ks by HTT status
p6 <- ggplot(te_ks_valid, aes(x = is_htt, y = ks, fill = is_htt)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_y_log10(labels = scales::scientific) +
  scale_x_discrete(labels = c("FALSE" = "Non-HTT", "TRUE" = "HTT")) +
  scale_fill_manual(
    values = c("FALSE" = "#AAAAAA", "TRUE" = "#2A9D8F"),
    guide = "none"
  ) +
  labs(x = "", y = "Ks (log scale)") +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "htt_ks_boxplot.pdf"), p6, width = 5, height = 6)
ggsave(file.path(results_dir, "htt_ks_boxplot.png"), p6, width = 5, height = 6, dpi = 300)

# Plot 7: HTT chord diagram

in_file   <- file.path(data_dir, "htt_by_clade_pair.tsv")
clade_map <- file.path(data_dir, "clade_assignment.tsv")

out_pdf <- file.path(results_dir, "htt_chord_diagram.pdf")
out_png <- file.path(results_dir, "htt_chord_diagram.png")
out_tsv <- file.path(results_dir, "htt_by_clade_pair_between_only_undirected.tsv")

stopifnot(file.exists(in_file), file.exists(clade_map))

clade_assign <- read.delim(clade_map, stringsAsFactors = FALSE) %>%
  mutate(sample = as.character(sample), clade = as.character(clade))

# collapsed labels with member strains (MULTI-LINE)
collapsed_members <- clade_assign %>%
  filter(grepl("^collapsed_node_", clade)) %>%
  group_by(clade) %>%
  summarise(
    members = paste(sample, collapse = "\n"),
    .groups = "drop"
  )

pretty_label <- function(x) {
  if (grepl("^collapsed_node_", x)) {
    num <- str_match(x, "^collapsed_node_(\\d+)$")[,2]
    memb <- collapsed_members$members[match(x, collapsed_members$clade)]
    paste0("Clade ", num, ":\n", memb)
  } else {
    x
  }
}

links0 <- read.delim(in_file, stringsAsFactors = FALSE) %>%
  mutate(
    qclade = as.character(qclade),
    sclade = as.character(sclade),
    htt_count = as.numeric(htt_count)
  ) %>%
  filter(is.finite(htt_count) & htt_count > 0) %>%
  filter(qclade != sclade)

# undirected combine (A-B + B-A)
links_u <- links0 %>%
  mutate(a = pmin(qclade, sclade), b = pmax(qclade, sclade)) %>%
  group_by(a, b) %>%
  summarise(htt_count = sum(htt_count), .groups = "drop") %>%
  rename(from = a, to = b)

write.table(links_u, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

# suppress weak links
links_plot <- links_u %>% filter(htt_count >= 5)
if (nrow(links_plot) == 0)
  stop("No links left after filtering (htt_count >= 5). Lower the cutoff.")

# sector order
sectors <- sort(unique(c(links_plot$from, links_plot$to)))

collapsed <- sectors[grepl("^collapsed_node_", sectors)]
singletons <- setdiff(sectors, collapsed)

clade_num <- suppressWarnings(as.integer(str_match(collapsed, "^collapsed_node_(\\d+)$")[,2]))
collapsed <- collapsed[order(clade_num)]

sector_order <- c(collapsed, singletons)

# colors
clade_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(sector_order)),
  sector_order
)

# label map
sector_labels <- setNames(
  vapply(sector_order, pretty_label, character(1)),
  sector_order
)

draw_chord <- function() {
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 4,
    track.margin = c(0.01, 0.01),
    canvas.xlim = c(-1.35, 1.35),
    canvas.ylim = c(-1.35, 1.35)
  )
  
  df <- links_plot %>% select(from, to, htt_count)
  
  chordDiagram(
    x = df,
    order = sector_order,
    grid.col = clade_cols,
    transparency = 0.7,
    link.lwd = rescale(df$htt_count, to = c(0.5, 6)),
    annotationTrack = "grid",
    preAllocateTracks = 1
  )
  
  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      s <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(
        x = mean(xlim),
        y = ylim[1] + 0.15,
        labels = sector_labels[[s]],
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.6
      )
    }
  )
  
  circos.clear()
}

pdf(out_pdf, width = 12, height = 12)
draw_chord()
dev.off()

png(out_png, width = 3600, height = 3600, res = 300)
draw_chord()
dev.off()

cat("Saved:\n")
cat("  - ", out_pdf, "\n", sep = "")
cat("  - ", out_png, "\n", sep = "")
cat("  - ", out_tsv, "\n", sep = "")

# --- Summary Tables and Interpretation ----------------------------------------
summary_stats <- data.frame(
  metric = c(
    "Total TE pairs with Ks",
    "HTT candidates (0.5% recomputed)",
    "Non-HTT pairs",
    "HTT percentage",
    "HTT Ks median",
    "HTT Ks min",
    "HTT Ks max",
    "Non-HTT Ks median",
    "Min BUSCO 0.5% threshold",
    "Max BUSCO 0.5% threshold"
  ),
  value = c(
    nrow(te_ks_valid),
    sum(te_ks_valid$is_htt, na.rm = TRUE),
    sum(!te_ks_valid$is_htt, na.rm = TRUE),
    round(100 * sum(te_ks_valid$is_htt, na.rm = TRUE) / nrow(te_ks_valid), 4),
    ifelse(sum(te_ks_valid$is_htt, na.rm = TRUE) > 0,
           round(median(te_ks_valid$ks[te_ks_valid$is_htt == TRUE]), 6), NA),
    ifelse(sum(te_ks_valid$is_htt, na.rm = TRUE) > 0,
           round(min(te_ks_valid$ks[te_ks_valid$is_htt == TRUE]), 6), NA),
    ifelse(sum(te_ks_valid$is_htt, na.rm = TRUE) > 0,
           round(max(te_ks_valid$ks[te_ks_valid$is_htt == TRUE]), 6), NA),
    round(median(te_ks_valid$ks[te_ks_valid$is_htt == FALSE]), 6),
    round(min(busco_valid$ks_0.5_quantile, na.rm = TRUE), 6),
    round(max(busco_valid$ks_0.5_quantile, na.rm = TRUE), 6)
  )
)

write.table(summary_stats, file.path(results_dir, "htt_summary_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

baseline <- quantile_results$htt_count[quantile_results$quantile == "0.5%"]

interpretation <- quantile_results %>%
  mutate(
    fold_change = round(htt_count / baseline, 2),
    additional_htt = htt_count - baseline,
    interpretation = case_when(
      quantile == "0.5%" ~ "Baseline (Romeijn et al. 2025)",
      fold_change < 2 ~ "Minimal increase",
      fold_change < 5 ~ "Moderate increase",
      fold_change < 10 ~ "Substantial increase",
      TRUE ~ "Large increase"
    )
  )

write.table(interpretation %>% select(-quantile_col),
            file.path(results_dir, "htt_quantile_interpretation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

max_fc <- max(interpretation$fold_change)
max_q <- as.character(interpretation$quantile[which.max(interpretation$fold_change)])
max_htt <- max(interpretation$htt_count)
max_pct <- max(interpretation$htt_percentage)

# --- Session info -------------------------------------------------------------
sessionInfo()
