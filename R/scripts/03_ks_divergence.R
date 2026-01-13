# ==============================================================================
# 03_ks_divergence.R
# Ks divergence analysis visualization for Cenococcum geophilum
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(jsonlite)
library(ape)
library(ggtree)
library(ggpubr)

# --- Paths --------------------------------------------------------------------
data_dir    <- "data/03_ks_divergence"
tree_dir    <- "data/01_phylo_tree"
results_dir <- "results/03_ks_divergence"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Parameters ---------------------------------------------------------------
KS_THRESHOLD <- 0.01  # primary threshold for HTT detection

# --- Load data ----------------------------------------------------------------
ks_raw <- read.delim(file.path(data_dir, "ks_values_by_node.tsv"), 
                     stringsAsFactors = FALSE)
node_summary <- read.delim(file.path(data_dir, "node_ks_summary.tsv"), 
                           stringsAsFactors = FALSE)
nodes_info <- fromJSON(file.path(data_dir, "nodes_info.json"))
tree <- read.tree(file.path(tree_dir, "phylogeny_rooted.treefile"))

# --- Fix zero-length root branch for aesthetics -------------------------------
outgroup_full <- "PsefloM405_1_AssemblyScaffolds_2024-02-18"
root_node <- length(tree$tip.label) + 1
root_edges <- which(tree$edge[, 1] == root_node)

if (any(tree$edge.length[root_edges] == 0)) {
  outgroup_idx <- which(tree$tip.label == outgroup_full)
  
  for (i in root_edges) {
    child <- tree$edge[i, 2]
    if (child == outgroup_idx) {
      other_edge <- setdiff(root_edges, i)
      total_length <- tree$edge.length[i] + tree$edge.length[other_edge]
      tree$edge.length[i] <- total_length / 2
      tree$edge.length[other_edge] <- total_length / 2
      break
    }
  }
}

# Simplify tip labels - remove _1_AssemblyScaffolds... to get just Cenge1005
tree$tip.label <- gsub("_1_AssemblyScaffolds.*", "", tree$tip.label)

# Remove bootstrap labels
tree$node.label <- NULL

cat("Loaded", nrow(ks_raw), "Ks values across", length(unique(ks_raw$node)), "nodes\n")

cat("Node summary:", nrow(node_summary), "nodes with sufficient data\n")
cat("Using Ks 0.5% quantile threshold:", KS_THRESHOLD, "\n")

# --- Create node description --------------------------------------------------
# Use the same pattern as tree tip labels for consistency
simplify_name <- function(x) gsub("_1_AssemblyScaffolds.*", "", x)

nodes_info$clade1_samples <- sapply(nodes_info$daughter_clades, function(x) {
  paste(simplify_name(x[[1]]), collapse = ", ")
})
nodes_info$clade2_samples <- sapply(nodes_info$daughter_clades, function(x) {
  paste(simplify_name(x[[2]]), collapse = ", ")
})
nodes_info$n_samples <- sapply(nodes_info$all_leaves, length)
nodes_info$all_leaves_clean <- lapply(nodes_info$all_leaves, simplify_name)

nodes_info$description <- sapply(seq_len(nrow(nodes_info)), function(i) {
  c1 <- simplify_name(nodes_info$daughter_clades[[i]][[1]])
  c2 <- simplify_name(nodes_info$daughter_clades[[i]][[2]])
  paste0(c1[1], " vs ", c2[1], 
         ifelse(length(c1) > 1 | length(c2) > 1, " (+more)", ""))
})

node_summary_annotated <- node_summary %>%
  left_join(
    nodes_info %>% select(node_name, description, n_samples, clade1_samples, clade2_samples, all_leaves_clean),
    by = c("node" = "node_name")
  ) %>%
  mutate(pass_threshold = ks_0.5_quantile > KS_THRESHOLD)

# --- Collapse nodes logic -----------------------------------------------------
# Identify failing nodes and determine collapsed clades

failing_nodes <- node_summary_annotated %>%
  filter(!pass_threshold) %>%
  arrange(n_samples)  # process smallest nodes first

# Get all sample names
all_samples <- unique(tree$tip.label)

# Initialize: each sample starts as its own clade
sample_to_clade <- setNames(all_samples, all_samples)
clade_members <- as.list(setNames(all_samples, all_samples))

# Process failing nodes from smallest to largest
# Samples within a failing node get merged into one clade
for (i in seq_len(nrow(failing_nodes))) {
  node_row <- failing_nodes[i, ]
  node_samples <- node_row$all_leaves_clean[[1]]
  
  # Find current clade assignments for these samples
  current_clades <- unique(sample_to_clade[node_samples])
  
  # Create new merged clade name
  new_clade_name <- paste0("collapsed_", node_row$node)
  
  # Collect all samples from all clades being merged
  samples_to_merge <- unique(unlist(lapply(current_clades, function(c) {
    if (c %in% names(clade_members)) {
      clade_members[[c]]
    } else {
      c
    }
  })))
  
  # Update mappings
  for (s in samples_to_merge) {
    sample_to_clade[s] <- new_clade_name
  }
  
  # Remove old clades and add new one
  clade_members[current_clades] <- NULL
  clade_members[[new_clade_name]] <- samples_to_merge
}

# Create final clade assignments
clade_assignments <- data.frame(
  sample = names(sample_to_clade),
  clade = unname(sample_to_clade),
  stringsAsFactors = FALSE
) %>%
  arrange(clade, sample)

# Summarize clades
clade_summary <- clade_assignments %>%
  group_by(clade) %>%
  summarise(
    n_samples = n(),
    samples = paste(sample, collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    is_collapsed = grepl("^collapsed_", clade),
    clade_type = ifelse(is_collapsed, "Collapsed", "Single sample")
  ) %>%
  arrange(desc(n_samples), clade)

# --- Print collapsing results -------------------------------------------------
cat("\n")
cat("=== Node Collapsing Results ===\n")
cat("Threshold: Ks 0.5% quantile >", KS_THRESHOLD, "\n\n")

cat("Failing nodes (will be collapsed):\n")
for (i in seq_len(nrow(failing_nodes))) {
  row <- failing_nodes[i, ]
  cat("  ", row$node, ": ", row$n_samples, " samples, Ks 0.5% = ", 
      round(row$ks_0.5_quantile, 4), "\n", sep = "")
}

cat("\n")
cat("Total samples:", length(all_samples), "\n")
cat("Total clades after collapsing:", nrow(clade_summary), "\n")
cat("  - Collapsed clades:", sum(clade_summary$is_collapsed), "\n")
cat("  - Single-sample clades:", sum(!clade_summary$is_collapsed), "\n")
cat("\n")

cat("Collapsed clades:\n")
collapsed_clades <- clade_summary %>% filter(is_collapsed)
for (i in seq_len(nrow(collapsed_clades))) {
  row <- collapsed_clades[i, ]
  cat("  ", row$clade, " (", row$n_samples, " samples): ", row$samples, "\n", sep = "")
}

# --- Print node mapping -------------------------------------------------------
cat("\n=== Node to Clade Mapping ===\n\n")
for (i in seq_len(nrow(node_summary_annotated))) {
  row <- node_summary_annotated[i, ]
  cat(row$node, " (", row$n_samples, " samples, Ks 0.5% = ", 
      round(row$ks_0.5_quantile, 4), "):\n", sep = "")
  cat("  Clade 1:", substr(row$clade1_samples, 1, 80), 
      ifelse(nchar(row$clade1_samples) > 80, "...", ""), "\n")
  cat("  Clade 2:", substr(row$clade2_samples, 1, 80),
      ifelse(nchar(row$clade2_samples) > 80, "...", ""), "\n")
  cat("  Pass threshold:", row$pass_threshold, "\n\n")
}

# --- Filter valid Ks values ---------------------------------------------------
ks_valid <- ks_raw %>%
  filter(is.finite(ks) & ks > 0 & ks < 10)

cat("Valid Ks values:", nrow(ks_valid), "/", nrow(ks_raw), "\n")

# --- Plots --------------------------------------------------------------------
# Plot 1: Overall Ks distribution
p1 <- ggplot(ks_valid, aes(x = ks)) +
  geom_histogram(bins = 50, fill = "#2E86AB", color = "white", alpha = 0.8) +
  geom_vline(xintercept = KS_THRESHOLD, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
  scale_x_log10(
    breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1),
    labels = c("0.001", "0.005", "0.01", "0.05", "0.1", "0.5", "1")
  ) +
  labs(
    x = "Ks (synonymous substitutions per site)",
    y = "Count"
  ) +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "ks_distribution_overall.pdf"), p1, width = 8, height = 5)
ggsave(file.path(results_dir, "ks_distribution_overall.png"), p1, width = 8, height = 5, dpi = 300)

# Plot 2: Ks boxplot per node
node_order <- ks_valid %>%
  group_by(node) %>%
  summarise(median_ks = median(ks)) %>%
  arrange(median_ks) %>%
  pull(node)

node_labels <- node_summary_annotated %>%
  mutate(label = paste0(node, "\n(", n_samples, " samples)")) %>%
  select(node, label, pass_threshold)

ks_valid_labeled <- ks_valid %>%
  left_join(node_labels, by = "node") %>%
  mutate(label = factor(label, levels = node_labels$label[match(node_order, node_labels$node)]))

p2 <- ggplot(ks_valid_labeled, aes(x = label, y = ks, fill = pass_threshold)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  geom_hline(yintercept = KS_THRESHOLD, linetype = "dashed", color = "#E63946", linewidth = 0.6) +
  scale_y_log10(
    breaks = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1),
    labels = c("0.001", "0.005", "0.01", "0.05", "0.1", "0.5", "1")
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#2E86AB", "FALSE" = "#CCCCCC"),
    name = paste0("Pass Ks > ", KS_THRESHOLD),
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  labs(
    x = "Node",
    y = "Ks (log scale)"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "ks_boxplot_per_node.pdf"), p2, width = 12, height = 6)
ggsave(file.path(results_dir, "ks_boxplot_per_node.png"), p2, width = 12, height = 6, dpi = 300)

# Plot 3: Node summary - 0.5% quantile
node_summary_ordered <- node_summary_annotated %>%
  arrange(ks_0.5_quantile) %>%
  mutate(node = factor(node, levels = node))

p3 <- ggplot(node_summary_ordered, aes(x = node, y = ks_0.5_quantile, fill = pass_threshold)) +
  geom_col(alpha = 0.9) +
  geom_hline(yintercept = KS_THRESHOLD, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
  scale_fill_manual(
    values = c("TRUE" = "#2E86AB", "FALSE" = "#CCCCCC"),
    name = paste0("Pass Ks > ", KS_THRESHOLD),
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  labs(
    x = "Node",
    y = "Ks 0.5% quantile"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "ks_quantile_per_node.pdf"), p3, width = 10, height = 6)
ggsave(file.path(results_dir, "ks_quantile_per_node.png"), p3, width = 10, height = 6, dpi = 300)

# Plot 4: Ka/Ks distribution
ks_valid_kaks <- ks_valid %>%
  filter(is.finite(ka) & ka > 0) %>%
  mutate(ka_ks = ka / ks)

p4 <- ggplot(ks_valid_kaks, aes(x = ka_ks)) +
  geom_histogram(bins = 50, fill = "#2A9D8F", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
  scale_x_log10() +
  labs(
    x = "Ka/Ks (log scale)",
    y = "Count"
  ) +
  theme_pubr() +
  theme(
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "ka_ks_distribution.pdf"), p4, width = 8, height = 5)
ggsave(file.path(results_dir, "ka_ks_distribution.png"), p4, width = 8, height = 5, dpi = 300)

# Plot 5: Tree with Ks values at nodes
tree_data <- fortify(tree)

find_mrca_node <- function(tree, samples) {
  samples_clean <- simplify_name(samples)
  tip_indices <- which(tree$tip.label %in% samples_clean)
  if (length(tip_indices) < 2) return(NA)
  mrca_node <- getMRCA(tree, tip_indices)
  return(mrca_node)
}

node_mapping <- data.frame(
  analysis_node = nodes_info$node_name,
  stringsAsFactors = FALSE
)
node_mapping$tree_node <- sapply(seq_len(nrow(nodes_info)), function(i) {
  find_mrca_node(tree, nodes_info$all_leaves[[i]])
})

node_mapping <- node_mapping %>%
  left_join(node_summary_annotated %>% select(node, ks_0.5_quantile, pass_threshold), 
            by = c("analysis_node" = "node")) %>%
  filter(!is.na(tree_node))

# Create node annotation data frame
node_annot_df <- data.frame(
  node = node_mapping$tree_node,
  ks_quantile = node_mapping$ks_0.5_quantile,
  pass = node_mapping$pass_threshold,
  node_name = node_mapping$analysis_node
)

# Merge with tree_data for plotting
tree_data <- left_join(tree_data, node_annot_df, by = "node")

# Get internal node data
internal_data <- tree_data %>% filter(!isTip & !is.na(ks_quantile))

x_max <- max(tree_data$x)

p5 <- ggtree(tree, ladderize = TRUE) +
  geom_tiplab(size = 3, fontface = "italic") +
  # colored circles at internal nodes
  geom_point(
    data = internal_data,
    aes(x = x, y = y, color = ks_quantile, shape = pass),
    size = 3,
    inherit.aes = FALSE
  ) +
  # node names next to the circles
  geom_text(
    data = internal_data,
    aes(x = x, y = y, label = node_name),
    hjust = -0.2, vjust = 0.5, size = 2.5,
    inherit.aes = FALSE
  ) +
  scale_color_viridis_c(
    name = "Ks 0.5%\nquantile",
    option = "plasma",
    trans = "log10"
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    name = paste0("Pass ", KS_THRESHOLD),
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  geom_treescale(x = 0, y = -1, width = 0.01, fontsize = 3) +
  xlim(0, x_max * 1.6) +
  theme_tree2() +
  theme(
    plot.margin = margin(10, 40, 10, 10),
    axis.text.x = element_text(size = 8),
    legend.position = c(0.85, 0.2)
  )

ggsave(file.path(results_dir, "tree_with_ks_nodes.pdf"), p5, width = 8, height = 10)
ggsave(file.path(results_dir, "tree_with_ks_nodes.png"), p5, width = 8, height = 10, dpi = 300)

# Plot 6: Tree with collapsed clades
# Add clade assignment to tip data
tip_clade_data <- clade_assignments %>%
  mutate(is_collapsed = grepl("^collapsed_", clade))

p6 <- ggtree(tree, ladderize = TRUE) %<+% tip_clade_data +
  geom_tiplab(aes(color = is_collapsed), size = 2.5, fontface = "italic") +
  geom_tippoint(aes(color = is_collapsed), size = 2) +
  scale_color_manual(
    values = c("TRUE" = "#E63946", "FALSE" = "#2E86AB"),
    name = "Collapsed",
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  geom_treescale(x = 0, y = -1, width = 0.01, fontsize = 2.5) +
  xlim(0, max(tree_data$x) * 1.5) +
  theme_tree2() +
  theme(
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(results_dir, "tree_collapsed_clades.pdf"), p6, width = 12, height = 10)
ggsave(file.path(results_dir, "tree_collapsed_clades.png"), p6, width = 12, height = 10, dpi = 300)

# Plot 7: Summary statistics comparison
node_summary_long <- node_summary %>%
  select(node, ks_min, ks_0.5_quantile, ks_5_quantile, ks_median, ks_mean) %>%
  pivot_longer(cols = -node, names_to = "statistic", values_to = "value") %>%
  mutate(
    statistic = factor(statistic, 
                       levels = c("ks_min", "ks_0.5_quantile", "ks_5_quantile", "ks_median", "ks_mean"),
                       labels = c("Min", "0.5% quantile", "5% quantile", "Median", "Mean"))
  )

node_order_summary <- node_summary %>%
  arrange(ks_0.5_quantile) %>%
  pull(node)
node_summary_long$node <- factor(node_summary_long$node, levels = node_order_summary)

p7 <- ggplot(node_summary_long, aes(x = node, y = value, color = statistic, group = statistic)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = KS_THRESHOLD, linetype = "dashed", color = "#E63946", linewidth = 0.5) +
  scale_y_log10() +
  scale_color_manual(values = c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51")) +
  labs(
    x = "Node",
    y = "Ks (log scale)",
    color = "Statistic"
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(results_dir, "ks_summary_statistics.pdf"), p7, width = 12, height = 6)
ggsave(file.path(results_dir, "ks_summary_statistics.png"), p7, width = 12, height = 6, dpi = 300)

# --- Summary ------------------------------------------------------------------
cat("\n")
cat("=== Summary ===\n")
cat("Threshold used: Ks 0.5% quantile >", KS_THRESHOLD, "\n")
cat("Nodes passing threshold:", sum(node_summary_annotated$pass_threshold), "/", 
    nrow(node_summary_annotated), "\n\n")

cat("Nodes with sufficient divergence (Ks 0.5% >", KS_THRESHOLD, "):\n")
passing_nodes <- node_summary_annotated %>%
  filter(pass_threshold) %>%
  arrange(desc(ks_0.5_quantile)) %>%
  select(node, n_samples, ks_0.5_quantile, description) %>%
  as.data.frame()
print(passing_nodes)

cat("\nNodes with insufficient divergence (Ks 0.5% <=", KS_THRESHOLD, "):\n")
failing_nodes_print <- node_summary_annotated %>%
  filter(!pass_threshold) %>%
  arrange(ks_0.5_quantile) %>%
  select(node, n_samples, ks_0.5_quantile, description) %>%
  as.data.frame()
print(failing_nodes_print)

cat("\n=== Clade Assignments for HTT Analysis ===\n")
cat("Total clades:", nrow(clade_summary), "\n")
cat("  - Collapsed clades:", sum(clade_summary$is_collapsed), 
    "(", sum(clade_summary$n_samples[clade_summary$is_collapsed]), "samples )\n")
cat("  - Single-sample clades:", sum(!clade_summary$is_collapsed), "\n\n")

# --- Save outputs -------------------------------------------------------------
write.table(node_summary_annotated %>% select(-all_leaves_clean), 
            file.path(results_dir, "node_summary_annotated.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(clade_assignments, 
            file.path(results_dir, "clade_assignments.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(clade_summary,
            file.path(results_dir, "clade_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Session info -------------------------------------------------------------
sessionInfo()
