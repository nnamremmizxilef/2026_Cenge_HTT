# ==============================================================================
# 07_08_htt_clustering_events.R
# HTT clustering and minimal events visualization for Cenococcum geophilum
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(ape)
library(ggtree)
library(ggpubr)
library(jsonlite)

# --- Paths --------------------------------------------------------------------
data_dir    <- "data/07_08_clustering_htt_events"
tree_dir    <- "data/01_phylo_tree"
results_dir <- "results/07_08_clustering_htt_events"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ----------------------------------------------------------------
tree <- read.tree(file.path(tree_dir, "phylogeny_rooted.treefile"))

events_by_node <- read.delim(
  file.path(data_dir, "minimal_events_by_node.tsv"),
  stringsAsFactors = FALSE
)

events_by_component <- read.delim(
  file.path(data_dir, "minimal_events_by_component.tsv"),
  stringsAsFactors = FALSE
)

htt_pairs <- read.delim(
  file.path(data_dir, "htt_pairs_with_community.tsv"),
  stringsAsFactors = FALSE
)

cluster_summary <- read.delim(
  file.path(data_dir, "htt_cluster_summary.tsv"),
  stringsAsFactors = FALSE
)

# Load nodes_info.json to get correct node-to-leaves mapping
nodes_info <- fromJSON(file.path(data_dir, "nodes_info.json"))

cat("Loaded", nrow(htt_pairs), "HTT candidate pairs\n")
cat("Communities:", length(unique(htt_pairs$community_id)), "\n")
cat("Components:", nrow(events_by_component), "\n")
cat("Minimal HTT events:", sum(events_by_node$min_events_total), "\n")

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

# Simplify tip labels
tree$tip.label <- gsub("_1_AssemblyScaffolds.*", "", tree$tip.label)
tree$node.label <- NULL

# --- Map node names to tree nodes using MRCA ----------------------------------
# For each node in nodes_info, find the corresponding tree node by computing
# the MRCA of its leaf set

simplify_leaf_name <- function(leaf_name) {
  # Remove _1_AssemblyScaffolds... suffix
  simple <- gsub("_1_AssemblyScaffolds.*", "", leaf_name)
  return(simple)
}

node_mapping <- data.frame(
  node_name = character(),
  tree_node = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(nodes_info))) {
  node_name <- nodes_info$node_name[i]
  all_leaves <- nodes_info$all_leaves[[i]]
  
  # Simplify leaf names to match tree tip labels
  simple_leaves <- sapply(all_leaves, simplify_leaf_name)
  
  # Find tip indices in tree
  tip_indices <- match(simple_leaves, tree$tip.label)
  tip_indices <- tip_indices[!is.na(tip_indices)]
  
  if (length(tip_indices) >= 2) {
    # Find MRCA of these tips
    tree_node <- getMRCA(tree, tip_indices)
  } else if (length(tip_indices) == 1) {
    tree_node <- tip_indices[1]
  } else {
    tree_node <- NA
  }
  
  node_mapping <- rbind(node_mapping, data.frame(
    node_name = node_name,
    tree_node = tree_node,
    stringsAsFactors = FALSE
  ))
}

cat("\nNode mapping (nodes_info -> tree):\n")
print(node_mapping)

# Join with events data
node_mapping <- node_mapping %>%
  left_join(events_by_node, by = c("node_name" = "node")) %>%
  mutate(min_events_total = ifelse(is.na(min_events_total), 0, min_events_total))

# --- Prepare tree data --------------------------------------------------------
tree_data <- fortify(tree)

node_annot_df <- node_mapping %>%
  filter(!is.na(tree_node)) %>%
  select(tree_node, min_events_total, node_name) %>%
  rename(node = tree_node, min_events = min_events_total)

tree_data <- left_join(tree_data, node_annot_df, by = "node")

internal_data <- tree_data %>% 
  filter(!isTip & !is.na(min_events) & min_events > 0)

cat("\nNodes with HTT events:\n")
print(internal_data %>% select(node, label, min_events, node_name))

x_max <- max(tree_data$x)

# --- Plot 1: Tree with HTT events at nodes ------------------------------------
p1 <- ggtree(tree, ladderize = TRUE) +
  geom_tiplab(size = 3, fontface = "italic") +
  geom_point(
    data = internal_data,
    aes(x = x, y = y, size = min_events),
    color = "#E63946",
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = internal_data,
    aes(x = x, y = y, label = min_events),
    hjust = 0.5, vjust = 0.5, size = 3, fontface = "bold",
    color = "white",
    inherit.aes = FALSE
  ) +
  scale_size_continuous(
    range = c(5, 12),
    breaks = c(1, 3, 5, 7),
    name = "HTT events"
  ) +
  geom_treescale(x = 0, y = -1, width = 0.01, fontsize = 3) +
  xlim(0, x_max * 1.6) +
  theme_tree2() +
  theme(
    plot.margin = margin(10, 40, 10, 10),
    axis.text.x = element_text(size = 8),
    legend.position = c(0.85, 0.2)
  )

ggsave(file.path(results_dir, "htt_events_tree.pdf"), p1, width = 8, height = 10)
ggsave(file.path(results_dir, "htt_events_tree.png"), p1, width = 8, height = 10, dpi = 300)

# --- Summary ------------------------------------------------------------------
cat("\n")
cat("=== HTT Summary ===\n")
cat("Total HTT candidate pairs:", nrow(htt_pairs), "\n")
cat("Total communities:", length(unique(htt_pairs$community_id)), "\n")
cat("Total components:", nrow(events_by_component), "\n")
cat("Minimal HTT events:", sum(events_by_node$min_events_total), "\n")
cat("\n")

cat("Events by node:\n")
print(events_by_node %>% arrange(desc(min_events_total)))

cat("\n")
cat("Top clade pairs by HTT candidates:\n")
print(cluster_summary %>% 
        arrange(desc(n_candidates)) %>% 
        head(10) %>%
        select(cladeA, cladeB, n_candidates, n_clusters))

# --- Save outputs -------------------------------------------------------------
write.table(events_by_node,
            file.path(results_dir, "events_by_node.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(events_by_component,
            file.path(results_dir, "events_by_component.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Output -------------------------------------------------------------------
cat("\n")
cat("Plot saved to:", results_dir, "\n")
cat("  - htt_events_tree.pdf/png\n")
cat("\n")
cat("Tables saved to:", results_dir, "\n")
cat("  - events_by_node.tsv\n")
cat("  - events_by_component.tsv\n")

# --- Session info -------------------------------------------------------------
sessionInfo()

