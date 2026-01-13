# ==============================================================================
# 02_div_times.R
# Phylogenetic tree visualization with RED values
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(dplyr)

# --- Paths --------------------------------------------------------------------
data_dir    <- "data/02_div_times"
results_dir <- "results/02_div_times"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ----------------------------------------------------------------
# Decorated tree (bootstraps + taxonomy, no RED inside)
tree  <- read.tree(file.path(data_dir, "decorated_tree"))

# PhyloRank stats table (just for sanity / logging)
stats <- read.delim(file.path(data_dir, "decorated_tree-table"), header = TRUE)

cat("Tree loaded:", Ntip(tree), "tips,", Nnode(tree), "internal nodes\n")
cat("Statistics table:\n")
print(stats)

# --- Fix zero-length root branch for aesthetics -------------------------------
outgroup <- "PsefloM405_1_AssemblyScaffolds_2024-02-18"
root_node <- length(tree$tip.label) + 1
root_edges <- which(tree$edge[, 1] == root_node)

if (any(tree$edge.length[root_edges] == 0)) {
  outgroup_idx <- which(tree$tip.label == outgroup)
  
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

# --- Compute RED values on this tree (Parks et al. formula) -------------------
# Node indices in ape: 1..Ntip are tips; (Ntip+1)..(Ntip+Nnode) internal nodes
n_tips  <- Ntip(tree)
n_nodes <- tree$Nnode
n_total <- n_tips + n_nodes

edges    <- tree$edge
edge_len <- tree$edge.length
depth    <- node.depth.edgelength(tree)

# Root = node that never appears as a child
root <- setdiff(edges[, 1], edges[, 2])[1]

# Map: parent -> children
children_of <- split(edges[, 2], edges[, 1])

# RED vector
red <- rep(NA_real_, n_total)
red[root] <- 0  # RED(root) = 0

# helper: mean distance from parent to all descendant tips of a child
mean_parent_to_tips <- function(parent, child, depth, tree) {
  tips <- phangorn::Descendants(tree, child, type = "tips")[[1]]
  if (length(tips) == 0) return(0)
  mean(depth[tips] - depth[parent])
}

# recursive traversal from root downwards
compute_red <- function(node) {
  if (!as.character(node) %in% names(children_of)) return()
  for (child in children_of[[as.character(node)]]) {
    idx <- which(edges[, 1] == node & edges[, 2] == child)
    d   <- edge_len[idx]
    u   <- mean_parent_to_tips(node, child, depth, tree)
    if (u <= 0 || !is.finite(u)) {
      red[child] <<- red[node]
    } else {
      red[child] <<- red[node] + (d / u) * (1 - red[node])
    }
    compute_red(child)
  }
}
compute_red(root)

cat("RED sanity check (tip range): ",
    sprintf("%.4f - %.4f", min(red[1:n_tips]), max(red[1:n_tips])), "\n")

# Save RED table (optional but handy)
red_table <- data.frame(
  node   = 1:n_total,
  label  = c(tree$tip.label,
             rep(NA_character_, n_nodes)),
  red    = red,
  is_tip = c(rep(1L, n_tips), rep(0L, n_nodes))
)
write.table(
  red_table,
  file = file.path(results_dir, "red_values_from_R.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# --- Prepare tree / metadata --------------------------------------------------
# Keep original full tip labels to detect genus
orig_tip_labels <- tree$tip.label

# Simplify labels for plotting (remove big suffix)
tree$tip.label <- gsub("_AssemblyScaffolds.*", "", tree$tip.label)
simple_labels  <- tree$tip.label

# Tip metadata: genus (for color) + tip RED (will be ~1)
tip_meta <- data.frame(
  label     = simple_labels,
  full_name = orig_tip_labels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    genus   = ifelse(grepl("^Pseflo", full_name),
                     "Pseudocenococcum",
                     "Cenococcum"),
    tip_red = red[1:n_tips]
  )

# Base ggtree object for coordinates
base_tree <- ggtree(tree, ladderize = TRUE)
tree_data <- base_tree$data

# Add bootstrap and RED to tree_data (joined by node index)
tree_data$bootstrap <- suppressWarnings(as.numeric(tree_data$label))
red_df <- data.frame(node = 1:length(red), red = red)
tree_data <- left_join(tree_data, red_df, by = "node")

x_max <- max(tree_data$x)

# --- Plots --------------------------------------------------------------------
# Plot 1: Basic tree with bootstrap support
p1 <- ggtree(tree, ladderize = TRUE) +
  geom_tiplab(size = 3, fontface = "italic") +
  geom_nodepoint(
    aes(subset = !isTip & as.numeric(label) >= 95),
    size = 2, color = "black", shape = 16
  ) +
  geom_nodepoint(
    aes(subset = !isTip & as.numeric(label) >= 70 & as.numeric(label) < 95),
    size = 2, color = "grey50", shape = 16
  ) +
  geom_treescale(x = 0, y = -1, width = 0.01, fontsize = 2.5) +
  xlim(0, x_max * 1.4) +
  theme_tree2() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    axis.text.x = element_text(size = 8)
  )

ggsave(
  file.path(results_dir, "phylo_tree_bootstrap.pdf"),
  plot = p1, width = 8, height = 10, units = "in"
)
ggsave(
  file.path(results_dir, "phylo_tree_bootstrap.png"),
  plot = p1, width = 8, height = 10, units = "in", dpi = 300
)

# Plot 2: Tree colored by genus
p2 <- ggtree(tree, ladderize = TRUE) %<+% tip_meta +
  geom_tiplab(aes(color = genus), size = 3, fontface = "italic") +
  geom_tippoint(aes(color = genus), size = 1.5) +
  geom_nodepoint(
    aes(subset = !isTip & as.numeric(label) >= 95),
    size = 2, color = "black", shape = 16
  ) +
  scale_color_manual(
    values = c("Cenococcum" = "#2E86AB", "Pseudocenococcum" = "#A23B72"),
    name = "Genus"
  ) +
  geom_treescale(x = 0, y = -1, width = 0.01, fontsize = 2.5) +
  xlim(0, x_max * 1.4) +
  theme_tree2() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    axis.text.x = element_text(size = 8),
    legend.position = c(0.15, 0.9)
  )

ggsave(
  file.path(results_dir, "phylo_tree_genus.pdf"),
  plot = p2, width = 8, height = 10, units = "in"
)
ggsave(
  file.path(results_dir, "phylo_tree_genus.png"),
  plot = p2, width = 8, height = 10, units = "in", dpi = 300
)

# Plot 3: Tree with internal nodes colored & labeled by RED
internal_data <- tree_data %>% filter(!isTip)

p3 <- ggtree(tree, ladderize = TRUE) +
  geom_tiplab(size = 3, fontface = "italic") +
  # colored circles at internal nodes
  geom_point(
    data = internal_data,
    aes(x = x, y = y, color = red),
    size = 2,
    inherit.aes = FALSE
  ) +
  # numeric RED values next to the circles
  geom_text(
    data = internal_data,
    aes(x = x, y = y, label = sprintf("%.2f", red)),
    hjust = -0.2, vjust = 0.5, size = 2.5,
    inherit.aes = FALSE
  ) +
  scale_color_viridis_c(name = "RED (internal)", limits = c(0, 1)) +
  geom_treescale(x = 0, y = -1, width = 0.01, fontsize = 3) +
  xlim(0, x_max * 1.6) +
  theme_tree2() +
  theme(
    plot.margin = margin(10, 40, 10, 10),
    axis.text.x = element_text(size = 8),
    legend.position = c(0.85, 0.2)
  )

ggsave(
  file.path(results_dir, "phylo_tree_red_internal.pdf"),
  plot = p3, width = 8, height = 10, units = "in"
)
ggsave(
  file.path(results_dir, "phylo_tree_red_internal.png"),
  plot = p3, width = 8, height = 10, units = "in", dpi = 300
)

# --- Session info -------------------------------------------------------------
sessionInfo()
