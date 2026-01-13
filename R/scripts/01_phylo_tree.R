# ==============================================================================
# 01_phylo_tree.R
# Phylogenetic tree visualization for Cenococcum geophilum
# ==============================================================================

# --- Libraries ----------------------------------------------------------------
library(ape)
library(ggtree)
library(ggplot2)

# --- Paths --------------------------------------------------------------------
data_dir    <- "data/01_phylo_tree"
results_dir <- "results/01_phylo_tree"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load tree ----------------------------------------------------------------
tree <- read.tree(file.path(data_dir, "phylogeny_rooted.treefile"))

# Simplify tip labels (extract strain ID)
tree$tip.label <- gsub("_AssemblyScaffolds.*", "", tree$tip.label)

# Identify outgroup
outgroup <- "PsefloM405_1"

# --- Fix zero-length root branch for aesthetics -------------------------------
root_node <- length(tree$tip.label) + 1
root_edges <- which(tree$edge[, 1] == root_node)

if (any(tree$edge.length[root_edges] == 0)) {
  # Find outgroup and ingroup edges from root
  outgroup_idx <- which(tree$tip.label == outgroup)
  
  for (i in root_edges) {
    child <- tree$edge[i, 2]
    # Check if this edge leads to outgroup (directly or indirectly)
    if (child == outgroup_idx) {
      # Direct edge to outgroup - split in half
      other_edge <- setdiff(root_edges, i)
      total_length <- tree$edge.length[i] + tree$edge.length[other_edge]
      tree$edge.length[i] <- total_length / 2
      tree$edge.length[other_edge] <- total_length / 2
      break
    }
  }
}

# --- Basic tree plot ----------------------------------------------------------
p <- ggtree(tree, ladderize = TRUE) +
  geom_tiplab(size = 3, fontface = "italic") +
  geom_nodepoint(
    aes(subset = !isTip & as.numeric(label) >= 95),
    size = 2, color = "black", shape = 16
  ) +
  geom_treescale(x = 0, y = -1, width = 0.01, fontsize = 2.5) +
  xlim(0, max(fortify(tree)$x) * 1.4) +
  theme_tree2() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    axis.text.x = element_text(size = 8)
  )

# --- Save plots ---------------------------------------------------------------
# PDF for publication
ggsave(
  file.path(results_dir, "phylo_tree.pdf"),
  plot = p,
  width = 8,
  height = 10,
  units = "in"
)

# PNG for quick viewing
ggsave(
  file.path(results_dir, "phylo_tree.png"),
  plot = p,
  width = 8,
  height = 10,
  units = "in",
  dpi = 300
)

cat("Tree plots saved to:", results_dir, "\n")

# --- Session info -------------------------------------------------------------
sessionInfo()
