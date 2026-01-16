# ==============================================================================
# HTT contribution to mobilome size analysis
# ==============================================================================
# Quantifies what fraction of each genome's TE content is HTT-derived
# and tests correlation with genome size
# Uses FILTERED HTT set (0% flank synteny, >75% TE identity)
#
# Updates:
#  - Genome sizes loaded from data/genome_sizes.tsv (NOT from highLevelCount files)
#  - Exclude any genomes assigned to collapsed_node_* in data/clade_assignments.tsv
#  - Still excludes OUTGROUP (PsefloM405)
#
# NOTE:
#  - Avoid using a variable named `data` (conflicts with base::data()).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
data <- "data"
out_dir  <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Top-level data folder (for genome_sizes.tsv + clade_assignments.tsv)
data_root <- "data"

# IMPORTANT: use strictly filtered HTT set
htt_file <- file.path(data_dir, "htt_candidates_flank_0pct_teid_gt75.tsv")
stopifnot(file.exists(htt_file))

# Outgroup to exclude
OUTGROUP <- "PsefloM405"

# ------------------------------------------------------------------------------
# 0. Load clade assignments and define genomes to exclude (collapsed nodes)
# ------------------------------------------------------------------------------
cat("Loading clade assignments...\n")
clade_file <- file.path(data_root, "clade_assignment.tsv")
stopifnot(file.exists(clade_file))

clades <- read_tsv(clade_file, show_col_types = FALSE)

need_clade_cols <- c("sample", "clade")
miss_clade <- setdiff(need_clade_cols, colnames(clades))
if (length(miss_clade) > 0) {
  stop("Clade file is missing columns: ", paste(miss_clade, collapse = ", "))
}

collapsed_genomes <- clades %>%
  filter(str_detect(clade, "^collapsed_node_")) %>%
  pull(sample) %>%
  unique()

exclude_genomes <- sort(unique(c(OUTGROUP, collapsed_genomes)))

cat("  Outgroup excluded:", OUTGROUP, "\n")
cat("  Collapsed-node genomes excluded:", length(collapsed_genomes), "\n")
if (length(collapsed_genomes) > 0) {
  cat("  Collapsed-node genomes:", paste(collapsed_genomes, collapse = ", "), "\n")
}
cat("  Total excluded genomes:", length(exclude_genomes), "\n\n")

# ------------------------------------------------------------------------------
# 1. Load HTT candidates (filtered)
# ------------------------------------------------------------------------------
cat("Loading FILTERED HTT candidates...\n")
cat("  File:", htt_file, "\n")
htt <- read_tsv(htt_file, show_col_types = FALSE)

# Check required columns
need_cols <- c("qseqid", "sseqid", "qgenome", "sgenome")
miss <- setdiff(need_cols, colnames(htt))
if (length(miss) > 0) stop("Filtered HTT file is missing columns: ", paste(miss, collapse = ", "))

# Exclude outgroup + collapsed-node genomes
htt <- htt %>%
  filter(
    !qgenome %in% exclude_genomes,
    !sgenome %in% exclude_genomes
  )

cat("  Excluded genomes:", paste(exclude_genomes, collapse = ", "), "\n")
cat("  HTT pairs after filtering:", nrow(htt), "\n")

# Extract unique TE copies involved in HTT (both query and subject)
htt_tes <- bind_rows(
  htt %>% select(te_id = qseqid, genome = qgenome),
  htt %>% select(te_id = sseqid, genome = sgenome)
) %>%
  distinct() %>%
  filter(!genome %in% exclude_genomes)

cat("  Unique TE copies involved in filtered HTT:", nrow(htt_tes), "\n")

# ------------------------------------------------------------------------------
# 2. Parse TE coordinates from IDs to get bp per HTT TE
# ------------------------------------------------------------------------------
# Format: genome|family::contig:start-end::contig:start-end(strand)

parse_te_bp <- function(te_id) {
  match <- str_match(te_id, "::([^:]+):(\\d+)-(\\d+)")
  if (is.na(match[1])) return(NA_integer_)
  start <- as.integer(match[3])
  end <- as.integer(match[4])
  end - start + 1
}

htt_tes <- htt_tes %>%
  mutate(te_bp = map_int(te_id, parse_te_bp)) %>%
  filter(!is.na(te_bp))

# Summarize HTT bp per genome
htt_per_genome <- htt_tes %>%
  group_by(genome) %>%
  summarise(
    htt_te_count = n(),
    htt_te_bp = sum(te_bp),
    .groups = "drop"
  )

cat("  HTT TEs with valid coordinates:", nrow(htt_tes), "\n")

# ------------------------------------------------------------------------------
# 3. Load total TE content per genome from BED files
# ------------------------------------------------------------------------------
cat("\nLoading TE BED files...\n")

bed_files <- list.files(data_dir, pattern = "\\.filteredRepeats\\.bed$", full.names = TRUE)

# Exclude outgroup + collapsed-node genomes by filename-derived genome ID
bed_genomes <- str_remove(basename(bed_files), "\\.filteredRepeats\\.bed$")
bed_files <- bed_files[!bed_genomes %in% exclude_genomes]

cat("  Found", length(bed_files), "BED files (excluding outgroup + collapsed nodes)\n")

parse_te_bed <- function(filepath) {
  genome <- str_remove(basename(filepath), "\\.filteredRepeats\\.bed$")
  
  bed <- tryCatch({
    read_tsv(
      filepath,
      col_names = c("chrom", "start", "end", "name", "score", "strand"),
      show_col_types = FALSE,
      col_types = "ciicdc"
    )
  }, error = function(e) {
    tryCatch({
      read_tsv(
        filepath,
        col_names = c("chrom", "start", "end", "name"),
        show_col_types = FALSE,
        col_types = "ciic"
      )
    }, error = function(e2) NULL)
  })
  
  if (is.null(bed) || nrow(bed) == 0) return(NULL)
  
  bed %>%
    mutate(
      genome = genome,
      te_bp = end - start
    )
}

all_tes <- map_dfr(bed_files, parse_te_bed)

if (nrow(all_tes) == 0) {
  stop("No TE data loaded. Check BED file paths and format.")
}

# Summarize total mobilome per genome
mobilome_per_genome <- all_tes %>%
  group_by(genome) %>%
  summarise(
    total_te_count = n(),
    total_te_bp = sum(te_bp),
    .groups = "drop"
  )

cat("  Loaded TEs from", n_distinct(all_tes$genome), "genomes\n")
cat("  Total TE copies:", nrow(all_tes), "\n")

# ------------------------------------------------------------------------------
# 4. Load genome sizes from data/genome_sizes.tsv
# ------------------------------------------------------------------------------
cat("\nLoading genome sizes from genome_sizes.tsv...\n")

genome_sizes_file <- file.path(data_root, "genome_sizes.tsv")
stopifnot(file.exists(genome_sizes_file))

genome_sizes <- read_tsv(genome_sizes_file, show_col_types = FALSE) %>%
  rename(genome_size = any_of(c("genome_bp", "genome_size"))) %>%
  select(genome, genome_size) %>%
  distinct() %>%
  filter(!genome %in% exclude_genomes)

# sanity checks
if (!all(c("genome", "genome_size") %in% colnames(genome_sizes))) {
  stop("genome_sizes.tsv must contain columns: genome and genome_bp (or genome_size)")
}
if (any(is.na(genome_sizes$genome_size))) stop("Missing genome_size values in genome_sizes.tsv")
if (any(genome_sizes$genome_size <= 0)) stop("Non-positive genome_size values in genome_sizes.tsv")

cat("  Genome sizes loaded for", nrow(genome_sizes), "genomes (excluding outgroup + collapsed nodes)\n")

# ------------------------------------------------------------------------------
# 5. Merge and calculate HTT contribution
# ------------------------------------------------------------------------------
cat("\nCalculating HTT contributions...\n")

results <- mobilome_per_genome %>%
  left_join(htt_per_genome, by = "genome") %>%
  left_join(genome_sizes, by = "genome") %>%
  mutate(
    htt_te_count = replace_na(htt_te_count, 0),
    htt_te_bp = replace_na(htt_te_bp, 0),
    htt_prop_count = htt_te_count / total_te_count,
    htt_prop_bp = htt_te_bp / total_te_bp,
    te_prop_genome = total_te_bp / genome_size
  )

# fail fast if any genome size is missing after join
missing_sizes <- results %>% filter(is.na(genome_size)) %>% pull(genome)
if (length(missing_sizes) > 0) {
  stop(
    "Missing genome sizes for: ", paste(missing_sizes, collapse = ", "),
    "\nCheck that genome IDs match between BED filenames and genome_sizes.tsv."
  )
}

cat("  Genomes in final dataset:", nrow(results), "\n")

# ------------------------------------------------------------------------------
# 6. Summary statistics
# ------------------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("SUMMARY (using filtered HTT set)\n")
cat(strrep("=", 60), "\n\n", sep = "")

cat("Genomes analyzed:", nrow(results), "\n\n")

cat("HTT-derived TE content (% of mobilome):\n")
cat("  Mean:   ", round(mean(results$htt_prop_bp) * 100, 2), "%\n")
cat("  Median: ", round(median(results$htt_prop_bp) * 100, 2), "%\n")
cat("  Range:  ", round(min(results$htt_prop_bp) * 100, 2), " - ",
    round(max(results$htt_prop_bp) * 100, 2), "%\n\n")

cat("HTT-derived TE content (Mb):\n")
cat("  Mean:   ", round(mean(results$htt_te_bp) / 1e6, 2), " Mb\n")
cat("  Median: ", round(median(results$htt_te_bp) / 1e6, 2), " Mb\n")
cat("  Range:  ", round(min(results$htt_te_bp) / 1e6, 2), " - ",
    round(max(results$htt_te_bp) / 1e6, 2), " Mb\n")

# ------------------------------------------------------------------------------
# 7. Statistical tests
# ------------------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("STATISTICAL TESTS\n")
cat(strrep("=", 60), "\n\n", sep = "")

# Correlation: HTT bp vs genome size
cor_genome <- cor.test(results$htt_te_bp, results$genome_size, method = "spearman")
cat("HTT bp vs genome size (Spearman):\n")
cat("  rho =", round(cor_genome$estimate, 3), "\n")
cat("  p   =", format.pval(cor_genome$p.value), "\n\n")

# Correlation: HTT bp vs total mobilome bp
cor_mob <- cor.test(results$htt_te_bp, results$total_te_bp, method = "spearman")
cat("HTT bp vs total mobilome bp (Spearman):\n")
cat("  rho =", round(cor_mob$estimate, 3), "\n")
cat("  p   =", format.pval(cor_mob$p.value), "\n\n")

# Correlation: HTT proportion vs TE content
cor_te <- cor.test(results$htt_prop_bp, results$te_prop_genome, method = "spearman")
cat("HTT proportion vs TE content (Spearman):\n")
cat("  rho =", round(cor_te$estimate, 3), "\n")
cat("  p   =", format.pval(cor_te$p.value), "\n")

# ------------------------------------------------------------------------------
# 8. Plots
# ------------------------------------------------------------------------------
cat("\nGenerating plots...\n")

cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")

# Plot 1: HTT proportion per genome (bar plot)
p1 <- ggplot(results, aes(x = reorder(genome, -htt_prop_bp), y = htt_prop_bp * 100)) +
  geom_col(fill = cols[1], alpha = 0.8) +
  geom_hline(yintercept = mean(results$htt_prop_bp) * 100,
             linetype = "dashed", color = "gray40") +
  labs(
    x = NULL,
    y = "HTT-derived TE content (%)",
    title = "Proportion of mobilome derived from horizontal transfer"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(out_dir, "htt_proportion_barplot.png"), p1,
       width = 12, height = 6, dpi = 300)
ggsave(file.path(out_dir, "htt_proportion_barplot.pdf"), p1,
       width = 12, height = 6)

# Plot 2: HTT bp vs genome size
p2 <- ggplot(results, aes(x = genome_size / 1e6, y = htt_te_bp / 1e6)) +
  geom_point(size = 3, alpha = 0.7, color = cols[1]) +
  geom_smooth(method = "lm", se = TRUE, color = cols[4], fill = cols[4], alpha = 0.2) +
  geom_text_repel(aes(label = genome), size = 3, max.overlaps = 15) +
  labs(
    x = "Genome size (Mb)",
    y = "HTT-derived TE content (Mb)",
    title = "HTT contribution vs genome size",
    subtitle = sprintf("Spearman ρ = %.2f, p = %s",
                       cor_genome$estimate, format.pval(cor_genome$p.value))
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "htt_vs_genome_size.png"), p2,
       width = 8, height = 7, dpi = 300)
ggsave(file.path(out_dir, "htt_vs_genome_size.pdf"), p2,
       width = 8, height = 7)

# Plot 3: HTT bp vs total mobilome
p3 <- ggplot(results, aes(x = total_te_bp / 1e6, y = htt_te_bp / 1e6)) +
  geom_point(size = 3, alpha = 0.7, color = cols[2]) +
  geom_smooth(method = "lm", se = TRUE, color = cols[4], fill = cols[4], alpha = 0.2) +
  geom_text_repel(aes(label = genome), size = 3, max.overlaps = 15) +
  labs(
    x = "Total mobilome (Mb)",
    y = "HTT-derived TE content (Mb)",
    title = "HTT contribution vs total TE content",
    subtitle = sprintf("Spearman ρ = %.2f, p = %s",
                       cor_mob$estimate, format.pval(cor_mob$p.value))
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "htt_vs_mobilome.png"), p3,
       width = 8, height = 7, dpi = 300)
ggsave(file.path(out_dir, "htt_vs_mobilome.pdf"), p3,
       width = 8, height = 7)

# Plot 4: HTT proportion vs TE content
p4 <- ggplot(results, aes(x = te_prop_genome * 100, y = htt_prop_bp * 100)) +
  geom_point(size = 3, alpha = 0.7, color = cols[3]) +
  geom_smooth(method = "lm", se = TRUE, color = cols[4], fill = cols[4], alpha = 0.2) +
  geom_text_repel(aes(label = genome), size = 3, max.overlaps = 15) +
  labs(
    x = "TE content (% of genome)",
    y = "HTT-derived (% of mobilome)",
    title = "HTT proportion vs overall TE content",
    subtitle = sprintf("Spearman ρ = %.2f, p = %s",
                       cor_te$estimate, format.pval(cor_te$p.value))
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "htt_prop_vs_te_content.png"), p4,
       width = 8, height = 7, dpi = 300)
ggsave(file.path(out_dir, "htt_prop_vs_te_content.pdf"), p4,
       width = 8, height = 7)

# ------------------------------------------------------------------------------
# 9. Save results
# ------------------------------------------------------------------------------
results_out <- results %>%
  mutate(
    htt_prop_bp_pct = round(htt_prop_bp * 100, 2),
    htt_prop_count_pct = round(htt_prop_count * 100, 2),
    te_prop_genome_pct = round(te_prop_genome * 100, 2),
    genome_size_mb = round(genome_size / 1e6, 2),
    total_te_mb = round(total_te_bp / 1e6, 2),
    htt_te_mb = round(htt_te_bp / 1e6, 2)
  ) %>%
  select(
    genome,
    genome_size_mb,
    total_te_count,
    total_te_mb,
    te_prop_genome_pct,
    htt_te_count,
    htt_te_mb,
    htt_prop_bp_pct,
    htt_prop_count_pct
  ) %>%
  arrange(desc(htt_prop_bp_pct))

write_tsv(results_out, file.path(out_dir, "htt_mobilome_contribution.tsv"))

# Save summary stats
sink(file.path(out_dir, "htt_summary_stats.txt"))
cat("HTT Contribution to Mobilome - Summary Statistics\n")
cat("=================================================\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("HTT input file:", basename(htt_file), "\n")
cat("  (Filtered: 0% flank synteny, >75% TE identity)\n\n")
cat("Outgroup excluded:", OUTGROUP, "\n")
cat("Collapsed-node genomes excluded:", length(collapsed_genomes), "\n")
if (length(collapsed_genomes) > 0) {
  cat("  ", paste(collapsed_genomes, collapse = ", "), "\n", sep = "")
}
cat("\nGenomes analyzed:", nrow(results), "\n\n")
cat("HTT-derived TE content (% of mobilome):\n")
cat("  Mean:   ", round(mean(results$htt_prop_bp) * 100, 2), "%\n")
cat("  Median: ", round(median(results$htt_prop_bp) * 100, 2), "%\n")
cat("  SD:     ", round(sd(results$htt_prop_bp) * 100, 2), "%\n")
cat("  Range:  ", round(min(results$htt_prop_bp) * 100, 2), " - ",
    round(max(results$htt_prop_bp) * 100, 2), "%\n\n")
cat("Correlations (Spearman):\n")
cat("  HTT bp vs genome size:    rho =", round(cor_genome$estimate, 3),
    ", p =", format.pval(cor_genome$p.value), "\n")
cat("  HTT bp vs mobilome:       rho =", round(cor_mob$estimate, 3),
    ", p =", format.pval(cor_mob$p.value), "\n")
cat("  HTT prop vs TE content:   rho =", round(cor_te$estimate, 3),
    ", p =", format.pval(cor_te$p.value), "\n")
sink()

cat("\n", strrep("=", 60), "\n", sep = "")
cat("OUTPUT FILES\n")
cat(strrep("=", 60), "\n\n", sep = "")
cat(file.path(out_dir, "htt_mobilome_contribution.tsv"), "\n")
cat(file.path(out_dir, "htt_summary_stats.txt"), "\n")
cat(file.path(out_dir, "htt_proportion_barplot.png/pdf"), "\n")
cat(file.path(out_dir, "htt_vs_genome_size.png/pdf"), "\n")
cat(file.path(out_dir, "htt_vs_mobilome.png/pdf"), "\n")
cat(file.path(out_dir, "htt_prop_vs_te_content.png/pdf"), "\n")

