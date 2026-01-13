# ==============================================================================
# 11_enrichment_analysis_v4.R
# HTT enrichment analysis - proportion plots + triangular heatmap
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggpubr)
library(forcats)

# --- Paths --------------------------------------------------------------------
data_dir    <- "data/X_enrichment_analysis"
results_dir <- "results/X_enrichment_analysis"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- Theme --------------------------------------------------------------------
theme_pub <- function() {
  theme_pubr() +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      legend.title     = element_text(size = 10),
      legend.text      = element_text(size = 9),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

# --- Load data ----------------------------------------------------------------
cat("Loading data...\n")

te_info        <- read.delim(file.path(data_dir, "all_tes_info.tsv"), stringsAsFactors = FALSE)
htt_candidates <- read.delim(file.path(data_dir, "htt_candidates.tsv"), stringsAsFactors = FALSE)
clade_assign   <- read.delim(file.path(data_dir, "clade_assignment.tsv"), stringsAsFactors = FALSE)

cat("Total TEs: ", nrow(te_info), "\n")
cat("HTT pairs: ", nrow(htt_candidates), "\n")

# --- Parse TE classification --------------------------------------------------
te_info <- te_info %>%
  mutate(
    te_annot = sub("::.*$", "", original_id),
    te_class = case_when(
      grepl("^LTR",      te_annot) ~ "LTR",
      grepl("^LINE",     te_annot) ~ "LINE",
      grepl("^SINE",     te_annot) ~ "SINE",
      grepl("^DNA",      te_annot) ~ "DNA",
      grepl("^RC|Helitron", te_annot) ~ "RC/Helitron",
      grepl("^PLE",      te_annot) ~ "PLE",
      TRUE ~ "Unclassified"
    ),
    te_superfamily_raw = sub("^[^/]*/", "", sub("::.*", "", te_annot)),
    te_superfamily_raw = ifelse(te_superfamily_raw == te_annot, NA_character_, te_superfamily_raw),
    superfamily = case_when(
      grepl("Copia",    te_annot, ignore.case = TRUE) ~ "LTR/Copia",
      grepl("Gypsy",    te_annot, ignore.case = TRUE) ~ "LTR/Gypsy",
      grepl("TcMar|Tc1|Mariner|Fot1", te_annot, ignore.case = TRUE) ~ "DNA/Tc1-Mariner",
      grepl("hAT",      te_annot, ignore.case = TRUE) ~ "DNA/hAT",
      grepl("Mutator|MULE", te_annot, ignore.case = TRUE) ~ "DNA/Mutator",
      grepl("PIF|Harbinger", te_annot, ignore.case = TRUE) ~ "DNA/PIF-Harbinger",
      grepl("CACTA|CMC",te_annot, ignore.case = TRUE) ~ "DNA/CACTA",
      grepl("Helitron", te_annot, ignore.case = TRUE) ~ "RC/Helitron",
      grepl("PiggyBac", te_annot, ignore.case = TRUE) ~ "DNA/PiggyBac",
      grepl("Kolobok",  te_annot, ignore.case = TRUE) ~ "DNA/Kolobok",
      grepl("Zisupton|Zator", te_annot, ignore.case = TRUE) ~ "DNA/Zisupton",
      grepl("L1",       te_annot, ignore.case = TRUE) ~ "LINE/L1",
      grepl("R2",       te_annot, ignore.case = TRUE) ~ "LINE/R2",
      grepl("RTE",      te_annot, ignore.case = TRUE) ~ "LINE/RTE",
      grepl("Tad1",     te_annot, ignore.case = TRUE) ~ "LINE/Tad1",
      grepl("I-Jockey|I/Jockey|Jockey", te_annot, ignore.case = TRUE) ~ "LINE/I-Jockey",
      grepl("LINE/I$|LINE/I::", te_annot) ~ "LINE/I",
      !is.na(te_superfamily_raw) ~ paste0(te_class, "/", te_superfamily_raw),
      TRUE ~ te_class
    ),
    te_type = case_when(
      te_class %in% c("LTR", "LINE", "SINE", "PLE") ~ "RNA",
      te_class %in% c("DNA", "RC/Helitron") ~ "DNA",
      TRUE ~ "Other"
    )
  ) %>%
  select(-te_annot, -te_superfamily_raw)

te_info <- te_info %>%
  rename(genome = strain) %>%
  left_join(clade_assign, by = c("genome" = "sample"))

htt_tes <- unique(c(htt_candidates$qseqid, htt_candidates$sseqid))
te_info <- te_info %>%
  mutate(in_htt = te_id %in% htt_tes)

cat("Unique TEs in HTT: ", sum(te_info$in_htt), "\n")
cat("TEs not in HTT:    ", sum(!te_info$in_htt), "\n\n")

# ==============================================================================
# ANALYSIS A: Superfamily enrichment - proportion plot
# ==============================================================================
cat("=== Analysis A: Superfamily enrichment ===\n")

sf_counts <- te_info %>%
  group_by(superfamily) %>%
  summarise(
    total   = n(),
    in_htt  = sum(in_htt),
    not_htt = sum(!in_htt),
    .groups = "drop"
  )

total_htt     <- sum(sf_counts$in_htt)
total_not_htt <- sum(sf_counts$not_htt)

sf_enrichment <- sf_counts %>%
  rowwise() %>%
  mutate(
    c_other_htt     = total_htt - in_htt,
    d_other_not_htt = total_not_htt - not_htt,
    fisher_p   = fisher.test(matrix(c(in_htt, not_htt, c_other_htt, d_other_not_htt), 
                                    nrow = 2), alternative = "greater")$p.value,
    odds_ratio = as.numeric(fisher.test(matrix(c(in_htt, not_htt, c_other_htt, d_other_not_htt), 
                                               nrow = 2))$estimate)
  ) %>%
  ungroup() %>%
  mutate(
    p_adjusted  = p.adjust(fisher_p, method = "BH"),
    significant = p_adjusted < 0.05,
    pct_htt     = 100 * in_htt / total,
    pct_of_all_htt = 100 * in_htt / total_htt
  ) %>%
  arrange(desc(in_htt))

write.table(sf_enrichment, file.path(results_dir, "table_A_superfamily_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nSuperfamily enrichment:\n")
print(sf_enrichment %>% select(superfamily, total, in_htt, pct_htt, odds_ratio, p_adjusted, significant))

# Prepare data for proportion plot
sf_for_plot <- sf_enrichment %>%
  filter(in_htt > 0) %>%
  mutate(superfamily = fct_reorder(superfamily, in_htt, .desc = TRUE))

# Color palette
te_colors <- c(
  "LINE/Tad1" = "#1f78b4",
  "LINE/I-Jockey" = "#a6cee3", 
  "LINE/I" = "#6baed6",
  "LINE" = "#c6dbef",
  "LTR/Gypsy" = "#33a02c",
  "LTR/Copia" = "#b2df8a",
  "LTR/ERV1" = "#74c476",
  "LTR/Ngaro" = "#c7e9c0",
  "LTR" = "#e5f5e0",
  "DNA/Tc1-Mariner" = "#e31a1c",
  "DNA/PiggyBac" = "#fb9a99",
  "DNA/PIF-Harbinger" = "#fc9272",
  "DNA/Kolobok" = "#fcbba1",
  "DNA" = "#fee0d2",
  "RC/Helitron" = "#ff7f00",
  "Unclassified" = "#999999"
)

sf_in_data <- unique(as.character(sf_for_plot$superfamily))
te_colors_use <- te_colors[names(te_colors) %in% sf_in_data]
missing_sf <- setdiff(sf_in_data, names(te_colors_use))
if (length(missing_sf) > 0) {
  extra_colors <- scales::hue_pal()(length(missing_sf))
  names(extra_colors) <- missing_sf
  te_colors_use <- c(te_colors_use, extra_colors)
}

# Proportion data
sf_props <- sf_for_plot %>%
  mutate(
    pct_of_total = 100 * total / sum(total),
    pct_of_htt   = 100 * in_htt / sum(in_htt)
  ) %>%
  select(superfamily, pct_of_total, pct_of_htt, significant) %>%
  pivot_longer(cols = c(pct_of_total, pct_of_htt), 
               names_to = "category", values_to = "pct") %>%
  mutate(category = factor(category, levels = c("pct_of_total", "pct_of_htt"),
                           labels = c("% of all TEs", "% of HTT TEs")))

# Calculate positions for asterisks on "% of HTT TEs" bar
sf_cumsum <- sf_for_plot %>%
  arrange(desc(as.numeric(superfamily))) %>%  # Match stacking order
  mutate(
    pct_of_htt = 100 * in_htt / sum(in_htt),
    cumsum_pct = cumsum(pct_of_htt),
    ypos = cumsum_pct - pct_of_htt / 2
  )

sig_positions <- sf_cumsum %>%
  filter(significant) %>%
  mutate(category = "% of HTT TEs") %>%
  select(superfamily, ypos, category)

p_a <- ggplot(sf_props, aes(x = category, y = pct, fill = superfamily)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = te_colors_use, name = "TE superfamily") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_pub() +
  theme(legend.position = "right")

# Add asterisks for significant enrichment
if (nrow(sig_positions) > 0) {
  p_a <- p_a +
    geom_text(data = sig_positions,
              aes(x = category, y = ypos, label = "*"),
              inherit.aes = FALSE, size = 8, color = "white", fontface = "bold")
}

ggsave(file.path(results_dir, "A_superfamily_enrichment.pdf"), p_a, width = 8, height = 6)
ggsave(file.path(results_dir, "A_superfamily_enrichment.png"), p_a, width = 8, height = 6, dpi = 300)

# ==============================================================================
# ANALYSIS B: Clade enrichment - proportion plot
# ==============================================================================
cat("\n=== Analysis B: Clade enrichment ===\n")

clade_counts <- te_info %>%
  filter(!is.na(clade)) %>%
  group_by(clade) %>%
  summarise(
    n_genomes = n_distinct(genome),
    total     = n(),
    in_htt    = sum(in_htt),
    not_htt   = sum(!in_htt),
    .groups   = "drop"
  )

total_htt_clade     <- sum(clade_counts$in_htt)
total_not_htt_clade <- sum(clade_counts$not_htt)

clade_enrichment <- clade_counts %>%
  rowwise() %>%
  mutate(
    c_other_htt     = total_htt_clade - in_htt,
    d_other_not_htt = total_not_htt_clade - not_htt,
    fisher_p   = fisher.test(matrix(c(in_htt, not_htt, c_other_htt, d_other_not_htt),
                                    nrow = 2), alternative = "greater")$p.value,
    odds_ratio = as.numeric(fisher.test(matrix(c(in_htt, not_htt, c_other_htt, d_other_not_htt),
                                               nrow = 2))$estimate)
  ) %>%
  ungroup() %>%
  mutate(
    p_adjusted  = p.adjust(fisher_p, method = "BH"),
    significant = p_adjusted < 0.05,
    pct_htt     = 100 * in_htt / total
  ) %>%
  arrange(desc(in_htt))

write.table(clade_enrichment, file.path(results_dir, "table_B_clade_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nClade enrichment:\n")
print(clade_enrichment %>% select(clade, n_genomes, total, in_htt, pct_htt, odds_ratio, p_adjusted, significant))

# Proportion plot for clades
clade_for_plot <- clade_enrichment %>%
  mutate(clade = fct_reorder(clade, in_htt, .desc = TRUE))

n_clades <- n_distinct(clade_for_plot$clade)
clade_colors <- scales::hue_pal()(n_clades)
names(clade_colors) <- levels(clade_for_plot$clade)

clade_props <- clade_for_plot %>%
  mutate(
    pct_of_total = 100 * total / sum(total),
    pct_of_htt   = 100 * in_htt / sum(in_htt)
  ) %>%
  select(clade, pct_of_total, pct_of_htt, significant) %>%
  pivot_longer(cols = c(pct_of_total, pct_of_htt),
               names_to = "category", values_to = "pct") %>%
  mutate(category = factor(category, levels = c("pct_of_total", "pct_of_htt"),
                           labels = c("% of all TEs", "% of HTT TEs")))

# Calculate positions for asterisks
clade_cumsum <- clade_for_plot %>%
  arrange(desc(as.numeric(clade))) %>%
  mutate(
    pct_of_htt = 100 * in_htt / sum(in_htt),
    cumsum_pct = cumsum(pct_of_htt),
    ypos = cumsum_pct - pct_of_htt / 2
  )

sig_clade_positions <- clade_cumsum %>%
  filter(significant) %>%
  mutate(category = "% of HTT TEs") %>%
  select(clade, ypos, category)

p_b <- ggplot(clade_props, aes(x = category, y = pct, fill = clade)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = clade_colors, name = "Clade") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_pub() +
  theme(legend.position = "right")

if (nrow(sig_clade_positions) > 0) {
  p_b <- p_b +
    geom_text(data = sig_clade_positions,
              aes(x = category, y = ypos, label = "*"),
              inherit.aes = FALSE, size = 8, color = "white", fontface = "bold")
}

ggsave(file.path(results_dir, "B_clade_enrichment.pdf"), p_b, width = 8, height = 6)
ggsave(file.path(results_dir, "B_clade_enrichment.png"), p_b, width = 8, height = 6, dpi = 300)

# ==============================================================================
# ANALYSIS C: Genome × Superfamily enrichment (bubble plot)
# ==============================================================================
cat("\n=== Analysis C: Genome × Superfamily enrichment ===\n")

genome_sf_counts <- te_info %>%
  filter(!is.na(clade)) %>%
  group_by(genome, clade, superfamily) %>%
  summarise(
    total   = n(),
    in_htt  = sum(in_htt),
    not_htt = sum(!in_htt),
    .groups = "drop"
  )

genome_sf_enrichment <- genome_sf_counts %>%
  group_by(genome) %>%
  mutate(
    genome_total_htt     = sum(in_htt),
    genome_total_not_htt = sum(not_htt)
  ) %>%
  ungroup() %>%
  filter(in_htt > 0) %>%
  rowwise() %>%
  mutate(
    c_other_htt     = genome_total_htt - in_htt,
    d_other_not_htt = genome_total_not_htt - not_htt,
    fisher_p = tryCatch(
      fisher.test(matrix(c(in_htt, not_htt, c_other_htt, d_other_not_htt), nrow = 2),
                  alternative = "greater")$p.value,
      error = function(e) NA_real_
    ),
    odds_ratio = tryCatch(
      as.numeric(fisher.test(matrix(c(in_htt, not_htt, c_other_htt, d_other_not_htt), 
                                    nrow = 2))$estimate),
      error = function(e) NA_real_
    )
  ) %>%
  ungroup() %>%
  mutate(
    p_adjusted = p.adjust(fisher_p, method = "BH"),
    enriched   = !is.na(p_adjusted) & p_adjusted < 0.05
  )

write.table(genome_sf_enrichment, file.path(results_dir, "table_C_genome_superfamily_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

genomes_per_clade <- te_info %>%
  filter(!is.na(clade)) %>%
  distinct(genome, clade) %>%
  count(clade, name = "n_genomes_total")

clade_sf_summary <- genome_sf_enrichment %>%
  group_by(clade, superfamily) %>%
  summarise(
    htt_involved_tes     = sum(in_htt),
    n_genomes_with_htt   = n_distinct(genome),
    n_genomes_enriched   = sum(enriched, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(genomes_per_clade, by = "clade") %>%
  mutate(
    frac_enriched = n_genomes_enriched / n_genomes_total,
    any_enriched  = n_genomes_enriched > 0
  )

top_sf <- clade_sf_summary %>%
  group_by(superfamily) %>%
  summarise(total = sum(htt_involved_tes)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 12) %>%
  pull(superfamily)

plot_c <- clade_sf_summary %>%
  filter(superfamily %in% top_sf) %>%
  mutate(
    superfamily = fct_reorder(superfamily, htt_involved_tes, .fun = sum, .desc = TRUE),
    clade       = fct_reorder(clade, htt_involved_tes, .fun = sum, .desc = TRUE),
    enrich_label = ifelse(n_genomes_enriched > 0, 
                          paste0(n_genomes_enriched, "/", n_genomes_total), "")
  )

p_c <- ggplot(plot_c, aes(x = superfamily, y = clade)) +
  geom_point(aes(size = htt_involved_tes, color = any_enriched), alpha = 0.85) +
  geom_text(aes(label = enrich_label), size = 2.5, vjust = 0.4) +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#2A9D8F"),
                     labels = c("No genome enriched", "≥1 genome enriched")) +
  scale_size_continuous(range = c(2, 16), breaks = pretty_breaks(4), labels = comma) +
  labs(x = "TE superfamily", y = "Clade", size = "TEs in HTT", color = NULL) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "C_clade_superfamily_bubble.pdf"), p_c, width = 10, height = 6)
ggsave(file.path(results_dir, "C_clade_superfamily_bubble.png"), p_c, width = 10, height = 6, dpi = 300)

# ==============================================================================
# ANALYSIS D: Clade-pair heatmap - triangular with DNA/RNA split
# ==============================================================================
cat("\n=== Analysis D: Clade-pair heatmap ===\n")

htt_with_type <- htt_candidates %>%
  mutate(
    te_type = case_when(
      grepl("\\|LTR|\\|LINE|\\|SINE|\\|PLE", qseqid) ~ "RNA",
      grepl("\\|DNA|\\|RC|\\|Helitron", qseqid) ~ "DNA",
      TRUE ~ "Other"
    )
  )

# Aggregate by clade pair and TE type (use ordered pairs for symmetry)
clade_pair_counts <- htt_with_type %>%
  filter(te_type != "Other") %>%
  mutate(
    clade_a = pmin(qclade, sclade),
    clade_b = pmax(qclade, sclade)
  ) %>%
  group_by(clade_a, clade_b, te_type) %>%
  summarise(n_htt = n(), .groups = "drop")

write.table(clade_pair_counts, file.path(results_dir, "table_D_clade_pair_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Order clades by total HTT involvement
clade_totals <- bind_rows(
  htt_with_type %>% count(clade = qclade, name = "n"),
  htt_with_type %>% count(clade = sclade, name = "n")
) %>%
  group_by(clade) %>%
  summarise(total = sum(n)) %>%
  arrange(desc(total))

clade_order <- clade_totals$clade
n_clades <- length(clade_order)

# Build plot data directly
# Upper triangle (idx1 < idx2) = DNA
# Lower triangle (idx1 > idx2) = RNA
# But y-axis is reversed, so we need to swap the assignments

# RNA goes to upper-left visually (which is lower triangle in reversed y coords)
# DNA goes to lower-right visually (which is upper triangle in reversed y coords)

rna_data <- clade_pair_counts %>%
  filter(te_type == "RNA") %>%
  mutate(
    # Upper-left visually: smaller x index, larger y index (in reversed coords = smaller original)
    clade1 = factor(clade_a, levels = clade_order),
    clade2 = factor(clade_b, levels = rev(clade_order))
  ) %>%
  select(clade1, clade2, n_htt)

dna_data <- clade_pair_counts %>%
  filter(te_type == "DNA") %>%
  mutate(
    # Lower-right visually: swap coordinates
    clade1 = factor(clade_b, levels = clade_order),
    clade2 = factor(clade_a, levels = rev(clade_order))
  ) %>%
  select(clade1, clade2, n_htt)

plot_data <- bind_rows(
  dna_data %>% mutate(triangle = "DNA"),
  rna_data %>% mutate(triangle = "RNA")
)

# Create the heatmap
p_d <- ggplot(plot_data, aes(x = clade1, y = clade2, fill = n_htt)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n_htt > 0, n_htt, "")), size = 2.5) +
  # Add diagonal line from top-left to bottom-right
  geom_segment(aes(x = 0.5, y = n_clades + 0.5, xend = n_clades + 0.5, yend = 0.5),
               linewidth = 1, color = "black", inherit.aes = FALSE) +
  scale_fill_gradient(low = "white", high = "#E76F51", trans = "sqrt",
                      breaks = c(10, 100, 500, 1000), limits = c(0, NA),
                      na.value = "white") +
  # Add labels - RNA upper-left, DNA lower-right
  annotate("text", x = 3, y = n_clades - 1.5, label = "RNA", 
           fontface = "bold", size = 5, color = "grey30") +
  annotate("text", x = n_clades - 2, y = 2.5, label = "DNA", 
           fontface = "bold", size = 5, color = "grey30") +
  labs(x = NULL, y = NULL, fill = "HTT events") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  ) +
  coord_fixed()

ggsave(file.path(results_dir, "D_clade_pair_heatmap.pdf"), p_d, width = 10, height = 9)
ggsave(file.path(results_dir, "D_clade_pair_heatmap.png"), p_d, width = 10, height = 9, dpi = 300)

# ==============================================================================
# Summary
# ==============================================================================
cat("\n=== Summary ===\n")
cat("Total TEs:", nrow(te_info), "\n")
cat("TEs in HTT:", sum(te_info$in_htt), sprintf("(%.1f%%)\n", 100*sum(te_info$in_htt)/nrow(te_info)))

cat("\nHTT by TE type:\n")
print(htt_with_type %>% count(te_type) %>% mutate(pct = round(100 * n / sum(n), 1)))

cat("\nSignificantly enriched superfamilies:\n")
print(sf_enrichment %>% filter(significant) %>% select(superfamily, in_htt, pct_htt, odds_ratio, p_adjusted))

cat("\nTop clade pairs (RNA):\n")
print(clade_pair_counts %>% filter(te_type == "RNA") %>% arrange(desc(n_htt)) %>% head(5))

cat("\nTop clade pairs (DNA):\n")
print(clade_pair_counts %>% filter(te_type == "DNA") %>% arrange(desc(n_htt)) %>% head(5))

cat("\nOutput files in:", results_dir, "\n")
sessionInfo()