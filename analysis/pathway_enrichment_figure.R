#!/usr/bin/env Rscript
# ============================================================================
# Extended Data Figure: Pathway enrichment analysis of MyoScore dimensions
# Multi-panel: A) Dot plot  B) Heatmap  C) Direction bars  D) Gene convergence
# ============================================================================

library(grDevices)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(cowplot)
library(patchwork)
library(ggrepel)

# -- Project colors (from CLAUDE.md) --
# Strength, Mass, LeanMuscle, Youth, Resilience
dim_colors <- c(
  "Strength"   = "#50327b",
  "Mass"       = "#46508b",
  "LeanMuscle" = "#f4e030",
  "Youth"      = "#72c95e",
  "Resilience" = "#31848f"
)

# Direction colors
col_healthy  <- "#f4e030"
col_disease  <- "#50327b"

# Publication theme base
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      text = element_text(family = "sans", color = "grey10"),
      plot.title = element_text(face = "bold", size = base_size + 1,
                                margin = margin(b = 8)),
      axis.text = element_text(color = "grey20"),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# -- Load data --
df <- read.csv("MyoScore/data/pathway_enrichment_by_dimension.csv",
               stringsAsFactors = FALSE)

dim_levels <- c("Strength", "Mass", "LeanMuscle", "Youth", "Resilience")

sig <- df %>%
  filter(padj < 0.05) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    term_short = str_remove(term, " \\(GO:\\d+\\)"),
    term_short = str_wrap(term_short, width = 32),
    dimension = factor(dimension, levels = dim_levels),
    bio_theme = case_when(
      str_detect(term, "thiol|CoA") ~ "Acetyl-CoA metabolism",
      str_detect(term, "thanol") ~ "Ethanol metabolism",
      str_detect(term, "icrotubule|ubulin") ~ "Microtubule/cytoskeleton",
      str_detect(term, "FGF|ibroblast") ~ "FGF signaling",
      str_detect(term, "ropanoate") ~ "Propanoate metabolism",
      str_detect(term, "erroptosis") ~ "Ferroptosis",
      str_detect(term, "MHC|HLA") ~ "MHC/immune surveillance",
      str_detect(term, "ethyl") ~ "mRNA methylation (m6A)",
      str_detect(term, "entrosome") ~ "Centrosome/cell division",
      str_detect(term, "asopressin") ~ "Vesicular transport",
      str_detect(term, "almonella") ~ "Infection response",
      str_detect(term, "olgi") ~ "Golgi-PM transport",
      str_detect(term, "ron.*ulfur|Fe-S") ~ "Iron-sulfur cluster",
      str_detect(term, "ulfat") ~ "Sulfation",
      str_detect(term, "mine") ~ "Amine metabolism",
      str_detect(term, "lcohol") ~ "Alcohol catabolism",
      str_detect(term, "itochondrion.*ransport") ~ "Mitochondrial transport",
      str_detect(term, "ubiquitin") ~ "Ubiquitin binding",
      str_detect(term, "Nucleoside") ~ "Kinase activity",
      str_detect(term, "ransmembrane.*yrosine") ~ "RTK activity",
      TRUE ~ term_short
    )
  )

# ============================================================================
# Panel A: Dot plot - top enrichments per dimension
# ============================================================================

top_terms <- sig %>%
  group_by(dimension, bio_theme) %>%
  slice_min(padj, n = 1) %>%
  ungroup() %>%
  group_by(dimension, bio_theme) %>%
  arrange(match(direction_group, c("all", "healthy_up", "unhealthy_up"))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(dimension) %>%
  slice_min(padj, n = 4) %>%
  ungroup() %>%
  mutate(bio_theme_wrap = str_wrap(bio_theme, width = 28))

p_a <- ggplot(top_terms,
              aes(x = dimension,
                  y = reorder(bio_theme_wrap, neg_log10_padj))) +
  geom_point(aes(size = overlap_count, fill = dimension),
             shape = 21, color = "white", stroke = 0.8, alpha = 0.9) +
  scale_fill_manual(values = dim_colors, guide = "none") +
  scale_size_continuous(range = c(4, 11), name = "Gene count",
                        breaks = c(2, 3, 5)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL,
       title = "Pathway enrichments by MyoScore dimension") +
  theme_pub(10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 10,
                               color = dim_colors[dim_levels]),
    axis.text.y = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = alpha("white", 0.8),
                                     color = NA)
  )

# ============================================================================
# Panel B: Heatmap - biological themes across dimensions
# ============================================================================

summary_data <- sig %>%
  filter(!is.na(bio_theme)) %>%
  group_by(dimension, bio_theme) %>%
  summarise(
    best_padj = min(padj),
    neg_log10 = -log10(min(padj)),
    gene_count = max(overlap_count),
    .groups = "drop"
  )

# Order bio_themes by total significance
theme_order <- summary_data %>%
  group_by(bio_theme) %>%
  summarise(total_sig = sum(neg_log10), .groups = "drop") %>%
  arrange(total_sig) %>%
  pull(bio_theme)

summary_data$bio_theme <- factor(summary_data$bio_theme, levels = theme_order)

p_b <- ggplot(summary_data, aes(x = dimension, y = bio_theme)) +
  geom_tile(aes(fill = neg_log10), color = "white", linewidth = 1.2) +
  geom_text(aes(label = gene_count), size = 3.2, fontface = "bold",
            color = "white") +
  scale_fill_gradient(low = alpha("#50327b", 0.15), high = "#50327b",
                      name = expression(-log[10]~italic(P)[adj]),
                      limits = c(1, NA),
                      breaks = c(1, 1.5, 2, 2.5)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Biological themes across dimensions") +
  theme_pub(10) +
  theme(
    axis.text.x.top = element_text(face = "bold", size = 10,
                                   color = dim_colors[dim_levels]),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm")
  )

# ============================================================================
# Panel C: Direction-stratified bar plot
# ============================================================================

dir_terms <- sig %>%
  filter(direction_group != "all") %>%
  group_by(dimension, bio_theme) %>%
  slice_min(padj, n = 1) %>%
  ungroup() %>%
  mutate(
    neg_log10_padj_signed = ifelse(direction_group == "unhealthy_up",
                                   -neg_log10_padj, neg_log10_padj),
    bio_theme_wrap = str_wrap(bio_theme, width = 22),
    dimension = factor(dimension, levels = dim_levels)
  )

p_c <- ggplot(dir_terms,
              aes(x = neg_log10_padj_signed,
                  y = reorder(bio_theme_wrap, neg_log10_padj_signed))) +
  geom_col(aes(fill = direction_group), width = 0.65) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  facet_grid(dimension ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_fill_manual(
    values = c("healthy_up" = col_healthy, "unhealthy_up" = col_disease),
    labels = c("Healthy direction", "Disease direction"),
    name = NULL
  ) +
  labs(x = expression(-log[10]~italic(P)[adj]),
       y = NULL,
       title = "Direction-stratified pathway enrichment") +
  theme_pub(10) +
  theme(
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 9,
                                     color = "grey20"),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey96", color = NA),
    axis.text.y = element_text(size = 7.5),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "bottom",
    legend.margin = margin(t = -2)
  )

# Add direction labels
p_c <- p_c +
  annotate("text", x = Inf, y = -Inf, label = "Healthy >>",
           hjust = 1.05, vjust = -0.3,
           size = 2.8, color = "#9E8A00", fontface = "bold") +
  annotate("text", x = -Inf, y = -Inf, label = "<< Disease",
           hjust = -0.05, vjust = -0.3,
           size = 2.8, color = col_disease, fontface = "bold")

# ============================================================================
# Panel D: Gene-pathway convergence tile (Strength dimension)
# ============================================================================

strength_sig <- sig %>%
  filter(dimension == "Strength",
         direction_group %in% c("all", "healthy_up")) %>%
  group_by(bio_theme) %>%
  slice_min(padj, n = 1) %>%
  ungroup()

gene_pathway <- strength_sig %>%
  select(bio_theme, genes, neg_log10_padj) %>%
  separate_rows(genes, sep = ";") %>%
  mutate(genes = trimws(genes)) %>%
  distinct()

pathway_order <- strength_sig %>%
  arrange(padj) %>%
  pull(bio_theme) %>%
  unique()
gene_pathway$bio_theme <- factor(gene_pathway$bio_theme,
                                  levels = rev(pathway_order))

gene_counts <- gene_pathway %>%
  group_by(genes) %>%
  summarise(n_pathways = n(), .groups = "drop") %>%
  arrange(desc(n_pathways))

gene_pathway$genes <- factor(gene_pathway$genes, levels = gene_counts$genes)
gene_pathway <- gene_pathway %>%
  left_join(gene_counts, by = "genes")

p_d <- ggplot(gene_pathway, aes(x = genes, y = bio_theme)) +
  geom_tile(aes(fill = n_pathways), color = "white", linewidth = 1) +
  geom_point(color = "white", size = 1.2, alpha = 0.6) +
  scale_fill_gradient(low = alpha("#50327b", 0.2), high = "#50327b",
                      name = "Shared\npathways",
                      breaks = c(1, 2, 3, 4)) +
  labs(x = NULL, y = NULL,
       title = "Strength dimension: gene-pathway convergence") +
  theme_pub(10) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, face = "italic",
                               size = 8.5),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.7, "cm")
  )

# ============================================================================
# Assemble: 2x2 layout
# ============================================================================

top_row    <- p_a + p_b + plot_layout(widths = c(1.1, 1))
bottom_row <- p_c + p_d + plot_layout(widths = c(1.1, 1))

final <- (top_row / bottom_row) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.margin = margin(5, 5, 5, 5)
    )
  ) &
  theme(
    plot.tag = element_text(size = 16, face = "bold", family = "sans")
  )

# Save
ggsave("MyoScore/figures/ExtData_pathway_enrichment_4panel.pdf",
       final, width = 15, height = 13.5, device = cairo_pdf)
ggsave("MyoScore/figures/ExtData_pathway_enrichment_4panel.png",
       final, width = 15, height = 13.5, dpi = 300)

cat("Saved: ExtData_pathway_enrichment_4panel.pdf/png\n")
