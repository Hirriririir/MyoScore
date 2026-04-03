#!/usr/bin/env Rscript
# ── QC Paired Plots: 4-panel UMAP + stats (A–D) ──
# A: Sequencing platform (DNBSEQ-T7 vs NovaSeq 6000, n=3)
# B: Library preparation (polyA vs riboD, n=12)
# C: Within-individual muscle (RF vs VL, n=14)
# D: GTEx ischemia time (n=803)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(cowplot)

# ── Paths ──
# setwd("path/to/Myopathy_Spectrum_Multiomics")  # adjust to your project root
output_dir <- "MyoScore/figures"

# ── Colors ──
col_purple <- "#50327b"
col_yellow <- "#f4e030"
col_red    <- "#e63946"
col_gray   <- "#d3d3d3"

# ── Load data ──
umap_df <- read.csv("MyoScore/data/Meta_umap_1722_tuned.csv", row.names = 1)
umap_df$Sample <- rownames(umap_df)

full_meta <- read.csv("Meta/RNA-seq integration meta (MultiOmics).csv", row.names = 1)
scores <- read.csv("GWAS_TWAS/gmhs_dimensions/dimension_scores_v33.csv", row.names = 1)
gtex_clin <- read.csv(
    "GWAS_TWAS/gmhs_dimensions/gtex_clinical_validation/gtex_clinical_merged.csv",
    row.names = 1)

# Merge scores + meta
umap_df$GMHS_v33 <- scores[umap_df$Sample, "GMHS_v33"]
umap_df$Geo_accession <- full_meta[umap_df$Sample, "Geo_accession"]
umap_df$BiopsySite <- full_meta[umap_df$Sample, "BiopsySite"]

# ── UMAP base function ──
make_umap_base <- function() {
    list(
        geom_point(data = umap_df, aes(UMAP_X, UMAP_Y),
                   color = col_gray, size = 0.8, alpha = 0.5),
        theme_minimal(base_size = 10),
        theme(panel.grid = element_blank(),
              axis.title = element_text(size = 9),
              plot.title = element_text(face = "bold", size = 9, hjust = 0),
              legend.position = "top", legend.title = element_blank(),
              aspect.ratio = 1),
        labs(x = "UMAP 1", y = "UMAP 2")
    )
}

# ══════════════════════════════════════════════════
# A: Sequencing Platform (HuashanMuscle, n=3)
# ══════════════════════════════════════════════════
hm_pairs <- data.frame(
    dnb = c("M1600_lb", "M1942_lb", "M2968_lb"),
    nova = c("M1600_nh", "M1942_nh", "M2968_nh"),
    stringsAsFactors = FALSE
)

hm_pts <- bind_rows(
    umap_df %>% filter(Sample %in% hm_pairs$dnb) %>% mutate(Group = "DNBSEQ-T7"),
    umap_df %>% filter(Sample %in% hm_pairs$nova) %>% mutate(Group = "NovaSeq 6000")
)
hm_lines <- hm_pairs %>%
    inner_join(umap_df[, c("Sample", "UMAP_X", "UMAP_Y")], by = c("dnb" = "Sample")) %>%
    rename(x1 = UMAP_X, y1 = UMAP_Y) %>%
    inner_join(umap_df[, c("Sample", "UMAP_X", "UMAP_Y")], by = c("nova" = "Sample")) %>%
    rename(x2 = UMAP_X, y2 = UMAP_Y)

p_a1 <- ggplot() + make_umap_base() +
    geom_segment(data = hm_lines, aes(x = x1, y = y1, xend = x2, yend = y2),
                 color = "#333", linewidth = 0.5, alpha = 0.6) +
    geom_point(data = hm_pts, aes(UMAP_X, UMAP_Y, color = Group, shape = Group),
               size = 4, stroke = 0.6) +
    scale_color_manual(values = c("DNBSEQ-T7" = col_purple, "NovaSeq 6000" = col_yellow)) +
    scale_shape_manual(values = c("DNBSEQ-T7" = 16, "NovaSeq 6000" = 17)) +
    ggtitle("Sequencing Platform Comparison\nHuashanMuscle: paired samples (n = 3)")

# Stats
ms_dnb <- umap_df[hm_pairs$dnb, "GMHS_v33"]
ms_nova <- umap_df[hm_pairs$nova, "GMHS_v33"]
pval_a <- tryCatch(wilcox.test(ms_dnb, ms_nova, paired = TRUE)$p.value,
                   error = function(e) NA)

hm_long <- data.frame(
    pair = rep(1:3, 2),
    Group = rep(c("DNBSEQ-T7", "NovaSeq 6000"), each = 3),
    MyoScore = c(ms_dnb, ms_nova)
)
p_a2 <- ggplot(hm_long, aes(Group, MyoScore)) +
    geom_boxplot(aes(fill = Group), alpha = 0.4, width = 0.5, outlier.shape = NA) +
    geom_line(aes(group = pair), color = "#666", alpha = 0.6) +
    geom_point(aes(color = Group), size = 3) +
    scale_fill_manual(values = c("DNBSEQ-T7" = col_purple, "NovaSeq 6000" = col_yellow)) +
    scale_color_manual(values = c("DNBSEQ-T7" = col_purple, "NovaSeq 6000" = col_yellow)) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(),
          aspect.ratio = 1) +
    labs(x = NULL, y = "MyoScore",
         title = sprintf("p = %.2f", ifelse(is.na(pval_a), 1, pval_a)))

# ══════════════════════════════════════════════════
# B: Library Prep (Myofin, polyA vs riboD, n=12)
# ══════════════════════════════════════════════════
ribo_ids <- grep("_riboD$", umap_df$Sample, value = TRUE)
mf_pairs <- data.frame(
    polya = sub("_riboD$", "_polyA", ribo_ids),
    ribo = ribo_ids, stringsAsFactors = FALSE
) %>% filter(polya %in% umap_df$Sample, ribo %in% umap_df$Sample)

mf_pts <- bind_rows(
    umap_df %>% filter(Sample %in% mf_pairs$polya) %>% mutate(Group = "polyA"),
    umap_df %>% filter(Sample %in% mf_pairs$ribo) %>% mutate(Group = "riboD")
)
mf_lines <- mf_pairs %>%
    inner_join(umap_df[, c("Sample", "UMAP_X", "UMAP_Y")], by = c("polya" = "Sample")) %>%
    rename(x1 = UMAP_X, y1 = UMAP_Y) %>%
    inner_join(umap_df[, c("Sample", "UMAP_X", "UMAP_Y")], by = c("ribo" = "Sample")) %>%
    rename(x2 = UMAP_X, y2 = UMAP_Y)

p_b1 <- ggplot() + make_umap_base() +
    geom_segment(data = mf_lines, aes(x = x1, y = y1, xend = x2, yend = y2),
                 color = "#333", linewidth = 0.3, alpha = 0.4) +
    geom_point(data = mf_pts, aes(UMAP_X, UMAP_Y, color = Group, shape = Group),
               size = 4, stroke = 0.6) +
    scale_color_manual(values = c("polyA" = col_purple, "riboD" = col_yellow)) +
    scale_shape_manual(values = c("polyA" = 16, "riboD" = 17)) +
    ggtitle(sprintf("Library Preparation Comparison\nMyofin: paired samples (n = %d)",
                    nrow(mf_pairs)))

ms_polya <- umap_df[mf_pairs$polya, "GMHS_v33"]
ms_ribo <- umap_df[mf_pairs$ribo, "GMHS_v33"]
pval_b <- tryCatch(wilcox.test(ms_polya, ms_ribo, paired = TRUE)$p.value,
                   error = function(e) NA)

mf_long <- data.frame(
    pair = rep(seq_len(nrow(mf_pairs)), 2),
    Group = rep(c("polyA", "riboD"), each = nrow(mf_pairs)),
    MyoScore = c(ms_polya, ms_ribo)
)
p_b2 <- ggplot(mf_long, aes(Group, MyoScore)) +
    geom_boxplot(aes(fill = Group), alpha = 0.4, width = 0.5, outlier.shape = NA) +
    geom_line(aes(group = pair), color = "#666", alpha = 0.4) +
    geom_point(aes(color = Group), size = 2) +
    scale_fill_manual(values = c("polyA" = col_purple, "riboD" = col_yellow)) +
    scale_color_manual(values = c("polyA" = col_purple, "riboD" = col_yellow)) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(),
          aspect.ratio = 1) +
    labs(x = NULL, y = "MyoScore",
         title = sprintf("p = %.2f", ifelse(is.na(pval_b), 1, pval_b)))

# ══════════════════════════════════════════════════
# C: Within-Individual (GSE202745: RF vs VL, n=14)
# ══════════════════════════════════════════════════
g202 <- umap_df %>% filter(Geo_accession == "GSE202745")
rf <- g202 %>% filter(grepl("Rectus Femoris|rectus femoris", BiopsySite, ignore.case = TRUE))
vl <- g202 %>% filter(grepl("Vastus Lateralis|vastus lateralis", BiopsySite, ignore.case = TRUE))

rf_nums <- setNames(rf$Sample, as.integer(gsub("\\D", "", rf$Sample)))
vl_nums <- setNames(vl$Sample, as.integer(gsub("\\D", "", vl$Sample)))

rv_pairs <- data.frame(vl = character(), rf = character(), stringsAsFactors = FALSE)
for (vn in sort(as.integer(names(vl_nums)))) {
    vs <- vl_nums[as.character(vn)]
    prefix <- substr(vs, 1, 1)
    rn <- as.character(vn + 1)
    if (rn %in% names(rf_nums) && substr(rf_nums[rn], 1, 1) == prefix) {
        rv_pairs <- rbind(rv_pairs,
                          data.frame(vl = vs, rf = rf_nums[rn], stringsAsFactors = FALSE))
    }
}

rv_pts <- bind_rows(
    umap_df %>% filter(Sample %in% rv_pairs$rf) %>% mutate(Group = "RF"),
    umap_df %>% filter(Sample %in% rv_pairs$vl) %>% mutate(Group = "VL")
)
rv_lines <- rv_pairs %>%
    inner_join(umap_df[, c("Sample", "UMAP_X", "UMAP_Y")], by = c("vl" = "Sample")) %>%
    rename(x1 = UMAP_X, y1 = UMAP_Y) %>%
    inner_join(umap_df[, c("Sample", "UMAP_X", "UMAP_Y")], by = c("rf" = "Sample")) %>%
    rename(x2 = UMAP_X, y2 = UMAP_Y)

p_c1 <- ggplot() + make_umap_base() +
    geom_segment(data = rv_lines, aes(x = x1, y = y1, xend = x2, yend = y2),
                 color = "#333", linewidth = 0.4, alpha = 0.5) +
    geom_point(data = rv_pts, aes(UMAP_X, UMAP_Y, color = Group, shape = Group),
               size = 4, stroke = 0.6) +
    scale_color_manual(values = c("RF" = col_purple, "VL" = col_yellow)) +
    scale_shape_manual(values = c("RF" = 16, "VL" = 17)) +
    ggtitle(sprintf("Within-Individual Muscle Correlation\nGSE202745: RF vs VL (n = %d)",
                    nrow(rv_pairs)))

ms_vl <- umap_df[rv_pairs$vl, "GMHS_v33"]
ms_rf <- umap_df[rv_pairs$rf, "GMHS_v33"]
cor_rv <- cor.test(ms_vl, ms_rf)

p_c2 <- ggplot(data.frame(VL = ms_vl, RF = ms_rf), aes(VL, RF)) +
    geom_smooth(method = "lm", se = TRUE, color = col_red, linewidth = 0.8, alpha = 0.15) +
    geom_point(color = col_purple, size = 2.5, stroke = 0.3) +
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1) +
    labs(x = "MyoScore (VL)", y = "MyoScore (RF)",
         title = sprintf("r = %.2f, p = %.3f", cor_rv$estimate, cor_rv$p.value))

# ══════════════════════════════════════════════════
# D: GTEx Ischemia Time
# ══════════════════════════════════════════════════
gtex_sub <- umap_df %>% filter(Data_source == "GTEx")
gtex_sub$SMTSISCH <- gtex_clin[gtex_sub$Sample, "SMTSISCH"]
gtex_sub$DTHHRDY <- gtex_clin[gtex_sub$Sample, "DTHHRDY"]
gtex_sub <- gtex_sub %>% filter(!is.na(SMTSISCH))
gtex_sub$isch_hr <- gtex_sub$SMTSISCH / 60

death_labels <- c("0" = "Ventilator", "1" = "Fast-Natural", "2" = "Fast-Violent",
                  "3" = "Intermediate", "4" = "Slow")
death_colors <- c("Ventilator" = "#3b3f8a", "Fast-Natural" = "#f4e030",
                  "Fast-Violent" = "#c8d62b", "Intermediate" = "#72c95e",
                  "Slow" = "#50327b")
gtex_sub$DeathType <- death_labels[as.character(gtex_sub$DTHHRDY)]

p_d1 <- ggplot() +
    geom_point(data = umap_df, aes(UMAP_X, UMAP_Y),
               color = col_gray, size = 0.8, alpha = 0.5) +
    geom_point(data = gtex_sub, aes(UMAP_X, UMAP_Y, color = SMTSISCH),
               size = 1, alpha = 0.8) +
    scale_color_viridis_c(name = "Ischemia\nTime (min)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), aspect.ratio = 1,
          plot.title = element_text(face = "bold", size = 9)) +
    labs(x = "UMAP 1", y = "UMAP 2",
         title = sprintf("GTEx Ischemia Time\nGTEx samples (n = %d)", nrow(gtex_sub)))

cor_isch <- cor.test(gtex_sub$isch_hr, gtex_sub$GMHS_v33)

# Partial correlation controlling death type
isch_resid <- residuals(lm(rank(SMTSISCH) ~ factor(DTHHRDY), data = gtex_sub))
score_resid <- residuals(lm(rank(GMHS_v33) ~ factor(DTHHRDY), data = gtex_sub))
partial_r <- cor(isch_resid, score_resid)

p_d2 <- ggplot(gtex_sub %>% filter(!is.na(DeathType)),
               aes(isch_hr, GMHS_v33, color = DeathType)) +
    geom_point(size = 1, alpha = 0.5) +
    scale_color_manual(values = death_colors, name = "Death Type") +
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, legend.position = "right",
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6)) +
    labs(x = "Ischemia Time (hours)", y = "MyoScore",
         title = sprintf("r = %.2f -> %.2f\npost-mortem interval",
                         cor_isch$estimate, partial_r)) +
    annotate("text", x = max(gtex_sub$isch_hr) * 0.7,
             y = max(gtex_sub$GMHS_v33) - 1,
             label = sprintf("Partial r = %.2f\n(after death type adj.)", partial_r),
             size = 3, hjust = 0)

# ══════════════════════════════════════════════════
# Combine: 2 rows × 4 cols
# ══════════════════════════════════════════════════
top_row <- plot_grid(p_a1, p_b1, p_c1, p_d1,
                     nrow = 1, labels = c("A", "B", "C", "D"),
                     label_size = 14)
bot_row <- plot_grid(p_a2, p_b2, p_c2, p_d2,
                     nrow = 1)
final <- plot_grid(top_row, bot_row, nrow = 2, rel_heights = c(1, 0.85))

ggsave(file.path(output_dir, "QC_UMAP_4panel.pdf"), final,
       width = 22, height = 10, dpi = 300)
ggsave(file.path(output_dir, "QC_UMAP_4panel.png"), final,
       width = 22, height = 10, dpi = 300)

# ── Individual panels (UMAP + stats combined per panel) ──
p_A <- plot_grid(p_a1, p_a2, nrow = 2, rel_heights = c(1, 0.85))
p_B <- plot_grid(p_b1, p_b2, nrow = 2, rel_heights = c(1, 0.85))
p_C <- plot_grid(p_c1, p_c2, nrow = 2, rel_heights = c(1, 0.85))
p_D <- plot_grid(p_d1, p_d2, nrow = 2, rel_heights = c(1, 0.85))

ggsave(file.path(output_dir, "QC_UMAP_A_platform.pdf"), p_A, width = 6, height = 10)
ggsave(file.path(output_dir, "QC_UMAP_B_library.pdf"), p_B, width = 6, height = 10)
ggsave(file.path(output_dir, "QC_UMAP_C_muscle.pdf"), p_C, width = 6, height = 10)
ggsave(file.path(output_dir, "QC_UMAP_D_ischemia.pdf"), p_D, width = 7, height = 10)

# ── Print summary ──
cat("=== QC Summary ===\n")
cat(sprintf("  A: Platform (n=3), p = %.2f\n", ifelse(is.na(pval_a), 1, pval_a)))
cat(sprintf("  B: Library (n=%d), p = %.2f\n", nrow(mf_pairs), ifelse(is.na(pval_b), 1, pval_b)))
cat(sprintf("  C: RF vs VL (n=%d), r = %.2f, p = %.3f\n",
            nrow(rv_pairs), cor_rv$estimate, cor_rv$p.value))
cat(sprintf("  D: Ischemia (n=%d), r = %.2f, partial r = %.2f\n",
            nrow(gtex_sub), cor_isch$estimate, partial_r))
cat("Saved: QC_UMAP_4panel.pdf/png\n")
