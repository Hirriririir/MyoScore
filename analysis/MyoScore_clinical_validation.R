# MyoScore_clinical_validation.R
# MyoScore clinical validation - FSHD, CDM, LGMDR12, Overweight Pre/Post
# Color scheme: green (#72c95e)

library(tidyverse)
library(readxl)
library(patchwork)

# setwd("path/to/Myopathy_Spectrum_Multiomics")  # adjust to your project root

# Green color
green_color <- '#72c95e'

# Read MyoScore data
scores <- read.csv('GWAS_TWAS/gmhs_dimensions/dimension_scores_v33.csv', row.names = 1)
scores$Sample_id <- rownames(scores)

# ========== Read validation data (skipped if file not present) ==========
plots <- list()

# Plot A: FSHD Inflammation Score
if (file.exists('Validation/FSHD_meta.xlsx')) {
  fshd_meta <- read_excel('Validation/FSHD_meta.xlsx')
  fshd_scores <- merge(fshd_meta, scores, by.x = 'Sample_ID', by.y = 'Sample_id')
  fshd_scores <- fshd_scores %>% filter(!is.na(Inflammation_Score))
  cor_fshd_inflam <- cor.test(fshd_scores$GMHS_v33, fshd_scores$Inflammation_Score, method = 'spearman')
  cat('FSHD Inflammation Score: rho =', round(cor_fshd_inflam$estimate, 2),
      ', p =', round(cor_fshd_inflam$p.value, 3), ', n =', nrow(fshd_scores), '\n')
  plots[['A']] <- ggplot(fshd_scores, aes(x = GMHS_v33, y = Inflammation_Score)) +
    geom_point(color = green_color, size = 2.5, alpha = 0.7) +
    geom_smooth(method = 'lm', color = green_color, fill = green_color, alpha = 0.2) +
    labs(x = 'MyoScore (GMHS v3.3)', y = 'FSHD Inflammation Score', title = 'A') +
    annotate('text', x = mean(fshd_scores$GMHS_v33),
             y = max(fshd_scores$Inflammation_Score) * 0.95,
             label = sprintf('ρ = %.2f, p = %.3f (n=%d)',
                            cor_fshd_inflam$estimate, cor_fshd_inflam$p.value, nrow(fshd_scores)),
             hjust = 0.5, size = 3.5) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0, face = 'bold', size = 14), aspect.ratio = 1)
} else {
  cat('Skipping FSHD (Validation/FSHD_meta.xlsx not found)\n')
}

# Plot B: CDM CTG Expansion
if (file.exists('Validation/CDM_meta.xlsx')) {
  cdm_meta <- read_excel('Validation/CDM_meta.xlsx')
  cdm_scores <- merge(cdm_meta, scores, by.x = 'Sample_ID', by.y = 'Sample_id')
  cdm_scores <- cdm_scores %>% filter(!is.na(CTG_Expansion))
  cor_cdm_ctg <- cor.test(cdm_scores$GMHS_v33, cdm_scores$CTG_Expansion, method = 'spearman')
  cat('CDM CTG Expansion: rho =', round(cor_cdm_ctg$estimate, 2),
      ', p =', round(cor_cdm_ctg$p.value, 3), ', n =', nrow(cdm_scores), '\n')
  plots[['B']] <- ggplot(cdm_scores, aes(x = GMHS_v33, y = CTG_Expansion)) +
    geom_point(color = green_color, size = 2.5, alpha = 0.7) +
    geom_smooth(method = 'lm', color = green_color, fill = green_color, alpha = 0.2) +
    labs(x = 'MyoScore (GMHS v3.3)', y = 'CDM CTG Expansion', title = 'B') +
    annotate('text', x = mean(cdm_scores$GMHS_v33),
             y = max(cdm_scores$CTG_Expansion) * 0.95,
             label = sprintf('ρ = %.2f, p = %.3f (n=%d)',
                            cor_cdm_ctg$estimate, cor_cdm_ctg$p.value, nrow(cdm_scores)),
             hjust = 0.5, size = 3.5) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0, face = 'bold', size = 14), aspect.ratio = 1)
} else {
  cat('Skipping CDM (Validation/CDM_meta.xlsx not found)\n')
}

# Plot C: LGMDR12 Mercuri Score
if (file.exists('Validation/LGMDR12_meta.xlsx')) {
  lgmdr12_meta <- read_excel('Validation/LGMDR12_meta.xlsx')
  lgmdr12_scores <- merge(lgmdr12_meta, scores, by.x = 'Sample_ID', by.y = 'Sample_id')
  lgmdr12_scores <- lgmdr12_scores %>% filter(!is.na(Mercuri_Score))
  cor_lgmdr12_mercuri <- cor.test(lgmdr12_scores$GMHS_v33, lgmdr12_scores$Mercuri_Score, method = 'spearman')
  cat('LGMDR12 Mercuri Score: rho =', round(cor_lgmdr12_mercuri$estimate, 2),
      ', p =', round(cor_lgmdr12_mercuri$p.value, 3), ', n =', nrow(lgmdr12_scores), '\n')
  plots[['C']] <- ggplot(lgmdr12_scores, aes(x = GMHS_v33, y = Mercuri_Score)) +
    geom_point(color = green_color, size = 2.5, alpha = 0.7) +
    geom_smooth(method = 'lm', color = green_color, fill = green_color, alpha = 0.2) +
    labs(x = 'MyoScore (GMHS v3.3)', y = 'LGMDR12 Mercuri Score', title = 'C') +
    annotate('text', x = mean(lgmdr12_scores$GMHS_v33),
             y = max(lgmdr12_scores$Mercuri_Score) * 0.95,
             label = sprintf('ρ = %.2f, p < 0.001 (n=%d)',
                            cor_lgmdr12_mercuri$estimate, nrow(lgmdr12_scores)),
             hjust = 0.5, size = 3.5) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0, face = 'bold', size = 14), aspect.ratio = 1)
} else {
  cat('Skipping LGMDR12 (Validation/LGMDR12_meta.xlsx not found)\n')
}

# ========== Plot D: Overweight Pre vs Post Lifestyle Intervention (GSE282733) ==========

meta <- read.csv('Meta/RNA-seq integration meta (MultiOmics).csv', row.names = 1)
ow_meta <- meta %>% filter(grepl('Overweight', Phenotype_2))

# Merge with scores
ow_scores <- merge(ow_meta, scores, by.x = 'row.names', by.y = 'Sample_id')
ow_scores <- ow_scores %>%
  mutate(
    timepoint = ifelse(grepl('Pre', Phenotype_2), 'Pre', 'Post'),
    subject_id = sub('_Pre$|_Post$', '', Row.names)
  )

# Build paired data
ow_pre <- ow_scores %>% filter(timepoint == 'Pre') %>%
  select(subject_id, GMHS_v33_pre = GMHS_v33)
ow_post <- ow_scores %>% filter(timepoint == 'Post') %>%
  select(subject_id, GMHS_v33_post = GMHS_v33)
ow_paired <- inner_join(ow_pre, ow_post, by = 'subject_id')

# Wilcoxon signed-rank test (paired)
wilcox_ow <- wilcox.test(ow_paired$GMHS_v33_pre, ow_paired$GMHS_v33_post, paired = TRUE)
cat('\nOverweight Pre vs Post (Wilcoxon paired, n =', nrow(ow_paired), 'pairs):\n')
cat('  Pre:  mean =', round(mean(ow_paired$GMHS_v33_pre), 2),
    '±', round(sd(ow_paired$GMHS_v33_pre), 2), '\n')
cat('  Post: mean =', round(mean(ow_paired$GMHS_v33_post), 2),
    '±', round(sd(ow_paired$GMHS_v33_post), 2), '\n')
cat('  Δ =', round(mean(ow_paired$GMHS_v33_post - ow_paired$GMHS_v33_pre), 2),
    ', p =', round(wilcox_ow$p.value, 4), '\n')

# Paired lines Pre→Post
ow_long <- ow_paired %>%
  pivot_longer(cols = c(GMHS_v33_pre, GMHS_v33_post),
               names_to = 'time', values_to = 'GMHS_v33') %>%
  mutate(time = factor(ifelse(time == 'GMHS_v33_pre', 'Pre', 'Post'),
                       levels = c('Pre', 'Post')))

plots[['D']] <- ggplot(ow_long, aes(x = time, y = GMHS_v33)) +
  geom_line(aes(group = subject_id), color = '#cccccc', linewidth = 0.4) +
  geom_point(aes(color = time), size = 2, alpha = 0.7) +
  geom_boxplot(aes(fill = time), width = 0.3, alpha = 0.3, outlier.shape = NA) +
  scale_color_manual(values = c('Pre' = '#50327b', 'Post' = '#f4e030')) +
  scale_fill_manual(values = c('Pre' = '#50327b', 'Post' = '#f4e030')) +
  labs(
    x = 'Lifestyle Intervention (16 weeks)',
    y = 'MyoScore (GMHS v3.3)',
    title = paste0(LETTERS[length(plots) + 1])
  ) +
  annotate('text', x = 1.5, y = max(ow_long$GMHS_v33) + 0.5,
           label = sprintf('Wilcoxon paired p = %.3f\nn = %d pairs, Δ = %+.2f (ns)',
                          wilcox_ow$p.value, nrow(ow_paired),
                          mean(ow_paired$GMHS_v33_post - ow_paired$GMHS_v33_pre)),
           hjust = 0.5, size = 3) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 14),
    legend.position = 'none',
    aspect.ratio = 1
  )

# ========== Combine and save ==========
n_plots <- length(plots)
combined <- wrap_plots(plots, nrow = 1)
w <- 4 * n_plots

ggsave('MyoScore/figures/MyoScore_clinical_validation.pdf', combined,
       width = w, height = 4, device = cairo_pdf)
ggsave('MyoScore/figures/MyoScore_clinical_validation.png', combined,
       width = w, height = 4, dpi = 300)

cat('\nPlot saved (', n_plots, 'panels) to MyoScore/figures/MyoScore_clinical_validation.pdf\n')
