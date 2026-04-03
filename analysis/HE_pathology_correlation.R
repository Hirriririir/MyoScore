# HE_pathology_correlation_full.R
# GTEx and HuashanMuscle HE pathology quantification vs MyoScore dimensions (4x4 full grid)
# Pathology variables: Fat, Fibrosis, Fiber CV, Nuclear Centralization
# GTEx: yellow (#f4e030)
# HuashanMuscle: green (#72c95e)

library(tidyverse)
library(patchwork)

# setwd("path/to/Myopathy_Spectrum_Multiomics")  # adjust to your project root

# Read HE data
he_data <- read.csv('HE slide results/HE slide results [Shanghai+GTEx] 20260103.csv')

# Read MyoScore data
scores <- read.csv('GWAS_TWAS/gmhs_dimensions/dimension_scores_v33.csv', row.names = 1)
scores$Sample_id <- rownames(scores)

# ========== Merge data ==========
# HuashanMuscle (Shanghai) samples
shanghai_he <- he_data %>% filter(grepl('^M[0-9]', sample_id))
shanghai_merged <- merge(shanghai_he, scores, by.x = 'sample_id', by.y = 'Sample_id')
shanghai_merged$Source <- 'HuashanMuscle'

# GTEx samples
gtex_he <- he_data %>% filter(grepl('^GTEX', sample_id))
gtex_scores <- scores %>% filter(grepl('^GTEX', Sample_id))
gtex_scores$he_key <- gsub('-SM-.*', '', gtex_scores$Sample_id)
gtex_merged <- merge(gtex_he, gtex_scores, by.x = 'sample_id', by.y = 'he_key')
gtex_merged$Source <- 'GTEx'

# Combine
all_data <- bind_rows(shanghai_merged, gtex_merged)

# Remove outlier: GTEX-117YX-2526 has fiber_cv = 15.7 (far exceeds median 0.46)
all_data <- all_data %>% filter(fiber_cv < 10)

# Colors
yellow_color <- '#f4e030'
purple_color <- '#50327b'
green_color <- '#72c95e'

# ========== Create correlation plot function ==========
create_cor_plot <- function(data, x_var, y_var, x_label, y_label, color, title_prefix = '') {
  cor_test <- cor.test(data[[x_var]], data[[y_var]], method = 'spearman')

  p_text <- if(cor_test$p.value < 0.001) 'p < 0.001' else sprintf('p = %.3f', cor_test$p.value)

  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(alpha = 0.5, color = color, size = 1.5) +
    geom_smooth(method = 'lm', color = purple_color, fill = purple_color, alpha = 0.2) +
    labs(x = x_label, y = y_label, title = title_prefix) +
    annotate('text', x = Inf, y = Inf,
             label = sprintf('ρ = %.2f, %s\n(n=%d)', cor_test$estimate, p_text, nrow(data)),
             hjust = 1.1, vjust = 1.5, size = 3) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0, face = 'bold', size = 12), aspect.ratio = 1)
}

# ========== GTEx plots - Full 4x4 grid ==========
gtex_data <- all_data %>% filter(Source == 'GTEx')

# Row 1: Fat
p1 <- create_cor_plot(gtex_data, 'fat_area_pct', 'GMHS_v33', 'Fat (%)', 'MyoScore', yellow_color, 'A')
p2 <- create_cor_plot(gtex_data, 'fat_area_pct', 'LeanMuscle_score', 'Fat (%)', 'LeanMuscle', yellow_color, 'B')
p3 <- create_cor_plot(gtex_data, 'fat_area_pct', 'Youth_score', 'Fat (%)', 'Youth', yellow_color, 'C')
p4 <- create_cor_plot(gtex_data, 'fat_area_pct', 'Resilience_score', 'Fat (%)', 'Resilience', yellow_color, 'D')

# Row 2: Fibrosis
p5 <- create_cor_plot(gtex_data, 'fibrosis_pct', 'GMHS_v33', 'Fibrosis (%)', 'MyoScore', yellow_color, 'E')
p6 <- create_cor_plot(gtex_data, 'fibrosis_pct', 'LeanMuscle_score', 'Fibrosis (%)', 'LeanMuscle', yellow_color, 'F')
p7 <- create_cor_plot(gtex_data, 'fibrosis_pct', 'Youth_score', 'Fibrosis (%)', 'Youth', yellow_color, 'G')
p8 <- create_cor_plot(gtex_data, 'fibrosis_pct', 'Resilience_score', 'Fibrosis (%)', 'Resilience', yellow_color, 'H')

# Row 3: Fiber CV
p9 <- create_cor_plot(gtex_data, 'fiber_cv', 'GMHS_v33', 'Fiber CV', 'MyoScore', yellow_color, 'I')
p10 <- create_cor_plot(gtex_data, 'fiber_cv', 'LeanMuscle_score', 'Fiber CV', 'LeanMuscle', yellow_color, 'J')
p11 <- create_cor_plot(gtex_data, 'fiber_cv', 'Youth_score', 'Fiber CV', 'Youth', yellow_color, 'K')
p12 <- create_cor_plot(gtex_data, 'fiber_cv', 'Resilience_score', 'Fiber CV', 'Resilience', yellow_color, 'L')

# Row 4: Nuclear Centralization
p13 <- create_cor_plot(gtex_data, 'nuclear_centralization_index', 'GMHS_v33', 'Nuclear Central.', 'MyoScore', yellow_color, 'M')
p14 <- create_cor_plot(gtex_data, 'nuclear_centralization_index', 'LeanMuscle_score', 'Nuclear Central.', 'LeanMuscle', yellow_color, 'N')
p15 <- create_cor_plot(gtex_data, 'nuclear_centralization_index', 'Youth_score', 'Nuclear Central.', 'Youth', yellow_color, 'O')
p16 <- create_cor_plot(gtex_data, 'nuclear_centralization_index', 'Resilience_score', 'Nuclear Central.', 'Resilience', yellow_color, 'P')

# Combine GTEx plots (4x4 grid)
gtex_combined <- (p1 | p2 | p3 | p4) / (p5 | p6 | p7 | p8) / (p9 | p10 | p11 | p12) / (p13 | p14 | p15 | p16) +
  plot_annotation(title = 'GTEx HE Pathology vs MyoScore Dimensions (n=399)',
                  theme = theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 14)))

ggsave('MyoScore/figures/GTEx_HE_pathology_correlation_full.pdf', gtex_combined, width = 14, height = 14, device = cairo_pdf)
ggsave('MyoScore/figures/GTEx_HE_pathology_correlation_full.png', gtex_combined, width = 14, height = 14, dpi = 300)
cat('GTEx full plot saved\n')

# ========== HuashanMuscle plots - Full 4x4 grid ==========
hs_data <- all_data %>% filter(Source == 'HuashanMuscle')

# Row 1: Fat
q1 <- create_cor_plot(hs_data, 'fat_area_pct', 'GMHS_v33', 'Fat (%)', 'MyoScore', green_color, 'A')
q2 <- create_cor_plot(hs_data, 'fat_area_pct', 'LeanMuscle_score', 'Fat (%)', 'LeanMuscle', green_color, 'B')
q3 <- create_cor_plot(hs_data, 'fat_area_pct', 'Youth_score', 'Fat (%)', 'Youth', green_color, 'C')
q4 <- create_cor_plot(hs_data, 'fat_area_pct', 'Resilience_score', 'Fat (%)', 'Resilience', green_color, 'D')

# Row 2: Fibrosis
q5 <- create_cor_plot(hs_data, 'fibrosis_pct', 'GMHS_v33', 'Fibrosis (%)', 'MyoScore', green_color, 'E')
q6 <- create_cor_plot(hs_data, 'fibrosis_pct', 'LeanMuscle_score', 'Fibrosis (%)', 'LeanMuscle', green_color, 'F')
q7 <- create_cor_plot(hs_data, 'fibrosis_pct', 'Youth_score', 'Fibrosis (%)', 'Youth', green_color, 'G')
q8 <- create_cor_plot(hs_data, 'fibrosis_pct', 'Resilience_score', 'Fibrosis (%)', 'Resilience', green_color, 'H')

# Row 3: Fiber CV
q9 <- create_cor_plot(hs_data, 'fiber_cv', 'GMHS_v33', 'Fiber CV', 'MyoScore', green_color, 'I')
q10 <- create_cor_plot(hs_data, 'fiber_cv', 'LeanMuscle_score', 'Fiber CV', 'LeanMuscle', green_color, 'J')
q11 <- create_cor_plot(hs_data, 'fiber_cv', 'Youth_score', 'Fiber CV', 'Youth', green_color, 'K')
q12 <- create_cor_plot(hs_data, 'fiber_cv', 'Resilience_score', 'Fiber CV', 'Resilience', green_color, 'L')

# Row 4: Nuclear Centralization
q13 <- create_cor_plot(hs_data, 'nuclear_centralization_index', 'GMHS_v33', 'Nuclear Central.', 'MyoScore', green_color, 'M')
q14 <- create_cor_plot(hs_data, 'nuclear_centralization_index', 'LeanMuscle_score', 'Nuclear Central.', 'LeanMuscle', green_color, 'N')
q15 <- create_cor_plot(hs_data, 'nuclear_centralization_index', 'Youth_score', 'Nuclear Central.', 'Youth', green_color, 'O')
q16 <- create_cor_plot(hs_data, 'nuclear_centralization_index', 'Resilience_score', 'Nuclear Central.', 'Resilience', green_color, 'P')

# Combine HuashanMuscle plots (4x4 grid)
hs_combined <- (q1 | q2 | q3 | q4) / (q5 | q6 | q7 | q8) / (q9 | q10 | q11 | q12) / (q13 | q14 | q15 | q16) +
  plot_annotation(title = 'HuashanMuscle HE Pathology vs MyoScore Dimensions (n=74)',
                  theme = theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 14)))

ggsave('MyoScore/figures/HuashanMuscle_HE_pathology_correlation_full.pdf', hs_combined, width = 14, height = 14, device = cairo_pdf)
ggsave('MyoScore/figures/HuashanMuscle_HE_pathology_correlation_full.png', hs_combined, width = 14, height = 14, dpi = 300)
cat('HuashanMuscle full plot saved\n')

cat('\nAll plots saved to MyoScore/figures/\n')
