# Age_Youth_correlation_combined.R
# Age vs Youth Score combined figure
# A: Scatter correlation (all samples with age data)
# B: Boxplot (GTEx healthy samples by age group)

library(tidyverse)
library(clinfun)
library(patchwork)

# setwd("path/to/Myopathy_Spectrum_Multiomics") 

# ========== Plot A: Scatter plot - All samples ==========
data_all <- read.csv('MyoScore/data/Age_MyoScore_data.csv')
cat('All samples:', nrow(data_all), '\n')

# Spearman correlation
cor_test <- cor.test(data_all$Age_numeric, data_all$Youth_score, method = 'spearman')
cat('Spearman rho:', round(cor_test$estimate, 2), '\n')
cat('p-value:', cor_test$p.value, '\n')

# Plot A: Scatter correlation
p1 <- ggplot(data_all, aes(x = Age_numeric, y = Youth_score)) +
  geom_point(alpha = 0.4, color = '#f4e030', size = 1.5) +
  geom_smooth(method = 'lm', color = '#50327b', se = TRUE, fill = '#50327b', alpha = 0.2) +
  labs(
    x = 'Age (years)',
    y = 'Youth Score',
    title = 'A'
  ) +
  annotate('text', x = 55, y = 85,
           label = sprintf('All Samples (n=%d)\nρ = %.2f, p = %.3f',
                          nrow(data_all), cor_test$estimate, cor_test$p.value),
           hjust = 0.5, size = 3.5) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 14),
    aspect.ratio = 1
  ) +
  coord_cartesian(xlim = c(10, 80), ylim = c(0, 95))

# ========== Plot B: Boxplot - GTEx healthy samples (fast/unexpected death) ==========
gtex <- read.csv('GWAS_TWAS/gmhs_dimensions/gtex_clinical_validation/gtex_clinical_merged.csv')

# Filter for healthy muscle samples: fast/unexpected deaths (DTHHRDY = 1 or 2)
gtex_healthy <- gtex %>%
  filter(DTHHRDY %in% c(1, 2))

# Create ordered age factor
gtex_healthy$Age_ordered <- factor(gtex_healthy$AGE,
                           levels = c('20-29', '30-39', '40-49', '50-59', '60-69'),
                           ordered = TRUE)

# Remove NA age ranges
gtex_healthy <- gtex_healthy %>% filter(!is.na(Age_ordered))

cat('GTEx healthy samples:', nrow(gtex_healthy), '\n')

# J-T test for Youth Score
jt_youth <- jonckheere.test(gtex_healthy$Youth_score, as.numeric(gtex_healthy$Age_ordered),
                             alternative = 'two.sided')
cat('J-T test p-value:', jt_youth$p.value, '\n')

# Plot B: Boxplot GTEx healthy
p2 <- ggplot(gtex_healthy, aes(x = Age_ordered, y = Youth_score)) +
  geom_boxplot(fill = '#f4e030', alpha = 0.7, outlier.shape = NA, color = 'black') +
  geom_jitter(width = 0.2, alpha = 0.5, color = '#f4e030', size = 1.5) +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 4, color = '#50327b') +
  labs(
    x = 'Age Range',
    y = 'Youth Score',
    title = 'B'
  ) +
  annotate('text', x = 3, y = 85,
           label = sprintf('GTEx Healthy (n=%d)\nJonckheere-Terpstra Test, p = %.3f',
                          nrow(gtex_healthy), jt_youth$p.value),
           hjust = 0.5, size = 3.5) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 14),
    aspect.ratio = 1
  ) +
  coord_cartesian(ylim = c(15, 95))

# Combine plots
combined <- p1 + p2

# Save
ggsave('MyoScore/figures/Age_Youth_correlation_combined.pdf', combined, width = 10, height = 5,
       device = cairo_pdf)
ggsave('MyoScore/figures/Age_Youth_correlation_combined.png', combined, width = 10, height = 5, dpi = 300)

cat('\nPlot saved to MyoScore/figures/Age_Youth_correlation_combined.pdf\n')
