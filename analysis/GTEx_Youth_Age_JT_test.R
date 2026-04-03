# GTEx_Youth_Age_JT_test.R
# GTEx healthy muscle samples: Youth Score vs Age Range
# Jonckheere-Terpstra trend test

library(tidyverse)
library(clinfun)

# setwd("path/to/Myopathy_Spectrum_Multiomics")  # adjust to your project root

# Read GTEx clinical data with scores
gtex <- read.csv('GWAS_TWAS/gmhs_dimensions/gtex_clinical_validation/gtex_clinical_merged.csv')

# Filter for healthy muscle samples: fast/unexpected deaths (DTHHRDY = 1 or 2)
# DTHHRDY mapping:
# 0 = Ventilator Case
# 1 = Fast death of violent causes
# 2 = Fast death of natural causes
# 3 = Intermediate death
# 4 = Slow death
gtex_healthy <- gtex %>%
  filter(DTHHRDY %in% c(1, 2))

cat('Total GTEx samples:', nrow(gtex), '\n')
cat('Healthy (fast/unexpected death) samples:', nrow(gtex_healthy), '\n\n')

# Create ordered age factor
gtex_healthy$Age_ordered <- factor(gtex_healthy$AGE,
                           levels = c('20-29', '30-39', '40-49', '50-59', '60-69'),
                           ordered = TRUE)

# Remove NA age ranges
gtex_healthy <- gtex_healthy %>% filter(!is.na(Age_ordered))

cat('Samples with valid age range:', nrow(gtex_healthy), '\n')
cat('Age range distribution:\n')
print(table(gtex_healthy$Age_ordered))

# Summary stats
summary_stats <- gtex_healthy %>%
  group_by(Age_ordered) %>%
  summarise(
    n = n(),
    mean_Youth = mean(Youth_score),
    sd_Youth = sd(Youth_score),
    se_Youth = sd(Youth_score) / sqrt(n())
  )
print(summary_stats)

# Jonckheere-Terpstra Test
jt_result <- jonckheere.test(gtex_healthy$Youth_score, as.numeric(gtex_healthy$Age_ordered),
                              alternative = 'two.sided')

cat('\nJonckheere-Terpstra Test Results (Healthy GTEx):\n')
print(jt_result)

# Create plot
p <- ggplot(gtex_healthy, aes(x = Age_ordered, y = Youth_score)) +
  geom_boxplot(fill = '#46508b', alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, color = '#50327b', size = 1) +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 4, color = '#f4e030') +
  labs(
    x = 'Age Range',
    y = 'Youth Score',
    title = 'GTEx Healthy Muscle: Youth Score by Age'
  ) +
  annotate('text', x = 3, y = 85,
           label = paste0('Jonckheere-Terpstra Test\np = ',
                         format(round(jt_result$p.value, 3), nsmall = 3)),
           hjust = 0.5, size = 3.5) +
  annotate('text', x = 3, y = 78,
           label = paste0('n = ', nrow(gtex_healthy), ' (fast/unexpected death)'),
           hjust = 0.5, size = 3, color = 'gray40') +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    aspect.ratio = 1,
    text = element_text(family = 'sans')
  )

# Save
ggsave('MyoScore/figures/GTEx_Youth_Age_JT_test.pdf', p, width = 5, height = 5,
       device = cairo_pdf)
ggsave('MyoScore/figures/GTEx_Youth_Age_JT_test.png', p, width = 5, height = 5, dpi = 300)

cat('\nPlot saved to MyoScore/figures/GTEx_Youth_Age_JT_test.pdf\n')
