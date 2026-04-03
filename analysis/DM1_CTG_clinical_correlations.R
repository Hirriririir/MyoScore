# DM1_CTG_clinical_correlations.R
# DM1 CTG Repeats vs MyoScore, 10MWT, Grip Strength correlations
# DM1 clinical correlations (n=27)

library(tidyverse)

# setwd("path/to/Myopathy_Spectrum_Multiomics")  

# Read DM1 clinical data
dm1_data <- read.csv('MyoScore/data/HuashanMuscle_DM1_clinical_meta.csv')

cat('DM1 samples:', nrow(dm1_data), '\n')

# Calculate correlations
cor_ctg_myoscore <- cor.test(dm1_data$GMHS_v33, dm1_data$CTG_numeric, method = 'spearman')
cor_myoscore_10mwt <- cor.test(dm1_data$GMHS_v33, dm1_data$`X10MWT_sec`, method = 'spearman')
cor_myoscore_grip <- cor.test(dm1_data$GMHS_v33, dm1_data$Grip_strength_kg, method = 'spearman')

cat('\nCorrelations:\n')
cat('MyoScore vs CTG: rho =', round(cor_ctg_myoscore$estimate, 2),
    ', p =', round(cor_ctg_myoscore$p.value, 3), '\n')
cat('MyoScore vs 10MWT: rho =', round(cor_myoscore_10mwt$estimate, 2),
    ', p =', round(cor_myoscore_10mwt$p.value, 3), '\n')
cat('MyoScore vs Grip: rho =', round(cor_myoscore_grip$estimate, 2),
    ', p =', round(cor_myoscore_grip$p.value, 3), '\n')

# Purple color for all points
purple_color <- '#50327b'

# Plot A: MyoScore vs CTG Repeats
p1 <- ggplot(dm1_data, aes(x = GMHS_v33, y = CTG_numeric)) +
  geom_point(color = purple_color, size = 2.5, alpha = 0.7) +
  geom_smooth(method = 'lm', color = purple_color, fill = purple_color, alpha = 0.2) +
  labs(
    x = 'MyoScore (GMHS v3.3)',
    y = 'CTG Repeats',
    title = 'A'
  ) +
  annotate('text', x = mean(dm1_data$GMHS_v33), y = max(dm1_data$CTG_numeric) * 0.95,
           label = sprintf('ρ = %.2f, p = %.3f',
                          cor_ctg_myoscore$estimate, cor_ctg_myoscore$p.value),
           hjust = 0.5, size = 3.5) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 14),
    aspect.ratio = 1
  )

# Plot B: MyoScore vs 10MWT
p2 <- ggplot(dm1_data, aes(x = GMHS_v33, y = `X10MWT_sec`)) +
  geom_point(color = purple_color, size = 2.5, alpha = 0.7) +
  geom_smooth(method = 'lm', color = purple_color, fill = purple_color, alpha = 0.2) +
  labs(
    x = 'MyoScore (GMHS v3.3)',
    y = '10-Meter Walk Test (sec)',
    title = 'B'
  ) +
  annotate('text', x = mean(dm1_data$GMHS_v33), y = max(dm1_data$`X10MWT_sec`) * 0.95,
           label = sprintf('ρ = %.2f, p = %.3f',
                          cor_myoscore_10mwt$estimate, cor_myoscore_10mwt$p.value),
           hjust = 0.5, size = 3.5) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 14),
    aspect.ratio = 1
  )

# Plot C: MyoScore vs Grip Strength
p3 <- ggplot(dm1_data, aes(x = GMHS_v33, y = Grip_strength_kg)) +
  geom_point(color = purple_color, size = 2.5, alpha = 0.7) +
  geom_smooth(method = 'lm', color = purple_color, fill = purple_color, alpha = 0.2) +
  labs(
    x = 'MyoScore (GMHS v3.3)',
    y = 'Grip Strength (kg)',
    title = 'C'
  ) +
  annotate('text', x = mean(dm1_data$GMHS_v33), y = max(dm1_data$Grip_strength_kg) * 0.95,
           label = sprintf('ρ = %.2f, p = %.3f',
                          cor_myoscore_grip$estimate, cor_myoscore_grip$p.value),
           hjust = 0.5, size = 3.5) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = 'bold', size = 14),
    aspect.ratio = 1
  )

# Combine plots
library(patchwork)
combined <- p1 + p2 + p3

# Save
ggsave('MyoScore/figures/DM1_CTG_clinical_correlations.pdf', combined,
       width = 12, height = 4, device = cairo_pdf)
ggsave('MyoScore/figures/DM1_CTG_clinical_correlations.png', combined,
       width = 12, height = 4, dpi = 300)

cat('\nPlot saved to MyoScore/figures/DM1_CTG_clinical_correlations.pdf\n')
