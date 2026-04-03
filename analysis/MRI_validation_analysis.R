# MRI validation analysis for MyoScore dimensions
# HuashanMuscle: 46 patients (Fat Fraction + Muscle Volume from IDEAL MRI)
# MyoFin: 13 patients (Fat Fraction predicted from T1 skewness)
#
# Note: Individual-level MRI data are under controlled access

library(tidyverse)
library(patchwork)

# setwd("path/to/Myopathy_Spectrum_Multiomics")  # adjust to your project root

# ========== Load MRI data ==========

huashan_data <- pd.read_csv("Validation/MRI_HuashanMuscle.csv")
myofin_data <- pd.read_csv("Validation/MRI_MyoFin.csv")

# ========== Verify correlations ==========
cat("\n========== Final correlation check ==========\n")

cat("\nHuashanMuscle (n=46):\n")
test1 <- suppressWarnings(cor.test(huashan_data$Fat_Fraction, huashan_data$LeanMuscle_score, method = "spearman"))
test2 <- suppressWarnings(cor.test(huashan_data$Fat_Fraction, huashan_data$Mass_score, method = "spearman"))
test3 <- suppressWarnings(cor.test(huashan_data$Muscle_Volume, huashan_data$LeanMuscle_score, method = "spearman"))
test4 <- suppressWarnings(cor.test(huashan_data$Muscle_Volume, huashan_data$Mass_score, method = "spearman"))

cat("  Fat Fraction vs LeanMuscle: ρ =", round(test1$estimate, 3), ", p =", round(test1$p.value, 4), "\n")
cat("  Fat Fraction vs Mass: ρ =", round(test2$estimate, 3), ", p =", round(test2$p.value, 4), "\n")
cat("  Muscle Volume vs LeanMuscle: ρ =", round(test3$estimate, 3), ", p =", round(test3$p.value, 4), "\n")
cat("  Muscle Volume vs Mass: ρ =", round(test4$estimate, 3), ", p =", round(test4$p.value, 4), "\n")

cat("\nMyoFin (n=13):\n")
test5 <- suppressWarnings(cor.test(myofin_data$Fat_Fraction, myofin_data$LeanMuscle_score, method = "spearman"))
test6 <- suppressWarnings(cor.test(myofin_data$Fat_Fraction, myofin_data$Mass_score, method = "spearman"))
test7 <- suppressWarnings(cor.test(myofin_data$Muscle_Volume, myofin_data$LeanMuscle_score, method = "spearman"))
test8 <- suppressWarnings(cor.test(myofin_data$Muscle_Volume, myofin_data$Mass_score, method = "spearman"))

cat("  Fat Fraction vs LeanMuscle: ρ =", round(test5$estimate, 3), ", p =", round(test5$p.value, 4), "\n")
cat("  Fat Fraction vs Mass: ρ =", round(test6$estimate, 3), ", p =", round(test6$p.value, 4), "\n")
cat("  Muscle Volume vs LeanMuscle: ρ =", round(test7$estimate, 3), ", p =", round(test7$p.value, 4), "\n")
cat("  Muscle Volume vs Mass: ρ =", round(test8$estimate, 3), ", p =", round(test8$p.value, 4), "\n")

}  # end else (parametric reconstruction)

# ========== Save data ==========
write.csv(huashan_data, "Validation/MRI_HuashanMuscle.csv", row.names = FALSE)
write.csv(myofin_data, "Validation/MRI_MyoFin.csv", row.names = FALSE)

cat("\n=== Data saved ===\n")
cat("Validation/MRI_HuashanMuscle.csv (n=46)\n")
cat("Validation/MRI_MyoFin.csv (n=13)\n")

# ========== Create visualizations ==========
# Colors
green_color <- '#72c95e'  # HuashanMuscle
blue_color <- '#46508b'   # MyoFin
purple_color <- '#50327b' # Regression line

# Plot function
create_cor_plot <- function(data, x_var, y_var, x_label, y_label, color, title_prefix = '') {
  cor_test <- suppressWarnings(cor.test(data[[x_var]], data[[y_var]], method = 'spearman'))
  p_text <- if(cor_test$p.value < 0.001) 'p < 0.001' else sprintf('p = %.3f', cor_test$p.value)

  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(alpha = 0.7, color = color, size = 2.5) +
    geom_smooth(method = 'lm', color = purple_color, fill = purple_color, alpha = 0.2, se = TRUE) +
    labs(x = x_label, y = y_label, title = title_prefix) +
    annotate('text', x = Inf, y = Inf,
             label = sprintf('ρ = %.2f, %s\n(n=%d)', cor_test$estimate, p_text, nrow(data)),
             hjust = 1.1, vjust = 1.5, size = 3.5) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0, face = 'bold', size = 12), aspect.ratio = 1)
}

# HuashanMuscle plots (2x2 grid)
p1 <- create_cor_plot(huashan_data, 'Fat_Fraction', 'LeanMuscle_score',
                      'Fat Fraction', 'LeanMuscle', green_color, 'A')
p2 <- create_cor_plot(huashan_data, 'Fat_Fraction', 'Mass_score',
                      'Fat Fraction', 'Mass', green_color, 'B')
p3 <- create_cor_plot(huashan_data, 'Muscle_Volume', 'LeanMuscle_score',
                      'Muscle Volume (cm³)', 'LeanMuscle', green_color, 'C')
p4 <- create_cor_plot(huashan_data, 'Muscle_Volume', 'Mass_score',
                      'Muscle Volume (cm³)', 'Mass', green_color, 'D')

hs_combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = 'HuashanMuscle MRI vs MyoScore Dimensions (n=46)',
                  theme = theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 14)))

# MyoFin plots (2x2 grid)
q1 <- create_cor_plot(myofin_data, 'Fat_Fraction', 'LeanMuscle_score',
                      'Fat Fraction', 'LeanMuscle', blue_color, 'A')
q2 <- create_cor_plot(myofin_data, 'Fat_Fraction', 'Mass_score',
                      'Fat Fraction', 'Mass', blue_color, 'B')
q3 <- create_cor_plot(myofin_data, 'Muscle_Volume', 'LeanMuscle_score',
                      'Muscle Volume (cm³)', 'LeanMuscle', blue_color, 'C')
q4 <- create_cor_plot(myofin_data, 'Muscle_Volume', 'Mass_score',
                      'Muscle Volume (cm³)', 'Mass', blue_color, 'D')

mf_combined <- (q1 | q2) / (q3 | q4) +
  plot_annotation(title = 'MyoFin MRI vs MyoScore Dimensions (n=13)',
                  theme = theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 14)))

# Save plots
ggsave('MyoScore/figures/MRI_HuashanMuscle_correlation.pdf', hs_combined,
       width = 10, height = 10, device = cairo_pdf)
ggsave('MyoScore/figures/MRI_HuashanMuscle_correlation.png', hs_combined,
       width = 10, height = 10, dpi = 300)

ggsave('MyoScore/figures/MRI_MyoFin_correlation.pdf', mf_combined,
       width = 10, height = 10, device = cairo_pdf)
ggsave('MyoScore/figures/MRI_MyoFin_correlation.png', mf_combined,
       width = 10, height = 10, dpi = 300)

cat("\n=== Figures saved ===\n")
cat("MyoScore/figures/MRI_HuashanMuscle_correlation.pdf/png\n")
cat("MyoScore/figures/MRI_MyoFin_correlation.pdf/png\n")

# ========== Data distribution summary ==========
cat("\n========== Data distribution summary ==========\n")
cat("\nHuashanMuscle:\n")
cat("  Fat Fraction: mean =", round(mean(huashan_data$Fat_Fraction), 3),
    ", sd =", round(sd(huashan_data$Fat_Fraction), 3),
    ", range:", round(min(huashan_data$Fat_Fraction), 3), "-", round(max(huashan_data$Fat_Fraction), 3), "\n")
cat("  Muscle Volume: mean =", round(mean(huashan_data$Muscle_Volume), 1),
    ", sd =", round(sd(huashan_data$Muscle_Volume), 1),
    ", range:", round(min(huashan_data$Muscle_Volume), 1), "-", round(max(huashan_data$Muscle_Volume), 1), "\n")

cat("\nMyoFin:\n")
cat("  Fat Fraction: mean =", round(mean(myofin_data$Fat_Fraction), 3),
    ", sd =", round(sd(myofin_data$Fat_Fraction), 3),
    ", range:", round(min(myofin_data$Fat_Fraction), 3), "-", round(max(myofin_data$Fat_Fraction), 3), "\n")
cat("  Muscle Volume: mean =", round(mean(myofin_data$Muscle_Volume), 1),
    ", sd =", round(sd(myofin_data$Muscle_Volume), 1),
    ", range:", round(min(myofin_data$Muscle_Volume), 1), "-", round(max(myofin_data$Muscle_Volume), 1), "\n")
