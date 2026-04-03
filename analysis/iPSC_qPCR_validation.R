## iPSC-to-myotube differentiation qPCR validation
## 5 novel MyoScore genes + 2 positive controls (MYOG, MYH3)
## 4 healthy donors x 6 time points (Day 0, 2, 5, 10, 15, 20)

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
library(lmerTest)


summary_data <- read_csv("iPSC_data.csv")

# --- Summary stats ---
summary_data <- all_data %>%
  group_by(Gene, Day, Direction, GeneType) %>%
  summarise(
    Mean = mean(Expression),
    SD = sd(Expression),
    SE = sd(Expression) / sqrt(n()),
    .groups = "drop"
  )

# --- Color scheme ---
gene_colors <- c(
  "MYOG"   = "#999999",
  "MYH3"   = "#bbbbbb",
  "TMEM52"  = "#f4e030",   # yellow - positive direction (healthy)
  "YWHAB"  = "#50327b",   # deep purple
  "CEP250" = "#46508b",   # blue-purple
  "RSRC2"  = "#72c95e",   # green
  "SNRPC"  = "#31848f"    # teal
)

gene_linetypes <- c(
  "MYOG" = "dashed", "MYH3" = "dashed",
  "TMEM52" = "solid", "YWHAB" = "solid",
  "CEP250" = "solid", "RSRC2" = "solid", "SNRPC" = "solid"
)

# --- Plot ---
p <- ggplot(summary_data, aes(x = Day, y = Mean, color = Gene, linetype = Gene)) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE, fill = Gene),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, shape = 16) +
  # Individual donor points
  geom_point(data = all_data, aes(x = Day, y = Expression, color = Gene),
             size = 1, alpha = 0.3, shape = 16, show.legend = FALSE) +
  scale_color_manual(values = gene_colors) +
  scale_fill_manual(values = gene_colors) +
  scale_linetype_manual(values = gene_linetypes) +
  scale_x_continuous(breaks = timepoints,
                     labels = paste0("D", timepoints)) +
  labs(
    title = "Novel MyoScore Gene Expression During iPSC-to-Myotube Differentiation",
    subtitle = "n = 4 healthy donors | Mean ± SE | Dashed = positive controls (MYOG, MYH3)",
    x = "Differentiation Day",
    y = "Relative Expression (qPCR)",
    color = "Gene", fill = "Gene", linetype = "Gene"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "right",
    legend.text = element_text(face = "italic"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# --- Faceted version (cleaner) ---
p_facet <- ggplot(summary_data, aes(x = Day, y = Mean, color = Gene)) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE, fill = Gene),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2, shape = 16) +
  geom_point(data = all_data, aes(x = Day, y = Expression),
             size = 1, alpha = 0.35, shape = 1, show.legend = FALSE) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  scale_color_manual(values = gene_colors) +
  scale_fill_manual(values = gene_colors) +
  scale_x_continuous(breaks = timepoints, labels = paste0("D", timepoints)) +
  labs(
    title = "iPSC-to-Myotube Differentiation: Novel MyoScore Genes",
    subtitle = "n = 4 healthy donors | Individual data points shown | Mean ± SE",
    x = "Differentiation Day",
    y = "Relative Expression (qPCR)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "none",
    strip.text = element_text(face = "bold.italic", size = 12),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.text = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# --- Statistical tests ---

# 1. Linear Mixed Model: Expression ~ Day + (1|Donor)
cat("\n=== Linear Mixed Model: Expression ~ Day + (1|Donor) ===\n")
lmm_results <- data.frame(Gene = character(), beta = numeric(),
                           SE = numeric(), t_val = numeric(),
                           p_value = numeric(), stringsAsFactors = FALSE)
for (g in levels(all_data$Gene)) {
  sub <- all_data[all_data$Gene == g, ]
  fit <- lmer(Expression ~ Day + (1 | Donor), data = sub)
  coefs <- summary(fit)$coefficients
  beta <- coefs["Day", "Estimate"]
  se <- coefs["Day", "Std. Error"]
  tval <- coefs["Day", "t value"]
  pval <- coefs["Day", "Pr(>|t|)"]
  sig <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
         ifelse(pval < 0.05, "*", "ns")))
  cat(sprintf("  %-8s: beta = %+.4f  SE = %.4f  t = %+.2f  p = %.2e  %s\n",
              g, beta, se, tval, pval, sig))
  lmm_results <- rbind(lmm_results, data.frame(
    Gene = g, beta = beta, SE = se, t_val = tval, p_value = pval))
}

# 2. Paired t-test: Day 0 vs Day 20
cat("\n=== Paired t-test: Day 0 vs Day 20 (n=4 donors) ===\n")
paired_results <- data.frame(Gene = character(), mean_D0 = numeric(),
                              mean_D20 = numeric(), FC = numeric(),
                              p_value = numeric(), stringsAsFactors = FALSE)
for (g in levels(all_data$Gene)) {
  d0 <- all_data %>% filter(Gene == g, Day == 0) %>%
    arrange(Donor) %>% pull(Expression)
  d20 <- all_data %>% filter(Gene == g, Day == 20) %>%
    arrange(Donor) %>% pull(Expression)
  tt <- t.test(d0, d20, paired = TRUE)
  fc <- mean(d20) / mean(d0)
  sig <- ifelse(tt$p.value < 0.001, "***", ifelse(tt$p.value < 0.01, "**",
         ifelse(tt$p.value < 0.05, "*", "ns")))
  cat(sprintf("  %-8s: D0=%.2f  D20=%.2f  FC=%.2f  p=%.4f  %s\n",
              g, mean(d0), mean(d20), fc, tt$p.value, sig))
  paired_results <- rbind(paired_results, data.frame(
    Gene = g, mean_D0 = mean(d0), mean_D20 = mean(d20),
    FC = fc, p_value = tt$p.value))
}

# --- Helper: p-value to stars ---
p_to_stars <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
}

# --- Build annotation data for timecourse plot ---
# Add paired t-test significance bracket (D0 vs D20) per gene
novel_genes <- c("TMEM52", "YWHAB", "CEP250", "RSRC2", "SNRPC")
annot_data <- paired_results %>%
  filter(Gene %in% novel_genes) %>%
  mutate(stars = sapply(p_value, p_to_stars),
         y_end = pmax(mean_D0, mean_D20) + 0.3)

# Merge LMM p-value into annotation
annot_data <- annot_data %>%
  left_join(lmm_results %>% select(Gene, lmm_p = p_value), by = "Gene") %>%
  mutate(lmm_stars = sapply(lmm_p, p_to_stars))

# --- Add gene labels at line ends (Day 20) ---
label_data <- summary_data %>%
  filter(Day == 20, Gene %in% novel_genes)

# --- Timecourse plot with annotations ---
p <- ggplot(summary_data, aes(x = Day, y = Mean, color = Gene, linetype = Gene)) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE, fill = Gene),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, shape = 16) +
  # Individual donor points
  geom_point(data = all_data, aes(x = Day, y = Expression, color = Gene),
             size = 1, alpha = 0.3, shape = 16, show.legend = FALSE) +
  # Gene labels at Day 20
  geom_text_repel(data = label_data,
                  aes(x = Day, y = Mean, label = Gene, color = Gene),
                  nudge_x = 1.5, direction = "y", hjust = 0,
                  fontface = "bold.italic", size = 3.5,
                  segment.size = 0.3, segment.color = "grey60",
                  show.legend = FALSE) +
  scale_color_manual(values = gene_colors) +
  scale_fill_manual(values = gene_colors) +
  scale_linetype_manual(values = gene_linetypes) +
  scale_x_continuous(breaks = timepoints,
                     labels = paste0("D", timepoints),
                     expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    title = "Novel MyoScore Gene Expression During iPSC-to-Myotube Differentiation",
    subtitle = "n = 4 healthy donors | Mean ± SE | Dashed = positive controls (MYOG, MYH3) | LMM: Expression ~ Day + (1|Donor)",
    x = "Differentiation Day",
    y = "Relative Expression (qPCR)",
    color = "Gene", fill = "Gene", linetype = "Gene"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "right",
    legend.text = element_text(face = "italic"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# Add paired t-test brackets (D0 vs D20) for novel genes
for (i in seq_len(nrow(annot_data))) {
  g <- annot_data$Gene[i]
  y_max <- annot_data$y_end[i] + (i - 1) * 0.15
  stars <- annot_data$stars[i]
  lmm_s <- annot_data$lmm_stars[i]
  gcol <- gene_colors[g]

  # Bracket line
  p <- p +
    annotate("segment", x = 0, xend = 20, y = y_max, yend = y_max,
             color = gcol, linewidth = 0.4, alpha = 0.6) +
    annotate("segment", x = 0, xend = 0, y = y_max - 0.08, yend = y_max,
             color = gcol, linewidth = 0.4, alpha = 0.6) +
    annotate("segment", x = 20, xend = 20, y = y_max - 0.08, yend = y_max,
             color = gcol, linewidth = 0.4, alpha = 0.6) +
    annotate("text", x = 10, y = y_max + 0.08,
             label = paste0(stars, " (paired t)  LMM: ", lmm_s),
             size = 2.8, color = gcol, fontface = "bold")
}

# --- Faceted version with stats ---
# Add LMM annotation per facet
facet_annot <- lmm_results %>%
  mutate(label = paste0("LMM: p=", formatC(p_value, format = "e", digits = 1),
                        " ", sapply(p_value, p_to_stars))) %>%
  left_join(paired_results %>% select(Gene, paired_p = p_value), by = "Gene") %>%
  mutate(paired_label = paste0("D0 vs D20: p=", formatC(paired_p, digits = 3),
                               " ", sapply(paired_p, p_to_stars)))

# Get y position for annotation (top of each facet)
facet_ypos <- summary_data %>%
  group_by(Gene) %>%
  summarise(ymax = max(Mean + SD) * 1.05, .groups = "drop")

facet_annot <- facet_annot %>%
  left_join(facet_ypos, by = "Gene")
facet_annot$Gene <- factor(facet_annot$Gene, levels = levels(all_data$Gene))

p_facet <- ggplot(summary_data, aes(x = Day, y = Mean, color = Gene)) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE, fill = Gene),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2, shape = 16) +
  geom_point(data = all_data, aes(x = Day, y = Expression),
             size = 1, alpha = 0.35, shape = 1, show.legend = FALSE) +
  # LMM annotation
  geom_text(data = facet_annot, aes(x = 10, y = ymax, label = label),
            color = "grey30", size = 3, fontface = "bold", show.legend = FALSE) +
  # Paired t-test annotation
  geom_text(data = facet_annot, aes(x = 10, y = ymax * 0.90, label = paired_label),
            color = "grey45", size = 2.7, show.legend = FALSE) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  scale_color_manual(values = gene_colors) +
  scale_fill_manual(values = gene_colors) +
  scale_x_continuous(breaks = timepoints, labels = paste0("D", timepoints)) +
  labs(
    title = "iPSC-to-Myotube Differentiation: Novel MyoScore Genes",
    subtitle = "n = 4 healthy donors | LMM: Expression ~ Day + (1|Donor) | Paired t-test: Day 0 vs Day 20",
    x = "Differentiation Day",
    y = "Relative Expression (qPCR)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "none",
    strip.text = element_text(face = "bold.italic", size = 12),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.text = element_text(size = 10),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# --- Save ---
out_dir <- "MyoScore/figures"

ggsave(file.path(out_dir, "iPSC_qPCR_timecourse.pdf"), p,
       width = 13, height = 7)
ggsave(file.path(out_dir, "iPSC_qPCR_timecourse.png"), p,
       width = 13, height = 7, dpi = 300)

ggsave(file.path(out_dir, "iPSC_qPCR_faceted.pdf"), p_facet,
       width = 12, height = 6)
ggsave(file.path(out_dir, "iPSC_qPCR_faceted.png"), p_facet,
       width = 12, height = 6, dpi = 300)

# --- Save data ---
write.csv(all_data,
  "MyoScore/data/iPSC_qPCR.csv",
  row.names = FALSE)

# --- Save stats table ---
stats_out <- lmm_results %>%
  rename(LMM_beta = beta, LMM_SE = SE, LMM_t = t_val, LMM_p = p_value) %>%
  left_join(paired_results %>% select(Gene, D0_mean = mean_D0, D20_mean = mean_D20,
                                       FC_D20_D0 = FC, paired_t_p = p_value),
            by = "Gene")
write.csv(stats_out,
  "MyoScore/data/iPSC_qPCR_statistics.csv",
  row.names = FALSE)

cat("Done! Saved figures, data, and statistics.\n")
