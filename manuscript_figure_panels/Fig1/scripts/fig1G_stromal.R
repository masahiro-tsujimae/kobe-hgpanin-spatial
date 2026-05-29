# ==============================================================================
# Fig1 G 再出力: Visium 間質割合 boxplot (N/LG/HG, Jonckheere-Terpstra)
#   入力: Visium_analysis/output/verification/pathology_annotation_summary.csv
#         （非-Other 分母 = Ductal+Acinar+Stromal+Endo+Fat+Lymph の Stromal_pct）
#   検証: Jonckheere P=0.0357 ≈ 公開 0.036。
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(clinfun) })
MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))   # theme_cd, colors_diagnosis
OUT <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig1", "panels"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  cat(sprintf("  [panel] %-26s %.1f x %.1f mm\n", name, w, h))
}
d <- read.csv(file.path(MOONSHOT, "Visium_analysis/output/verification/pathology_annotation_summary.csv"))
d$diagnosis_category <- factor(d$diagnosis_category, levels = c("Normal", "LG-PanIN", "HG-PanIN"))
jp <- suppressWarnings(jonckheere.test(d$Stromal_pct, as.numeric(d$diagnosis_category), alternative = "increasing")$p.value)
cat(sprintf("\n================ Fig1 G ================\nJonckheere P=%.4f (pub 0.036)\n", jp))
ymax <- max(d$Stromal_pct)
fig1g <- ggplot(d, aes(diagnosis_category, Stromal_pct, fill = diagnosis_category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
  scale_fill_manual(values = colors_diagnosis, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22)), limits = c(0, NA)) +
  annotate("text", x = 2, y = ymax * 1.14, label = sprintf("Jonckheere-Terpstra\nP=%.3f", jp), size = 2.2, color = "gray30") +
  labs(x = NULL, y = "Stromal proportion (%)") +
  theme_cd(base_size = 8) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_panel(fig1g, "Fig1G_stromal_box", 53, 55)
cat("\nFig1 G 完了\n")
