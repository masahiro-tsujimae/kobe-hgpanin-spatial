# ==============================================================================
# Figure 1 C / D 再出力 (Masson trichrome 定量, IHC_analysis.R と同一データ・ロジック)
#   C: Acinar percentage boxplot   D: Stromal percentage boxplot
# データは IHC_analysis.R にインラインで埋め込まれた値（外部ファイル不要）。
# 出力を配置実寸・大きめフォントに変更（解析・p値は不変）。
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(ggsignif) })
MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))   # theme_cd, colors_diagnosis
OUT <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig1", "panels"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  cat(sprintf("  [panel] %-26s %.1f x %.1f mm\n", name, w, h))
}
sig_sym <- function(a, b) { p <- wilcox.test(a, b)$p.value
  if (p < 1e-4) "****" else if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < 0.05) "*" else "ns" }
fmt_kw <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

stroma_data <- data.frame(
  HE_max = c(rep("HG-PanIN", 27), rep("LG-PanIN", 15), rep("Normal", 23)),
  Stroma_percent = c(
    66.5335,91.4714,57.5820,55.9365,65.9865,46.8543,89.2963,74.2051,53.6875,47.3259,
    68.6643,41.1064,82.2487,66.7571,49.6285,46.4688,68.1005,47.0395,57.3515,95.3974,
    46.8215,64.4020,65.0704,54.2374,49.0806,65.1374,47.7961,
    7.4284,26.6127,23.7258,44.1570,42.6979,5.9338,24.9566,7.9202,31.8691,29.5809,
    31.9292,41.4097,36.7947,55.3543,21.3633,
    32.3551,21.8185,63.4235,39.1845,8.9912,12.6651,5.5916,8.4104,6.1173,6.8347,
    17.8506,17.2497,24.0248,16.8840,20.1510,17.4100,13.0546,30.0977,16.7446,10.4846,
    13.6512,18.5027,23.5777))
stroma_data$HE_max <- factor(stroma_data$HE_max, levels = c("Normal", "LG-PanIN", "HG-PanIN"))
stroma_data$Acinar_percent <- 100 - stroma_data$Stroma_percent

mk_box <- function(yvar, ylab) {
  ymax <- max(stroma_data[[yvar]])
  g <- stroma_data; lv <- levels(g$HE_max)
  s12 <- sig_sym(g[[yvar]][g$HE_max==lv[1]], g[[yvar]][g$HE_max==lv[2]])
  s23 <- sig_sym(g[[yvar]][g$HE_max==lv[2]], g[[yvar]][g$HE_max==lv[3]])
  s13 <- sig_sym(g[[yvar]][g$HE_max==lv[1]], g[[yvar]][g$HE_max==lv[3]])
  kw <- kruskal.test(g[[yvar]] ~ g$HE_max)$p.value
  ggplot(g, aes(HE_max, .data[[yvar]], fill = HE_max)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.1, shape = 16) +
    scale_fill_manual(values = colors_diagnosis, guide = "none") +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    coord_cartesian(ylim = c(0, ymax * 1.72)) +
    labs(x = NULL, y = ylab) +
    geom_signif(comparisons = list(lv[1:2]), annotations = s12, y_position = ymax*1.05, tip_length = 0.02, textsize = 2.3, vjust = -0.2) +
    geom_signif(comparisons = list(lv[2:3]), annotations = s23, y_position = ymax*1.17, tip_length = 0.02, textsize = 2.3, vjust = -0.2) +
    geom_signif(comparisons = list(lv[c(1,3)]), annotations = s13, y_position = ymax*1.30, tip_length = 0.02, textsize = 2.3, vjust = -0.2) +
    annotate("text", x = 2, y = ymax*1.62, label = paste("Kruskal-Wallis,", fmt_kw(kw)), size = 2.1, fontface = "italic") +
    theme_cd(base_size = 8) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}
cat("\n================ Fig1 C/D (IHC Masson) ================\n")
cat(sprintf("Acinar KW p = %.3g ; Stroma KW p = %.3g\n",
            kruskal.test(Acinar_percent ~ HE_max, stroma_data)$p.value,
            kruskal.test(Stroma_percent ~ HE_max, stroma_data)$p.value))
save_panel(mk_box("Acinar_percent", "Acinar percentage (%)"), "Fig1C_acinar_box", 48, 55)
save_panel(mk_box("Stroma_percent", "Stromal percentage (%)"), "Fig1D_stroma_box", 48, 55)
cat("\nFig1 C/D 完了\n")
