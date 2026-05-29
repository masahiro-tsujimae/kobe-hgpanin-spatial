# ==============================================================================
# Figure 3 D / G 再出力
#   D: HD epithelial subtype composition stacked bar  [Visium_HD/data/fig3_subtype_counts.csv]
#   G: TP63 / KRT5 / GATA6 IHC boxplots               [Fig3 Excel; IHC_analysis.R logic]
#   (E density は HD オブジェクトが必要なため別スクリプト)
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(ggsignif); library(readxl) })
MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))   # theme_cd, colors_stage, colors_epi
OUT <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig3", "panels"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  cat(sprintf("  [panel] %-26s %.1f x %.1f mm\n", name, w, h))
}
sig_sym <- function(a, b) { p <- wilcox.test(a, b)$p.value
  if (p < 1e-4) "****" else if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < 0.05) "*" else "ns" }
fmt_kw <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
ihc_dir <- file.path(MOONSHOT, "IHC_analysis")  # de-identified per-duct IHC scores (columns: duct, score)

cat("\n================ Fig3 D/G ================\n")

# ---- G: TP63 / KRT5 / GATA6 boxplots ----
mk_ihc_box <- function(csv, ylab, name, gata = FALSE) {
  d <- read.csv(file.path(ihc_dir, csv))
  df <- data.frame(duct = d$duct, val = d$score)
  df <- df %>% filter(!is.na(duct), !is.na(val)) %>% filter(duct %in% c("N", "LG", "HG")) %>%
    mutate(duct = factor(duct, levels = c("N", "LG", "HG")))
  if (gata) { mx <- max(df$val, na.rm = TRUE); df$val <- mx - df$val }   # GATA6 expression score = max - negative
  ymax <- max(df$val, na.rm = TRUE)
  s12 <- sig_sym(df$val[df$duct=="N"], df$val[df$duct=="LG"])
  s23 <- sig_sym(df$val[df$duct=="LG"], df$val[df$duct=="HG"])
  s13 <- sig_sym(df$val[df$duct=="N"], df$val[df$duct=="HG"])
  kw <- kruskal.test(val ~ duct, df)$p.value
  cat(sprintf("  %s KW p = %.3g\n", name, kw))
  ggplot(df, aes(duct, val, fill = duct)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, linewidth = 0.4) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.0, shape = 16) +
    scale_fill_manual(values = colors_stage, guide = "none") +
    coord_cartesian(ylim = c(0, ymax * (if (gata) 1.92 else 1.72))) +
    labs(x = NULL, y = ylab) +
    geom_signif(comparisons = list(c("N","LG")), annotations = s12, y_position = ymax*1.05, tip_length = 0.02, textsize = 2.2, vjust = -0.2) +
    geom_signif(comparisons = list(c("LG","HG")), annotations = s23, y_position = ymax*1.17, tip_length = 0.02, textsize = 2.2, vjust = -0.2) +
    geom_signif(comparisons = list(c("N","HG")), annotations = s13, y_position = ymax*1.30, tip_length = 0.02, textsize = 2.2, vjust = -0.2) +
    annotate("text", x = 2, y = ymax * (if (gata) 1.78 else 1.62), label = paste("KW,", fmt_kw(kw)), size = 2.0, fontface = "italic") +
    theme_cd(base_size = 8)
}
save_panel(mk_ihc_box("fig3_TP63.csv", "TP63+ cells/HPF", "Fig3G_TP63_box"), "Fig3G_TP63_box", 38, 52)
save_panel(mk_ihc_box("fig3_KRT5.csv", "KRT5+ cells/HPF", "Fig3G_KRT5_box"), "Fig3G_KRT5_box", 38, 52)
save_panel(mk_ihc_box("fig3_GATA6.csv", "GATA6 expression score", "Fig3G_GATA6_box", gata = TRUE), "Fig3G_GATA6_box", 38, 52)

# ---- D: HD subtype composition stacked bar (from CSV) ----
sc <- read.csv(file.path(MOONSHOT, "Visium_HD/data/fig3_subtype_counts.csv"))
sample_ids <- c("7531-22H03674_11", "7532-23H11594_3", "K22H06985_7", "K23H01635_6")
sc$Epi_subtype <- factor(sc$Epi_subtype, levels = c("Classical", "Intermediate", "Basal-like"))
sc$sample_id <- factor(sc$sample_id, levels = sample_ids)
sc$diagnosis <- factor(sc$diagnosis, levels = c("LG-PanIN", "HG-PanIN"))
sc$lab <- factor(sc$sample_id, levels = sample_ids, labels = c("LG-1", "LG-2", "HG-1", "HG-2"))
fig3d <- ggplot(sc, aes(lab, pct, fill = Epi_subtype)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.3) +
  facet_grid(~ diagnosis, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colors_epi, name = "Subtype") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_cd(base_size = 8) +
  theme(legend.position = "right", legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6.5), legend.title = element_text(size = 7),
        strip.text = element_text(size = 7, face = "bold"))
save_panel(fig3d, "Fig3D_subtype_composition", 70, 45)

cat("\nFig3 D/G 完了\n")
