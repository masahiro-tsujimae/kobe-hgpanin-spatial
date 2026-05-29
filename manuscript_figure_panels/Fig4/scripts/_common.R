# ==============================================================================
# Fig4 再出力 共通設定
# 配置実寸(mm)で書き出し、合体後もラベルが読めるようにする。
# 解析は再実行せず、既存の precomputed L-R スコア CSV から作図のみ再出力。
#   A-F: Visium_analysis/output/LR_pathology/*.csv  (fig4_AF.R と同一ロジック)
#   G/H/J: Visium_HD/data/lr_scores_HD.csv, lr_scores_4combination_HD.csv (fig4_GK_LR / 4combo)
#   I (network): CellChat 由来。CellChat 未インストールのため別扱い。
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(ComplexHeatmap); library(circlize); library(ggpubr)
})

MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))  # theme_cd, colors_*
OUT <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig4", "panels"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PANEL_DPI <- 600
mm2in <- function(mm) mm / 25.4

# L-R ペア表示用: アンダースコア -> ハイフン(en-dash)。小サイズで _ が打ち消し線に
# 見える問題を解消し、ligand-receptor として読みやすくする。全パネル共通で使用。
lr_label <- function(x) gsub("_", "–", x)  # – = en-dash

save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = PANEL_DPI, bg = "white")
  cat(sprintf("  [panel]   %-30s %.1f x %.1f mm\n", name, w, h))
}
save_heatmap <- function(ht, name, w, h, ...) {
  cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = mm2in(w), height = mm2in(h)); draw(ht, ...); dev.off()
  png(file.path(OUT, paste0(name, ".png")), width = mm2in(w), height = mm2in(h), units = "in", res = PANEL_DPI); draw(ht, ...); dev.off()
  cat(sprintf("  [heatmap] %-30s %.1f x %.1f mm\n", name, w, h))
}

# pathway 色（master の colors_pathway を使い、欠けはパレットで補完）
pathway_palette <- function(pathways) {
  pw <- unique(pathways)
  cols <- colors_pathway[pw]
  miss <- pw[is.na(cols)]
  if (length(miss) > 0) {
    extra <- setNames(scales::hue_pal()(length(miss)), miss)
    cols[miss] <- extra
  }
  setNames(as.character(cols), pw)
}
