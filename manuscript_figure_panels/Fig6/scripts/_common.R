# ==============================================================================
# Fig6 / Sup Fig 10 再出力 共通設定
# ------------------------------------------------------------------------------
# 目的: 各パネルを「最終図中での配置実寸(mm)」で書き出し、合体後もラベルが
#       意図した pt (6-8pt) を保つようにする。
#       解析の数値は公開版と同一（変更するのは出力寸法とフォントのみ）。
#
# 元スクリプト:
#   RNA_seq/result_RNAseq_organoids3_PSC/DEG_Pathway_Analysis/PSC_organoids3_DEG_pathway.R
#   RNA_seq/result_RNAseq_organoids4_PSC/DEG_Pathway_Analysis/PSC_conditioned_medium_DEG_pathway.R
#   RNA_seq/result_RNAseq_2d_PSC/DEG_Pathway_Analysis/PSC_2d_cellline_CM_DEG_pathway.R
# ==============================================================================

# --- 必要パッケージ (pacman は使わず直接ロード / pathway 系は不要) ---
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

# --- パス ---
MOONSHOT   <- "."
RNASEQ_DIR <- file.path(MOONSHOT, "RNA_seq")
OUTPUT_DIR <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig6", "panels")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 統一テーマ・カラー ---
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))
theme_nature <- theme_cd

# --- 出力解像度 ---
PANEL_DPI <- 600

mm2in <- function(mm) mm / 25.4

# --- ggplot パネル保存 (mm 実寸; PDF=ベクター + PNG) ---
save_panel <- function(plot, name, width_mm, height_mm) {
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".pdf")), plot,
         width = width_mm, height = height_mm, units = "mm",
         device = cairo_pdf, bg = "white")
  ggsave(file.path(OUTPUT_DIR, paste0(name, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm",
         dpi = PANEL_DPI, bg = "white")
  cat(sprintf("  [panel] %-28s %.1f x %.1f mm\n", name, width_mm, height_mm))
}

# --- ComplexHeatmap 保存 (mm 実寸; PDF + PNG)。... は draw() に渡す ---
save_heatmap <- function(ht, name, width_mm, height_mm, ...) {
  w_in <- mm2in(width_mm); h_in <- mm2in(height_mm)
  cairo_pdf(file.path(OUTPUT_DIR, paste0(name, ".pdf")), width = w_in, height = h_in)
  draw(ht, ...); dev.off()
  png(file.path(OUTPUT_DIR, paste0(name, ".png")),
      width = w_in, height = h_in, units = "in", res = PANEL_DPI)
  draw(ht, ...); dev.off()
  cat(sprintf("  [heatmap] %-26s %.1f x %.1f mm\n", name, width_mm, height_mm))
}

# --- スコア検定の検証用 (t.test と Wilcoxon を両方算出して比較) ---
verify_pvals <- function(score_df, value_col, group_col, comparisons, label) {
  cat(sprintf("\n--- p-value 検証: %s ---\n", label))
  for (cmp in comparisons) {
    g1 <- cmp[1]; g2 <- cmp[2]
    v1 <- score_df[[value_col]][score_df[[group_col]] == g1]
    v2 <- score_df[[value_col]][score_df[[group_col]] == g2]
    pt <- tryCatch(t.test(v1, v2)$p.value, error = function(e) NA)
    pw <- tryCatch(suppressWarnings(wilcox.test(v1, v2)$p.value), error = function(e) NA)
    cat(sprintf("  %-22s t.test=%.5g  wilcox=%.5g\n",
                paste(g1, "vs", g2), pt, pw))
  }
}

cat("Fig6 _common.R loaded.\n")
