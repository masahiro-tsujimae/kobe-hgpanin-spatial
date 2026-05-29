# ==============================================================================
# Fig4 全パネル再出力 (一括実行)  ※ panel I (CellChat network) は別扱い・未生成
# 実行: & "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla run_all.R
# ==============================================================================
sdir <- "./manuscript_figure_panels/Fig4/scripts"
source(file.path(sdir, "fig4_AF_panels.R"))   # A B C D E F (standard Visium)
source(file.path(sdir, "fig4_GHJ_panels.R"))  # G H J (Visium HD)
cat("\n==== Fig4 A-H, J 再出力完了 (I=CellChat network は別扱い) ====\n")
