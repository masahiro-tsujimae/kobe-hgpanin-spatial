# ==============================================================================
# Fig6 / Sup Fig 10 全パネル再出力 (一括実行)
# 実行: & "C:/Program Files/R/R-4.4.2/bin/Rscript.exe" --vanilla run_all.R
# ==============================================================================
sdir <- "./manuscript_figure_panels/Fig6/scripts"
source(file.path(sdir, "fig6_org3_panels.R"))   # Fig6 B / C / D
source(file.path(sdir, "fig6_org4_panels.R"))   # Fig6 F / G / H
source(file.path(sdir, "supfig10_2d_panels.R")) # SupFig10 C / D / E
cat("\n==== すべてのパネル再出力が完了しました ====\n")
