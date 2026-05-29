# ============================================================
# MASTER: Unified Theme & Color Palettes
# Cancer Discovery / Nature Medicine style
# E:\moonshot\common\theme_cancer_discovery.R
#
# Usage: source("E:/moonshot/common/theme_cancer_discovery.R")
#
# ---- フォントサイズについて ----
# base_size = 8pt は Cancer Discovery の1段組(88mm)に
# 実際の印刷サイズで出力したときに読めるサイズ。
# ggsave の width/height は必ず最終印刷寸法(mm)に合わせること。
#   1段組: width = 88mm
#   1.5段組: width = 120mm
#   2段組: width = 183mm
# この原則を守れば、文字は常に base_size pt で印刷される。
# ============================================================

library(ggplot2)
library(scales)

# ============================================================
# テーマ
# ============================================================
theme_cd <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text            = element_text(family = "sans", color = "black"),
      # タイトル: base + 1pt, センタリング
      plot.title      = element_text(size = base_size + 1, face = "bold",
                                     hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle   = element_text(size = base_size - 0.5, color = "gray30",
                                     hjust = 0.5, margin = margin(b = 3)),
      # 軸
      axis.title      = element_text(size = base_size, face = "bold"),
      axis.text       = element_text(size = base_size - 0.5, color = "black"),
      axis.line       = element_line(color = "black", linewidth = 0.4),
      axis.ticks      = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(1.5, "mm"),
      # 凡例
      legend.title    = element_text(size = base_size - 0.5, face = "bold"),
      legend.text     = element_text(size = base_size - 1),
      legend.key.size = unit(3.5, "mm"),
      legend.background = element_blank(),
      legend.key      = element_blank(),
      # パネル
      panel.background = element_blank(),
      panel.grid       = element_blank(),
      panel.border     = element_blank(),
      # ファセット
      strip.background = element_blank(),
      strip.text       = element_text(size = base_size - 0.5, face = "bold"),
      # マージン・背景
      plot.margin     = margin(4, 4, 4, 4),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# ============================================================
# カラーパレット（全解析共通・色覚バリアフリー対応）
# ============================================================

# --- 疾患ステージ (全解析共通) ---
# Normal=グレー / LG-PanIN=青 / HG-PanIN=赤 の軸で統一
colors_diagnosis <- c(
  "Normal"   = "#808080",   # neutral gray
  "LG-PanIN" = "#3182BD",   # RdYlBu blue
  "HG-PanIN" = "#CB181D"    # RdYlBu red
)
# IHC_analysis の略称にも対応
colors_stage <- c(
  "N"  = "#808080",
  "LG" = "#3182BD",
  "HG" = "#CB181D"
)

# --- 上皮サブタイプ ---
# RdYlBu パレット由来: blue=Classical, orange=Intermediate, red=Basal
colors_subtype <- c(
  "Classical-like" = "#4575B4",
  "Intermediate"   = "#D9D9D9",
  "Basal-like"     = "#D73027"
)
colors_epi <- c(
  "Classical"    = "#4575B4",
  "Intermediate" = "#D9D9D9",
  "Basal-like"   = "#D73027"
)
# Visium スペクトラム用 (連続グラデーション端点)
colors_spectrum <- c(
  "Classical" = "#4575B4",
  "Neutral"   = "#F7F7F7",
  "Basal"     = "#D73027"
)

# --- CAF サブタイプ ---
colors_caf <- c(
  "myCAF" = "#2171B5",   # blue
  "iCAF"  = "#CB181D",   # red
  "apCAF" = "#238B45",   # green
  "Mixed" = "#BDBDBD"    # gray
)

# --- 疾患ステージ 2群 ---
colors_disease <- c(
  "HG-PanIN" = "#3182BD",
  "PDAC"     = "#CB181D"
)

# --- サンプル個別 (進行グラジエント) ---
colors_sample <- c(
  "KKp048" = "#74ADD1",
  "KKp055" = "#4575B4",
  "KKp011" = "#F46D43",
  "KKp060" = "#D73027"
)

# --- RNA-seq グループカラー (PSC CM 実験) ---
cols_group_organoids4 <- c(
  "Control"          = "#808080",
  "KKp048_Early"     = "#3182BD",
  "KKp011_Advanced"  = "#CB181D",
  "KKp060_Advanced"  = "#F46D43"
)
cols_group_2d <- c(
  "Control" = "#808080",
  "Panc1"   = "#CB181D",
  "BxPC3"   = "#3182BD"
)
cols_group_organoids3 <- c(
  "Control"     = "#808080",
  "KKp48_Early" = "#3182BD"
)

# --- 方向性 (L-R, DEG) ---
colors_direction <- c(
  "Increased" = "#CB181D",
  "Decreased" = "#3182BD",
  "Up"        = "#CB181D",
  "Down"      = "#3182BD"
)

# --- スコアグレード (IHC) ---
colors_score <- c(
  "0" = "#F7F7F7",
  "1" = "#FDD49E",
  "2" = "#FC8D59",
  "3" = "#CB181D"
)

# --- ヒートマップ ---
heatmap_colors_diverging  <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
heatmap_colors_sequential <- colorRampPalette(c("#F7F7F7", "#FDD49E", "#B2182B"))(100)

# --- パスウェイ注釈 ---
colors_pathway <- c(
  "TGFb"          = "#CB181D",
  "PDGF"          = "#3182BD",
  "HGF"           = "#238B45",
  "FGF"           = "#3C5488",
  "VEGF"          = "#00B4D8",
  "IGF"           = "#45B7D1",
  "EGF"           = "#F39B7F",
  "BMP"           = "#8491B4",
  "Chemokine"     = "#91D1C2",
  "IL6"           = "#808080",
  "IL6_family"    = "#808080",
  "IL1"           = "#DC0000",
  "WNT"           = "#FDAE61",
  "Hedgehog"      = "#7570B3",
  "Integrin"      = "#B09C85",
  "ECM_Integrin"  = "#B09C85",
  "Matricellular" = "#7E6148"
)

# --- 細胞種 (Visium spot) ---
colors_celltype <- c(
  "Ductal"  = "#FDB462",
  "Acinar"  = "#80B1D3",
  "Stromal" = "#FB8072"
)

# ============================================================
# 出力ヘルパー
# ============================================================
OUTPUT_DPI <- 600

# ggplot 用 (PNG + PDF)
save_figure_base <- function(plot, filename, width, height, units = "mm",
                              output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  ggsave(file.path(output_dir, paste0(filename, ".png")),
         plot, width = width, height = height, units = units,
         dpi = OUTPUT_DPI, bg = "white")

  ggsave(file.path(output_dir, paste0(filename, ".pdf")),
         plot, width = width, height = height, units = units,
         device = cairo_pdf, bg = "white")

  cat(sprintf("  Saved: %s (.png + .pdf) [%g x %g %s]\n",
              filename, width, height, units))
}

# ComplexHeatmap 用
save_heatmap_base <- function(draw_expr, filename, width, height, units = "in",
                               output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  png(file.path(output_dir, paste0(filename, ".png")),
      width = width, height = height, units = units, res = OUTPUT_DPI)
  print(draw_expr)
  dev.off()

  pdf(file.path(output_dir, paste0(filename, ".pdf")),
      width = width, height = height)
  print(draw_expr)
  dev.off()

  cat(sprintf("  Saved: %s (.png + .pdf) [%g x %g %s]\n",
              filename, width, height, units))
}

cat("theme_cancer_discovery (master) loaded.\n")
