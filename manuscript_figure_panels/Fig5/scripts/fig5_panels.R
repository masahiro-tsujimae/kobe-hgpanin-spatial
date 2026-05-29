# ==============================================================================
# Figure 5  C / D / E / F / G / H  再出力 (scRNAseq organoid)
# ------------------------------------------------------------------------------
# 元データ(正): past_analysis/Fig_organoid_singlecell_20251130/merged_subtype_annotated.rds
#   → 本文 legend の組成・細胞数を完全再現する“正”のオブジェクト。
#   （scRNAseq/data/02_*.rds は別系統の新しい再解析でズレるため不使用）
# 元ロジック: 同フォルダ 3_Fig.R (C/D/E/F) と scRNAseq/04 + 7_ (G/H heatmap)
# 出力: 各パネルを Figure5_ver4.pptx 計測の配置実寸で、ラベルを大きく再出力。
#       Intermediate は公開版どおり黄 (#FEE090)。
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(ggridges); library(pheatmap); library(scales); library(grid)
})

MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))  # theme_cd
OUT <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig5", "panels"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PANEL_DPI <- 600
mm2in <- function(mm) mm / 25.4

save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  # PNG via ragg (agg_png): much better anti-aliasing for dense, small,
  # semi-transparent scatter points so the raster no longer looks faded vs PDF.
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = PANEL_DPI, bg = "white", device = ragg::agg_png)
  cat(sprintf("  [panel]   %-26s %.1f x %.1f mm\n", name, w, h))
}
save_ph <- function(ph, name, w, h) {
  png(file.path(OUT, paste0(name, ".png")), width = w, height = h, units = "mm", res = PANEL_DPI)
  grid.newpage(); grid.draw(ph$gtable); dev.off()
  cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = mm2in(w), height = mm2in(h))
  grid.newpage(); grid.draw(ph$gtable); dev.off()
  cat(sprintf("  [heatmap] %-26s %.1f x %.1f mm\n", name, w, h))
}

# --- palettes (公開版どおり) ---
cols_subtype <- c("Classical-like" = "#4393C3", "Intermediate" = "#FEE090", "Basal-like" = "#D6604D")
cols_sample  <- c("KKp048" = "#92C5DE", "KKp055" = "#4393C3", "KKp011" = "#F4A582", "KKp060" = "#B2182B")
sample_order <- c("KKp048", "KKp055", "KKp011", "KKp060")
labs_sample  <- c("KKp048" = "KKp048\n(HG-PanIN)", "KKp055" = "KKp055\n(HG-PanIN)",
                  "KKp011" = "KKp011\n(Inv.PDAC)", "KKp060" = "KKp060\n(Inv.PDAC)")

# --- load ---
cat("Loading Nov30 object...\n")
m <- readRDS(file.path(MOONSHOT, "past_analysis/Fig_organoid_singlecell_20251130/merged_subtype_annotated.rds"))
m$orig.ident <- factor(m$orig.ident, levels = sample_order)
m$Subtype_3levels <- factor(m$Subtype_3levels, levels = c("Classical-like", "Intermediate", "Basal-like"))
m$Disease_Type <- factor(m$Disease_Type, levels = c("HG-PanIN", "Invasive PDAC"))
cat(sprintf("cells: %d\n", ncol(m)))

# --- composition (verification + panel D) ---
comp_sample <- m@meta.data %>% group_by(orig.ident, Disease_Type, Subtype_3levels) %>%
  summarise(n = n(), .groups = "drop") %>% group_by(orig.ident) %>%
  mutate(pct = 100 * n / sum(n)) %>% ungroup()
comp_dis <- m@meta.data %>% group_by(Disease_Type, Subtype_3levels) %>%
  summarise(n = n(), .groups = "drop") %>% group_by(Disease_Type) %>%
  mutate(pct = 100 * n / sum(n)) %>% ungroup()
cat("\n=== per-sample composition (legend照合) ===\n")
print(comp_sample %>% transmute(orig.ident, Subtype_3levels, pct = round(pct, 1)) %>% as.data.frame())

# ==============================================================================
# Panel C-left: UMAP by sample  (53.5 x 45.8 mm)
# ==============================================================================
umap_df <- data.frame(
  UMAP1 = Embeddings(m, "umap")[, 1], UMAP2 = Embeddings(m, "umap")[, 2],
  Sample = m$orig.ident, Subtype = m$Subtype_3levels,
  Classical = m$Classical_Score1, Basal = m$Basal_Score1, ScoreDiff = m$Score_Diff
)
set.seed(1); umap_df <- umap_df[sample(nrow(umap_df)), ]

p_c_sample <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Sample)) +
  geom_point(size = 0.28, alpha = 0.6, stroke = 0) +  # was 0.12/0.5 — too faint after rasterization
  scale_color_manual(values = cols_sample, name = NULL) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_cd(base_size = 8) +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 6),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA)) +
  guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1)))
save_panel(p_c_sample, "Fig5C_UMAP_sample", 53.5, 45.8)

# ==============================================================================
# Panel C-right: UMAP by subtype  (51.5 x 43.5 mm)
# ==============================================================================
umap_sub <- umap_df %>% arrange(factor(Subtype, levels = c("Classical-like", "Intermediate", "Basal-like")))
p_c_sub <- ggplot(umap_sub, aes(UMAP1, UMAP2, color = Subtype)) +
  geom_point(size = 0.28, alpha = 0.6, stroke = 0) +  # was 0.12/0.5 — too faint after rasterization
  scale_color_manual(values = cols_subtype, name = NULL) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_cd(base_size = 8) +
  theme(legend.position = c(0.02, 0.02), legend.justification = c(0, 0),
        legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 6),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA)) +
  guides(color = guide_legend(override.aes = list(size = 1.6, alpha = 1)))
save_panel(p_c_sub, "Fig5C_UMAP_subtype", 51.5, 43.5)

# ==============================================================================
# Panel D-left: disease-level composition (35 x 56 mm; 公開22mmは狭すぎるため拡幅)
# ==============================================================================
lab_dis <- comp_dis %>% filter(pct >= 3) %>% group_by(Disease_Type) %>%
  arrange(desc(Subtype_3levels)) %>% mutate(cum = cumsum(pct), pos = cum - pct / 2)
p_d_dis <- ggplot(comp_dis, aes(Disease_Type, pct, fill = Subtype_3levels)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(data = lab_dis, aes(y = pos, label = sprintf("%.0f%%", pct)),
            color = "black", size = 2.4, fontface = "bold") +
  scale_fill_manual(values = cols_subtype, name = NULL) +
  scale_x_discrete(labels = c("HG-PanIN" = "HG-PanIN", "Invasive PDAC" = "Inv. PDAC")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_cd(base_size = 8) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6.5))
save_panel(p_d_dis, "Fig5D_composition_disease", 42, 56)

# ==============================================================================
# Panel D-right: per-sample composition (62.4 x 38.3 mm)
# ==============================================================================
lab_smp <- comp_sample %>% filter(pct >= 4) %>% group_by(orig.ident) %>%
  arrange(desc(Subtype_3levels)) %>% mutate(cum = cumsum(pct), pos = cum - pct / 2)
p_d_smp <- ggplot(comp_sample, aes(orig.ident, pct, fill = Subtype_3levels)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.25) +
  geom_text(data = lab_smp, aes(y = pos, label = sprintf("%.0f%%", pct)),
            color = "black", size = 1.9, fontface = "bold") +
  facet_grid(~ Disease_Type, scales = "free_x", space = "free_x",
             labeller = labeller(Disease_Type = c("HG-PanIN" = "HG-PanIN",
                                                  "Invasive PDAC" = "Inv. PDAC"))) +
  scale_fill_manual(values = cols_subtype, name = NULL) +
  scale_x_discrete(labels = c("KKp048" = "KKp048", "KKp055" = "KKp055",
                              "KKp011" = "KKp011", "KKp060" = "KKp060")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_cd(base_size = 8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5),
        legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 6),
        strip.text = element_text(size = 6.5, face = "bold"))
save_panel(p_d_smp, "Fig5D_composition_sample", 62.4, 38.3)

# ==============================================================================
# Panel E: ridge of Score_Diff by sample (63.5 x 46.1 mm)
# ==============================================================================
ridge_df <- data.frame(ScoreDiff = m$Score_Diff,
                       Sample = factor(m$orig.ident, levels = rev(sample_order)))
p_e <- ggplot(ridge_df, aes(ScoreDiff, Sample, fill = Sample)) +
  geom_density_ridges(alpha = 0.85, scale = 1.4, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray30", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray50", linewidth = 0.3) +
  scale_fill_manual(values = cols_sample) +
  scale_y_discrete(labels = c("KKp048" = "KKp048", "KKp055" = "KKp055",
                              "KKp011" = "KKp011", "KKp060" = "KKp060"), expand = c(0.01, 0.4)) +
  scale_x_continuous(limits = c(-1.2, 1.4), breaks = seq(-1, 1, 0.5)) +
  labs(x = "Subtype Score (Classical − Basal)", y = NULL) +
  theme_cd(base_size = 8) +
  theme(legend.position = "none", axis.text.y = element_text(face = "bold"))
save_panel(p_e, "Fig5E_ridge", 63.5, 46.1)

# ==============================================================================
# Panel F: Classical vs Basal scatter, facet by sample (123.4 x 33.7 mm)
# ==============================================================================
umap_df$Sample_lab <- factor(umap_df$Sample, levels = sample_order,
  labels = c("KKp048 (HG-PanIN)", "KKp055 (HG-PanIN)", "KKp011 (Inv. PDAC)", "KKp060 (Inv. PDAC)"))
p_f <- ggplot(umap_df, aes(Classical, Basal, color = Subtype)) +
  # size/alpha bumped (was 0.1/0.3): the original tiny, faint points rendered
  # fine in the vector PDF but washed out in the rasterized PNG, where sub-pixel
  # points get averaged toward white on display down-scaling.
  geom_point(size = 0.25, alpha = 0.5, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.25) +
  scale_color_manual(values = cols_subtype, name = NULL) +
  facet_wrap(~ Sample_lab, ncol = 4) +
  labs(x = "Classical Score", y = "Basal Score") +
  theme_cd(base_size = 8) +
  theme(legend.position = "right", strip.text = element_text(size = 6.5),
        legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 6.5)) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1)))
save_panel(p_f, "Fig5F_score_scatter", 123.4, 33.7)

# ==============================================================================
# Panels G / H: CAF-activating ligand heatmaps
# ==============================================================================
caf_ligands <- c("BMP4","BMP7","GDF15","SHH","IHH","WNT2B","WNT7B","PDGFB","CXCL1","CXCL3",
                 "CXCL5","CXCL8","CCL2","CCL20","AREG","EREG","HBEGF","VEGFA","MMP7","MMP14",
                 "LOXL2","LOXL4","TIMP1","TIMP2","IL1B","THBS2","SPP1")
caf_cat <- c(BMP4="TGFb_family",BMP7="TGFb_family",GDF15="TGFb_family",SHH="Hedgehog",IHH="Hedgehog",
             WNT2B="Wnt",WNT7B="Wnt",PDGFB="PDGF_family",CXCL1="Chemokines",CXCL3="Chemokines",
             CXCL5="Chemokines",CXCL8="Chemokines",CCL2="Chemokines",CCL20="Chemokines",
             AREG="Growth_factors",EREG="Growth_factors",HBEGF="Growth_factors",VEGFA="Growth_factors",
             MMP7="ECM_modifiers",MMP14="ECM_modifiers",LOXL2="ECM_modifiers",LOXL4="ECM_modifiers",
             TIMP1="ECM_modifiers",TIMP2="ECM_modifiers",IL1B="Inflammatory",
             THBS2="Other_CAF_activators",SPP1="Other_CAF_activators")
cat_colors <- c(TGFb_family="#7FBC41", Hedgehog="#F0E442", Wnt="#CC79A7", PDGF_family="#FFC0CB",
                Chemokines="#56B4E9", Growth_factors="#00CED1", ECM_modifiers="#E69F00",
                Inflammatory="#FFB6C1", Other_CAF_activators="#999999")
hm_cols <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

m <- tryCatch(JoinLayers(m), error = function(e) m)
m$Stage_Subtype <- paste0(ifelse(m$Disease_Type == "HG-PanIN", "HG-", "Inv-"),
                          ifelse(m$Subtype_3levels == "Basal-like", "Basal",
                                 ifelse(m$Subtype_3levels == "Classical-like", "Classical", "Intermediate")))

# ---- Panel G: HG-PanIN by 3 subtypes (87.1 x 97.1 mm) ----
hg <- subset(m, Disease_Type == "HG-PanIN")
Idents(hg) <- factor(hg$Subtype_3levels, levels = c("Classical-like", "Intermediate", "Basal-like"))
avg_g <- as.matrix(AverageExpression(hg, features = caf_ligands, assays = "RNA", layer = "data")$RNA)
colnames(avg_g) <- gsub("-like", "", colnames(avg_g))
avg_g <- avg_g[, c("Classical", "Intermediate", "Basal")]
zg <- t(scale(t(avg_g)))
direction <- ifelse(avg_g[, "Basal"] > avg_g[, "Classical"], "Basal-up", "Classical-up")
row_ann <- data.frame(Direction = direction, Category = caf_cat[rownames(zg)], row.names = rownames(zg))
col_ann_g <- data.frame(Subtype = colnames(zg), row.names = colnames(zg))
ann_cols_g <- list(Subtype = c(Classical = "#4393C3", Intermediate = "#FEE090", Basal = "#D6604D"),
                   Direction = c("Basal-up" = "#D6604D", "Classical-up" = "#4393C3"),
                   Category = cat_colors)
ph_g <- pheatmap(zg, cluster_rows = TRUE, cluster_cols = FALSE,
                 annotation_row = row_ann, annotation_col = col_ann_g, annotation_colors = ann_cols_g,
                 show_colnames = FALSE, color = hm_cols, border_color = NA,
                 fontsize = 6, fontsize_row = 7, treeheight_row = 12, silent = TRUE,
                 main = "")
save_ph(ph_g, "Fig5G_ligand_heatmap_HGPanIN", 87.1, 97.1)

# ---- Panel H: disease progression (94 x 99.3 mm) ----
keep_h <- m$Stage_Subtype %in% c("HG-Classical","HG-Basal","Inv-Classical","Inv-Basal")
mh <- subset(m, cells = colnames(m)[keep_h])
Idents(mh) <- mh$Stage_Subtype
avg_h <- as.matrix(AverageExpression(mh, features = caf_ligands, assays = "RNA", layer = "data")$RNA)
ord_h <- c("HG-Classical","HG-Basal","Inv-Classical","Inv-Basal")
avg_h <- avg_h[, ord_h]
zh <- t(scale(t(avg_h)))
col_ann_h <- data.frame(Stage = c("HG-PanIN","HG-PanIN","InvPDAC","InvPDAC"),
                        Subtype = c("Classical","Basal","Classical","Basal"), row.names = ord_h)
ann_cols_h <- list(Stage = c("HG-PanIN" = "#009E73", "InvPDAC" = "#D55E00"),
                   Subtype = c(Classical = "#4393C3", Basal = "#D6604D"),
                   Direction = c("Basal-up" = "#D6604D", "Classical-up" = "#4393C3"),
                   Category = cat_colors)
row_ann_h <- data.frame(Direction = direction[rownames(zh)], Category = caf_cat[rownames(zh)], row.names = rownames(zh))
ph_h <- pheatmap(zh, cluster_rows = TRUE, cluster_cols = FALSE, gaps_col = 2,
                 annotation_row = row_ann_h, annotation_col = col_ann_h, annotation_colors = ann_cols_h,
                 show_colnames = FALSE, color = hm_cols, border_color = NA,
                 fontsize = 6, fontsize_row = 7, treeheight_row = 12, silent = TRUE,
                 main = "")
save_ph(ph_h, "Fig5H_ligand_heatmap_progression", 94, 99.3)

cat("\nFig5 C-H 完了\n")
