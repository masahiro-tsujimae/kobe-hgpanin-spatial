# ==============================================================================
# Figure 6  F / G / H  再出力  (organoids4: Control / KKp048 / KKp011 / KKp060)
#   F = PCA            (配置実寸 47.5 x 37.7 mm)
#   G = myCAF / iCAF score boxplot  (各 41 x 36 mm)
#   H = iCAF/myCAF marker heatmap + 3x log2FC bar  (102.8 x 91.4 mm)
# 解析は元スクリプト PSC_conditioned_medium_DEG_pathway.R と同一。出力寸法/フォントのみ変更。
# ==============================================================================
source(file.path("./manuscript_figure_panels/Fig6/scripts", "_common.R"))

cols_group <- cols_group_organoids4
input_file <- file.path(RNASEQ_DIR, "result_RNAseq_organoids4_PSC",
                        "Expression_profile", "StringTie",
                        "Expression_Profile.GRCh38.gene.xlsx")

cat("\n================ Fig6 F/G/H (organoids4) ================\n")

# --- Step 1: データ読み込み (KKp114 除外) ---
data <- read_excel(input_file, sheet = 1)
count_cols <- grep("Read_Count$", colnames(data), value = TRUE)
count_cols <- count_cols[!grepl("KKp114", count_cols)]
count_matrix <- as.matrix(data[, count_cols])
rownames(count_matrix) <- make.unique(as.character(data$Gene_Symbol))
colnames(count_matrix) <- gsub("_Read_Count", "", colnames(count_matrix))
storage.mode(count_matrix) <- "integer"
keep <- rowSums(count_matrix >= 10) >= 3
count_matrix_filtered <- count_matrix[keep, ]
colnames(count_matrix_filtered) <- gsub("-", ".", colnames(count_matrix_filtered))
cat(sprintf("samples: %s\n", paste(colnames(count_matrix_filtered), collapse = ", ")))
cat(sprintf("genes after filter: %d\n", nrow(count_matrix_filtered)))

# --- Step 2: DESeq2 ---
lvls <- c("Control", "KKp048_Early", "KKp011_Advanced", "KKp060_Advanced")
sample_info <- data.frame(
  sample = colnames(count_matrix_filtered),
  condition = case_when(
    grepl("Control", colnames(count_matrix_filtered)) ~ "Control",
    grepl("KKp048",  colnames(count_matrix_filtered)) ~ "KKp048_Early",
    grepl("KKp011",  colnames(count_matrix_filtered)) ~ "KKp011_Advanced",
    grepl("KKp060",  colnames(count_matrix_filtered)) ~ "KKp060_Advanced"
  ),
  row.names = colnames(count_matrix_filtered)
)
sample_info$condition <- factor(sample_info$condition, levels = lvls)
dds <- DESeqDataSetFromMatrix(count_matrix_filtered, sample_info, design = ~ condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# --- DEG (log2FC バー用、各群 vs Control) ---
res_048 <- results(dds, contrast = c("condition", "KKp048_Early",    "Control")); res_048 <- res_048[!is.na(res_048$padj), ]
res_011 <- results(dds, contrast = c("condition", "KKp011_Advanced", "Control")); res_011 <- res_011[!is.na(res_011$padj), ]
res_060 <- results(dds, contrast = c("condition", "KKp060_Advanced", "Control")); res_060 <- res_060[!is.na(res_060$padj), ]

# ==============================================================================
# Panel F: PCA  (47.5 x 37.7 mm)
# ==============================================================================
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
cat(sprintf("PCA variance: PC1=%d%%  PC2=%d%%  (legend: PC1 53%%, PC2 24%%)\n",
            percentVar[1], percentVar[2]))

disp_lab <- c("Control" = "Control", "KKp048_Early" = "KKp048",
              "KKp011_Advanced" = "KKp011", "KKp060_Advanced" = "KKp060")
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 2.0, stroke = 0.5) +
  scale_color_manual(values = cols_group, name = NULL, labels = disp_lab) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL, labels = disp_lab) +
  xlab(paste0("PC1: ", percentVar[1], "% var")) +
  ylab(paste0("PC2: ", percentVar[2], "% var")) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_nature(base_size = 8) +
  theme(legend.position = "bottom", legend.box.spacing = unit(1, "mm"),
        legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0))
save_panel(p_pca, "Fig6F_PCA", 47.5, 37.7)

# ==============================================================================
# Step 4: マーカー / スコア
# ==============================================================================
myCAF_genes <- c("ACTA2","CNN1","TAGLN","MYL9","TPM1","TPM2","ACTG2",
                 "MYLK","CALD1","LMOD1","CSRP1","CSRP2","FLNA","VIM","DES")
iCAF_genes  <- c("IL6","IL1B","CXCL1","CXCL2","CXCL8","CXCL12","CCL2",
                 "CCL5","PTGS2","NFKB1","STAT3","CEBPB","LIF","SAA1",
                 "ICAM1","VCAM1","SERPINE1","MMP1","MMP3","TIMP1")
myCAF_present <- myCAF_genes[myCAF_genes %in% rownames(vsd_mat)]
iCAF_present  <- iCAF_genes[iCAF_genes  %in% rownames(vsd_mat)]

score_df <- data.frame(
  sample    = colnames(vsd_mat),
  condition = sample_info$condition,
  iCAF_score  = colMeans(vsd_mat[iCAF_present, , drop = FALSE]),
  myCAF_score = colMeans(vsd_mat[myCAF_present, , drop = FALSE])
)
cmp_list <- list(c("Control","KKp048_Early"),
                 c("Control","KKp011_Advanced"),
                 c("Control","KKp060_Advanced"))
verify_pvals(score_df, "myCAF_score", "condition", cmp_list, "Fig6G myCAF")
verify_pvals(score_df, "iCAF_score",  "condition", cmp_list, "Fig6G iCAF")

# ==============================================================================
# Panel G: score boxplots (各 41 x 36 mm)  myCAF=左, iCAF=右
# ==============================================================================
xlabs <- c("Control" = "Ctrl", "KKp048_Early" = "KKp048",
           "KKp011_Advanced" = "KKp011", "KKp060_Advanced" = "KKp060")
box_panel <- function(yvar, ylab_txt) {
  ggplot(score_df, aes(x = condition, y = .data[[yvar]], fill = condition)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 0.9, alpha = 0.9) +
    scale_fill_manual(values = cols_group) +
    stat_compare_means(comparisons = cmp_list, method = "t.test",
                       size = 1.9, tip.length = 0.01) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    scale_x_discrete(labels = xlabs) +
    labs(x = NULL, y = ylab_txt) +
    theme_nature(base_size = 8) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5))
}
save_panel(box_panel("myCAF_score", "myCAF score"), "Fig6G_myCAF_box", 41, 36)
save_panel(box_panel("iCAF_score",  "iCAF score"),  "Fig6G_iCAF_box",  41, 36)

# ==============================================================================
# Panel H: iCAF/myCAF marker heatmap + 3x log2FC bar (102.8 x 91.4 mm)
# ==============================================================================
all_caf <- c(iCAF_present, myCAF_present)
col_order <- order(factor(sample_info$condition, levels = lvls))
hm <- t(scale(t(vsd_mat[all_caf, col_order])))
caf_subtype <- ifelse(all_caf %in% iCAF_genes, "iCAF", "myCAF")

getfc <- function(res) { v <- setNames(rep(NA_real_, length(all_caf)), all_caf)
  for (g in all_caf) if (g %in% rownames(res)) v[g] <- res[g, "log2FoldChange"]; v }
fc048 <- getfc(res_048); fc011 <- getfc(res_011); fc060 <- getfc(res_060)

cond_disp <- factor(disp_lab[as.character(sample_info$condition[col_order])],
                    levels = c("Control","KKp048","KKp011","KKp060"))
cols_disp <- setNames(unname(cols_group), c("Control","KKp048","KKp011","KKp060"))

ha_row <- rowAnnotation(
  CAF = caf_subtype,
  KKp048 = anno_barplot(fc048, gp = gpar(fill = "#3182BD"), width = unit(7, "mm")),
  KKp011 = anno_barplot(fc011, gp = gpar(fill = "#CB181D"), width = unit(7, "mm")),
  KKp060 = anno_barplot(fc060, gp = gpar(fill = "#F46D43"), width = unit(7, "mm")),
  col = list(CAF = c("iCAF" = "#E64B35", "myCAF" = "#4393C3")),
  annotation_name_gp = gpar(fontsize = 5.5), simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(CAF = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"), labels_gp = gpar(fontsize = 6),
    grid_height = unit(3, "mm"), grid_width = unit(3, "mm")))
)
ha_top <- HeatmapAnnotation(
  Group = cond_disp, col = list(Group = cols_disp),
  annotation_name_gp = gpar(fontsize = 6), simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(Group = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"), labels_gp = gpar(fontsize = 6),
    grid_height = unit(3, "mm"), grid_width = unit(3, "mm")))
)
ht <- Heatmap(
  hm, name = "Z-score", top_annotation = ha_top, right_annotation = ha_row,
  col = colorRamp2(c(-2, 0, 2), c("#4393C3", "white", "#D6604D")),
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 6.5, fontface = "italic"),
  show_column_names = FALSE, row_dend_width = unit(5, "mm"),
  cluster_columns = FALSE, row_split = caf_subtype,
  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  column_title = NULL,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7, fontface = "bold"),
                              labels_gp = gpar(fontsize = 6),
                              grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm"))
)
save_heatmap(ht, "Fig6H_heatmap", 102.8, 91.4,
             heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
             merge_legends = TRUE)

cat("\nFig6 F/G/H 完了\n")
