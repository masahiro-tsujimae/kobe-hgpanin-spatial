# ==============================================================================
# Figure 6  B / C / D  再出力  (organoids3: Control vs KKp48 = KKp048 Early)
#   B = PCA            (配置実寸 61.3 x 49.0 mm)
#   C = myCAF / iCAF score boxplot  (各 30 x 35 mm)
#   D = iCAF/myCAF marker heatmap + log2FC bar  (90.1 x 103.0 mm)
# 解析は元スクリプト PSC_organoids3_DEG_pathway.R と同一。出力寸法/フォントのみ変更。
# ==============================================================================
source(file.path("./manuscript_figure_panels/Fig6/scripts", "_common.R"))

cols_group <- cols_group_organoids3   # Control=gray, KKp48_Early=blue
input_file <- file.path(RNASEQ_DIR, "result_RNAseq_organoids3_PSC",
                        "Expression_Profile.GRCh38.gene.xlsx")

cat("\n================ Fig6 B/C/D (organoids3) ================\n")

# --- Step 1: データ読み込み (KKp60 除外) ---
data <- read_excel(input_file, sheet = 1)
count_cols <- grep("Read_Count$", colnames(data), value = TRUE)
count_cols <- count_cols[!grepl("KKp60", count_cols)]
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
sample_info <- data.frame(
  sample = colnames(count_matrix_filtered),
  condition = case_when(
    grepl("Control", colnames(count_matrix_filtered)) ~ "Control",
    grepl("KKp48",   colnames(count_matrix_filtered)) ~ "KKp48_Early"
  ),
  row.names = colnames(count_matrix_filtered)
)
sample_info$condition <- factor(sample_info$condition, levels = c("Control", "KKp48_Early"))
dds <- DESeqDataSetFromMatrix(count_matrix_filtered, sample_info, design = ~ condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# --- DEG (log2FC バー用) ---
res <- results(dds, contrast = c("condition", "KKp48_Early", "Control"))
res <- res[!is.na(res$padj), ]

# ==============================================================================
# Panel B: PCA  (61.3 x 49.0 mm)
# ==============================================================================
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
cat(sprintf("PCA variance: PC1=%d%%  PC2=%d%%  (legend: PC1 94%%)\n",
            percentVar[1], percentVar[2]))

disp_lab <- c("Control" = "Control", "KKp48_Early" = "KKp048")
p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 2.4, stroke = 0.5) +
  scale_color_manual(values = cols_group, name = NULL, labels = disp_lab) +
  scale_shape_manual(values = c(16, 17), name = NULL, labels = disp_lab) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_nature(base_size = 8) +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.key.size = unit(3, "mm"))
save_panel(p_pca, "Fig6B_PCA", 61.3, 49.0)

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

# p値検証 (legend: myCAF P=0.0028, iCAF P=0.0014)
verify_pvals(score_df, "myCAF_score", "condition",
             list(c("Control","KKp48_Early")), "Fig6C myCAF")
verify_pvals(score_df, "iCAF_score", "condition",
             list(c("Control","KKp48_Early")), "Fig6C iCAF")

# ==============================================================================
# Panel C: score boxplots (各 30 x 35 mm)  myCAF=左, iCAF=右
# ==============================================================================
box_panel <- function(yvar, ylab_txt) {
  ggplot(score_df, aes(x = condition, y = .data[[yvar]], fill = condition)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, linewidth = 0.4) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.9) +
    scale_fill_manual(values = cols_group) +
    stat_compare_means(comparisons = list(c("Control","KKp48_Early")),
                       method = "t.test", size = 2.3, tip.length = 0.01) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    scale_x_discrete(labels = c("Control" = "Ctrl", "KKp48_Early" = "KKp048")) +
    labs(x = NULL, y = ylab_txt) +
    theme_nature(base_size = 8) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}
save_panel(box_panel("myCAF_score", "myCAF score"), "Fig6C_myCAF_box", 30, 35)
save_panel(box_panel("iCAF_score",  "iCAF score"),  "Fig6C_iCAF_box",  30, 35)

# ==============================================================================
# Panel D: iCAF/myCAF marker heatmap + log2FC bar (90.1 x 103.0 mm)
# ==============================================================================
all_caf <- c(iCAF_present, myCAF_present)
hm <- t(scale(t(vsd_mat[all_caf, ])))
caf_subtype <- ifelse(all_caf %in% iCAF_genes, "iCAF", "myCAF")
log2fc <- setNames(rep(NA_real_, length(all_caf)), all_caf)
for (g in all_caf) if (g %in% rownames(res)) log2fc[g] <- res[g, "log2FoldChange"]

cond_disp <- factor(ifelse(sample_info$condition == "Control", "Control", "KKp048"),
                    levels = c("Control", "KKp048"))
cols_disp <- c("Control" = unname(cols_group["Control"]),
               "KKp048"  = unname(cols_group["KKp48_Early"]))

ha_row <- rowAnnotation(
  CAF = caf_subtype,
  log2FC = anno_barplot(log2fc, gp = gpar(fill = "#4DBBD5"), width = unit(9, "mm")),
  col = list(CAF = c("iCAF" = "#E64B35", "myCAF" = "#4393C3")),
  annotation_name_gp = gpar(fontsize = 6),
  simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(CAF = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"), labels_gp = gpar(fontsize = 6),
    grid_height = unit(3, "mm"), grid_width = unit(3, "mm")))
)
ha_top <- HeatmapAnnotation(
  Group = cond_disp,
  col = list(Group = cols_disp),
  annotation_name_gp = gpar(fontsize = 6),
  simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(Group = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"), labels_gp = gpar(fontsize = 6),
    grid_height = unit(3, "mm"), grid_width = unit(3, "mm")))
)
ht <- Heatmap(
  hm, name = "Z-score", top_annotation = ha_top, right_annotation = ha_row,
  col = colorRamp2(c(-2, 0, 2), c("#4393C3", "white", "#D6604D")),
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 7, fontface = "italic"),
  show_column_names = FALSE, row_dend_width = unit(6, "mm"),
  cluster_columns = FALSE, row_split = caf_subtype,
  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  column_title = NULL,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 7, fontface = "bold"),
                              labels_gp = gpar(fontsize = 6),
                              grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm"),
                              legend_height = unit(14, "mm"))
)
save_heatmap(ht, "Fig6D_heatmap", 90.1, 103.0,
             heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
             merge_legends = TRUE)

cat("\nFig6 B/C/D 完了\n")
