# ==============================================================================
# Supplemental Fig 10  C / D / E  再出力  (2D cell line CM: Control / Panc1 / BxPC3)
#   C = PCA                         (配置実寸 62.0 x 49.6 mm)
#   D = iCAF score boxplot          (66.5 x 58.2 mm)
#   E = iCAF core marker heatmap    (114.3 x 114.3 mm)
# 解析は元スクリプト PSC_2d_cellline_CM_DEG_pathway.R と同一。出力寸法/フォントのみ変更。
# ==============================================================================
source(file.path("./manuscript_figure_panels/Fig6/scripts", "_common.R"))

cols_group <- cols_group_2d   # Control=gray, Panc1=red, BxPC3=blue
input_file <- file.path(RNASEQ_DIR, "result_RNAseq_2d_PSC",
                        "Expression_Profile.GRCh38.2dcells_PSC.gene.xlsx")

cat("\n================ SupFig10 C/D/E (2d cell line CM) ================\n")

# --- Step 1: データ読み込み (AxPC 除外) ---
data <- read_excel(input_file, sheet = 1)
count_cols <- grep("Read_Count$", colnames(data), value = TRUE)
count_cols <- count_cols[!grepl("AxPC", count_cols)]
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
lvls <- c("Control", "Panc1", "BxPC3")
sample_info <- data.frame(
  sample = colnames(count_matrix_filtered),
  condition = case_when(
    grepl("Control", colnames(count_matrix_filtered)) ~ "Control",
    grepl("BxPC3",   colnames(count_matrix_filtered)) ~ "BxPC3",
    grepl("Panc1",   colnames(count_matrix_filtered)) ~ "Panc1"
  ),
  row.names = colnames(count_matrix_filtered)
)
sample_info$condition <- factor(sample_info$condition, levels = lvls)
dds <- DESeqDataSetFromMatrix(count_matrix_filtered, sample_info, design = ~ condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# ==============================================================================
# Panel C: PCA  (62.0 x 49.6 mm)
# ==============================================================================
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
cat(sprintf("PCA variance: PC1=%d%%  PC2=%d%%  (legend: PC1 77%%, PC2 21%%)\n",
            percentVar[1], percentVar[2]))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_point(size = 2.2, stroke = 0.5) +
  scale_color_manual(values = cols_group, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_nature(base_size = 8) +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.key.size = unit(3, "mm"))
save_panel(p_pca, "SupFig10C_PCA", 62.0, 49.6)

# ==============================================================================
# Step 4: マーカー / スコア
# ==============================================================================
iCAF_genes <- c("IL6","IL1B","CXCL1","CXCL2","CXCL8","CXCL12","CCL2",
                "CCL5","PTGS2","NFKB1","STAT3","CEBPB","LIF","SAA1",
                "ICAM1","VCAM1","SERPINE1","MMP1","MMP3","TIMP1")
iCAF_core  <- c("CFD","DPT","PDGFRA","CXCL12","LIF","CCL2","IL1A",
                "CXCL2","HAS1","CXCL1","CXCL8","IL6","CXCL3")
iCAF_present      <- iCAF_genes[iCAF_genes %in% rownames(vsd_mat)]
iCAF_core_present <- iCAF_core[iCAF_core   %in% rownames(vsd_mat)]

score_df <- data.frame(
  sample    = colnames(vsd_mat),
  condition = sample_info$condition,
  iCAF_score = colMeans(vsd_mat[iCAF_present, , drop = FALSE])
)
cmp_list <- list(c("Control","Panc1"), c("Control","BxPC3"), c("Panc1","BxPC3"))
# legend: Panc1 vs Ctrl 0.0012, BxPC3 vs Ctrl 0.00011, BxPC3 vs Panc1 7e-5
verify_pvals(score_df, "iCAF_score", "condition", cmp_list, "SupFig10D iCAF")

# ==============================================================================
# Panel D: iCAF score boxplot (66.5 x 58.2 mm)
# ==============================================================================
p_box <- ggplot(score_df, aes(x = condition, y = iCAF_score, fill = condition)) +
  geom_violin(alpha = 0.3, linewidth = 0.3) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.08, size = 1.3, alpha = 0.9) +
  scale_fill_manual(values = cols_group) +
  stat_compare_means(comparisons = cmp_list, method = "t.test",
                     size = 2.3, tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = "iCAF score") +
  theme_nature(base_size = 8) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
save_panel(p_box, "SupFig10D_iCAF_box", 66.5, 58.2)

# ==============================================================================
# Panel E: iCAF core marker heatmap (114.3 x 114.3 mm)
# ==============================================================================
col_order <- order(factor(sample_info$condition, levels = lvls))
hm <- t(scale(t(vsd_mat[iCAF_core_present, col_order])))
cond_disp <- sample_info$condition[col_order]

ha_top <- HeatmapAnnotation(
  Group = cond_disp, col = list(Group = cols_group),
  annotation_name_gp = gpar(fontsize = 7), simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(Group = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7),
    grid_height = unit(3.5, "mm"), grid_width = unit(3.5, "mm")))
)
ht <- Heatmap(
  hm, name = "Z-score", top_annotation = ha_top,
  col = colorRamp2(c(-2, 0, 2), c("#4393C3", "white", "#D6604D")),
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  row_names_side = "left", show_column_names = FALSE,
  cluster_columns = FALSE, cluster_rows = TRUE, row_dend_width = unit(8, "mm"),
  column_title = NULL,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "bold"),
                              labels_gp = gpar(fontsize = 7),
                              grid_height = unit(3, "mm"), grid_width = unit(3, "mm"))
)
save_heatmap(ht, "SupFig10E_heatmap", 114.3, 114.3)

cat("\nSupFig10 C/D/E 完了\n")
