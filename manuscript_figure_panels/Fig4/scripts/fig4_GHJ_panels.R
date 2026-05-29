# ==============================================================================
# Figure 4 G / H / J 再出力 (Visium HD L-R)  — 配置実寸で出力
#   G: Epi->Str L-R heatmap, subtype x stage (row z-score)  48.6 x 57.4 mm
#   H: Basal-like Epi->Str log2FC bar (HG vs LG)            54.4 x 57.4 mm
#   J: Comprehensive heatmap, sender x receiver x stage     58.4 x 112.1 mm  [ComplexHeatmap]
#   (I = CellChat network: 別扱い。CellChat 未インストール)
# 入力: Visium_HD/data/lr_scores_HD.csv, lr_scores_4combination_HD.csv
# ==============================================================================
source(file.path("./manuscript_figure_panels/Fig4/scripts", "_common.R"))

hd_dir <- file.path(MOONSHOT, "Visium_HD/data")
lr_hd <- read.csv(file.path(hd_dir, "lr_scores_HD.csv"))
lr_4c <- read.csv(file.path(hd_dir, "lr_scores_4combination_HD.csv"))

cat("\n================ Fig4 G/H/J ================\n")

# ---- G: subtype x stage heatmap (row z-score) ----
g_df <- lr_hd %>% filter(epi_subtype %in% c("Classical", "Basal-like")) %>%
  group_by(LR_pair, pathway, diagnosis, epi_subtype) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(diag_short = ifelse(diagnosis == "HG-PanIN", "HG", "LG"),
         col = paste0(diag_short, "-", epi_subtype))
col_levels_g <- c("LG-Classical", "LG-Basal-like", "HG-Classical", "HG-Basal-like")
g_wide <- g_df %>% select(LR_pair, col, mean_score) %>%
  pivot_wider(names_from = col, values_from = mean_score)
g_mat <- as.matrix(g_wide[, col_levels_g]); rownames(g_mat) <- g_wide$LR_pair
g_mat[is.na(g_mat)] <- 0
g_z <- t(scale(t(g_mat))); g_z[is.na(g_z)] <- 0
row_ord_g <- rownames(g_z)[hclust(dist(g_z))$order]
g_long <- as.data.frame(g_z) %>% tibble::rownames_to_column("LR_pair") %>%
  pivot_longer(-LR_pair, names_to = "col", values_to = "z") %>%
  mutate(LR_pair = factor(LR_pair, levels = row_ord_g),
         col = factor(col, levels = col_levels_g))
fig4g <- ggplot(g_long, aes(col, LR_pair, fill = z)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(colors = c("#2166AC", "#F7F7F7", "#B2182B"), limits = c(-1.5, 1.5),
                       name = "Row\nz-score", oob = squish) +
  scale_x_discrete(labels = c("LG-Classical" = "LG\nCl", "LG-Basal-like" = "LG\nBa",
                              "HG-Classical" = "HG\nCl", "HG-Basal-like" = "HG\nBa")) +
  scale_y_discrete(labels = lr_label) +
  labs(title = "Epi-to-Str (HD)", x = NULL, y = NULL) +
  theme_cd(base_size = 7) +
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 5.5),
        legend.title = element_text(size = 6))
save_panel(fig4g, "Fig4G_HD_subtype_stage_heatmap", 48.6, 57.4)

# ---- H: Basal-like log2FC bar (HG vs LG) ----
h_df <- lr_hd %>% filter(epi_subtype == "Basal-like", diagnosis %in% c("LG-PanIN", "HG-PanIN")) %>%
  group_by(LR_pair, pathway) %>%
  summarise(LG = mean(score[diagnosis == "LG-PanIN"], na.rm = TRUE),
            HG = mean(score[diagnosis == "HG-PanIN"], na.rm = TRUE),
            log2FC = log2((HG + 0.001) / (LG + 0.001)),
            p = tryCatch(wilcox.test(score[diagnosis == "HG-PanIN"], score[diagnosis == "LG-PanIN"])$p.value,
                         error = function(e) NA), .groups = "drop") %>%
  arrange(desc(log2FC)) %>%
  mutate(LR_pair = factor(LR_pair, levels = LR_pair),
         sig = case_when(is.na(p) ~ "", p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ "+", TRUE ~ ""),
         dir = ifelse(log2FC > 0, "Up", "Down"))
fig4h <- ggplot(h_df, aes(LR_pair, log2FC, fill = dir)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sig, y = log2FC + 0.18 * sign(log2FC)), size = 2.6, fontface = "bold") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = colors_direction, guide = "none") +
  scale_x_discrete(labels = lr_label) + coord_flip() +
  # 上端に少し余白（バーが軸線に張り付かないように）
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(title = "Basal-like Epi→Str (HD)", subtitle = "log2FC (HG vs LG)", x = NULL, y = "log2 fold change") +
  theme_cd(base_size = 7) +
  theme(axis.text.y = element_text(size = 7), plot.subtitle = element_text(size = 6.5),
        # タイトルをパネルではなく図全体の幅基準でセンタリング → 長いタイトルが
        # 右端で切れる問題を解消。サイズも 8->7.2pt に下げ右余白も確保。
        plot.title = element_text(size = 7.2, face = "bold", hjust = 0.5),
        plot.title.position = "plot",
        plot.margin = margin(4, 6, 4, 4))
save_panel(fig4h, "Fig4H_HD_basal_log2FC", 54.4, 57.4)

# ---- J: comprehensive heatmap (sender x receiver x stage), raw L-R score ----
j_df <- lr_4c %>%
  filter(sender_subtype %in% c("Classical", "Basal-like"),
         receiver_subtype %in% c("iCAF", "myCAF"),
         diagnosis %in% c("LG-PanIN", "HG-PanIN")) %>%
  group_by(LR_pair, pathway, diagnosis, sender_subtype, receiver_subtype) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(diag_short = ifelse(diagnosis == "HG-PanIN", "HG", "LG"),
         send_short = ifelse(sender_subtype == "Basal-like", "Basal", "Classical"),
         col = paste(diag_short, send_short, receiver_subtype, sep = "_"))
col_levels_j <- c("LG_Classical_iCAF", "LG_Classical_myCAF", "LG_Basal_iCAF", "LG_Basal_myCAF",
                  "HG_Classical_iCAF", "HG_Classical_myCAF", "HG_Basal_iCAF", "HG_Basal_myCAF")
col_levels_j <- col_levels_j[col_levels_j %in% unique(j_df$col)]
j_wide <- j_df %>% select(LR_pair, pathway, col, mean_score) %>%
  pivot_wider(names_from = col, values_from = mean_score, values_fn = mean)
j_pw <- j_wide %>% select(LR_pair, pathway) %>% distinct()
# order rows by pathway then LR_pair
j_wide <- j_wide %>% arrange(pathway, LR_pair)
j_mat <- as.matrix(j_wide[, col_levels_j]); rownames(j_mat) <- lr_label(j_wide$LR_pair); j_mat[is.na(j_mat)] <- 0
pw_info_j <- j_wide$pathway

parse_col <- function(x) do.call(rbind, strsplit(x, "_"))
cc <- parse_col(col_levels_j)
col_ann_j <- HeatmapAnnotation(
  Diagnosis = cc[, 1], Sender = cc[, 2], Receiver = cc[, 3],
  col = list(Diagnosis = c("LG" = colors_diagnosis[["LG-PanIN"]], "HG" = colors_diagnosis[["HG-PanIN"]]),
             Sender = c("Classical" = "#4575B4", "Basal" = "#D73027"),
             Receiver = c("iCAF" = "#CB181D", "myCAF" = "#2171B5")),
  simple_anno_size = unit(2.5, "mm"), annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(
    Diagnosis = list(title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 5.5), grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm")),
    Sender    = list(title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 5.5), grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm")),
    Receiver  = list(title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 5.5), grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm"))))
row_ann_j <- rowAnnotation(Pathway = pw_info_j, col = list(Pathway = pathway_palette(pw_info_j)),
  simple_anno_size = unit(2.5, "mm"), annotation_name_gp = gpar(fontsize = 6),
  annotation_legend_param = list(Pathway = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5), grid_height = unit(2.3, "mm"), grid_width = unit(2.3, "mm"))))
col_fun_j <- colorRamp2(c(0, 0.15, 0.3), c("#FFF7EC", "#FDBB84", "#B30000"))
ht_j <- Heatmap(j_mat, name = "L-R score", col = col_fun_j,
  top_annotation = col_ann_j, left_annotation = row_ann_j,
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
  column_split = factor(cc[, 1], levels = c("LG", "HG")), column_gap = unit(1.2, "mm"),
  column_title = c("LG", "HG"), column_title_gp = gpar(fontsize = 6, fontface = "bold"),
  row_names_gp = gpar(fontsize = 5.5),
  width = unit(20, "mm"), height = unit(92, "mm"),
  rect_gp = gpar(col = "grey90", lwd = 0.2),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5), grid_height = unit(2.2, "mm"), grid_width = unit(2.2, "mm")))
# Canvas widened (58.4 -> 92) so the right-hand legends (incl. "L-R score") are
# no longer clipped, and top padding added so the LG/HG column titles are not cut
# off at the upper edge. Heatmap body stays fixed at 20 x 92 mm.
save_heatmap(ht_j, "Fig4J_HD_comprehensive_heatmap", 92, 120,
             heatmap_legend_side = "right", annotation_legend_side = "right", merge_legends = TRUE,
             padding = unit(c(2, 2, 5, 2), "mm"))

cat("\nFig4 G/H/J 完了\n")
