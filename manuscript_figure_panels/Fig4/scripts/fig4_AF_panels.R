# ==============================================================================
# Figure 4 A-F 再出力 (standard Visium L-R)  — fig4_AF.R と同一ロジック、配置実寸で出力
#   A: Epi->Str heatmap (Z-score)      56.0 x 56.0 mm   [ComplexHeatmap]
#   B: Str->Epi heatmap (Z-score)      46.6 x 51.7 mm   [ComplexHeatmap]
#   C: Top-20 changed L-R bar          51.7 x 52.9 mm
#   D: 9 key-pair boxplots (LG vs HG)  57.4 x 49.8 mm
#   E: mean L-R heatmap (9 pairs)      48.5 x 42.9 mm
#   F: log2FC bar (Basal-like, HG/LG)  48.5 x 51.4 mm
# 入力: Visium_analysis/output/LR_pathology/{sample_level_LR_scores,LR_scores_all}.csv
# ==============================================================================
source(file.path("./manuscript_figure_panels/Fig4/scripts", "_common.R"))

lr_dir <- file.path(MOONSHOT, "Visium_analysis/output/LR_pathology")
lr_scores_all <- read.csv(file.path(lr_dir, "LR_scores_all.csv"))
sample_lr <- read.csv(file.path(lr_dir, "sample_level_LR_scores.csv"))

epi_to_str <- sample_lr %>% filter(direction == "Epi_to_Str") %>%
  mutate(diagnosis = factor(diagnosis, levels = c("LG-PanIN", "HG-PanIN")))
str_to_epi <- sample_lr %>% filter(direction == "Str_to_Epi") %>%
  mutate(diagnosis = factor(diagnosis, levels = c("LG-PanIN", "HG-PanIN")))

diag2 <- colors_diagnosis[c("LG-PanIN", "HG-PanIN")]
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))

# 凡例ブロック周りの余白を詰める（既定2mm）。密ヒートマップ(A/B)で本体↔凡例の
# 隙間を減らし、配置幅(56/46.6mm)に確実に収める。
ht_opt$HEATMAP_LEGEND_PADDING    <- unit(0.5, "mm")
ht_opt$ANNOTATION_LEGEND_PADDING <- unit(0.5, "mm")

cat("\n================ Fig4 A-F ================\n")

# ---- helper: build per-sample Z-score heatmap (A / B) ----
make_lr_heatmap <- function(df, title, name, w, h) {
  # 各 L-R ペアを 1 行に集約（diagnosis を id 列から外す）。旧版は diagnosis を
  # 残していたため各ペアが LG行/HG行 の 2 行に重複展開され（片側が NA->0 埋め）、
  # 19 ペアが 38(->33) 行に膨らんでラベルが約4.5ptまで縮小していた。全 10 サンプル
  # 横断で Z-score 化すれば 1 ペア 1 行・人工ゼロなしで、行数が半減し可読になる。
  wide <- df %>% select(LR_pair, pathway, sample_id, score) %>%
    pivot_wider(names_from = sample_id, values_from = score)
  pw_lookup <- wide %>% select(LR_pair, pathway) %>% distinct()
  mat <- wide %>% select(-LR_pair, -pathway) %>% as.matrix()
  rownames(mat) <- wide$LR_pair; mat[is.na(mat)] <- 0
  mat <- mat[apply(mat, 1, var) > 1e-10, , drop = FALSE]
  pw_info <- pw_lookup$pathway[match(rownames(mat), pw_lookup$LR_pair)]
  sdiag <- df %>% select(sample_id, diagnosis) %>% distinct()
  dvec <- setNames(as.character(sdiag$diagnosis), sdiag$sample_id)
  ord <- names(sort(factor(dvec[colnames(mat)], levels = c("LG-PanIN", "HG-PanIN"))))
  mat <- mat[, ord, drop = FALSE]
  matz <- t(scale(t(mat))); matz[is.na(matz)] <- 0
  rownames(matz) <- lr_label(rownames(matz))   # _ -> en-dash (ligand–receptor)
  # 本体サイズと行ラベルを配置寸法(w,h mm)に自動適応。A(19行/56mm)/B(25行/46.6mm)で
  # 行数も枠も違うため固定値だと一方が見切れる。高さは枠内に収め、行ラベルは行高に
  # 合わせて可読サイズへ（A≈6.5pt / B≈4.4pt）。
  nr <- nrow(matz); nc <- ncol(matz)
  body_h <- min(nr * 2.4, h - 10)
  body_w <- min(nc * 1.4, w - 30)                       # 1列凡例なので本体を広く取れる
  rfont  <- max(4.0, min(6.5, (body_h / nr) / 0.38))
  # 注釈は色のみ表示し、凡例は下で手動生成（自動だと Pathway と Z-score が右に2列に
  # 並んで横幅を食い、狭い B では端が見切れるため）。Diagnosis 凡例は列見出し LG/HG と
  # 重複するので出さない。
  ha_top <- HeatmapAnnotation(Diagnosis = dvec[colnames(mat)],
    col = list(Diagnosis = diag2), simple_anno_size = unit(2, "mm"),
    show_annotation_name = FALSE, show_legend = FALSE)
  ha_left <- rowAnnotation(Pathway = pw_info, col = list(Pathway = pathway_palette(pw_info)),
    simple_anno_size = unit(2, "mm"), show_annotation_name = FALSE, show_legend = FALSE)
  ht <- Heatmap(matz, name = "Z-score", col = col_fun, top_annotation = ha_top, left_annotation = ha_left,
    cluster_columns = FALSE, show_column_names = FALSE, show_row_dend = FALSE, show_heatmap_legend = FALSE,
    width = unit(body_w, "mm"), height = unit(body_h, "mm"),
    column_split = factor(dvec[colnames(mat)], levels = c("LG-PanIN", "HG-PanIN")),
    column_title = c("LG", "HG"), row_names_gp = gpar(fontsize = rfont),
    column_gap = unit(1, "mm"), column_title_gp = gpar(fontsize = 7, fontface = "bold"))
  # 凡例を手動で「縦1列」に固定（max_height を大きく取り packLegend が折り返さない）。
  pw_lv  <- sort(unique(pw_info))
  pw_col <- pathway_palette(pw_info)[pw_lv]
  lgd_pw <- Legend(labels = pw_lv, title = "Pathway", legend_gp = gpar(fill = pw_col),
                   title_gp = gpar(fontsize = 5.5, fontface = "bold"), labels_gp = gpar(fontsize = 5),
                   grid_height = unit(1.4, "mm"), grid_width = unit(1.8, "mm"))
  lgd_z  <- Legend(col_fun = col_fun, title = "Z-score", at = c(-2, -1, 0, 1, 2),
                   title_gp = gpar(fontsize = 5.5, fontface = "bold"), labels_gp = gpar(fontsize = 4.5),
                   grid_width = unit(1.4, "mm"), legend_height = unit(9, "mm"),
                   title_position = "leftcenter-rot")
  pd <- packLegend(lgd_pw, lgd_z, direction = "vertical", gap = unit(2, "mm"),
                   max_height = unit(h, "mm"))
  save_heatmap(ht, name, w, h, annotation_legend_list = pd,
               annotation_legend_side = "right", padding = unit(c(1, 1, 1, 1), "mm"))
}
make_lr_heatmap(epi_to_str, "Epithelial -> Stromal", "Fig4A_heatmap_epi_to_str", 56.0, 56.0)
make_lr_heatmap(str_to_epi, "Stromal -> Epithelial", "Fig4B_heatmap_str_to_epi", 46.6, 51.7)

# ---- C: Top-20 changed L-R pairs (HG vs LG) ----
change_data <- lr_scores_all %>% filter(diagnosis %in% c("LG-PanIN", "HG-PanIN")) %>%
  group_by(LR_pair) %>%
  summarise(LG = mean(score[diagnosis == "LG-PanIN"], na.rm = TRUE),
            HG = mean(score[diagnosis == "HG-PanIN"], na.rm = TRUE), .groups = "drop") %>%
  mutate(change = HG - LG, dir = ifelse(change > 0, "Increased", "Decreased")) %>%
  distinct(LR_pair, .keep_all = TRUE) %>% arrange(desc(abs(change))) %>% head(20) %>% arrange(change)
change_data$LR_pair <- factor(change_data$LR_pair, levels = change_data$LR_pair)
fig4c <- ggplot(change_data, aes(LR_pair, change, fill = dir)) +
  geom_col(width = 0.72) + geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = colors_direction, name = NULL) +
  scale_x_discrete(labels = lr_label) + coord_flip() +
  labs(title = "Top changed L-R pairs", subtitle = "HG-PanIN vs LG-PanIN", x = NULL, y = "Score change (HG - LG)") +
  theme_cd(base_size = 7) +
  theme(legend.position = c(0.78, 0.18), legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6), plot.subtitle = element_text(size = 6.5),
        axis.text.y = element_text(size = 6))
save_panel(fig4c, "Fig4C_top_changed_LR", 51.7, 52.9)

# ---- D: 9 key-pair boxplots (LG vs HG, Wilcoxon) ----
key_pairs <- c("FGF2_FGFR1","PDGFA_PDGFRA","PDGFB_PDGFRB","SHH_PTCH1","TGFB1_TGFBR1",
               "TGFB2_TGFBR2","VEGFA_KDR","WNT7A_FZD1","WNT7B_FZD7")
plot_d <- epi_to_str %>% filter(LR_pair %in% key_pairs) %>%
  mutate(LR_pair = factor(LR_pair, levels = key_pairs))
fig4d <- ggplot(plot_d, aes(diagnosis, score, fill = diagnosis)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.5, alpha = 0.7) +
  facet_wrap(~ LR_pair, scales = "free_y", ncol = 3, labeller = as_labeller(lr_label)) +
  scale_fill_manual(values = diag2, guide = "none") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 1.8, label.y.npc = 0.97) +
  scale_x_discrete(labels = c("LG-PanIN" = "LG", "HG-PanIN" = "HG")) +
  # p値が箱に重ならないよう上に余白を確保（free_y の各パネルに適用）
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = "L-R score") +
  theme_cd(base_size = 7) +
  theme(axis.text.x = element_text(size = 6), strip.text = element_text(size = 5, face = "bold"),
        strip.clip = "off",  # ストリップ見出し(遺伝子名)が枠で切れないように
        axis.text.y = element_text(size = 5.5), panel.spacing = unit(0.6, "mm"),
        plot.margin = margin(4, 7, 4, 4))  # 右端ファセットの見出しが切れない余白
save_panel(fig4d, "Fig4D_boxplots_epi_to_str", 57.4, 49.8)

# ---- E: mean L-R heatmap (9 pairs, LG vs HG) ----
hm_e <- epi_to_str %>% filter(LR_pair %in% key_pairs) %>%
  group_by(LR_pair, diagnosis) %>% summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
lr_ord <- hm_e %>% filter(diagnosis == "HG-PanIN") %>% arrange(desc(mean_score)) %>% pull(LR_pair)
hm_e <- hm_e %>% mutate(LR_pair = factor(LR_pair, levels = rev(lr_ord)))
fig4e <- ggplot(hm_e, aes(diagnosis, LR_pair, fill = mean_score)) +
  geom_tile(color = "white", linewidth = 0.7) +
  scale_fill_gradientn(colors = c("#2166AC", "#F7F7F7", "#B2182B"),
                       values = rescale(c(0, 0.12, 0.25)), limits = c(0, 0.3), name = "L-R\nscore", oob = squish) +
  scale_x_discrete(labels = c("LG-PanIN" = "LG", "HG-PanIN" = "HG")) +
  scale_y_discrete(labels = lr_label) +
  labs(title = "Epi-to-Str activity", x = NULL, y = NULL) +
  theme_cd(base_size = 8) +
  theme(axis.text.y = element_text(size = 7), legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6.5))
save_panel(fig4e, "Fig4E_heatmap_progression", 48.5, 42.9)

# ---- F: log2FC bar (9 pairs, HG vs LG) ----
stat_f <- epi_to_str %>% filter(LR_pair %in% key_pairs) %>% group_by(LR_pair) %>%
  summarise(LG = mean(score[diagnosis == "LG-PanIN"], na.rm = TRUE),
            HG = mean(score[diagnosis == "HG-PanIN"], na.rm = TRUE),
            log2FC = log2((HG + 0.001) / (LG + 0.001)),
            p = tryCatch(wilcox.test(score[diagnosis == "HG-PanIN"], score[diagnosis == "LG-PanIN"])$p.value,
                         error = function(e) NA), .groups = "drop") %>%
  arrange(desc(log2FC)) %>%
  mutate(LR_pair = factor(LR_pair, levels = LR_pair),
         sig = case_when(p < 0.01 ~ "**", p < 0.05 ~ "*", p < 0.1 ~ "+", TRUE ~ ""),
         dir = ifelse(log2FC > 0, "Up", "Down"))
fig4f <- ggplot(stat_f, aes(LR_pair, log2FC, fill = dir)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sig, y = log2FC + 0.12 * sign(log2FC)), size = 3, fontface = "bold") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = colors_direction, guide = "none") +
  scale_x_discrete(labels = lr_label) + coord_flip() +
  labs(title = "Epi-to-Str (Basal-like)", subtitle = "log2FC (HG vs LG)", x = NULL, y = "log2 fold change",
       caption = "Wilcoxon: **P<.01 *P<.05 +P<.1") +
  theme_cd(base_size = 7) +
  theme(axis.text.y = element_text(size = 6.5), plot.subtitle = element_text(size = 6.5),
        plot.caption = element_text(size = 5, color = "gray40", hjust = 0))
save_panel(fig4f, "Fig4F_log2FC_HGvsLG", 48.5, 51.4)

cat("\n=== Fig4D p-value 検証 (legend: TGFB1 0.032, WNT7B 0.0079) ===\n")
print(stat_f %>% transmute(LR_pair, log2FC = round(log2FC, 2), p = signif(p, 3)) %>% as.data.frame())
cat("\nFig4 A-F 完了\n")
