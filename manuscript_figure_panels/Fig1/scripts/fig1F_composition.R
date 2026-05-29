# ==============================================================================
# Fig1 F (cell type composition) + Fig2 A (CAF composition) 再出力
#   入力: Visium_analysis/output/verification/annotation_summary_per_sample.csv
#         （病理 cell-type 由来の per-sample 集計。total-spot 分母）
#   構成バーは total-spot 分母なので本CSVで再現可（boxplot G/C は非-Other 分母のため別途）。
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))   # theme_cd, colors_celltype, colors_caf
save_panel <- function(p, name, dir, w, h) {
  out <- file.path(MOONSHOT, "manuscript_figure_panels", dir, "panels"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  cat(sprintf("  [panel] %-26s %.1f x %.1f mm\n", name, w, h))
}
d <- read.csv(file.path(MOONSHOT, "Visium_analysis/output/verification/annotation_summary_per_sample.csv"))
d$diagnosis_category <- factor(d$diagnosis_category, levels = c("Normal", "LG-PanIN", "HG-PanIN"))

cat("\n================ Fig1 F + Fig2 A (composition, CSV) ================\n")

# ---- Fig1 F: cell type composition (Acinar/Ductal/Stromal/Other, total denom) ----
fF <- d %>% mutate(Other_pct = pmax(0, 100 - Ductal_pct - Acinar_pct - Stromal_pct)) %>%
  arrange(diagnosis_category, sample_id)
fF$sample_id <- factor(fF$sample_id, levels = unique(fF$sample_id))
fF_long <- fF %>% select(sample_id, diagnosis_category, Acinar_pct, Ductal_pct, Stromal_pct, Other_pct) %>%
  pivot_longer(-c(sample_id, diagnosis_category), names_to = "CellType", values_to = "pct") %>%
  mutate(CellType = factor(gsub("_pct", "", CellType), levels = c("Acinar", "Ductal", "Stromal", "Other")))
cols_ct <- c("Acinar" = unname(colors_celltype["Acinar"]), "Ductal" = unname(colors_celltype["Ductal"]),
             "Stromal" = unname(colors_celltype["Stromal"]), "Other" = "#BDBDBD")
fig1f <- ggplot(fF_long, aes(sample_id, pct, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.85, color = "black", linewidth = 0.25) +
  facet_grid(~ diagnosis_category, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cols_ct, name = "Cell type") +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Percentage (%)") +
  theme_cd(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4.5),
        legend.position = "right", legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6.5),
        strip.text = element_text(size = 6.5, face = "bold"), panel.spacing = unit(0.6, "mm"))
save_panel(fig1f, "Fig1F_celltype_composition", "Fig1", 110, 38)

# 注: Fig2 A (CAF構成) は本CSV(total分母)では公開 image3 と case順序・値が一致しないため
#     ここでは生成しない。Fig2 A は非-Other 分母 (= per-spot 病理データ) が必要でブロック中。
cat("\nFig1 F 完了\n")
