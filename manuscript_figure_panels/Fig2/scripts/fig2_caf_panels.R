# ==============================================================================
# Figure 2 A / C 再出力 (CAF subtype, fig2_AB.R ロジック)
#   A: 症例別 iCAF/myCAF 構成 stacked bar (Case 1-5, LG/HG)
#   C: myCAF / iCAF 割合 boxplot (N/LG/HG, Jonckheere-Terpstra)
# 入力: subtype_LR_analysis/visium_with_subtypes.rds (CAF_subtype/compartment 計算済み)
#   ※ 元の fig2_AB.R は F:/visium/cell_type の病理 cell-type で Stromal を定義するが、
#     その CSV が無いため、本オブジェクト既存の compartment=="Stromal" + CAF_subtype を使用。
#     Jonckheere p を公開値(myCAF 0.025, iCAF 0.233)と照合して整合を確認する。
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(tidyr); library(clinfun) })
MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))   # theme_cd, colors_diagnosis, colors_caf
OUT <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig2", "panels"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  cat(sprintf("  [panel] %-26s %.1f x %.1f mm\n", name, w, h))
}
cat("\n================ Fig2 A/C (CAF) ================\n")
o <- readRDS(file.path(MOONSHOT, "past_analysis/20251206_Fig4_normalVisium_LR/subtype_LR_analysis/visium_with_subtypes.rds"))
m <- o@meta.data
cat("CAF_subtype levels: "); print(table(m$CAF_subtype, useNA = "ifany"))

stats <- m %>%
  group_by(sample_id, diagnosis_category) %>%
  summarise(total = n(),
            stromal_n = sum(compartment == "Stromal", na.rm = TRUE),
            myCAF_n = sum(CAF_subtype == "myCAF", na.rm = TRUE),
            iCAF_n = sum(CAF_subtype == "iCAF", na.rm = TRUE), .groups = "drop") %>%
  mutate(panCAF_pct = stromal_n / total * 100,
         myCAF_pct = myCAF_n / total * 100,
         iCAF_pct = iCAF_n / total * 100,
         diagnosis_category = factor(diagnosis_category, levels = c("Normal", "LG-PanIN", "HG-PanIN")))

# --- verification: Jonckheere (Normal<LG<HG) ---
jt <- function(v) jonckheere.test(v, as.numeric(stats$diagnosis_category), alternative = "increasing")$p.value
cat(sprintf("Jonckheere  myCAF P=%.4f (pub 0.025) ; iCAF P=%.4f (pub 0.233) ; panCAF P=%.4f\n",
            jt(stats$myCAF_pct), jt(stats$iCAF_pct), jt(stats$panCAF_pct)))

# ==============================================================================
# Panel A: per-case CAF composition stacked bar (LG/HG, Case 1-5)
# ==============================================================================
case_map <- stats %>% filter(diagnosis_category != "Normal") %>%
  arrange(diagnosis_category, desc(panCAF_pct)) %>% group_by(diagnosis_category) %>%
  mutate(case_label = paste0("Case ", row_number())) %>% ungroup()
barA <- case_map %>%
  select(case_label, diagnosis_category, myCAF_pct, iCAF_pct) %>%
  pivot_longer(c(myCAF_pct, iCAF_pct), names_to = "CAF", values_to = "pct") %>%
  mutate(CAF = factor(ifelse(CAF == "myCAF_pct", "myCAF", "iCAF"), levels = c("iCAF", "myCAF")),
         case_label = factor(case_label, levels = paste0("Case ", 1:5)))
fig2a <- ggplot(barA, aes(case_label, pct, fill = CAF)) +
  geom_bar(stat = "identity", width = 0.78, color = "white", linewidth = 0.3) +
  facet_grid(~ diagnosis_category, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colors_caf, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Proportion (%)") +
  theme_cd(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
        legend.position = "top", legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 6.5), legend.margin = margin(0,0,-4,0),
        strip.text = element_text(size = 6.5, face = "bold"), panel.spacing = unit(0.6, "mm"))
save_panel(fig2a, "Fig2A_CAF_composition", 56, 33)

# ==============================================================================
# Panel C: myCAF / iCAF proportion boxplots (N/LG/HG, Jonckheere)
# ==============================================================================
mk_caf_box <- function(yvar, ylab, jp) {
  ymax <- max(stats[[yvar]])
  ggplot(stats, aes(diagnosis_category, .data[[yvar]], fill = diagnosis_category)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.8) +
    scale_fill_manual(values = colors_diagnosis, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22)), limits = c(0, NA)) +
    annotate("text", x = 2, y = ymax * 1.16, label = sprintf("Jonckheere\nP=%.3f", jp), size = 2.1, color = "gray30") +
    labs(x = NULL, y = ylab) +
    theme_cd(base_size = 8) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}
save_panel(mk_caf_box("myCAF_pct", "myCAF (%)", jt(stats$myCAF_pct)), "Fig2C_myCAF_box", 40, 52)
save_panel(mk_caf_box("iCAF_pct",  "iCAF (%)",  jt(stats$iCAF_pct)),  "Fig2C_iCAF_box",  40, 52)

cat("\nFig2 A/C 完了\n")
