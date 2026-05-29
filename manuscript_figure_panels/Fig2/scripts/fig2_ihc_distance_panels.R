# ==============================================================================
# Figure 2 D / F / H 再出力 (IHC_analysis.R + fig2_distance.R と同一データ・ロジック)
#   D: CAF-to-epithelium distance violin (iCAF vs myCAF, by diagnosis)  [CSV, cell-level]
#   F: alpha-SMA expression score stacked bar  [Fig1 Excel]
#   H: IL-6 expression score stacked bar       [Fig1 Excel]
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(ggsignif); library(readxl) })
MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))   # theme_cd, colors_*, colors_score, colors_caf
OUT <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig2", "panels"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  cat(sprintf("  [panel] %-26s %.1f x %.1f mm\n", name, w, h))
}
ihc_dir <- file.path(MOONSHOT, "IHC_analysis")  # de-identified per-duct IHC scores (columns: duct, score)
dunn_pairs <- function(x, g) {  # base-R Dunn-equivalent via pairwise wilcox + bonferroni
  lv <- levels(g); res <- list()
  for (cmb in combn(lv, 2, simplify = FALSE)) {
    p <- tryCatch(wilcox.test(x[g == cmb[1]], x[g == cmb[2]])$p.value, error = function(e) NA)
    res[[length(res)+1]] <- data.frame(g1 = cmb[1], g2 = cmb[2], p = p)
  }
  d <- do.call(rbind, res); d$padj <- p.adjust(d$p, method = "bonferroni"); d
}
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

cat("\n================ Fig2 D/F/H ================\n")

# ---- F / H: IHC score stacked bars ----
mk_stacked <- function(csv, name, w, h) {
  df <- read.csv(file.path(ihc_dir, csv)); names(df) <- c("duct", "score")
  df <- df %>% filter(!is.na(duct), !is.na(score)) %>%
    filter(duct %in% c("Normal", "LG", "HG")) %>%
    mutate(duct = recode(duct, "LG" = "LG-PanIN", "HG" = "HG-PanIN"),
           duct = factor(duct, levels = c("Normal", "LG-PanIN", "HG-PanIN")))
  cnt <- df %>% group_by(duct, score) %>% summarise(n = n(), .groups = "drop") %>%
    group_by(duct) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
  dn <- dunn_pairs(df$score, df$duct)
  comps <- lapply(seq_len(nrow(dn)), function(i) c(dn$g1[i], dn$g2[i]))
  labs_ <- star(dn$padj)
  ggplot(cnt, aes(duct, pct, fill = factor(score))) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.4) +
    # NOTE: labels must be a *named* vector keyed by score level. aSMA has no
    # score 0, so positional labels shifted everything down by one (score 3
    # "Strong" was mislabeled "2 Mod"). Named labels match by level regardless
    # of which scores are present.
    scale_fill_manual(values = colors_score, name = "Score",
                      labels = c("0" = "0 Neg", "1" = "1 Weak", "2" = "2 Mod", "3" = "3 Strong")) +
    scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 132), clip = "off") +
    geom_signif(comparisons = comps, annotations = labs_, y_position = c(104, 114, 124),
                tip_length = 0.01, textsize = 2.2, vjust = -0.3, size = 0.3) +
    labs(x = NULL, y = "Percentage (%)") +
    theme_cd(base_size = 8) +
    theme(legend.position = "right", legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6), legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 20, hjust = 1))
}
save_panel(mk_stacked("fig2_aSMA.csv", "Fig2F_aSMA_stacked", 55, 52), "Fig2F_aSMA_stacked", 55, 52)
save_panel(mk_stacked("fig2_IL6.csv", "Fig2H_IL6_stacked", 55, 52), "Fig2H_IL6_stacked", 55, 52)

# ---- D: distance violin ----
# NOTE: the spot-level Wilcoxon previously annotated here (e.g. P=2.2e-17 in HG)
# is pseudoreplicated — it pools thousands of spots from only n=5 biological
# samples, so the large spot count manufactures a tiny P from a negligible
# effect (Cliff's delta < 0.13). Statistics are now done at the biological-
# replicate level: per-sample medians (the proper unit, n=5/stage) compared by
# paired Wilcoxon. The spot violins/boxes remain for VISUALIZATION ONLY.
dist <- read.csv(file.path(MOONSHOT, "Visium_analysis/output/distance/caf_epithelial_distances.csv"))
dist$diagnosis <- factor(dist$diagnosis, levels = c("LG-PanIN", "HG-PanIN"))
dist$CAF_subtype <- factor(dist$CAF_subtype, levels = c("iCAF", "myCAF"))

per_samp <- dist %>% group_by(diagnosis, sample_id, CAF_subtype) %>%
  summarise(med = median(distance_um), .groups = "drop")
pw <- per_samp %>%
  pivot_wider(names_from = CAF_subtype, values_from = med) %>%
  filter(!is.na(iCAF), !is.na(myCAF)) %>%
  group_by(diagnosis) %>%
  summarise(n = dplyr::n(),
            p = tryCatch(wilcox.test(iCAF, myCAF, paired = TRUE)$p.value, error = function(e) NA),
            .groups = "drop") %>%
  mutate(label = sprintf("paired Wilcoxon\n(n=%d) P=%.2f", n, p),
         y = max(dist$distance_um) * 0.98)

fig2d <- ggplot(dist, aes(CAF_subtype, distance_um, fill = CAF_subtype)) +
  geom_violin(alpha = 0.45, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3, fill = "white", alpha = 0.8) +
  # per-sample medians (biological replicates), paired by sample
  geom_line(data = per_samp, aes(x = CAF_subtype, y = med, group = sample_id),
            inherit.aes = FALSE, color = "grey45", linewidth = 0.25, alpha = 0.7) +
  geom_point(data = per_samp, aes(x = CAF_subtype, y = med),
             inherit.aes = FALSE, shape = 21, fill = "white", color = "black",
             size = 1.1, stroke = 0.3) +
  facet_grid(~ diagnosis) +
  scale_fill_manual(values = c("iCAF" = colors_caf[["iCAF"]], "myCAF" = colors_caf[["myCAF"]]), guide = "none") +
  geom_text(data = pw, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 1.8, vjust = 1, lineheight = 0.9) +
  labs(x = NULL, y = expression("Distance to epithelium (" * mu * "m)")) +
  theme_cd(base_size = 8) +
  theme(strip.text = element_text(size = 7, face = "bold"))
save_panel(fig2d, "Fig2D_distance_violin", 60, 50)

cat("\nFig2 D/F/H 完了\n")
