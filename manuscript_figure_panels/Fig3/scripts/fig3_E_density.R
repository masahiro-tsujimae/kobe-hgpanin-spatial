# ==============================================================================
# Figure 3 E 再出力 (fig3_DE_epithelial.R fig3e ロジック)
#   E: Classical-to-Basal spectrum density (Subtype Index, LG vs HG), zone-colored
# 入力: Visium_HD/data/visium_hd_annotated.rds (1.8GB)。Ductal spot の Subtype_index を使用。
#   検証: KS p (pub 1.07e-2), Wilcoxon p (pub 3.73e-3)。
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })
MOONSHOT <- "."
source(file.path(MOONSHOT, "common", "theme_cancer_discovery.R"))   # theme_cd, colors_epi
OUT <- file.path(MOONSHOT, "manuscript_figure_panels", "Fig3", "panels"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
save_panel <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  cat(sprintf("  [panel] %-26s %.1f x %.1f mm\n", name, w, h))
}
cat("\n================ Fig3 E (HD density) ================\n")
cat("Loading HD object (1.8GB)...\n")
v <- readRDS(file.path(MOONSHOT, "Visium_HD/data/visium_hd_annotated.rds"))
epi <- v@meta.data %>% filter(compartment == "Ductal", !is.na(Subtype_index),
                              diagnosis %in% c("LG-PanIN", "HG-PanIN"))
rm(v); gc(verbose = FALSE)
epi$diagnosis <- factor(epi$diagnosis, levels = c("LG-PanIN", "HG-PanIN"))
cat(sprintf("Ductal spots: LG=%d HG=%d\n", sum(epi$diagnosis=="LG-PanIN"), sum(epi$diagnosis=="HG-PanIN")))

lg <- epi$Subtype_index[epi$diagnosis == "LG-PanIN"]
hg <- epi$Subtype_index[epi$diagnosis == "HG-PanIN"]
ks <- ks.test(lg, hg); wx <- suppressWarnings(wilcox.test(lg, hg))
cat(sprintf("KS p=%.2e (pub 1.07e-2) ; Wilcoxon p=%.2e (pub 3.73e-3)\n", ks$p.value, wx$p.value))

# zone-colored density per diagnosis
rng <- range(epi$Subtype_index)
dens <- do.call(rbind, lapply(levels(epi$diagnosis), function(dg) {
  d <- density(epi$Subtype_index[epi$diagnosis == dg], n = 1024, from = rng[1]-0.1, to = rng[2]+0.1)
  data.frame(x = d$x, y = d$y, diagnosis = dg)
}))
dens$diagnosis <- factor(dens$diagnosis, levels = c("LG-PanIN", "HG-PanIN"))
dens$zone <- factor(ifelse(dens$x < -0.2, "Basal-like", ifelse(dens$x > 0.2, "Classical", "Intermediate")),
                    levels = c("Basal-like", "Intermediate", "Classical"))
lab <- data.frame(x = c(-0.7, 0.7), y = Inf, label = c("Basal-like", "Classical"),
                  col = c(colors_epi[["Basal-like"]], colors_epi[["Classical"]]),
                  diagnosis = factor("LG-PanIN", levels = c("LG-PanIN","HG-PanIN")))
fig3e <- ggplot(dens, aes(x, y)) +
  geom_segment(aes(xend = x, yend = 0, color = zone), linewidth = 0.5, alpha = 0.85) +
  geom_line(linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "grey70", linewidth = 0.25) +
  scale_color_manual(values = c("Basal-like" = colors_epi[["Basal-like"]],
                                "Intermediate" = colors_epi[["Intermediate"]],
                                "Classical" = colors_epi[["Classical"]]), guide = "none") +
  facet_wrap(~ diagnosis, ncol = 1, scales = "free_y") +
  geom_text(data = lab, aes(x = x, y = y, label = label), color = lab$col, size = 2.2,
            fontface = "italic", vjust = 1.6, inherit.aes = FALSE) +
  labs(x = "Subtype Index (Classical − Basal)", y = "Density",
       caption = sprintf("KS p = %.2e, Wilcoxon p = %.2e", ks$p.value, wx$p.value)) +
  theme_cd(base_size = 8) +
  theme(strip.text = element_text(size = 7.5, face = "bold"),
        plot.caption = element_text(size = 5.5, color = "gray30", hjust = 0))
save_panel(fig3e, "Fig3E_subtype_index_density", 52, 62)
cat("\nFig3 E 完了\n")
