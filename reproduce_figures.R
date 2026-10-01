# =============================================================================
# Reproduce the data panels of Figures 1-6 and Supplementary Fig. 10.
#
# USAGE: run with the working directory set to this bundle root (the folder
# that contains this file). For example:
#     setwd("/path/to/code"); source("reproduce_figures.R")
#   or from a shell:
#     cd /path/to/code && Rscript reproduce_figures.R
#
# Panels are written to manuscript_figure_panels/<Fig>/panels/ as PDF + PNG.
# See environment.txt for required R packages and DATA.md for data inputs.
# =============================================================================

if (!file.exists("common/theme_cancer_discovery.R")) {
  stop("Run this from the bundle root (the folder containing this file). ",
       "Current working directory: ", getwd())
}

# Panels that require large processed Seurat objects NOT included in this
# bundle (see DATA.md). They are listed for completeness but skipped here.
needs_large_object <- c(
  "manuscript_figure_panels/Fig2/scripts/fig2_caf_panels.R",   # visium_with_subtypes.rds  (~0.9 GB)
  "manuscript_figure_panels/Fig3/scripts/fig3_E_density.R",    # visium_hd_annotated.rds   (~0.9 GB)
  "manuscript_figure_panels/Fig5/scripts/fig5_panels.R"        # merged_subtype_annotated.rds (~1.25 GB)
)

# Self-contained panels (run from the bundled processed data + de-identified CSVs).
runnable <- c(
  "manuscript_figure_panels/Fig1/scripts/fig1_ihc_panels.R",        # Fig 1C, 1D
  "manuscript_figure_panels/Fig1/scripts/fig1F_composition.R",      # Fig 1F
  "manuscript_figure_panels/Fig1/scripts/fig1G_stromal.R",          # Fig 1G
  "manuscript_figure_panels/Fig2/scripts/fig2_ihc_distance_panels.R", # Fig 2D, 2F, 2H
  "manuscript_figure_panels/Fig3/scripts/fig3_ihc_subtype_panels.R",  # Fig 3D, 3G
  "manuscript_figure_panels/Fig4/scripts/fig4_AF_panels.R",         # Fig 4A-F
  "manuscript_figure_panels/Fig4/scripts/fig4_GHJ_panels.R",        # Fig 4G, 4H, 4J
  "manuscript_figure_panels/Fig6/scripts/fig6_org3_panels.R",       # Fig 6B-D
  "manuscript_figure_panels/Fig6/scripts/fig6_org4_panels.R",       # Fig 6F-H
  "manuscript_figure_panels/Fig6/scripts/supfig10_2d_panels.R"      # Sup Fig 10C-E
)

results <- character(0)
for (s in runnable) {
  cat("\n>>>>>>>>>>>>>>>>>>>>", s, "\n")
  ok <- tryCatch({ source(s, local = new.env()); TRUE },
                 error = function(e) { cat("  !! ERROR:", conditionMessage(e), "\n"); FALSE })
  results[s] <- if (ok) "PASS" else "FAIL"
}

cat("\n================ SUMMARY ================\n")
for (s in names(results)) cat(sprintf("  [%s] %s\n", results[s], s))
cat("\nSkipped (require large objects not in this bundle; see DATA.md):\n")
for (s in needs_large_object) cat("  -", s, "\n")
