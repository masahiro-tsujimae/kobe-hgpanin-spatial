# Figure reproduction code

Analysis and figure-generation code for:

**Spatial epithelial and stromal remodeling in isolated human pancreatic intraepithelial neoplasia before invasion**

Tsujimae M, et al. Manuscript under peer review at *Nature Communications*
(NCOMMS-26-078461). A Zenodo DOI for the archived release will be added upon publication.

This bundle reproduces the **data panels** of Figures 1–6 and Supplementary
Fig. 10 from processed data. Sequencing data are in NCBI GEO SuperSeries
**GSE332553** (SubSeries GSE331451 / GSE331450 / GSE331447).

---

## Quick start

Requires **R 4.4.2** and the packages in [`environment.txt`](environment.txt).
Run with the working directory set to this bundle root:

```r
setwd("/path/to/code")
source("reproduce_figures.R")
```

or from a shell:

```sh
cd /path/to/code
Rscript reproduce_figures.R
```

Each panel is written as PDF (vector) + PNG (600 dpi) to
`manuscript_figure_panels/<Fig>/panels/`. `reproduce_figures.R` runs all
self-contained panels and reports PASS/FAIL; three panels that need large
objects not distributed here are skipped (see [`DATA.md`](DATA.md)).

Paths inside the scripts are relative to this bundle root (the project root
`MOONSHOT` is set to `.`), so no path editing is needed as long as you run from
here.

---

## Layout

```
code/
├─ README.md                     this file
├─ DATA.md                       data inputs, de-identification, GEO accessions
├─ environment.txt               R version + package versions
├─ reproduce_figures.R           top-level runner
├─ _assemble.py                  how this bundle was assembled (provenance)
├─ common/
│   └─ theme_cancer_discovery.R  shared ggplot theme + colour palettes
├─ manuscript_figure_panels/
│   ├─ Fig1..Fig6/scripts/*.R    panel scripts (write to ../panels/)
│   └─ Fig1..Fig6/panels/        reference outputs (PDF + PNG)
├─ Visium_analysis/output/...    processed CSVs (Visium)
├─ Visium_HD/data/...            processed CSVs (Visium HD)
├─ RNA_seq/...                   bulk RNA-seq expression tables (xlsx)
└─ IHC_analysis/*.csv            de-identified per-duct IHC scores
```

## Script → figure map

| Script | Panels | Inputs |
|---|---|---|
| `Fig1/scripts/fig1_ihc_panels.R` | 1C, 1D | inline (Masson quantification) |
| `Fig1/scripts/fig1F_composition.R` | 1F | annotation_summary_per_sample.csv |
| `Fig1/scripts/fig1G_stromal.R` | 1G | pathology_annotation_summary.csv |
| `Fig2/scripts/fig2_ihc_distance_panels.R` | 2D, 2F, 2H | caf_epithelial_distances.csv; IHC CSVs |
| `Fig2/scripts/fig2_caf_panels.R` | 2A, 2C | **needs large object** (visium_with_subtypes.rds) |
| `Fig3/scripts/fig3_ihc_subtype_panels.R` | 3D, 3G | fig3_subtype_counts.csv; IHC CSVs |
| `Fig3/scripts/fig3_E_density.R` | 3E | **needs large object** (visium_hd_annotated.rds) |
| `Fig4/scripts/fig4_AF_panels.R` | 4A–4F | LR_pathology CSVs |
| `Fig4/scripts/fig4_GHJ_panels.R` | 4G, 4H, 4J | lr_scores_HD.csv, lr_scores_4combination_HD.csv |
| `Fig5/scripts/fig5_panels.R` | 5C–5H | **needs large object** (merged_subtype_annotated.rds) |
| `Fig6/scripts/fig6_org3_panels.R` | 6B–6D | organoids3 expression xlsx |
| `Fig6/scripts/fig6_org4_panels.R` | 6F–6H | organoids4 expression xlsx |
| `Fig6/scripts/supfig10_2d_panels.R` | Supplementary Fig. 10c–e | 2D cell-line expression xlsx |

`Fig4/scripts/_common.R` and `Fig6/scripts/_common.R` hold shared helpers and
are sourced by the panel scripts. `Fig4` and `Fig6` also include a `run_all.R`.

## Statistics note (Fig. 2D)

Fig. 2D distances are summarised at the **biological-replicate level**:
per-sample medians compared by paired Wilcoxon test (n = 5 per stage; not
significant, P = 0.31). Spot-level violins are shown for visualisation only. An
earlier spot-pooled test (P = 2.2 × 10⁻¹⁷) was pseudoreplicated and is no longer
reported.

## What is NOT in this bundle

- **Photomicrographs, schematics, spatial coordinate maps** (Fig 1A/B, 2B/E/G,
  3A/B/C/F, 5A/B, 6A/E/I and the model schematic) — images/illustrations, not
  computed panels.
- **Fig 4I** (CellChat network) — requires CellChat; not bundled.
- **Three large processed objects** for Fig 2A/C, Fig 3E, Fig 5 — see
  [`DATA.md`](DATA.md).
- **Patient clinical database** — only de-identified IHC score columns are
  included (see [`DATA.md`](DATA.md)).

## Provenance

This bundle was assembled by `_assemble.py` from the project's
`manuscript_figure_panels/` directory. The only transformation applied to the
scripts is rewriting the absolute project root to `.` (and, for the two IHC
scripts, reading the de-identified CSVs instead of the clinical database).
Analysis logic and reported statistics are unchanged.

## License

MIT License (see [`LICENSE`](LICENSE)).

## Contact

Corresponding author: Atsuhiro Masuda — atmasuda@med.kobe-u.ac.jp
Study IRB: Kobe University B200341.
