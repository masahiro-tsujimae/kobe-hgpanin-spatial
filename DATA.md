# Data inputs and availability

All processed inputs needed to reproduce the **self-contained** panels are
included under this bundle at the relative paths the scripts expect. Raw
sequencing data are in GEO; very large processed objects and patient clinical
records are **not** distributed here.

## Bundled processed data (small)

| Path | Used by | Description |
|---|---|---|
| `Visium_analysis/output/verification/annotation_summary_per_sample.csv` | Fig 1F | per-sample cell-type composition |
| `Visium_analysis/output/verification/pathology_annotation_summary.csv` | Fig 1G | per-sample stromal proportion (non-Other denominator) |
| `Visium_analysis/output/distance/caf_epithelial_distances.csv` | Fig 2D | per-spot CAF-to-epithelium distances (spot-level; sample_id = specimen/slide code) |
| `Visium_analysis/output/LR_pathology/LR_scores_all.csv` | Fig 4A–F | standard-Visium ligand–receptor scores |
| `Visium_analysis/output/LR_pathology/sample_level_LR_scores.csv` | Fig 4A–F | sample-level L-R scores |
| `Visium_HD/data/lr_scores_HD.csv` | Fig 4G, 4H | Visium HD L-R scores (sample × pair × subtype) |
| `Visium_HD/data/lr_scores_4combination_HD.csv` | Fig 4J | Visium HD L-R, sender × receiver |
| `Visium_HD/data/fig3_subtype_counts.csv` | Fig 3D | Visium HD epithelial subtype composition |
| `RNA_seq/result_RNAseq_organoids3_PSC/Expression_Profile.GRCh38.gene.xlsx` | Fig 6B–D | bulk RNA-seq read counts (PSC + KKp048 CM) |
| `RNA_seq/result_RNAseq_organoids4_PSC/Expression_profile/StringTie/Expression_Profile.GRCh38.gene.xlsx` | Fig 6F–H | bulk RNA-seq read counts (PSC + 4 organoid CM) |
| `RNA_seq/result_RNAseq_2d_PSC/Expression_Profile.GRCh38.2dcells_PSC.gene.xlsx` | Sup Fig 10C–E | bulk RNA-seq read counts (PSC + 2D cell-line CM) |
| `IHC_analysis/fig2_aSMA.csv`, `fig2_IL6.csv`, `fig3_TP63.csv`, `fig3_KRT5.csv`, `fig3_GATA6.csv` | Fig 2F/H, Fig 3G | **de-identified** per-duct IHC scores (see below) |

Fig 1C / 1D (Masson's trichrome quantification) use values embedded directly in
`fig1_ihc_panels.R`; no external file is required.

## De-identification of IHC tables

The IHC panels were originally computed from a clinical pathology database that
contains patient identifiers (slide numbers, patient IDs). Those source files
are **not** included. Instead, only the two columns each figure uses — the duct
category (`Normal`/`LG`/`HG`) and the IHC score — were extracted into the
`IHC_analysis/*.csv` files above. No identifiers, dates, or other clinical
fields are present. (Extraction is performed by `_assemble.py`.)

## NOT included — large processed objects

Three panels require integrated Seurat objects that are too large to distribute
(~0.9–1.25 GB each). The scripts are included for inspection but will error
without these objects:

| Script | Required object | Source |
|---|---|---|
| `manuscript_figure_panels/Fig2/scripts/fig2_caf_panels.R` | `visium_with_subtypes.rds` (~0.9 GB) | derived from standard Visium (GSE331451) |
| `manuscript_figure_panels/Fig3/scripts/fig3_E_density.R` | `visium_hd_annotated.rds` (~0.9 GB) | derived from Visium HD (GSE331450) |
| `manuscript_figure_panels/Fig5/scripts/fig5_panels.R` | `merged_subtype_annotated.rds` (~1.25 GB) | derived from organoid scRNA-seq (GSE331447) |

## Raw data (GEO)

Raw spatial and single-cell sequencing data are deposited in NCBI GEO as
SuperSeries **GSE332553**:

- **GSE331451** — Visium FFPE spatial transcriptomics (12 specimens)
- **GSE331450** — Visium HD spatial transcriptomics (4 specimens)
- **GSE331447** — 10x Chromium single-cell RNA-seq (4 organoid lines)

The processed Seurat objects above can be regenerated from these raw data
(alignment/QC/integration as described in the manuscript Methods) or are
available from the corresponding author on reasonable request.
