"""Assemble the figure-reproduction code release bundle (`code/`).

This script gathers the curated figure-panel scripts, the shared plotting
theme, the small processed-data inputs, and DE-IDENTIFIED extracts of the
IHC scoring tables into a self-contained, path-portable bundle.

Design
------
* The bundle reproduces the original repo-relative directory layout, so each
  copied R script works after a single change: the hardcoded absolute project
  root ``C:/Users/kofau/OneDrive/Desktop/Moonshot`` is rewritten to ``.``
  (current working directory). Run any script with the working directory set
  to this bundle root (see README).
* Large processed Seurat objects (~0.9-1.25 GB each) are NOT bundled; the raw
  data are in GEO (GSE332553) and the processed objects are available on
  request. Scripts that need them (Fig2 CAF, Fig3E, Fig5) are included for
  inspection but will error without the objects (documented in data/README).
* Patient clinical databases are NOT copied. Only the columns the figures use
  (duct category + IHC score) are extracted into de-identified CSVs; all
  identifiers (slide number, patient ID, focus point) are dropped.

Run:  python code/_assemble.py
"""
from __future__ import annotations

import csv
import shutil
from pathlib import Path

import openpyxl

ROOT = Path(r"C:/Users/kofau/OneDrive/Desktop/Moonshot")
DST = ROOT / "code"
ABS = "C:/Users/kofau/OneDrive/Desktop/Moonshot"  # absolute root literal to neutralize

# --- figure-panel scripts to ship (dev/diagnostic helpers excluded) ---
FIG_SCRIPTS = {
    "Fig1": ["fig1_ihc_panels.R", "fig1F_composition.R", "fig1G_stromal.R"],
    "Fig2": ["fig2_caf_panels.R", "fig2_ihc_distance_panels.R"],
    "Fig3": ["fig3_E_density.R", "fig3_ihc_subtype_panels.R"],
    "Fig4": ["_common.R", "fig4_AF_panels.R", "fig4_GHJ_panels.R", "run_all.R"],
    "Fig5": ["fig5_panels.R"],
    "Fig6": ["_common.R", "fig6_org3_panels.R", "fig6_org4_panels.R",
             "supfig10_2d_panels.R", "run_all.R"],
}

# --- small processed inputs, mirrored at the same relative path scripts expect ---
DATA_FILES = [
    "Visium_analysis/output/verification/annotation_summary_per_sample.csv",
    "Visium_analysis/output/verification/pathology_annotation_summary.csv",
    "Visium_analysis/output/distance/caf_epithelial_distances.csv",
    "Visium_analysis/output/LR_pathology/LR_scores_all.csv",
    "Visium_analysis/output/LR_pathology/sample_level_LR_scores.csv",
    "Visium_HD/data/lr_scores_HD.csv",
    "Visium_HD/data/lr_scores_4combination_HD.csv",
    "Visium_HD/data/fig3_subtype_counts.csv",
    "RNA_seq/result_RNAseq_organoids3_PSC/Expression_Profile.GRCh38.gene.xlsx",
    "RNA_seq/result_RNAseq_organoids4_PSC/Expression_profile/StringTie/Expression_Profile.GRCh38.gene.xlsx",
    "RNA_seq/result_RNAseq_2d_PSC/Expression_Profile.GRCh38.2dcells_PSC.gene.xlsx",
]

# --- de-identified IHC extracts: (source xlsx, sheet, duct col, score col, out csv) ---
IHC_EXTRACTS = [
    ("past_analysis/Fig1_20251015/Database_earlyPDAC_20241127_kobepanc.xlsx",
     "aSMA", "aSMA_duct", "aSMA", "fig2_aSMA.csv"),
    ("past_analysis/Fig1_20251015/Database_earlyPDAC_20241127_kobepanc.xlsx",
     "IL6", "IL_6_duct", "IL_6", "fig2_IL6.csv"),
    ("past_analysis/Fig3_20251017/Database_earlyPDAC_20241127_kobepanc_fig3.xlsx",
     "TP63", "p63_1_組織", "p63(/HPF)_1カ所目", "fig3_TP63.csv"),
    ("past_analysis/Fig3_20251017/Database_earlyPDAC_20241127_kobepanc_fig3.xlsx",
     "KRT5", "KRT5_1_組織", "KRT5(/HPF)_1カ所目", "fig3_KRT5.csv"),
    ("past_analysis/Fig3_20251017/Database_earlyPDAC_20241127_kobepanc_fig3.xlsx",
     "GATA6", "GATA6_duct", "GATA6", "fig3_GATA6.csv"),
]


def neutralize(text: str) -> str:
    """Rewrite the absolute project root to '.' so paths resolve from the bundle."""
    return text.replace(ABS, ".")


def rewire_ihc(name: str, text: str) -> str:
    """Swap patient-Excel reads for de-identified CSV reads."""
    if name == "fig2_ihc_distance_panels.R":
        text = text.replace(
            'fig1_excel <- file.path(MOONSHOT, "past_analysis/Fig1_20251015/Database_earlyPDAC_20241127_kobepanc.xlsx")',
            'ihc_dir <- file.path(MOONSHOT, "IHC_analysis")  # de-identified per-duct IHC scores (columns: duct, score)')
        text = text.replace(
            'mk_stacked <- function(sheet, duct_col, score_col, name, w, h) {\n  d <- read_excel(fig1_excel, sheet = sheet)\n  df <- d[, c(duct_col, score_col)]; names(df) <- c("duct", "score")',
            'mk_stacked <- function(csv, name, w, h) {\n  df <- read.csv(file.path(ihc_dir, csv)); names(df) <- c("duct", "score")')
        text = text.replace(
            'save_panel(mk_stacked("aSMA", "aSMA_duct", "aSMA", "Fig2F_aSMA_stacked", 55, 52), "Fig2F_aSMA_stacked", 55, 52)',
            'save_panel(mk_stacked("fig2_aSMA.csv", "Fig2F_aSMA_stacked", 55, 52), "Fig2F_aSMA_stacked", 55, 52)')
        text = text.replace(
            'save_panel(mk_stacked("IL6", "IL_6_duct", "IL_6", "Fig2H_IL6_stacked", 55, 52), "Fig2H_IL6_stacked", 55, 52)',
            'save_panel(mk_stacked("fig2_IL6.csv", "Fig2H_IL6_stacked", 55, 52), "Fig2H_IL6_stacked", 55, 52)')
    elif name == "fig3_ihc_subtype_panels.R":
        text = text.replace(
            'fig3_excel <- file.path(MOONSHOT, "past_analysis/Fig3_20251017/Database_earlyPDAC_20241127_kobepanc_fig3.xlsx")',
            'ihc_dir <- file.path(MOONSHOT, "IHC_analysis")  # de-identified per-duct IHC scores (columns: duct, score)')
        text = text.replace(
            'mk_ihc_box <- function(sheet, duct_raw, val_raw, ylab, name, gata = FALSE) {\n  d <- read_excel(fig3_excel, sheet = sheet)\n  df <- data.frame(duct = d[[duct_raw]], val = d[[val_raw]])',
            'mk_ihc_box <- function(csv, ylab, name, gata = FALSE) {\n  d <- read.csv(file.path(ihc_dir, csv))\n  df <- data.frame(duct = d$duct, val = d$score)')
        text = text.replace(
            'cat(sprintf("  %s KW p = %.3g\\n", sheet, kw))',
            'cat(sprintf("  %s KW p = %.3g\\n", name, kw))')
        text = text.replace(
            'save_panel(mk_ihc_box("TP63", "p63_1_組織", "p63(/HPF)_1カ所目", "TP63+ cells/HPF", "Fig3G_TP63_box"), "Fig3G_TP63_box", 38, 52)',
            'save_panel(mk_ihc_box("fig3_TP63.csv", "TP63+ cells/HPF", "Fig3G_TP63_box"), "Fig3G_TP63_box", 38, 52)')
        text = text.replace(
            'save_panel(mk_ihc_box("KRT5", "KRT5_1_組織", "KRT5(/HPF)_1カ所目", "KRT5+ cells/HPF", "Fig3G_KRT5_box"), "Fig3G_KRT5_box", 38, 52)',
            'save_panel(mk_ihc_box("fig3_KRT5.csv", "KRT5+ cells/HPF", "Fig3G_KRT5_box"), "Fig3G_KRT5_box", 38, 52)')
        text = text.replace(
            'save_panel(mk_ihc_box("GATA6", "GATA6_duct", "GATA6", "GATA6 expression score", "Fig3G_GATA6_box", gata = TRUE), "Fig3G_GATA6_box", 38, 52)',
            'save_panel(mk_ihc_box("fig3_GATA6.csv", "GATA6 expression score", "Fig3G_GATA6_box", gata = TRUE), "Fig3G_GATA6_box", 38, 52)')
    return text


def main() -> None:
    report = {"scripts": 0, "data": 0, "ihc": 0, "missing": []}

    # 1) figure-panel scripts (+ output panel dirs)
    for fig, scripts in FIG_SCRIPTS.items():
        src_dir = ROOT / "manuscript_figure_panels" / fig / "scripts"
        dst_dir = DST / "manuscript_figure_panels" / fig / "scripts"
        dst_dir.mkdir(parents=True, exist_ok=True)
        (DST / "manuscript_figure_panels" / fig / "panels").mkdir(parents=True, exist_ok=True)
        for s in scripts:
            sp = src_dir / s
            if not sp.exists():
                report["missing"].append(str(sp)); continue
            text = neutralize(sp.read_text(encoding="utf-8"))
            text = rewire_ihc(s, text)
            (dst_dir / s).write_text(text, encoding="utf-8")
            report["scripts"] += 1

    # 2) shared theme
    theme_src = ROOT / "common" / "theme_cancer_discovery.R"
    (DST / "common").mkdir(parents=True, exist_ok=True)
    if theme_src.exists():
        (DST / "common" / "theme_cancer_discovery.R").write_text(
            neutralize(theme_src.read_text(encoding="utf-8")), encoding="utf-8")
    else:
        report["missing"].append(str(theme_src))

    # 3) small processed inputs (mirrored relative paths)
    for rel in DATA_FILES:
        sp = ROOT / rel
        dp = DST / rel
        if not sp.exists():
            report["missing"].append(str(sp)); continue
        dp.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(sp, dp)
        report["data"] += 1

    # 4) de-identified IHC extracts
    ihc_dir = DST / "IHC_analysis"
    ihc_dir.mkdir(parents=True, exist_ok=True)
    for rel, sheet, duct_col, score_col, out_csv in IHC_EXTRACTS:
        sp = ROOT / rel
        if not sp.exists():
            report["missing"].append(str(sp)); continue
        wb = openpyxl.load_workbook(sp, read_only=True, data_only=True)
        ws = wb[sheet]
        rows = list(ws.iter_rows(values_only=True))
        header = list(rows[0])
        di = header.index(duct_col); si = header.index(score_col)
        out_rows = []
        for r in rows[1:]:
            duct = r[di] if di < len(r) else None
            score = r[si] if si < len(r) else None
            if duct is None or score is None:
                continue
            out_rows.append((duct, score))
        wb.close()
        with open(ihc_dir / out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["duct", "score"]); w.writerows(out_rows)
        report["ihc"] += 1
        print(f"  IHC: {out_csv:18s} {len(out_rows):4d} rows  (from {sheet})")

    print("\n=== assembled ===")
    print(f"  scripts: {report['scripts']}")
    print(f"  data files: {report['data']}")
    print(f"  ihc csvs: {report['ihc']}")
    if report["missing"]:
        print("  MISSING:")
        for m in report["missing"]:
            print("   -", m)


if __name__ == "__main__":
    main()
