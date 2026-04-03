#!/usr/bin/env python3
"""Pathway enrichment analysis (GO + KEGG) for each GMHS V3.3 dimension.

Reads the 1,116-gene dimension table, filters non-coding RNA IDs, and runs
Enrichr-based enrichment via gseapy for each dimension (and per-direction
sub-groups).  Outputs combined CSV, summary CSV, and a publication-quality
dot plot.

Usage:
    python pathway_enrichment_by_dimension.py
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path

import gseapy as gp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Matplotlib settings for Illustrator-editable fonts
# ---------------------------------------------------------------------------
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
mpl.rcParams["font.size"] = 10

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)

# ---------------------------------------------------------------------------
# Paths (script assumes working directory is the project root)
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(".")
INPUT_CSV = PROJECT_ROOT / "GWAS_TWAS/gmhs_dimensions/dimension_genes_v33_telomere.csv"
OUT_DIR_DATA = PROJECT_ROOT / "MyoScore/data"
OUT_DIR_FIG = PROJECT_ROOT / "MyoScore/figures"
OUT_DIR_DATA.mkdir(parents=True, exist_ok=True)
OUT_DIR_FIG.mkdir(parents=True, exist_ok=True)

OUT_FULL = OUT_DIR_DATA / "pathway_enrichment_by_dimension.csv"
OUT_SUMMARY = OUT_DIR_DATA / "pathway_enrichment_summary.csv"
OUT_FIG_PDF = OUT_DIR_FIG / "pathway_enrichment_dotplot.pdf"
OUT_FIG_PNG = OUT_DIR_FIG / "pathway_enrichment_dotplot.png"

# ---------------------------------------------------------------------------
# Colour palette (per dimension)
# ---------------------------------------------------------------------------
DIM_COLORS: dict[str, str] = {
    "Strength": "#E64B35",
    "Mass": "#4DBBD5",
    "LeanMuscle": "#00A087",
    "Youth": "#F39B7F",
    "Resilience": "#8491B4",
}

DIMENSIONS = ["Strength", "Mass", "LeanMuscle", "Youth", "Resilience"]

# Enrichr gene-set libraries to query
LIBRARIES: dict[str, str] = {
    "GO_BP": "GO_Biological_Process_2023",
    "GO_MF": "GO_Molecular_Function_2023",
    "KEGG": "KEGG_2021_Human",
}

# Regex patterns for non-coding / non-standard gene symbols to exclude
_NCRNA_PATTERNS = re.compile(
    r"^("
    r"RP\d+-"           # RP11-, RP4-, RP1- clone-based IDs
    r"|AC\d+"           # AC034220.3 etc.
    r"|AL\d+"           # AL-prefixed Ensembl names
    r"|AP\d+"           # AP-prefixed
    r"|BX\d+"           # BX-prefixed
    r"|CR\d+"           # CR-prefixed
    r"|CT[BCD]-"        # CTB-, CTC-, CTD- clone IDs
    r"|MIR\d+"          # microRNAs
    r"|LINC\d+"         # long intergenic non-coding
    r"|LOC\d+"          # uncharacterised LOC genes
    r"|SNORD\d+"        # small nucleolar RNAs
    r"|SNHG\d+"         # small nucleolar RNA host genes
    r"|RNU\d+"          # small nuclear RNAs
    r"|SCARNA\d+"       # small Cajal body RNAs
    r"|LA16c"           # LA16c clone names
    r"|KB-"             # KB- clone names
    r"|XXbac-"          # XXbac clones
    r")",
    re.IGNORECASE,
)

PADJ_THRESHOLD = 0.05


def is_protein_coding_symbol(gene: str) -> bool:
    """Return True if the gene symbol looks like a standard protein-coding gene."""
    if _NCRNA_PATTERNS.match(gene):
        return False
    # Also drop IDs with dots (Ensembl-style versioned IDs like AC034220.3)
    if re.match(r"^[A-Z]{2}\d+\.\d+$", gene):
        return False
    return True


def filter_genes(genes: list[str]) -> list[str]:
    """Filter a gene list to keep only likely protein-coding symbols."""
    filtered = [g for g in genes if is_protein_coding_symbol(g)]
    return filtered


def run_enrichr(
    gene_list: list[str],
    library: str,
    description: str,
) -> pd.DataFrame:
    """Run Enrichr enrichment and return results DataFrame."""
    if len(gene_list) < 3:
        logger.warning(
            "Skipping enrichment for '%s': only %d genes after filtering.",
            description,
            len(gene_list),
        )
        return pd.DataFrame()

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=library,
            organism="human",
            outdir=None,  # do not write intermediate files
            no_plot=True,
            verbose=False,
        )
        res = enr.results.copy()
        return res
    except Exception as exc:
        logger.error("Enrichr failed for '%s' / %s: %s", description, library, exc)
        return pd.DataFrame()


def main() -> None:
    # ------------------------------------------------------------------
    # 1. Load genes
    # ------------------------------------------------------------------
    logger.info("Reading %s", INPUT_CSV)
    df = pd.read_csv(INPUT_CSV)
    logger.info("Total genes: %d across %d dimensions", len(df), df["dimension"].nunique())

    # ------------------------------------------------------------------
    # 2. Run enrichment per dimension (all + by direction)
    # ------------------------------------------------------------------
    all_results: list[pd.DataFrame] = []

    for dim in DIMENSIONS:
        dim_df = df[df["dimension"] == dim]
        all_genes = filter_genes(dim_df["ID"].tolist())
        logger.info(
            "Dimension %s: %d total genes, %d after ncRNA filter",
            dim, len(dim_df), len(all_genes),
        )

        # Groups: all, healthy (+1), unhealthy (-1)
        groups: dict[str, list[str]] = {
            "all": all_genes,
            "healthy_up": filter_genes(
                dim_df.loc[dim_df["direction_v3"] == 1.0, "ID"].tolist()
            ),
            "unhealthy_up": filter_genes(
                dim_df.loc[dim_df["direction_v3"] == -1.0, "ID"].tolist()
            ),
        }

        for group_name, gene_list in groups.items():
            for lib_short, lib_full in LIBRARIES.items():
                desc = f"{dim}_{group_name}_{lib_short}"
                logger.info("  Running %s (%d genes) ...", desc, len(gene_list))
                res = run_enrichr(gene_list, lib_full, desc)
                if res.empty:
                    continue
                res["dimension"] = dim
                res["direction_group"] = group_name
                res["database"] = lib_short
                all_results.append(res)

    if not all_results:
        logger.error("No enrichment results obtained. Exiting.")
        sys.exit(1)

    combined = pd.concat(all_results, ignore_index=True)
    logger.info("Total enrichment rows: %d", len(combined))

    # Rename Enrichr columns for clarity
    col_rename = {
        "Term": "term",
        "Overlap": "overlap",
        "P-value": "pvalue",
        "Adjusted P-value": "padj",
        "Odds Ratio": "odds_ratio",
        "Combined Score": "combined_score",
        "Genes": "genes",
        "Gene_set": "gene_set",
    }
    combined.rename(columns={k: v for k, v in col_rename.items() if k in combined.columns}, inplace=True)

    # Compute gene_ratio from overlap string (e.g. "5/300")
    if "overlap" in combined.columns:
        def parse_ratio(x: str) -> float:
            try:
                num, denom = x.split("/")
                return int(num) / int(denom) if int(denom) > 0 else 0.0
            except Exception:
                return 0.0

        combined["gene_ratio"] = combined["overlap"].apply(parse_ratio)
        combined["overlap_count"] = combined["overlap"].apply(
            lambda x: int(x.split("/")[0]) if "/" in str(x) else 0
        )

    # ------------------------------------------------------------------
    # 3. Save full results
    # ------------------------------------------------------------------
    combined.to_csv(OUT_FULL, index=False)
    logger.info("Full results saved to %s", OUT_FULL)

    # ------------------------------------------------------------------
    # 4. Build summary (top 5 per dimension per database, all-direction only)
    # ------------------------------------------------------------------
    sig = combined[
        (combined["padj"] < PADJ_THRESHOLD) & (combined["direction_group"] == "all")
    ].copy()

    if sig.empty:
        logger.warning("No significant results at padj < %.2f for 'all' group.", PADJ_THRESHOLD)
        # Relax: take top 5 regardless of significance for summary
        summary_parts = []
        for (dim, db), grp in combined[combined["direction_group"] == "all"].groupby(
            ["dimension", "database"]
        ):
            summary_parts.append(grp.nsmallest(5, "padj"))
        summary = pd.concat(summary_parts, ignore_index=True) if summary_parts else combined.head(0)
    else:
        summary_parts = []
        for (dim, db), grp in sig.groupby(["dimension", "database"]):
            summary_parts.append(grp.nsmallest(5, "padj"))
        summary = pd.concat(summary_parts, ignore_index=True)

    summary.to_csv(OUT_SUMMARY, index=False)
    logger.info("Summary saved to %s (%d rows)", OUT_SUMMARY, len(summary))

    # ------------------------------------------------------------------
    # 5. Dot plot
    # ------------------------------------------------------------------
    _make_dotplot(combined)
    logger.info("Done.")


def _make_dotplot(combined: pd.DataFrame) -> None:
    """Create a publication-quality dot plot of top pathways per dimension."""

    # Use 'all' direction group, significant results
    plot_df = combined[
        (combined["direction_group"] == "all") & (combined["padj"] < PADJ_THRESHOLD)
    ].copy()

    if plot_df.empty:
        logger.warning("No significant 'all'-group results for dot plot. Using top results by p-value.")
        plot_df = combined[combined["direction_group"] == "all"].copy()
        if plot_df.empty:
            logger.error("Cannot create dot plot: no results at all.")
            return

    # Select top 3 terms per dimension per database (by padj)
    top_parts = []
    for (dim, db), grp in plot_df.groupby(["dimension", "database"]):
        top_parts.append(grp.nsmallest(3, "padj"))
    top = pd.concat(top_parts, ignore_index=True)

    if top.empty:
        logger.error("No terms to plot.")
        return

    # Shorten term names (remove GO:xxxx or (GO:xxxx) suffix for readability)
    top["term_short"] = top["term"].apply(
        lambda t: re.sub(r"\s*\(GO:\d+\)$", "", str(t))[:60]
    )

    top["-log10(padj)"] = -np.log10(top["padj"].clip(lower=1e-30))

    # Order dimensions
    top["dim_order"] = top["dimension"].map({d: i for i, d in enumerate(DIMENSIONS)})
    top.sort_values(["dim_order", "database", "padj"], inplace=True)

    # Create a combined label for y-axis
    top["ylabel"] = top["term_short"] + "  [" + top["database"] + "]"

    # De-duplicate y-labels that appear for multiple dimensions: keep the best
    # (We want unique y entries; if same term appears in multiple dims, keep all)
    top = top.drop_duplicates(subset=["ylabel", "dimension"])

    # Limit total terms to avoid overcrowding (max ~50)
    if len(top) > 50:
        top = top.groupby("dimension").head(9)

    # Reverse so first dimension appears at top
    top = top.iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(12, max(10, len(top) * 0.28)))

    for dim in DIMENSIONS:
        subset = top[top["dimension"] == dim]
        if subset.empty:
            continue
        ax.scatter(
            subset["-log10(padj)"],
            subset["ylabel"],
            s=subset["gene_ratio"] * 800 + 20,  # scale dot size
            c=DIM_COLORS[dim],
            alpha=0.8,
            edgecolors="white",
            linewidths=0.5,
            label=dim,
            zorder=3,
        )

    ax.set_xlabel("-log$_{10}$(adjusted p-value)", fontsize=11)
    ax.set_ylabel("")
    ax.set_title("Pathway Enrichment by GMHS Dimension", fontsize=13, fontweight="bold")
    ax.axvline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8, zorder=1)
    ax.legend(
        title="Dimension",
        loc="lower right",
        framealpha=0.9,
        fontsize=9,
        title_fontsize=10,
    )

    # Add size legend
    for ratio_val, label in [(0.05, "5%"), (0.10, "10%"), (0.20, "20%")]:
        ax.scatter(
            [], [],
            s=ratio_val * 800 + 20,
            c="grey",
            alpha=0.5,
            label=f"Gene ratio {label}",
        )
    # Re-draw legend with both colour and size entries
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles, labels,
        loc="lower right",
        framealpha=0.9,
        fontsize=8,
        title="Dimension / Gene ratio",
        title_fontsize=9,
    )

    ax.tick_params(axis="y", labelsize=8)
    plt.tight_layout()

    fig.savefig(OUT_FIG_PDF, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_FIG_PNG, dpi=300, bbox_inches="tight")
    logger.info("Dot plot saved to %s and %s", OUT_FIG_PDF, OUT_FIG_PNG)
    plt.close(fig)


if __name__ == "__main__":
    main()
