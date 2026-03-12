"""GSEA prerank enrichment runner.

Wraps gseapy.prerank with configurable databases and parameters.
Provides ranked gene list creation, result parsing, and normalization.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Optional

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rnaseq_explorer.viz.theme import COLOR_UP, COLOR_DOWN


# ---------------------------------------------------------------------------
# GSEA column normalization
# ---------------------------------------------------------------------------


def normalize_gsea_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize gseapy result column names to canonical lowercase forms.

    Handles version differences in gseapy (uppercase vs lowercase,
    'FDR q-val' vs 'fdr' vs 'FDR', etc.).

    Parameters
    ----------
    df : pd.DataFrame
        Raw gseapy results DataFrame.

    Returns
    -------
    pd.DataFrame
        DataFrame with standardized column names.
    """
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip().replace(" ", "_").replace("-", "_")
        if cl in ("fdr", "fdr_q_val", "fdr_bh", "padj", "adj_p_value"):
            col_map[c] = "fdr"
        elif cl in ("nes", "normalized_enrichment_score"):
            col_map[c] = "nes"
        elif cl in ("term", "pathway", "gene_set"):
            col_map[c] = "Term"
        elif cl in ("geneset_size", "gene_set_size", "gs_size"):
            col_map[c] = "geneset_size"
        elif cl in ("lead_genes", "leading_edge", "lead_edge_genes"):
            col_map[c] = "lead_genes"
        elif cl in ("pvalue", "pval", "p_value", "nom_p_val", "nom_p_value"):
            col_map[c] = "pvalue"
        elif cl in ("es", "enrichment_score"):
            col_map[c] = "es"
        elif cl in ("tag_%", "tag_percent"):
            col_map[c] = "tag_pct"
        elif cl in ("gene_%", "gene_percent"):
            col_map[c] = "gene_pct"
    if col_map:
        df = df.rename(columns=col_map)

    # Derive geneset_size from tag_pct if missing (gseapy 1.1+)
    if "geneset_size" not in df.columns and "tag_pct" in df.columns:
        def _parse_gs_size(val):
            try:
                return int(str(val).split("/")[1])
            except (IndexError, ValueError):
                return 0
        df["geneset_size"] = df["tag_pct"].apply(_parse_gs_size)
    return df


# ---------------------------------------------------------------------------
# Ranked gene list creation
# ---------------------------------------------------------------------------


def create_ranked_list(
    deseq2_raw: pd.DataFrame,
    gene_name_col: str,
    ranking_col: str = "stat",
    log2fc_col: str = "log2FoldChange",
    stat_col: str = "stat",
    ranking_method: str = "stat",
) -> pd.Series:
    """Create a ranked gene list for GSEA prerank from full DESeq2 results.

    Parameters
    ----------
    deseq2_raw : pd.DataFrame
        Full unfiltered DESeq2 data (all genes).
    gene_name_col : str
        Column containing gene names/symbols.
    ranking_col : str
        Column to use for ranking (determined by ranking_method).
    log2fc_col : str
        Column name for log2 fold change.
    stat_col : str
        Column name for Wald test statistic.
    ranking_method : str
        'stat' (Wald statistic, preferred) or 'log2fc'.

    Returns
    -------
    pd.Series
        Ranked gene list (index=gene names, values=rank scores), sorted descending.
    """
    use_stat = (
        ranking_method == "stat"
        and stat_col in deseq2_raw.columns
        and deseq2_raw[stat_col].notna().sum() > 0
    )
    if use_stat:
        rank_col = stat_col
        print(f"    Ranking by Wald statistic (column: '{stat_col}') -- WSF method")
    else:
        rank_col = log2fc_col
        if ranking_method == "stat":
            print(f"    [INFO] Wald stat column not found -- falling back to log2FC ranking")
        else:
            print(f"    Ranking by log2FC (column: '{log2fc_col}')")

    ranked_genes = deseq2_raw[[gene_name_col, rank_col]].dropna()
    # Remove duplicate gene names (keep the one with largest absolute rank value)
    ranked_genes["abs_rank"] = ranked_genes[rank_col].abs()
    ranked_genes = (
        ranked_genes.sort_values("abs_rank", ascending=False)
        .drop_duplicates(subset=gene_name_col, keep="first")
        .drop(columns="abs_rank")
    )
    ranked_genes = ranked_genes.sort_values(rank_col, ascending=False)
    ranked_genes = ranked_genes.set_index(gene_name_col)[rank_col]
    return ranked_genes


# ---------------------------------------------------------------------------
# GSEA runner
# ---------------------------------------------------------------------------


def run_gsea_prerank(
    ranked_genes: pd.Series,
    databases: list[str],
    outdir: str | Path,
    min_size: int = 15,
    max_size: int = 500,
    permutations: int = 1000,
    fdr_threshold: float = 0.25,
    top_n: int = 5,
    seed: int = 42,
) -> dict[str, pd.DataFrame]:
    """Run GSEA prerank analysis across multiple databases.

    Parameters
    ----------
    ranked_genes : pd.Series
        Ranked gene list (index=genes, values=rank metric).
    databases : list[str]
        Enrichr-compatible database names.
    outdir : str or Path
        Directory to save gseapy output.
    min_size : int
        Minimum gene set size.
    max_size : int
        Maximum gene set size.
    permutations : int
        Number of permutations.
    fdr_threshold : float
        FDR threshold for significance.
    top_n : int
        Number of top pathways to return per database.
    seed : int
        Random seed.

    Returns
    -------
    dict
        Mapping of database name -> DataFrame of top significant pathways.
    """
    try:
        import gseapy as gp
    except ImportError:
        print("  WARNING: gseapy not installed, skipping GSEA analysis")
        print("  Install with: pip install gseapy")
        return {}

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    results: dict[str, pd.DataFrame] = {}

    for db in databases:
        try:
            print(f"    Running {db}...")
            prerank_res = gp.prerank(
                rnk=ranked_genes,
                gene_sets=db,
                outdir=str(outdir / db.replace(":", "_")),
                min_size=min_size,
                max_size=max_size,
                permutation_num=permutations,
                seed=seed,
                verbose=False,
            )

            res_df = normalize_gsea_cols(prerank_res.res2d)
            if len(res_df) > 0:
                sig_pathways = (
                    res_df[res_df["fdr"] < fdr_threshold]
                    .sort_values("nes", key=abs, ascending=False)
                    .head(top_n)
                )
                if len(sig_pathways) > 0:
                    results[db] = sig_pathways[
                        ["Term", "nes", "fdr", "geneset_size", "lead_genes"]
                    ]
                    print(f"      Found {len(sig_pathways)} significant pathways")
                else:
                    print(f"      No significant pathways (FDR < {fdr_threshold})")
            else:
                print(f"      No results from {db}")

        except Exception as e:
            print(f"      ERROR with {db}: {e}")
            continue

    return results


# ---------------------------------------------------------------------------
# Full GSEA enrichment pipeline
# ---------------------------------------------------------------------------


def run_gsea_enrichment(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    databases: list[str] | None = None,
    ranking_method: str = "stat",
    min_size: int = 15,
    max_size: int = 500,
    permutations: int = 1000,
) -> dict[str, dict[str, pd.DataFrame]]:
    """Run GSEA prerank enrichment analysis for each condition.

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Output directory (a gsea_results/ subdirectory will be created).
    cols : dict or None
        DESeq2 column name mapping.
    databases : list[str] or None
        GSEA databases to query.
    ranking_method : str
        'stat' or 'log2fc'.
    min_size : int
        Minimum gene set size.
    max_size : int
        Maximum gene set size.
    permutations : int
        Number of permutations.

    Returns
    -------
    dict
        gsea_results[cond_name] = {db_name: DataFrame}
    """
    from rnaseq_explorer.engine.deseq2 import DEFAULT_DESEQ2_COLS

    if cols is None:
        cols = DEFAULT_DESEQ2_COLS
    if databases is None:
        databases = [
            "GO_Biological_Process_2023",
            "GO_Cellular_Component_2023",
            "GO_Molecular_Function_2023",
            "KEGG_2021_Human",
            "Reactome_2022",
            "MSigDB_Hallmark_2020",
            "WikiPathway_2021_Human",
        ]

    outdir = Path(outdir)
    gsea_dir = outdir / "gsea_results"
    gsea_dir.mkdir(exist_ok=True)

    print("\n-- Running GSEA Enrichment --")

    gene_name_col = cols.get("gene_name", "gene_name")
    log2fc_col = cols.get("log2fc", "log2FoldChange")
    stat_col = cols.get("stat", "stat")

    gsea_results: dict[str, dict[str, pd.DataFrame]] = {}

    for cond_name, data in condition_results.items():
        cond_label = condition_labels.get(cond_name, cond_name)
        print(f"  Processing {cond_label}...")

        deseq2_full = data.get("deseq2_raw", pd.DataFrame())
        if len(deseq2_full) == 0:
            print(f"    No DESeq2 data for {cond_label}, skipping GSEA")
            continue

        ranked_genes = create_ranked_list(
            deseq2_full,
            gene_name_col=gene_name_col,
            log2fc_col=log2fc_col,
            stat_col=stat_col,
            ranking_method=ranking_method,
        )

        if len(ranked_genes) < 15:
            print(f"    Too few genes ({len(ranked_genes)}) for GSEA in {cond_label}")
            continue
        print(f"    Ranking {len(ranked_genes)} genes for GSEA")

        cond_gsea_dir = gsea_dir / cond_name
        cond_gsea_dir.mkdir(exist_ok=True)

        cond_results = run_gsea_prerank(
            ranked_genes,
            databases=databases,
            outdir=cond_gsea_dir,
            min_size=min_size,
            max_size=max_size,
            permutations=permutations,
        )

        gsea_results[cond_name] = cond_results

        # Export per-condition GSEA summary
        if cond_results:
            summary_rows = []
            for db, df in cond_results.items():
                for _, row in df.iterrows():
                    summary_rows.append({
                        "Database": db,
                        "Pathway": row["Term"],
                        "NES": row["nes"],
                        "FDR": row["fdr"],
                        "Genes": row["geneset_size"],
                        "LeadingEdge": row["lead_genes"],
                    })
            if summary_rows:
                summary_df = pd.DataFrame(summary_rows)
                summary_path = cond_gsea_dir / f"GSEA_summary_{cond_name}.xlsx"
                summary_df.to_excel(summary_path, index=False)
                print(f"    Saved: {summary_path.name}")

    print(f"  GSEA complete: {len(gsea_results)} conditions processed")
    return gsea_results


def parse_gsea_results(
    gsea_results: dict[str, dict[str, pd.DataFrame]],
    condition_labels: dict[str, str],
) -> pd.DataFrame:
    """Flatten GSEA results into a single summary DataFrame.

    Parameters
    ----------
    gsea_results : dict
        Output of run_gsea_enrichment().
    condition_labels : dict
        Maps condition name -> human-readable label.

    Returns
    -------
    pd.DataFrame
        Combined summary with columns: Condition, Database, Term, NES, FDR, geneset_size.
    """
    rows = []
    for cond_name, cond_data in gsea_results.items():
        cond_label = condition_labels.get(cond_name, cond_name)
        for db, df in cond_data.items():
            for _, row in df.iterrows():
                rows.append({
                    "Condition": cond_label,
                    "Database": db,
                    "Term": row.get("Term", ""),
                    "NES": row.get("nes", np.nan),
                    "FDR": row.get("fdr", np.nan),
                    "geneset_size": row.get("geneset_size", 0),
                })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


# ---------------------------------------------------------------------------
# GSEA visualization helpers
# ---------------------------------------------------------------------------

DB_SHORT: dict[str, str] = {
    "GO_Biological_Process_2023": "BP",
    "GO_Cellular_Component_2023": "CC",
    "GO_Molecular_Function_2023": "MF",
    "KEGG_2021_Human": "KEGG",
    "Reactome_2022": "Reactome",
    "MSigDB_Hallmark_2020": "Hallmark",
    "WikiPathway_2021_Human": "WikiPath",
}

# Reverse lookup: folder name on disk -> short label
# gseapy replaces colons with underscores in directory names
_DB_SHORT_FUZZY: dict[str, str] = {}
for _full, _short in DB_SHORT.items():
    _DB_SHORT_FUZZY[_full] = _short
    _DB_SHORT_FUZZY[_full.replace(":", "_")] = _short


def _db_short_label(db_name: str) -> str:
    """Return the short label for a database name, with fuzzy matching."""
    if db_name in _DB_SHORT_FUZZY:
        return _DB_SHORT_FUZZY[db_name]
    # Fallback: try substring matching
    for full, short in DB_SHORT.items():
        if full.replace(":", "_") in db_name or db_name in full:
            return short
    # Last resort: first 8 chars
    return db_name[:8]


def _collect_gsea_rows(
    gsea_results: dict | None,
    cond_name: str,
    gsea_dir: Path,
    top_n: int = 5,
    include_lead_genes: bool = False,
) -> list[dict]:
    """Collect GSEA pathway rows from in-memory results or disk CSV fallback.

    Returns a list of dicts with keys: Term, NES, FDR, Database,
    DatabaseShort, and optionally LeadingEdgeGenes.
    """
    all_rows: list[dict] = []

    # --- Try in-memory results first ---
    if gsea_results and cond_name in gsea_results:
        for db_name, pathways_df in gsea_results[cond_name].items():
            if pathways_df is None or len(pathways_df) == 0:
                continue
            df = pathways_df.copy()
            df["abs_nes"] = df["nes"].abs()
            df = df.sort_values("abs_nes", ascending=False).head(top_n)
            short = _db_short_label(db_name)
            for _, row in df.iterrows():
                entry: dict = {
                    "Term": str(row["Term"]),
                    "NES": float(row["nes"]),
                    "FDR": float(row["fdr"]),
                    "Database": db_name,
                    "DatabaseShort": short,
                }
                if include_lead_genes and "lead_genes" in row.index:
                    entry["LeadingEdgeGenes"] = (
                        str(row["lead_genes"]) if pd.notna(row["lead_genes"]) else ""
                    )
                all_rows.append(entry)

    # --- Fallback: scan disk CSV reports ---
    if not all_rows:
        cond_gsea_dir = gsea_dir / cond_name
        if cond_gsea_dir.is_dir():
            for db_subdir in sorted(cond_gsea_dir.iterdir()):
                if not db_subdir.is_dir():
                    continue
                report_csv = db_subdir / "gseapy.gene_set.prerank.report.csv"
                if not report_csv.exists():
                    continue
                try:
                    rpt = normalize_gsea_cols(pd.read_csv(report_csv))
                    if len(rpt) == 0:
                        continue
                    rpt["abs_nes"] = rpt["nes"].abs()
                    rpt = rpt.sort_values("abs_nes", ascending=False).head(top_n)
                    db_name_disk = db_subdir.name
                    short = _db_short_label(db_name_disk)
                    for _, row in rpt.iterrows():
                        entry = {
                            "Term": str(row["Term"]),
                            "NES": float(row["nes"]),
                            "FDR": float(row["fdr"]),
                            "Database": db_name_disk,
                            "DatabaseShort": short,
                        }
                        if include_lead_genes and "lead_genes" in row.index:
                            entry["LeadingEdgeGenes"] = (
                                str(row["lead_genes"]) if pd.notna(row["lead_genes"]) else ""
                            )
                        all_rows.append(entry)
                except Exception as e:
                    logging.debug("Skipping GSEA row: %s", e)
                    continue

    return all_rows


# ---------------------------------------------------------------------------
# GSEA visualization functions
# ---------------------------------------------------------------------------


def gsea_dotplot_legacy(
    gsea_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate horizontal dot plots of top GSEA pathways per condition.

    Falls back to reading gseapy CSV reports from disk when in-memory
    results are empty.

    NOTE: Legacy function, superseded by gsea_combined_plot().

    Parameters
    ----------
    gsea_results : dict
        Output of run_gsea_enrichment().
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Base output directory.
    fig_format : str
        Image format (e.g. 'png', 'svg').
    fig_dpi : int
        Resolution in DPI.
    """
    outdir = Path(outdir)
    gsea_dir = outdir / "gsea_results"
    gsea_dir.mkdir(exist_ok=True)

    print("\n-- Generating GSEA Dot Plots --")

    for cond_name, cond_label in condition_labels.items():
        all_rows: list[dict] = []

        if gsea_results and cond_name in gsea_results:
            for db_name, pathways_df in gsea_results[cond_name].items():
                if len(pathways_df) > 0:
                    for _, row in pathways_df.iterrows():
                        all_rows.append({
                            "Term": row["Term"],
                            "NES": float(row["nes"]),
                            "FDR": float(row["fdr"]),
                        })

        if not all_rows:
            # Fallback: scan disk CSV reports
            cond_gsea_dir = gsea_dir / cond_name
            if cond_gsea_dir.is_dir():
                for db_subdir in sorted(cond_gsea_dir.iterdir()):
                    if not db_subdir.is_dir():
                        continue
                    report_csv = db_subdir / "gseapy.gene_set.prerank.report.csv"
                    if not report_csv.exists():
                        continue
                    try:
                        rpt = normalize_gsea_cols(pd.read_csv(report_csv))
                        if len(rpt) == 0:
                            continue
                        rpt["abs_nes"] = rpt["nes"].abs()
                        rpt = rpt.sort_values("abs_nes", ascending=False).head(5)
                        for _, row in rpt.iterrows():
                            all_rows.append({
                                "Term": row["Term"],
                                "NES": float(row["nes"]),
                                "FDR": float(row["fdr"]),
                            })
                    except Exception as e:
                        logging.debug("Skipping GSEA fallback row: %s", e)
                        continue

        if not all_rows:
            print(f"  No GSEA data for {cond_label}, skipping dotplot")
            continue

        df = pd.DataFrame(all_rows)
        df["abs_NES"] = df["NES"].abs()
        df = df.sort_values("abs_NES", ascending=False).head(20)
        df = df.sort_values("NES", ascending=True)

        df["Term_short"] = df["Term"].str[:50]

        neg_log_fdr = -np.log10(df["FDR"].values + 1e-10)
        size_min, size_max = 20, 200
        if neg_log_fdr.max() > neg_log_fdr.min():
            sizes = size_min + (neg_log_fdr - neg_log_fdr.min()) / (
                neg_log_fdr.max() - neg_log_fdr.min()
            ) * (size_max - size_min)
        else:
            sizes = np.full_like(neg_log_fdr, (size_min + size_max) / 2)

        colors = [COLOR_UP if nes > 0 else COLOR_DOWN for nes in df["NES"].values]

        fig, ax = plt.subplots(figsize=(10, max(5, len(df) * 0.35)))
        ax.scatter(
            df["NES"].values, range(len(df)), s=sizes, c=colors,
            edgecolors="black", linewidths=0.5, zorder=3,
        )
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df["Term_short"].values, fontsize=9)
        ax.axvline(0, color="grey", linewidth=0.8, linestyle="--", zorder=1)
        ax.set_xlabel("NES (Normalized Enrichment Score)", fontsize=11)
        ax.set_title(
            f"GSEA Top Pathways \u2014 {cond_label}", fontsize=13, fontweight="bold"
        )

        for fdr_val, label in [(0.05, "FDR=0.05"), (0.25, "FDR=0.25")]:
            sz = size_min + (
                -np.log10(fdr_val + 1e-10) - neg_log_fdr.min()
            ) / max(neg_log_fdr.max() - neg_log_fdr.min(), 1e-10) * (size_max - size_min)
            sz = np.clip(sz, size_min, size_max)
            ax.scatter([], [], s=sz, c="grey", edgecolors="black",
                       linewidths=0.5, label=label)
        ax.legend(title="Dot size", loc="lower right", fontsize=8, title_fontsize=9)

        plt.tight_layout()
        outpath = gsea_dir / f"gsea_dotplot_{cond_name}.{fig_format}"
        fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {outpath.name}")


def gsea_combined_plot(
    gsea_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate combined GSEA dot plots showing pathways from ALL databases.

    Replaces the older gsea_dotplot() function. For each condition, creates
    one figure with database-tagged pathway names, NES on x-axis, dot size
    proportional to -log10(FDR), and color indicating enrichment direction.

    Falls back to reading gseapy CSV reports from disk when in-memory
    results are empty.

    Parameters
    ----------
    gsea_results : dict
        Output of run_gsea_enrichment().
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Base output directory.
    fig_format : str
        Image format.
    fig_dpi : int
        Resolution in DPI.
    """
    outdir = Path(outdir)
    gsea_dir = outdir / "gsea_results"
    gsea_dir.mkdir(exist_ok=True)

    print("\n-- Generating Combined GSEA Plots --")

    for cond_name, cond_label in condition_labels.items():
        rows = _collect_gsea_rows(gsea_results, cond_name, gsea_dir, top_n=5)

        if not rows:
            print(f"  No GSEA data for {cond_label}, skipping combined plot")
            continue

        df = pd.DataFrame(rows)

        # Cap at 35 pathways total (top 5 per database, 7 databases)
        df["abs_NES"] = df["NES"].abs()
        df = df.sort_values("abs_NES", ascending=False).head(35)

        # Sort by NES: positive at top, negative at bottom
        df = df.sort_values("NES", ascending=True).reset_index(drop=True)

        # Build tagged pathway names: "[KEGG] PI3K-Akt signaling"
        df["TaggedTerm"] = df.apply(
            lambda r: f"[{r['DatabaseShort']}] {r['Term'][:55]}", axis=1
        )

        # Deduplicate display names
        seen: dict[str, int] = {}
        unique_names: list[str] = []
        for name in df["TaggedTerm"]:
            if name in seen:
                seen[name] += 1
                unique_names.append(f"{name} ({seen[name]})")
            else:
                seen[name] = 0
                unique_names.append(name)
        df["TaggedTerm"] = unique_names

        n_pathways = len(df)
        n_databases = df["DatabaseShort"].nunique()

        # -log10(FDR) for dot size
        neg_log_fdr = -np.log10(df["FDR"].values + 1e-10)
        size_min, size_max = 30, 250
        if neg_log_fdr.max() > neg_log_fdr.min():
            sizes = size_min + (
                (neg_log_fdr - neg_log_fdr.min())
                / (neg_log_fdr.max() - neg_log_fdr.min())
                * (size_max - size_min)
            )
        else:
            sizes = np.full_like(neg_log_fdr, (size_min + size_max) / 2)

        # Dot color: orange for positive NES, blue for negative
        colors = [COLOR_UP if nes > 0 else COLOR_DOWN for nes in df["NES"].values]

        fig_h = max(8, n_pathways * 0.4)
        fig, ax = plt.subplots(figsize=(14, fig_h))

        ax.scatter(
            df["NES"].values, range(n_pathways), s=sizes, c=colors,
            edgecolors="black", linewidths=0.5, zorder=3, alpha=0.85,
        )

        ax.set_yticks(range(n_pathways))
        ax.set_yticklabels(df["TaggedTerm"].values, fontsize=9)
        ax.axvline(0, color="grey", linewidth=0.8, linestyle="--", zorder=1)
        ax.set_xlabel("NES (Normalized Enrichment Score)", fontsize=12)
        ax.set_title(
            f"GSEA Enrichment \u2014 {cond_label}", fontsize=14, fontweight="bold",
        )

        ax.set_axisbelow(True)
        ax.yaxis.grid(True, linestyle=":", alpha=0.3)

        # Legend: dot size scale
        legend_elements = []
        fdr_examples = [0.001, 0.01, 0.05, 0.25]
        for fdr_val in fdr_examples:
            nlf = -np.log10(fdr_val + 1e-10)
            if neg_log_fdr.max() > neg_log_fdr.min():
                sz = size_min + (nlf - neg_log_fdr.min()) / (
                    neg_log_fdr.max() - neg_log_fdr.min()
                ) * (size_max - size_min)
            else:
                sz = (size_min + size_max) / 2
            sz = np.clip(sz, size_min, size_max)
            legend_elements.append(
                ax.scatter(
                    [], [], s=sz, c="grey", edgecolors="black",
                    linewidths=0.5, label=f"FDR={fdr_val}",
                )
            )

        legend_elements.append(
            mpatches.Patch(
                facecolor=COLOR_UP, edgecolor="black",
                linewidth=0.5, label="Positive NES (up)",
            )
        )
        legend_elements.append(
            mpatches.Patch(
                facecolor=COLOR_DOWN, edgecolor="black",
                linewidth=0.5, label="Negative NES (down)",
            )
        )

        ax.legend(
            handles=legend_elements, title="Significance & Direction",
            loc="lower right", fontsize=8, title_fontsize=9, framealpha=0.9,
        )

        plt.tight_layout()
        outpath = gsea_dir / f"gsea_combined_{cond_name}.{fig_format}"
        fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
        plt.close(fig)
        print(
            f"  Saved: gsea_combined_{cond_name}.{fig_format} "
            f"({n_pathways} pathways from {n_databases} databases)"
        )


def gsea_enrichment_plots(
    gsea_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Collect or regenerate per-pathway enrichment score plots.

    Strategy:
      A) Look for pre-generated plots produced by gseapy.prerank() in the
         output directory structure.
      B) If none found, attempt to regenerate using gseapy's plotting API.

    Organizes all enrichment plots into:
        outdir/gsea_results/{cond_name}/enrichment_plots/

    Parameters
    ----------
    gsea_results : dict
        Output of run_gsea_enrichment().
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Base output directory.
    fig_format : str
        Image format.
    fig_dpi : int
        Resolution in DPI.
    """
    outdir = Path(outdir)
    gsea_dir = outdir / "gsea_results"
    gsea_dir.mkdir(exist_ok=True)

    print("\n-- Collecting GSEA Enrichment Plots --")

    for cond_name, cond_label in condition_labels.items():
        cond_gsea_dir = gsea_dir / cond_name
        if not cond_gsea_dir.is_dir():
            print(f"  No GSEA output directory for {cond_label}, skipping")
            continue

        plots_dest = cond_gsea_dir / "enrichment_plots"
        plots_dest.mkdir(exist_ok=True)

        collected = 0

        # --- Option A: scan for pre-generated plot files ---
        for db_subdir in sorted(cond_gsea_dir.iterdir()):
            if not db_subdir.is_dir() or db_subdir.name == "enrichment_plots":
                continue

            short_label = _db_short_label(db_subdir.name)

            plot_files = (
                list(db_subdir.glob("*.png"))
                + list(db_subdir.glob("*.pdf"))
                + list(db_subdir.glob("*.svg"))
            )

            plot_files = [
                f for f in plot_files
                if f.suffix.lower() in (".png", ".pdf", ".svg")
                and "report" not in f.stem.lower()
            ]

            for pf in plot_files:
                dest_name = f"{short_label}_{pf.name}"
                dest_path = plots_dest / dest_name
                try:
                    shutil.copy2(str(pf), str(dest_path))
                    collected += 1
                except Exception as e:
                    logging.debug("Could not copy %s: %s", pf, e)

        # --- Option B: try to regenerate if none found ---
        if collected == 0:
            try:
                from gseapy.plot import gseaplot
                _has_gseaplot = True
            except ImportError:
                _has_gseaplot = False

            if _has_gseaplot:
                for db_subdir in sorted(cond_gsea_dir.iterdir()):
                    if not db_subdir.is_dir() or db_subdir.name == "enrichment_plots":
                        continue

                    short_label = _db_short_label(db_subdir.name)
                    report_csv = db_subdir / "gseapy.gene_set.prerank.report.csv"
                    if not report_csv.exists():
                        continue

                    try:
                        rpt = normalize_gsea_cols(pd.read_csv(report_csv))
                        sig = rpt[rpt["fdr"] < 0.25].sort_values(
                            "nes", key=abs, ascending=False
                        )
                        if len(sig) == 0:
                            continue

                        for _, row in sig.head(10).iterrows():
                            term = str(row["Term"])
                            term_clean = (
                                term.replace("/", "_")
                                .replace("\\", "_")
                                .replace(":", "_")
                                .replace(" ", "_")
                            )
                            rank_file = db_subdir / term_clean / "gene_set_prerank.csv"
                            if not rank_file.exists():
                                rank_file = db_subdir / f"{term_clean}.csv"
                            if not rank_file.exists():
                                continue

                            try:
                                out_name = f"{short_label}_{term_clean[:60]}.{fig_format}"
                                out_path = plots_dest / out_name
                                gseaplot(
                                    rank_metric=rank_file,
                                    term=term,
                                    ofname=str(out_path),
                                )
                                collected += 1
                            except Exception as e:
                                logging.debug("Could not regenerate plot for %s: %s", term, e)
                                continue
                    except Exception as e:
                        logging.debug("Error processing GSEA database results: %s", e)
                        continue

        if collected > 0:
            print(
                f"  {cond_label}: {collected} enrichment plots -> "
                f"{plots_dest.relative_to(outdir)}/"
            )
        else:
            print(f"  {cond_label}: no enrichment plots found or generated")


def export_gsea_leading_edge(
    gsea_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
) -> None:
    """Export leading edge genes for significant GSEA pathways to Excel.

    For each condition and database, extracts pathways with FDR < 0.25 and
    their leading edge gene lists.  Saves to:
        outdir / "gsea_leading_edge_genes.xlsx"

    Falls back to reading disk CSV reports if in-memory results are empty.

    Parameters
    ----------
    gsea_results : dict
        Output of run_gsea_enrichment().
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Base output directory.
    """
    outdir = Path(outdir)
    gsea_dir = outdir / "gsea_results"
    gsea_dir.mkdir(exist_ok=True)

    print("\n-- Exporting GSEA Leading Edge Genes --")

    all_entries: list[dict] = []

    for cond_name, cond_label in condition_labels.items():
        rows = _collect_gsea_rows(
            gsea_results, cond_name, gsea_dir,
            top_n=999,
            include_lead_genes=True,
        )

        if not rows:
            print(f"  No GSEA data for {cond_label}, skipping")
            continue

        for r in rows:
            fdr = r.get("FDR", 1.0)
            if fdr >= 0.25:
                continue
            all_entries.append({
                "Condition": cond_label,
                "Database": r.get("Database", ""),
                "Pathway": r.get("Term", ""),
                "NES": r.get("NES", np.nan),
                "FDR": fdr,
                "LeadingEdgeGenes": r.get("LeadingEdgeGenes", ""),
            })

    if not all_entries:
        print("  No significant pathways (FDR < 0.25) found across any condition")
        return

    le_df = pd.DataFrame(all_entries)
    le_df = le_df.sort_values(["Condition", "FDR"], ascending=[True, True])

    le_path = outdir / "gsea_leading_edge_genes.xlsx"
    try:
        le_df.to_excel(le_path, index=False, sheet_name="Leading Edge Genes")
        print(
            f"  Saved: gsea_leading_edge_genes.xlsx "
            f"({len(le_df)} pathways across "
            f"{le_df['Condition'].nunique()} conditions)"
        )
    except Exception as e:
        print(f"  ERROR saving leading edge Excel: {e}")
        return

    # Try to append a "Leading Edge" sheet to the main cross-condition Excel
    main_xlsx = outdir / "cross_condition_summary.xlsx"
    if main_xlsx.exists():
        try:
            import openpyxl  # noqa: F401
            with pd.ExcelWriter(
                str(main_xlsx), engine="openpyxl", mode="a",
                if_sheet_exists="replace",
            ) as writer:
                le_df.to_excel(writer, index=False, sheet_name="Leading Edge")
            print(f"  Also added 'Leading Edge' sheet to {main_xlsx.name}")
        except ImportError:
            print("  openpyxl not installed; skipped appending to main Excel")
        except Exception as e:
            print(f"  Could not append to {main_xlsx.name}: {e}")
