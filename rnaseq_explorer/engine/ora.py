"""Over-Representation Analysis (ORA) runner.

Supports two backends:
- g:Profiler (gprofiler-official) with g:SCS FDR correction
- Enrichr (via gseapy) as legacy fallback

Both return DataFrames with a common schema: Term, Adjusted_P_value,
Overlap_count, Category.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional
from xml.dom import minidom

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rnaseq_explorer.engine.deseq2 import (
    DEFAULT_DESEQ2_COLS,
    best_gene_key,
)
from rnaseq_explorer.viz.theme import CATEGORY_COLORS

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_ORA_DATABASES: list[str] = [
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "GO_Molecular_Function_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
]

DB_TO_CATEGORY: dict[str, str] = {
    "GO_Biological_Process_2023": "BP",
    "GO_Cellular_Component_2023": "CC",
    "GO_Molecular_Function_2023": "MF",
    "KEGG_2021_Human": "KEGG",
    "Reactome_2022": "Reactome",
}

CATEGORY_COLORS: dict[str, str] = {
    "BP": "#0072B2",
    "CC": "#E69F00",
    "MF": "#009E73",
    "KEGG": "#CC79A7",
    "Reactome": "#56B4E9",
}

GPROFILER_SPECIES: dict[str, str] = {
    "human": "hsapiens",
    "mouse": "mmusculus",
    "rat": "rnorvegicus",
    "zebrafish": "drerio",
    "fly": "dmelanogaster",
    "worm": "celegans",
}


# ---------------------------------------------------------------------------
# g:Profiler ORA
# ---------------------------------------------------------------------------


def run_gprofiler_ora(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    species: str = "human",
    padj_threshold: float = 0.05,
    top_n_per_db: int = 10,
) -> dict[str, dict[str, pd.DataFrame]]:
    """Run g:Profiler over-representation analysis for up/down DEGs.

    Uses g:SCS FDR correction (hierarchy-aware, superior to BH for GO terms).

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Output directory.
    cols : dict or None
        DESeq2 column name mapping.
    species : str
        Species name.
    padj_threshold : float
        Adjusted p-value threshold for significance.
    top_n_per_db : int
        Number of top terms to keep per database.

    Returns
    -------
    dict
        go_results[cond_name] = {"up": DataFrame, "down": DataFrame}
    """
    try:
        from gprofiler import GProfiler
    except ImportError:
        print("  [INFO] gprofiler-official not installed -- falling back to Enrichr")
        print("  Install with: pip install gprofiler-official")
        return run_enrichr_ora(
            condition_results, condition_labels, outdir,
            cols=cols,
        )

    if cols is None:
        cols = DEFAULT_DESEQ2_COLS

    outdir = Path(outdir)
    ora_dir = outdir / "go_ora_gprofiler"
    ora_dir.mkdir(parents=True, exist_ok=True)

    print("\n-- Running g:Profiler ORA --")

    organism = GPROFILER_SPECIES.get(species.lower(), "hsapiens")
    print(f"  Organism: {organism}")

    gp = GProfiler(return_dataframe=True)

    source_to_category = {
        "GO:BP": "BP",
        "GO:MF": "MF",
        "GO:CC": "CC",
        "KEGG": "KEGG",
        "REAC": "Reactome",
    }
    sources = list(source_to_category.keys())

    go_results: dict[str, dict[str, pd.DataFrame]] = {}

    for cond_name, data in condition_results.items():
        cond_label = condition_labels.get(cond_name, cond_name)

        deg_df = data.get("deseq2_filtered", {}).get("all_genes", pd.DataFrame())
        if len(deg_df) == 0:
            print(f"  No DEGs for {cond_label}, skipping ORA")
            continue

        gene_col, gene_type = best_gene_key(deg_df, cols)
        if "Ensembl" in gene_type:
            name_col = cols.get("gene_name", "")
            if name_col and name_col in deg_df.columns:
                gene_col = name_col

        if "direction" not in deg_df.columns:
            print(f"  WARNING: no 'direction' column in {cond_label} DEGs, skipping")
            continue

        up_genes = deg_df.loc[deg_df["direction"] == "up", gene_col].dropna().astype(str).unique().tolist()
        down_genes = deg_df.loc[deg_df["direction"] == "down", gene_col].dropna().astype(str).unique().tolist()

        cond_go: dict[str, pd.DataFrame] = {}

        for direction, gene_list in [("up", up_genes), ("down", down_genes)]:
            n = len(gene_list)
            if n < 5:
                print(f"  {cond_label} ({direction}): only {n} genes -- too few for ORA")
                cond_go[direction] = pd.DataFrame()
                continue

            print(f"  Running g:Profiler ORA for {cond_label} ({direction}: {n} genes)...")

            try:
                result = gp.profile(
                    organism=organism,
                    query=gene_list,
                    sources=sources,
                    significance_threshold_method="g_SCS",
                )

                if result is None or len(result) == 0:
                    cond_go[direction] = pd.DataFrame()
                    continue

                if "p_value" in result.columns:
                    result = result[result["p_value"] < padj_threshold].copy()

                if len(result) == 0:
                    cond_go[direction] = pd.DataFrame()
                    continue

                rows = []
                for _, row in result.iterrows():
                    source = row.get("source", "")
                    category = source_to_category.get(source, source)
                    rows.append({
                        "Term": row.get("name", row.get("native", "")),
                        "Adjusted_P_value": row.get("p_value", 1.0),
                        "Overlap_count": int(row.get("intersection_size", 0)),
                        "Category": category,
                    })

                direction_df = pd.DataFrame(rows)

                top_rows = []
                for cat in direction_df["Category"].unique():
                    cat_df = direction_df[direction_df["Category"] == cat]
                    cat_df = cat_df.sort_values("Adjusted_P_value").head(top_n_per_db)
                    top_rows.append(cat_df)

                if top_rows:
                    direction_df = pd.concat(top_rows, ignore_index=True)
                    print(
                        f"    Found {len(direction_df)} enriched terms across "
                        f"{direction_df['Category'].nunique()} databases"
                    )
                else:
                    direction_df = pd.DataFrame()

                cond_go[direction] = direction_df

            except Exception as e:
                print(f"    ERROR running g:Profiler for {cond_label} ({direction}): {e}")
                cond_go[direction] = pd.DataFrame()

        go_results[cond_name] = cond_go

        for direction in ("up", "down"):
            df = cond_go.get(direction, pd.DataFrame())
            if len(df) > 0:
                fname = ora_dir / f"gprofiler_ora_{cond_name}_{direction}.csv"
                df.to_csv(fname, index=False)
                print(f"    Saved: {fname.name}")

    print(f"  g:Profiler ORA complete: {len(go_results)} conditions processed")
    return go_results


# ---------------------------------------------------------------------------
# Enrichr ORA (via gseapy)
# ---------------------------------------------------------------------------


def run_enrichr_ora(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    databases: list[str] | None = None,
    padj_threshold: float = 0.05,
    top_n_per_db: int = 10,
) -> dict[str, dict[str, pd.DataFrame]]:
    """Run Enrichr over-representation analysis for up/down DEGs.

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Output directory.
    cols : dict or None
        DESeq2 column name mapping.
    databases : list[str] or None
        Enrichr database names.
    padj_threshold : float
        Adjusted p-value threshold for significance.
    top_n_per_db : int
        Number of top terms per database.

    Returns
    -------
    dict
        go_results[cond_name] = {"up": DataFrame, "down": DataFrame}
    """
    try:
        import gseapy as gp
    except ImportError:
        print("  WARNING: gseapy not installed, skipping GO ORA")
        print("  Install with: pip install gseapy")
        return {}

    if cols is None:
        cols = DEFAULT_DESEQ2_COLS
    if databases is None:
        databases = DEFAULT_ORA_DATABASES

    outdir = Path(outdir)
    ora_dir = outdir / "go_ora_enrichr"
    ora_dir.mkdir(parents=True, exist_ok=True)

    print("\n-- Running GO / Pathway ORA (Enrichr) --")

    go_results: dict[str, dict[str, pd.DataFrame]] = {}

    for cond_name, data in condition_results.items():
        cond_label = condition_labels.get(cond_name, cond_name)

        deg_df = data.get("deseq2_filtered", {}).get("all_genes", pd.DataFrame())
        if len(deg_df) == 0:
            print(f"  No DEGs for {cond_label}, skipping ORA")
            continue

        gene_col, gene_type = best_gene_key(deg_df, cols)
        if "Ensembl" in gene_type:
            name_col = cols.get("gene_name", "")
            if name_col and name_col in deg_df.columns:
                gene_col = name_col

        if "direction" not in deg_df.columns:
            print(f"  WARNING: no 'direction' column in {cond_label} DEGs, skipping")
            continue

        up_genes = deg_df.loc[deg_df["direction"] == "up", gene_col].dropna().astype(str).unique().tolist()
        down_genes = deg_df.loc[deg_df["direction"] == "down", gene_col].dropna().astype(str).unique().tolist()

        cond_go: dict[str, pd.DataFrame] = {}

        for direction, gene_list in [("up", up_genes), ("down", down_genes)]:
            n = len(gene_list)
            if n < 5:
                print(f"  {cond_label} ({direction}): only {n} genes -- too few for ORA")
                cond_go[direction] = pd.DataFrame()
                continue

            print(f"  Running GO ORA for {cond_label} ({direction}: {n} genes)...")
            direction_rows = []

            for db in databases:
                category = DB_TO_CATEGORY.get(db, db)
                try:
                    enr = gp.enrich(
                        gene_list=gene_list,
                        gene_sets=db,
                        outdir=None,
                        no_plot=True,
                        verbose=False,
                    )
                    res_df = enr.results if hasattr(enr, "results") else enr.res2d
                    if res_df is None or len(res_df) == 0:
                        continue

                    if "Overlap" in res_df.columns:
                        res_df["Overlap_count"] = (
                            res_df["Overlap"].astype(str).str.split("/").str[0].astype(int)
                        )
                    else:
                        res_df["Overlap_count"] = 0

                    padj_col = None
                    for candidate in [
                        "Adjusted P-value", "Adjusted_P_value",
                        "FDR q-val", "padj", "fdr",
                    ]:
                        if candidate in res_df.columns:
                            padj_col = candidate
                            break
                    if padj_col is None:
                        for c in res_df.columns:
                            if "adjust" in c.lower() and "p" in c.lower():
                                padj_col = c
                                break
                    if padj_col is None:
                        continue

                    term_col = None
                    for candidate in ["Term", "term", "Pathway", "Gene_set"]:
                        if candidate in res_df.columns:
                            term_col = candidate
                            break
                    if term_col is None:
                        term_col = res_df.columns[0]

                    sig = res_df[res_df[padj_col].astype(float) < padj_threshold].copy()
                    if len(sig) == 0:
                        continue

                    sig = sig.sort_values(padj_col, ascending=True).head(top_n_per_db)

                    for _, row in sig.iterrows():
                        direction_rows.append({
                            "Term": str(row[term_col]),
                            "Adjusted_P_value": float(row[padj_col]),
                            "Overlap_count": int(row["Overlap_count"]),
                            "Category": category,
                        })

                except Exception as e:
                    print(f"    ERROR with {db}: {e}")
                    continue

            if direction_rows:
                result_df = pd.DataFrame(direction_rows)
                result_df = result_df.sort_values("Adjusted_P_value", ascending=True)
                cond_go[direction] = result_df

                xlsx_path = ora_dir / f"GO_ORA_{cond_name}_{direction}.xlsx"
                result_df.to_excel(xlsx_path, index=False)
                print(f"    Saved: {xlsx_path.name}")
            else:
                cond_go[direction] = pd.DataFrame()

        go_results[cond_name] = cond_go

    print(f"\n  GO ORA complete: {len(go_results)} conditions processed")
    return go_results


# ---------------------------------------------------------------------------
# Dual ORA runner
# ---------------------------------------------------------------------------


def run_dual_ora(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    species: str = "human",
    method: str = "both",
) -> tuple[dict, dict]:
    """Run ORA with configurable method: 'both', 'gprofiler', or 'enrichr'.

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Output directory.
    cols : dict or None
        DESeq2 column name mapping.
    species : str
        Species name.
    method : str
        'both' (run both side-by-side), 'gprofiler', or 'enrichr'.

    Returns
    -------
    tuple
        (enrichr_results, gprofiler_results) -- either may be empty dict.
    """
    enrichr_results: dict = {}
    gprofiler_results: dict = {}

    if method in ("both", "enrichr"):
        print("\n-- Running Enrichr ORA --")
        try:
            enrichr_results = run_enrichr_ora(
                condition_results, condition_labels, outdir, cols=cols
            )
        except Exception as e:
            print(f"  WARNING: Enrichr ORA failed: {e}")

    if method in ("both", "gprofiler"):
        print("\n-- Running g:Profiler ORA --")
        try:
            gprofiler_results = run_gprofiler_ora(
                condition_results, condition_labels, outdir,
                cols=cols, species=species,
            )
        except Exception as e:
            print(f"  WARNING: g:Profiler ORA failed: {e}")

    return enrichr_results, gprofiler_results


# ---------------------------------------------------------------------------
# ORA visualization functions
# ---------------------------------------------------------------------------


def go_enrichment_combined_plot(
    go_results: dict,
    condition_labels: dict[str, str],
    outdir: str | Path,
    fig_format: str = "png",
    fig_dpi: int = 300,
    filename_suffix: str = "",
) -> None:
    """Create a two-panel dot plot (up / down) for each condition.

    Parameters
    ----------
    go_results : dict
        Output of run_enrichr_ora() or run_gprofiler_ora().
    condition_labels : dict
        Maps condition name -> human-readable label.
    outdir : str or Path
        Output directory for figures.
    fig_format : str
        Image format (e.g. 'png', 'svg').
    fig_dpi : int
        Resolution in DPI.
    filename_suffix : str
        Optional suffix appended to output filenames (e.g. '_enrichr').
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for cond_name, directions in go_results.items():
        cond_label = condition_labels.get(cond_name, cond_name)
        up_df = directions.get("up", pd.DataFrame())
        down_df = directions.get("down", pd.DataFrame())

        has_up = len(up_df) > 0
        has_down = len(down_df) > 0

        if not has_up and not has_down:
            print(f"  No ORA results for {cond_label}, skipping plot")
            continue

        fig, axes = plt.subplots(1, 2, figsize=(20, 12))
        method_label = filename_suffix.replace("_", " ").strip().title()
        title_extra = f" ({method_label})" if method_label else ""
        fig.suptitle(
            f"GO/Pathway Enrichment \u2014 {cond_label}{title_extra}",
            fontsize=14, fontweight="bold", y=0.98,
        )

        for ax, df, title in [
            (axes[0], up_df, "Upregulated"),
            (axes[1], down_df, "Downregulated"),
        ]:
            if len(df) == 0:
                ax.set_title(title, fontsize=12, fontweight="bold")
                ax.text(
                    0.5, 0.5, "No significant enrichment",
                    ha="center", va="center", fontsize=11,
                    transform=ax.transAxes, color="#666666",
                )
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                continue

            # Top 10 per category, sorted by significance
            plot_rows = []
            for cat in ["BP", "CC", "MF", "KEGG", "Reactome"]:
                cat_df = df[df["Category"] == cat].copy()
                if len(cat_df) == 0:
                    continue
                cat_df = cat_df.sort_values("Adjusted_P_value").head(10)
                plot_rows.append(cat_df)

            if not plot_rows:
                ax.set_title(title, fontsize=12, fontweight="bold")
                ax.text(
                    0.5, 0.5, "No significant enrichment",
                    ha="center", va="center", fontsize=11,
                    transform=ax.transAxes, color="#666666",
                )
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                continue

            plot_df = pd.concat(plot_rows, ignore_index=True)
            plot_df = plot_df.sort_values("Adjusted_P_value", ascending=False)

            plot_df["neg_log10_padj"] = -np.log10(
                plot_df["Adjusted_P_value"].clip(lower=1e-300)
            )

            plot_df["label"] = plot_df["Category"] + ": " + plot_df["Term"]
            plot_df["label"] = plot_df["label"].apply(
                lambda x: x[:80] + "..." if len(x) > 83 else x
            )

            colors = [CATEGORY_COLORS.get(c, "#999999") for c in plot_df["Category"]]

            counts = plot_df["Overlap_count"].values
            min_size, max_size = 30, 300
            if counts.max() > counts.min():
                sizes = min_size + (counts - counts.min()) / (
                    counts.max() - counts.min()
                ) * (max_size - min_size)
            else:
                sizes = np.full(len(counts), (min_size + max_size) / 2)

            ax.scatter(
                plot_df["neg_log10_padj"], range(len(plot_df)),
                c=colors, s=sizes, edgecolors="white", linewidths=0.5, zorder=3,
            )
            ax.set_yticks(range(len(plot_df)))
            ax.set_yticklabels(plot_df["label"].tolist(), fontsize=9)
            ax.set_xlabel("-log$_{10}$(adjusted p-value)", fontsize=10)
            ax.set_title(title, fontsize=12, fontweight="bold")
            ax.grid(axis="x", alpha=0.3, linestyle="--")
            ax.set_axisbelow(True)

        # Legend for category colors
        cat_handles = [
            mlines.Line2D(
                [], [], marker="o", color="w",
                markerfacecolor=CATEGORY_COLORS[cat],
                markeredgecolor="white", markersize=9, label=cat,
            )
            for cat in ["BP", "CC", "MF", "KEGG", "Reactome"]
        ]

        # Legend for dot sizes
        all_counts: list = []
        for d in [up_df, down_df]:
            if len(d) > 0 and "Overlap_count" in d.columns:
                all_counts.extend(d["Overlap_count"].tolist())

        if all_counts:
            c_min, c_max = int(min(all_counts)), int(max(all_counts))
            if c_min == c_max:
                size_vals = [c_min]
            else:
                c_mid = (c_min + c_max) // 2
                size_vals = [c_min, c_mid, c_max]

            size_handles = []
            for sv in size_vals:
                if c_max > c_min:
                    s = min_size + (sv - c_min) / (c_max - c_min) * (max_size - min_size)
                else:
                    s = (min_size + max_size) / 2
                size_handles.append(
                    mlines.Line2D(
                        [], [], marker="o", color="w",
                        markerfacecolor="#AAAAAA",
                        markeredgecolor="gray",
                        markersize=np.sqrt(s) * 0.7,
                        label=f"{sv} genes",
                    )
                )
            all_handles = cat_handles + size_handles
        else:
            all_handles = cat_handles

        fig.legend(
            handles=all_handles, loc="lower center",
            ncol=len(all_handles), fontsize=9,
            frameon=True, fancybox=True, shadow=False,
            bbox_to_anchor=(0.5, 0.01),
        )

        plt.tight_layout(rect=[0, 0.06, 1, 0.95])

        fig_path = outdir / f"go_enrichment_combined_{cond_name}{filename_suffix}.{fig_format}"
        fig.savefig(fig_path, dpi=fig_dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {fig_path.name}")


def export_go_prism(
    go_results: dict,
    condition_labels: dict[str, str],
    prism_dir: str | Path,
    filename_suffix: str = "",
) -> None:
    """Export GO ORA results to GraphPad Prism .pzfx files.

    Parameters
    ----------
    go_results : dict
        Output of run_enrichr_ora() or run_gprofiler_ora().
    condition_labels : dict
        Maps condition name -> human-readable label.
    prism_dir : str or Path
        Directory for Prism files.
    filename_suffix : str
        Optional suffix for filenames.
    """
    prism_dir = Path(prism_dir)
    prism_dir.mkdir(parents=True, exist_ok=True)

    def _is_numeric_str(s):
        try:
            float(s)
            return True
        except (ValueError, TypeError):
            return False

    def _create_prism_xml():
        root = ET.Element("GraphPadPrismFile")
        root.set("PrismXMLVersion", "5.00")
        created = ET.SubElement(root, "Created")
        orig = ET.SubElement(created, "OriginalVersion")
        orig.set("CreatedByProgram", "GraphPad Prism")
        orig.set("CreatedByVersion", "6.0f.254")
        orig.set("Login", "")
        orig.set("DateTime", datetime.now().strftime("%Y-%m-%dT%H:%M:%S+00:00"))
        info_seq = ET.SubElement(root, "InfoSequence")
        ref = ET.SubElement(info_seq, "Ref")
        ref.set("ID", "Info0")
        ref.set("Selected", "1")
        info = ET.SubElement(root, "Info")
        info.set("ID", "Info0")
        ET.SubElement(info, "Title").text = "Project info 1"
        ET.SubElement(info, "Notes")
        table_seq = ET.SubElement(root, "TableSequence")
        table_seq.set("Selected", "1")
        root.set("_table_count", "0")
        return root

    def _add_table(root, table_name, columns_data):
        table_idx = int(root.get("_table_count", "0"))
        table_id = f"Table{table_idx}"
        root.set("_table_count", str(table_idx + 1))
        table_seq = root.find("TableSequence")
        ref = ET.SubElement(table_seq, "Ref")
        ref.set("ID", table_id)
        if table_idx == 0:
            ref.set("Selected", "1")
        table = ET.SubElement(root, "Table")
        table.set("ID", table_id)
        table.set("XFormat", "none")
        table.set("TableType", "OneWay")
        table.set("EVFormat", "AsteriskAfterNumber")
        ET.SubElement(table, "Title").text = table_name
        # Auto-detect row titles
        row_titles_col = None
        if columns_data:
            _, first_vals = columns_data[0]
            if all(isinstance(v, str) and not _is_numeric_str(v) for v in first_vals):
                row_titles_col = 0
        if row_titles_col is not None:
            _, rt_vals = columns_data[row_titles_col]
            rt_elem = ET.SubElement(table, "RowTitlesColumn")
            rt_elem.set("Width", "89")
            subcol = ET.SubElement(rt_elem, "Subcolumn")
            for val in rt_vals:
                ET.SubElement(subcol, "d").text = str(val)
            data_columns = [c for i, c in enumerate(columns_data) if i != row_titles_col]
        else:
            data_columns = list(columns_data)
        for col_name, values in data_columns:
            ycol = ET.SubElement(table, "YColumn")
            ycol.set("Width", "89")
            ycol.set("Decimals", "6")
            ycol.set("Subcolumns", "1")
            ET.SubElement(ycol, "Title").text = str(col_name)
            subcol = ET.SubElement(ycol, "Subcolumn")
            for val in values:
                d = ET.SubElement(subcol, "d")
                if val is None or (isinstance(val, float) and str(val) == "nan"):
                    d.text = ""
                else:
                    d.text = str(val)
        return table

    def _save_prism_xml(root, filepath):
        if "_table_count" in root.attrib:
            del root.attrib["_table_count"]
        xml_str = ET.tostring(root, encoding="unicode", xml_declaration=False)
        dom = minidom.parseString(xml_str)
        pretty_xml = dom.toprettyxml(indent="  ", encoding=None)
        if pretty_xml.startswith("<?xml"):
            first_newline = pretty_xml.index("\n")
            pretty_xml = pretty_xml[first_newline + 1:]
        with open(filepath, "w", encoding="utf-8") as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write(pretty_xml)

    # Generate one .pzfx per condition
    for cond_name, directions in go_results.items():
        cond_label = condition_labels.get(cond_name, cond_name)
        root = _create_prism_xml()
        tables_added = 0

        for direction in ("up", "down"):
            df = directions.get(direction, pd.DataFrame())
            if len(df) == 0:
                continue

            terms = df["Term"].tolist()
            categories = df["Category"].tolist()
            neg_log_pvals = (
                -np.log10(df["Adjusted_P_value"].clip(lower=1e-300))
            ).round(4).tolist()
            counts = df["Overlap_count"].astype(int).tolist()

            table_name = f"GO_ORA_{direction}_{cond_label}"
            _add_table(root, table_name, [
                ("Term", terms),
                ("Category", categories),
                ("-log10(padj)", neg_log_pvals),
                ("GeneCount", counts),
            ])
            tables_added += 1

        if tables_added > 0:
            pzfx_path = prism_dir / f"GO_ORA_{cond_name}{filename_suffix}.pzfx"
            _save_prism_xml(root, pzfx_path)
            print(f"  Saved Prism: {pzfx_path.name}")
        else:
            print(f"  No ORA data for {cond_label}, skipping Prism export")
