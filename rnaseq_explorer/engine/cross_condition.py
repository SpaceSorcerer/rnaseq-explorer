"""Cross-condition comparison analyses.

Provides Venn diagrams, UpSet plots, direction concordance heatmaps,
pairwise scatter plots, and multi-condition overlap analyses for both
DESeq2 (gene-level) and rMATS (event-level) results.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from rnaseq_explorer.engine.deseq2 import (
    DEFAULT_DESEQ2_COLS,
    best_gene_key,
    extract_gene_sets,
)
from rnaseq_explorer.engine.rmats import (
    DEFAULT_RMATS_COLS,
    RMATS_EVENT_TYPES,
    make_event_key,
)
import matplotlib.gridspec as gridspec

from rnaseq_explorer.viz.theme import (
    COLOR_DOWN,
    COLOR_NS,
    COLOR_UP,
    EVENT_COLORS,
    grid_dims,
    style_venn,
)

try:
    from matplotlib_venn import venn2, venn3
    _VENN_AVAILABLE = True
except ImportError:
    _VENN_AVAILABLE = False

try:
    from upsetplot import UpSet, from_memberships
    _UPSET_AVAILABLE = True
except ImportError:
    _UPSET_AVAILABLE = False

try:
    from scipy.stats import pearsonr, fisher_exact  # noqa: F401
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Venn diagram helpers
# ---------------------------------------------------------------------------


def compute_venn_data(
    gene_sets: dict[str, dict[str, set[str]]],
    condition_labels: dict[str, str],
) -> dict[str, dict]:
    """Compute Venn diagram intersection data for DESeq2 gene sets.

    Parameters
    ----------
    gene_sets : dict
        Output of extract_gene_sets().
    condition_labels : dict
        Maps condition name -> label.

    Returns
    -------
    dict
        Venn data for each subset type ('all', 'up', 'down').
    """
    names = list(gene_sets.keys())
    result: dict[str, dict] = {}

    for subset_key in ("all", "up", "down"):
        sets = {condition_labels.get(n, n): gene_sets[n].get(subset_key, set()) for n in names}
        result[subset_key] = sets

    return result


def deseq2_venn_diagrams(
    gene_sets: dict[str, dict[str, set[str]]],
    condition_labels: dict[str, str],
    outdir: str | Path,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate 3-way Venn diagrams for DESeq2 gene sets.

    Parameters
    ----------
    gene_sets : dict
        Output of extract_gene_sets().
    condition_labels : dict
        Maps condition name -> label.
    outdir : str or Path
        Directory to save plots.
    fig_format : str
        Figure format.
    fig_dpi : int
        Figure DPI.
    """
    if not _VENN_AVAILABLE:
        print("  Skipping Venn diagrams (matplotlib-venn not installed)")
        return

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    names = list(gene_sets.keys())
    if len(names) < 2:
        print("  Need at least 2 conditions for Venn diagrams")
        return

    for subset_key, title_suffix, fname_suffix in [
        ("all", "All Significant Genes", "all_sig_genes"),
        ("up", "Upregulated", "upregulated"),
        ("down", "Downregulated", "downregulated"),
    ]:
        sets_list = [gene_sets[n].get(subset_key, set()) for n in names]
        labels_list = [condition_labels.get(n, n) for n in names]

        fig, ax = plt.subplots(figsize=(8, 7))

        if len(names) == 2:
            v = venn2(sets_list[:2], set_labels=labels_list[:2], ax=ax)
            style_venn(v, 2)
        elif len(names) >= 3:
            v = venn3(sets_list[:3], set_labels=labels_list[:3], ax=ax)
            style_venn(v, 3)
        else:
            plt.close(fig)
            continue

        ax.set_title(f"Venn Diagram -- {title_suffix}")
        fname = outdir / f"venn_{fname_suffix}.{fig_format}"
        fig.savefig(fname, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {fname.name}")


# ---------------------------------------------------------------------------
# UpSet plots
# ---------------------------------------------------------------------------


def compute_upset_data(
    gene_sets: dict[str, dict[str, set[str]]],
    condition_labels: dict[str, str],
    subset_key: str = "all",
) -> Optional[pd.DataFrame]:
    """Compute UpSet plot membership data.

    Parameters
    ----------
    gene_sets : dict
        Output of extract_gene_sets().
    condition_labels : dict
        Maps condition name -> label.
    subset_key : str
        'all', 'up', or 'down'.

    Returns
    -------
    pd.DataFrame or None
        Membership data suitable for upsetplot.from_memberships().
    """
    if not _UPSET_AVAILABLE:
        return None

    names = list(gene_sets.keys())
    memberships = []
    for gene in set().union(*(gene_sets[n].get(subset_key, set()) for n in names)):
        member_of = [
            condition_labels.get(n, n)
            for n in names
            if gene in gene_sets[n].get(subset_key, set())
        ]
        memberships.append(member_of)

    if not memberships:
        return None

    return from_memberships(memberships)


def deseq2_upset_plot(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate UpSet plots for DESeq2 gene overlap.

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> label.
    outdir : str or Path
        Directory to save plots.
    cols : dict or None
        DESeq2 column name mapping.
    fig_format : str
        Figure format.
    fig_dpi : int
        Figure DPI.
    """
    if not _UPSET_AVAILABLE:
        print("  Skipping UpSet plots (upsetplot not installed)")
        return

    if cols is None:
        cols = DEFAULT_DESEQ2_COLS

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_sets = extract_gene_sets(condition_results, cols)

    for subset_key, fname_suffix in [
        ("all", "all_sig"),
        ("up", "up"),
        ("down", "down"),
    ]:
        data = compute_upset_data(gene_sets, condition_labels, subset_key)
        if data is None or len(data) == 0:
            continue

        try:
            upset = UpSet(data, show_counts=True, sort_by="cardinality")
            fig = plt.figure(figsize=(12, 6))
            upset.plot(fig=fig)
            fig.suptitle(f"UpSet -- DESeq2 {subset_key.title()} Genes", fontsize=12)
            fname = outdir / f"deseq2_upset_{fname_suffix}.{fig_format}"
            fig.savefig(fname, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {fname.name}")
        except Exception as e:
            print(f"  [WARNING] UpSet plot ({subset_key}) failed: {e}")


# ---------------------------------------------------------------------------
# Direction concordance
# ---------------------------------------------------------------------------


def compute_concordance(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    cols: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Compute pairwise direction concordance matrix for DESeq2 results.

    For each pair of conditions, counts shared DEGs that are up in both,
    down in both, or discordant (one up, one down).

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> label.
    cols : dict or None
        DESeq2 column name mapping.

    Returns
    -------
    pd.DataFrame
        Concordance matrix (n_conditions x n_conditions).
    """
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS

    names = list(condition_results.keys())
    labels = [condition_labels.get(n, n) for n in names]

    gene_directions: dict[str, dict[str, str]] = {}
    for cond_name, data in condition_results.items():
        deg_df = data.get("deseq2_filtered", {}).get("all_genes", pd.DataFrame())
        if len(deg_df) == 0:
            gene_directions[cond_name] = {}
            continue
        key_col, _ = best_gene_key(deg_df, cols)
        if key_col not in deg_df.columns or "direction" not in deg_df.columns:
            gene_directions[cond_name] = {}
            continue
        directions = {}
        for _, row in deg_df.iterrows():
            gene = str(row[key_col])
            directions[gene] = row["direction"]
        gene_directions[cond_name] = directions

    n = len(names)
    concordance = np.zeros((n, n))

    for i, j in combinations(range(n), 2):
        dir_i = gene_directions[names[i]]
        dir_j = gene_directions[names[j]]
        shared = set(dir_i.keys()) & set(dir_j.keys())
        if not shared:
            continue
        concordant = sum(1 for g in shared if dir_i[g] == dir_j[g])
        total = len(shared)
        rate = concordant / total if total > 0 else 0
        concordance[i, j] = rate
        concordance[j, i] = rate

    for i in range(n):
        concordance[i, i] = 1.0

    return pd.DataFrame(concordance, index=labels, columns=labels)


def compute_direction_heatmap(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> pd.DataFrame:
    """Generate and save direction concordance heatmap.

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> label.
    outdir : str or Path
        Directory to save plots.
    cols : dict or None
        DESeq2 column name mapping.
    fig_format : str
        Figure format.
    fig_dpi : int
        Figure DPI.

    Returns
    -------
    pd.DataFrame
        Concordance matrix.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    concordance_df = compute_concordance(condition_results, condition_labels, cols)

    if len(concordance_df) > 1:
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(
            concordance_df, annot=True, fmt=".2f", cmap="RdYlGn",
            vmin=0, vmax=1, ax=ax, square=True,
        )
        ax.set_title("DESeq2 Direction Concordance (shared DEGs)")
        fname = outdir / f"direction_concordance_heatmap.{fig_format}"
        fig.savefig(fname, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {fname.name}")

    return concordance_df


# ---------------------------------------------------------------------------
# Log2FC heatmap
# ---------------------------------------------------------------------------


def deseq2_log2fc_heatmap(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    max_genes: int = 100,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> pd.DataFrame:
    """Generate a clustered log2FC heatmap for DEGs across conditions.

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> label.
    outdir : str or Path
        Directory to save plots.
    cols : dict or None
        DESeq2 column name mapping.
    max_genes : int
        Maximum genes to show.
    fig_format : str
        Figure format.
    fig_dpi : int
        Figure DPI.

    Returns
    -------
    pd.DataFrame
        Log2FC matrix.
    """
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fc_col = cols["log2fc"]
    names = list(condition_results.keys())

    # Build matrix
    fc_dfs = []
    for cond_name in names:
        raw = condition_results[cond_name]["deseq2_raw"]
        key_col, _ = best_gene_key(raw, cols)
        if key_col not in raw.columns or fc_col not in raw.columns:
            continue

        # Get filtered DEGs for this condition
        filt = condition_results[cond_name]["deseq2_filtered"].get("all_genes", pd.DataFrame())
        if len(filt) == 0:
            continue
        deg_genes = set(filt[key_col].dropna().astype(str).unique()) if key_col in filt.columns else set()

        sub = raw[[key_col, fc_col]].dropna(subset=[key_col]).copy()
        sub[key_col] = sub[key_col].astype(str)
        sub = sub[sub[key_col].isin(deg_genes)]
        sub = sub.drop_duplicates(subset=[key_col]).set_index(key_col)
        fc_series = sub[fc_col]
        fc_series.name = condition_labels[cond_name]
        fc_dfs.append(fc_series)

    if not fc_dfs:
        return pd.DataFrame()

    matrix = pd.concat(fc_dfs, axis=1).dropna()
    if len(matrix) == 0:
        return pd.DataFrame()

    if len(matrix) > max_genes:
        matrix["max_abs_fc"] = matrix.abs().max(axis=1)
        matrix = matrix.nlargest(max_genes, "max_abs_fc").drop(columns="max_abs_fc")

    try:
        from rnaseq_explorer.viz.theme import diverging_cmap

        g = sns.clustermap(
            matrix, cmap=diverging_cmap(), center=0,
            figsize=(8, max(6, len(matrix) * 0.15)),
            row_cluster=True, col_cluster=False,
            yticklabels=(len(matrix) <= 60),
            linewidths=0.3, linecolor="white",
        )
        g.fig.suptitle(f"log2FC Heatmap (top {len(matrix)} DEGs)", y=1.02)
        fname = outdir / f"log2fc_heatmap.{fig_format}"
        g.savefig(fname, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
        plt.close(g.fig)
        print(f"  Saved: {fname.name}")
    except Exception as e:
        print(f"  [WARNING] log2FC heatmap failed: {e}")

    return matrix


# ---------------------------------------------------------------------------
# Pairwise log2FC scatter
# ---------------------------------------------------------------------------


def pairwise_log2fc_scatter(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate pairwise log2FC scatter plots for shared genes.

    Parameters
    ----------
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> label.
    outdir : str or Path
        Directory to save plots.
    cols : dict or None
        DESeq2 column name mapping.
    fig_format : str
        Figure format.
    fig_dpi : int
        Figure DPI.
    """
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fc_col = cols["log2fc"]
    names = list(condition_results.keys())
    pairs = list(combinations(names, 2))

    if not pairs:
        return

    nrows, ncols_grid = grid_dims(len(pairs))
    fig, axes = plt.subplots(nrows, ncols_grid, figsize=(ncols_grid * 6, nrows * 5), squeeze=False)
    axes_flat = [axes[r][c] for r in range(nrows) for c in range(ncols_grid)]

    for idx, (cA, cB) in enumerate(pairs):
        ax = axes_flat[idx]
        rawA = condition_results[cA]["deseq2_raw"]
        rawB = condition_results[cB]["deseq2_raw"]

        key_col, _ = best_gene_key(rawA, cols)
        if key_col not in rawA.columns or key_col not in rawB.columns:
            ax.set_visible(False)
            continue

        dfA = rawA[[key_col, fc_col]].dropna().drop_duplicates(subset=[key_col]).set_index(key_col)
        dfB = rawB[[key_col, fc_col]].dropna().drop_duplicates(subset=[key_col]).set_index(key_col)

        merged = dfA.join(dfB, lsuffix="_A", lsuffix_B="_B" if False else "", rsuffix="_B")
        if len(merged) == 0:
            ax.set_visible(False)
            continue

        colA = f"{fc_col}_A"
        colB = f"{fc_col}_B"

        ax.scatter(merged[colA], merged[colB], s=4, alpha=0.3, c=COLOR_NS, rasterized=True)
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.axvline(0, color="grey", linewidth=0.5)

        lbl_A = condition_labels.get(cA, cA)
        lbl_B = condition_labels.get(cB, cB)
        ax.set_xlabel(f"log2FC ({lbl_A})")
        ax.set_ylabel(f"log2FC ({lbl_B})")
        ax.set_title(f"{lbl_A} vs {lbl_B}", fontsize=10)

        if _SCIPY_AVAILABLE and len(merged) > 10:
            r, p = pearsonr(merged[colA], merged[colB])
            ax.text(0.05, 0.95, f"r = {r:.3f}", transform=ax.transAxes,
                    fontsize=9, va="top", fontweight="bold")

    for idx in range(len(pairs), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle("Pairwise log2FC Scatter", fontsize=12, y=1.01)
    plt.tight_layout()
    fname = outdir / f"pairwise_log2fc_scatter.{fig_format}"
    fig.savefig(fname, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fname.name}")


# ---------------------------------------------------------------------------
# DESeq2 DE counts chart
# ---------------------------------------------------------------------------


def deseq2_de_counts_chart(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Grouped bar chart: Up / Down / Total DE gene counts per condition."""
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    names = list(condition_results.keys())
    labels = [condition_labels[n] for n in names]
    ups, downs, totals = [], [], []
    for name in names:
        filt = condition_results[name]["deseq2_filtered"].get("all_genes", pd.DataFrame())
        if "direction" in filt.columns:
            u = int((filt["direction"] == "up").sum())
            d = int((filt["direction"] == "down").sum())
        else:
            u = d = 0
        ups.append(u)
        downs.append(d)
        totals.append(u + d)

    x = np.arange(len(names))
    w = 0.25
    fig, ax = plt.subplots(figsize=(max(8, len(names) * 2.5), 6))
    b1 = ax.bar(x - w, ups, w, label="Up", color=COLOR_UP)
    b2 = ax.bar(x, downs, w, label="Down", color=COLOR_DOWN)
    b3 = ax.bar(x + w, totals, w, label="Total", color="#888888")

    for bars in (b1, b2, b3):
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    h + 0.5,
                    str(int(h)),
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    fontweight="bold",
                )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Number of DE Genes")
    ax.set_title("DESeq2 \u2014 DE Gene Counts Overview Across Conditions")
    ax.legend()

    outpath = outdir / f"deseq2_de_counts_overview.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Pairwise DEG Venn diagrams
# ---------------------------------------------------------------------------


def pairwise_deg_venns(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Pairwise 5-panel Venn diagrams for differentially expressed genes."""
    if not _VENN_AVAILABLE:
        print("  Skipping pairwise DEG Venns (matplotlib-venn not installed)")
        return
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = list(condition_results.keys())

    for name_a, name_b in combinations(names, 2):
        label_a = condition_labels[name_a]
        label_b = condition_labels[name_b]

        filt_a = condition_results[name_a]["deseq2_filtered"]["all_genes"]
        filt_b = condition_results[name_b]["deseq2_filtered"]["all_genes"]

        key_col_a, _ = best_gene_key(filt_a, cols)
        key_col_b, _ = best_gene_key(filt_b, cols)
        if not key_col_a or key_col_a not in filt_a.columns:
            key_col_a = cols["gene_name"]
        if not key_col_b or key_col_b not in filt_b.columns:
            key_col_b = cols["gene_name"]

        all_a = set(filt_a[key_col_a].dropna().unique())
        all_b = set(filt_b[key_col_b].dropna().unique())
        up_a = set(filt_a.loc[filt_a["direction"] == "up", key_col_a].dropna().unique())
        up_b = set(filt_b.loc[filt_b["direction"] == "up", key_col_b].dropna().unique())
        down_a = set(filt_a.loc[filt_a["direction"] == "down", key_col_a].dropna().unique())
        down_b = set(filt_b.loc[filt_b["direction"] == "down", key_col_b].dropna().unique())

        panels = [
            ("Significant in Either Condition", all_a, all_b),
            ("Upregulated", up_a, up_b),
            ("Downregulated", down_a, down_b),
            (f"Up in {label_a} / Down in {label_b}", up_a, down_b),
            (f"Down in {label_a} / Up in {label_b}", down_a, up_b),
        ]

        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        axes_flat = axes.flatten()

        for idx, (panel_title, set_a, set_b) in enumerate(panels):
            ax = axes_flat[idx]
            v = venn2([set_a, set_b], set_labels=(label_a, label_b), ax=ax)
            style_venn(v, 2)
            ax.set_title(
                f"{panel_title}\n({label_a}: {len(set_a)} | {label_b}: {len(set_b)})",
                fontsize=11,
                fontweight="bold",
            )

        axes_flat[5].set_visible(False)
        fig.suptitle(
            f"Differentially Expressed Genes \u2014 {label_a} vs {label_b}",
            fontsize=14,
            fontweight="bold",
        )
        plt.tight_layout()
        outpath = outdir / f"venn_deg_{name_a}_vs_{name_b}.{fig_format}"
        fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {outpath.name}")


# ---------------------------------------------------------------------------
# rMATS cross-condition Venn diagrams
# ---------------------------------------------------------------------------


def rmats_cross_condition_venn(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
    match_by: str = "event",
) -> None:
    """Venn diagrams for significant splicing events across conditions."""
    if not _VENN_AVAILABLE:
        print("  Skipping rMATS Venn diagrams (matplotlib-venn not installed)")
        return

    cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)
    names = list(condition_results.keys())
    labels = [condition_labels[n] for n in names]
    gene_col = cols["gene_name"]
    id_col = cols["gene_id"]

    is_gene = match_by == "gene"
    level_label = "gene-level" if is_gene else "coordinate-level"
    fname_suffix = "_genelevel" if is_gene else ""

    if not is_gene:
        outdir = outdir / "event_level"
        outdir.mkdir(exist_ok=True, parents=True)

    # Combined Venn (all event types)
    event_sets = []
    for name in names:
        all_sig = set()
        for et, filt_df in condition_results[name]["rmats_filtered"].items():
            if is_gene:
                _match_col = id_col if id_col in filt_df.columns else gene_col
                if _match_col in filt_df.columns:
                    all_sig.update(filt_df[_match_col].dropna().unique())
            else:
                keys = make_event_key(filt_df, et)
                all_sig.update(keys[keys != ""].unique())
        event_sets.append(all_sig)

    fig, ax = plt.subplots(figsize=(8, 8))
    if len(event_sets) == 2:
        v = venn2(event_sets, set_labels=labels, ax=ax)
        style_venn(v, 2)
    elif len(event_sets) >= 3:
        v = venn3(event_sets[:3], set_labels=labels[:3], ax=ax)
        style_venn(v, 3)
    ax.set_title(
        f"rMATS \u2014 Significant Splicing Events ({level_label})",
        fontsize=13,
        fontweight="bold",
    )
    outpath = outdir / f"venn_rmats_events{fname_suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")

    # Per-event-type Venns
    n_types = len(RMATS_EVENT_TYPES)
    fig, axes = plt.subplots(1, n_types, figsize=(5 * n_types, 5))
    if n_types == 1:
        axes = [axes]

    for idx, et in enumerate(RMATS_EVENT_TYPES):
        ax = axes[idx]
        sets = []
        for name in names:
            if et in condition_results[name]["rmats_filtered"]:
                filt_df = condition_results[name]["rmats_filtered"][et]
                if is_gene:
                    _match_col = id_col if id_col in filt_df.columns else gene_col
                    items = (
                        set(filt_df[_match_col].dropna().unique())
                        if _match_col in filt_df.columns
                        else set()
                    )
                else:
                    keys = make_event_key(filt_df, et)
                    items = set(keys[keys != ""].unique())
            else:
                items = set()
            sets.append(items)

        if len(sets) == 2:
            v = venn2(sets, set_labels=labels, ax=ax)
            style_venn(v, 2)
        elif len(sets) >= 3:
            v = venn3(sets[:3], set_labels=labels[:3], ax=ax)
            style_venn(v, 3)
        ax.set_title(f"{et}", fontsize=12, fontweight="bold")

    fig.suptitle(
        f"rMATS \u2014 Splicing Event Overlap by Type ({level_label})",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout()
    outpath = outdir / f"venn_rmats_events_by_type{fname_suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# rMATS direction concordance
# ---------------------------------------------------------------------------


def rmats_direction_concordance(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> pd.DataFrame:
    """Compare dPSI direction for shared significant splicing genes."""
    cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    names = list(condition_results.keys())
    gene_col = cols["gene_name"]
    dpsi_col = cols["inclevel_diff"]

    concordance_rows = []

    for et in RMATS_EVENT_TYPES:
        dfs: dict[str, pd.Series] = {}
        for name in names:
            if et in condition_results[name]["rmats_filtered"]:
                df = condition_results[name]["rmats_filtered"][et]
                if len(df) > 0:
                    gene_dpsi = df.groupby(gene_col)[dpsi_col].mean()
                    dfs[name] = gene_dpsi

        if len(dfs) < 2:
            continue

        shared_genes = set.intersection(*[set(d.index) for d in dfs.values()])

        for gene in shared_genes:
            signs = {name: np.sign(dfs[name][gene]) for name in dfs}
            all_same = len(set(signs.values())) == 1
            row: dict = {"gene": gene, "event_type": et, "concordant": all_same}
            for n in dfs:
                row[f"dPSI_{condition_labels[n]}"] = dfs[n][gene]
            concordance_rows.append(row)

    if not concordance_rows:
        print("  No shared significant splicing events for concordance analysis")
        return pd.DataFrame()

    conc_df = pd.DataFrame(concordance_rows)

    event_order = [et for et in RMATS_EVENT_TYPES if et in conc_df["event_type"].values]
    concordant_counts = []
    discordant_counts = []
    for et in event_order:
        et_data = conc_df[conc_df["event_type"] == et]
        concordant_counts.append(int(et_data["concordant"].sum()))
        discordant_counts.append(int((~et_data["concordant"]).sum()))

    fig, ax = plt.subplots(figsize=(9, 5))
    x = range(len(event_order))
    width = 0.35
    ax.bar(
        [xi - width / 2 for xi in x],
        concordant_counts,
        width,
        label="Concordant (same dPSI sign)",
        color="#0072B2",
    )
    ax.bar(
        [xi + width / 2 for xi in x],
        discordant_counts,
        width,
        label="Discordant (opposite dPSI sign)",
        color="#E69F00",
    )
    ax.set_xticks(list(x))
    ax.set_xticklabels(event_order)
    ax.set_ylabel("Number of Shared Genes")
    ax.set_title("rMATS Direction Concordance - Shared Significant Splicing Genes")
    ax.legend()
    for xi, cc, dc in zip(x, concordant_counts, discordant_counts):
        ax.text(xi - width / 2, cc + 0.5, str(cc), ha="center", fontsize=9, fontweight="bold")
        ax.text(xi + width / 2, dc + 0.5, str(dc), ha="center", fontsize=9, fontweight="bold")

    outpath = outdir / f"rmats_direction_concordance.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")

    return conc_df


# ---------------------------------------------------------------------------
# Pairwise splicing Venns
# ---------------------------------------------------------------------------


def pairwise_splicing_venns(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    dpsi_cutoff: float = 0.1,
    fig_format: str = "png",
    fig_dpi: int = 300,
    match_by: str = "event",
) -> None:
    """Pairwise 5-panel Venn diagrams for splicing events."""
    if not _VENN_AVAILABLE:
        print("  Skipping pairwise splicing Venns (matplotlib-venn not installed)")
        return

    cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)
    names = list(condition_results.keys())
    dpsi_col = cols["inclevel_diff"]
    gene_col = cols["gene_name"]
    id_col = cols["gene_id"]

    is_gene = match_by == "gene"
    level_label = "gene-level" if is_gene else "event-level"
    fname_suffix = "_genelevel" if is_gene else ""

    if not is_gene:
        outdir = outdir / "event_level"
        outdir.mkdir(exist_ok=True, parents=True)

    for name_a, name_b in combinations(names, 2):
        label_a = condition_labels[name_a]
        label_b = condition_labels[name_b]

        for et in RMATS_EVENT_TYPES:
            df_a = condition_results[name_a]["rmats_filtered"].get(et)
            df_b = condition_results[name_b]["rmats_filtered"].get(et)

            if df_a is None or len(df_a) == 0 or df_b is None or len(df_b) == 0:
                continue

            if is_gene:
                _match_col = (
                    id_col
                    if id_col in df_a.columns and id_col in df_b.columns
                    else gene_col
                )
                if _match_col not in df_a.columns or _match_col not in df_b.columns:
                    continue

                events_all_a = set(df_a[_match_col].dropna().unique())
                events_all_b = set(df_b[_match_col].dropna().unique())
                events_inc_a = (
                    set(df_a.loc[df_a[dpsi_col] >= dpsi_cutoff, _match_col].dropna().unique())
                    if dpsi_col in df_a.columns
                    else set()
                )
                events_inc_b = (
                    set(df_b.loc[df_b[dpsi_col] >= dpsi_cutoff, _match_col].dropna().unique())
                    if dpsi_col in df_b.columns
                    else set()
                )
                events_exc_a = (
                    set(df_a.loc[df_a[dpsi_col] <= -dpsi_cutoff, _match_col].dropna().unique())
                    if dpsi_col in df_a.columns
                    else set()
                )
                events_exc_b = (
                    set(df_b.loc[df_b[dpsi_col] <= -dpsi_cutoff, _match_col].dropna().unique())
                    if dpsi_col in df_b.columns
                    else set()
                )
            else:
                key_a = make_event_key(df_a, et)
                key_b = make_event_key(df_b, et)
                if key_a.eq("").all() or key_b.eq("").all():
                    continue

                df_a = df_a.copy()
                df_b = df_b.copy()
                df_a["_ekey"] = key_a.values
                df_b["_ekey"] = key_b.values

                events_all_a = set(df_a["_ekey"].dropna().unique())
                events_all_b = set(df_b["_ekey"].dropna().unique())
                events_inc_a = set(
                    df_a.loc[df_a[dpsi_col] >= dpsi_cutoff, "_ekey"].dropna().unique()
                )
                events_inc_b = set(
                    df_b.loc[df_b[dpsi_col] >= dpsi_cutoff, "_ekey"].dropna().unique()
                )
                events_exc_a = set(
                    df_a.loc[df_a[dpsi_col] <= -dpsi_cutoff, "_ekey"].dropna().unique()
                )
                events_exc_b = set(
                    df_b.loc[df_b[dpsi_col] <= -dpsi_cutoff, "_ekey"].dropna().unique()
                )

            panels = [
                ("Significant in Either Condition", events_all_a, events_all_b),
                (f"Included in Both (dPSI \u2265 {dpsi_cutoff})", events_inc_a, events_inc_b),
                (f"Excluded in Both (dPSI \u2264 \u2212{dpsi_cutoff})", events_exc_a, events_exc_b),
                (f"Inc {label_a} / Exc {label_b}", events_inc_a, events_exc_b),
                (f"Exc {label_a} / Inc {label_b}", events_exc_a, events_inc_b),
            ]

            fig, axes = plt.subplots(2, 3, figsize=(20, 12))
            axes_flat = axes.flatten()

            for idx, (panel_title, set_a, set_b) in enumerate(panels):
                ax = axes_flat[idx]
                v = venn2([set_a, set_b], set_labels=(label_a, label_b), ax=ax)
                style_venn(v, 2)
                ax.set_title(
                    f"{panel_title}\n({label_a}: {len(set_a)} | {label_b}: {len(set_b)})",
                    fontsize=11,
                    fontweight="bold",
                )

            axes_flat[5].set_visible(False)
            fig.suptitle(
                f"{et} Splicing Events ({level_label}) \u2014 {label_a} vs {label_b}",
                fontsize=14,
                fontweight="bold",
            )
            plt.tight_layout()
            outpath = (
                outdir / f"venn_splicing_{et}_{name_a}_vs_{name_b}{fname_suffix}.{fig_format}"
            )
            fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved: {outpath.name}")


# ---------------------------------------------------------------------------
# rMATS event count comparison
# ---------------------------------------------------------------------------


def rmats_event_count_comparison(
    rmats_conditions: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Grouped bar chart: x=event type, groups=conditions, y=significant event count."""
    if len(rmats_conditions) == 0:
        return
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    names = list(rmats_conditions.keys())
    labels = [condition_labels[n] for n in names]
    palette = sns.color_palette("Set1", n_colors=len(names))

    x = np.arange(len(RMATS_EVENT_TYPES))
    w = 0.8 / max(len(names), 1)

    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (name, lbl, col) in enumerate(zip(names, labels, palette)):
        counts = []
        for et in RMATS_EVENT_TYPES:
            filt_df = rmats_conditions[name]["rmats_filtered"].get(et, pd.DataFrame())
            counts.append(len(filt_df))
        offset = (idx - len(names) / 2 + 0.5) * w
        bars = ax.bar(x + offset, counts, w, label=lbl, color=col, alpha=0.85)
        for bar, cnt in zip(bars, counts):
            if cnt > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    cnt + 0.5,
                    str(cnt),
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    fontweight="bold",
                )

    ax.set_xticks(x)
    ax.set_xticklabels(RMATS_EVENT_TYPES, fontsize=11)
    ax.set_ylabel("Significant Splicing Events")
    ax.set_title("rMATS \u2014 Significant Event Counts by Type Across Conditions")
    ax.legend(fontsize=9)

    outpath = outdir / f"rmats_event_count_comparison.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Pairwise dPSI scatter
# ---------------------------------------------------------------------------


def pairwise_dpsi_scatter(
    rmats_conditions: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
    match_by: str = "event",
) -> None:
    """Pairwise scatter of dPSI for shared splicing events."""
    cols = rmats_cols or DEFAULT_RMATS_COLS
    names = list(rmats_conditions.keys())
    if len(names) < 2:
        return

    dpsi_col = cols["inclevel_diff"]
    gene_col = cols["gene_name"]
    id_col = cols["gene_id"]
    pairs = list(combinations(names, 2))
    nrows, ncols_grid = grid_dims(len(pairs))

    is_gene = match_by == "gene"
    level_label = "gene-level" if is_gene else "coordinate-level"
    fname_suffix = "_genelevel" if is_gene else ""

    outdir = Path(outdir)
    if not is_gene:
        outdir = outdir / "event_level"
        outdir.mkdir(exist_ok=True, parents=True)

    fig, axes = plt.subplots(
        nrows, ncols_grid, figsize=(5 * ncols_grid, 5 * nrows), squeeze=False
    )
    axes_flat = [axes[r][c] for r in range(nrows) for c in range(ncols_grid)]

    for ax_idx, (nameA, nameB) in enumerate(pairs):
        ax = axes_flat[ax_idx]
        lblA = condition_labels[nameA]
        lblB = condition_labels[nameB]

        all_rows = []
        for et in RMATS_EVENT_TYPES:
            dfA = rmats_conditions[nameA]["rmats_filtered"].get(et, pd.DataFrame())
            dfB = rmats_conditions[nameB]["rmats_filtered"].get(et, pd.DataFrame())
            if len(dfA) == 0 or len(dfB) == 0:
                continue

            if is_gene:
                _match_col = (
                    id_col if id_col in dfA.columns and id_col in dfB.columns else gene_col
                )
                if _match_col not in dfA.columns or _match_col not in dfB.columns:
                    continue
                if dpsi_col not in dfA.columns or dpsi_col not in dfB.columns:
                    continue
                dfA = dfA.copy()
                dfB = dfB.copy()
                dfA["_abs_dpsi"] = dfA[dpsi_col].abs()
                dfB["_abs_dpsi"] = dfB[dpsi_col].abs()
                repA = dfA.loc[dfA.groupby(_match_col)["_abs_dpsi"].idxmax()]
                repB = dfB.loc[dfB.groupby(_match_col)["_abs_dpsi"].idxmax()]
                merged = repA[[_match_col, dpsi_col]].merge(
                    repB[[_match_col, dpsi_col]], on=_match_col, suffixes=("_A", "_B")
                )
                if len(merged) == 0:
                    continue
                merged["event_type"] = et
                all_rows.append(merged)
            else:
                keyA = make_event_key(dfA, et)
                keyB = make_event_key(dfB, et)
                if keyA.eq("").all() or keyB.eq("").all():
                    continue
                dfA = dfA.copy()
                dfB = dfB.copy()
                dfA["_ekey"] = keyA.values
                dfB["_ekey"] = keyB.values
                merged = dfA[["_ekey", dpsi_col]].merge(
                    dfB[["_ekey", dpsi_col]], on="_ekey", suffixes=("_A", "_B")
                )
                if len(merged) == 0:
                    continue
                merged["event_type"] = et
                all_rows.append(merged)

        if not all_rows:
            ax.set_visible(False)
            continue

        combined = pd.concat(all_rows, ignore_index=True)
        for et in RMATS_EVENT_TYPES:
            sub = combined[combined["event_type"] == et]
            if len(sub) == 0:
                continue
            ax.scatter(
                sub[f"{dpsi_col}_A"],
                sub[f"{dpsi_col}_B"],
                c=EVENT_COLORS.get(et, "#888888"),
                s=8,
                alpha=0.6,
                edgecolors="none",
                rasterized=True,
                label=f"{et} ({len(sub):,})",
            )

        lims = [
            min(combined[f"{dpsi_col}_A"].min(), combined[f"{dpsi_col}_B"].min()) - 0.05,
            max(combined[f"{dpsi_col}_A"].max(), combined[f"{dpsi_col}_B"].max()) + 0.05,
        ]
        ax.plot(lims, lims, "k--", lw=0.7, alpha=0.5)
        ax.axhline(0, color="black", lw=0.5)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        if _SCIPY_AVAILABLE and len(combined) >= 3:
            r, _ = pearsonr(combined[f"{dpsi_col}_A"], combined[f"{dpsi_col}_B"])
            x_arr = combined[f"{dpsi_col}_A"].values
            y_arr = combined[f"{dpsi_col}_B"].values
            slope, intercept = np.polyfit(x_arr, y_arr, 1)
            x_fit = np.linspace(x_arr.min(), x_arr.max(), 100)
            ax.plot(x_fit, slope * x_fit + intercept, color="#E69F00", lw=1.5, alpha=0.8)
            ax.annotate(
                f"R\u00b2 = {r**2:.3f}",
                xy=(0, 1),
                xycoords="axes fraction",
                xytext=(4, 4),
                textcoords="offset points",
                ha="left",
                va="bottom",
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="grey", alpha=0.9),
            )

        ax.set_xlabel(f"dPSI  {lblA}", fontsize=9)
        ax.set_ylabel(f"dPSI  {lblB}", fontsize=9)
        ax.set_title(f"{lblA} vs {lblB}", fontsize=10)
        ax.legend(fontsize=7, markerscale=2, loc="lower right")

    for i in range(len(pairs), nrows * ncols_grid):
        axes_flat[i].set_visible(False)

    fig.suptitle(
        f"Pairwise dPSI Comparison (Shared Splicing Events, {level_label})",
        fontsize=13,
        fontweight="bold",
    )
    plt.tight_layout()
    outpath = outdir / f"pairwise_dpsi_scatter{fname_suffix}.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# rMATS UpSet plot
# ---------------------------------------------------------------------------


def rmats_upset_plot(
    rmats_conditions: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
    match_by: str = "event",
) -> None:
    """UpSet plot for rMATS splicing event sets across 3+ conditions."""
    if not _UPSET_AVAILABLE:
        print("  Skipping rMATS UpSet: 'upsetplot' not installed")
        return
    if len(rmats_conditions) < 3:
        print("  Skipping rMATS UpSet: requires 3+ conditions")
        return

    cols = rmats_cols or DEFAULT_RMATS_COLS
    names = list(rmats_conditions.keys())
    labels = [condition_labels[n] for n in names]
    gene_col = cols["gene_name"]
    id_col = cols["gene_id"]

    is_gene = match_by == "gene"
    level_label = "gene-level" if is_gene else "coordinate-level"
    fname_suffix = "_genelevel" if is_gene else ""

    outdir = Path(outdir)
    if not is_gene:
        outdir = outdir / "event_level"
        outdir.mkdir(exist_ok=True, parents=True)

    # Combined upset
    event_sets: dict[str, set] = {}
    for name, lbl in zip(names, labels):
        items: set = set()
        for et, et_df in rmats_conditions[name]["rmats_filtered"].items():
            if len(et_df) > 0:
                if is_gene:
                    _match_col = id_col if id_col in et_df.columns else gene_col
                    if _match_col in et_df.columns:
                        items.update(et_df[_match_col].dropna().unique())
                else:
                    keys = make_event_key(et_df, et)
                    items.update(f"{et}|{k}" for k in keys[keys != ""].unique())
        event_sets[lbl] = items

    all_events = set.union(*event_sets.values()) if event_sets else set()
    if all_events:
        memberships = [
            tuple(lbl for lbl in labels if ev in event_sets[lbl]) for ev in all_events
        ]
        try:
            upset_data = from_memberships(memberships)
            upset_data = upset_data.groupby(
                level=list(range(upset_data.index.nlevels))
            ).sum()
            upset = UpSet(upset_data, show_counts=True, sort_by="cardinality")
            upset.plot()
            plt.suptitle(
                f"rMATS UpSet \u2014 Splicing Event Overlap ({level_label})",
                y=1.02,
                fontsize=12,
                fontweight="bold",
            )
            outpath = outdir / f"rmats_upset_events{fname_suffix}.{fig_format}"
            plt.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
            plt.close("all")
            print(f"  Saved: {outpath}")
        except Exception as e:
            print(f"  rMATS UpSet (events) failed: {e}")


# ---------------------------------------------------------------------------
# rMATS event heatmap
# ---------------------------------------------------------------------------


def rmats_event_heatmap(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    event_type: str,
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    dpsi_cutoff: float = 0.1,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Clustered heatmap of dPSI values across conditions for one event type."""
    from matplotlib.colors import LinearSegmentedColormap

    cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = list(condition_results.keys())
    labels = [condition_labels[n] for n in names]
    dpsi_col = cols["inclevel_diff"]
    gene_col = cols["gene_name"]

    event_data: dict[str, pd.Series] = {}
    gene_lookup: dict[str, str] = {}

    for name, lbl in zip(names, labels):
        filt_df = condition_results[name]["rmats_filtered"].get(event_type, pd.DataFrame())
        if len(filt_df) == 0:
            continue
        df = filt_df.copy()
        df["_ekey"] = make_event_key(df, event_type).values
        df = df[df["_ekey"] != ""]
        if len(df) == 0:
            continue

        grouped = df.groupby("_ekey")[dpsi_col].mean()
        event_data[lbl] = grouped

        if gene_col in df.columns:
            for ekey, gene in df.groupby("_ekey")[gene_col].first().items():
                if ekey not in gene_lookup and pd.notna(gene):
                    gene_lookup[ekey] = gene

    if len(event_data) < 2:
        print(f"  Skipping {event_type} heatmap (insufficient conditions)")
        return

    dpsi_df = pd.DataFrame(event_data).dropna(how="all")
    if len(dpsi_df) == 0:
        return

    max_events = 80
    if len(dpsi_df) > max_events:
        max_abs = dpsi_df.abs().max(axis=1)
        dpsi_df = dpsi_df.loc[max_abs.nlargest(max_events).index]

    row_labels = []
    for ekey in dpsi_df.index:
        gene = gene_lookup.get(ekey, "")
        if gene:
            coord_short = ekey.split(":")[0] + ":" + ekey.split(":")[-1]
            row_labels.append(f"{gene} ({coord_short})")
        else:
            row_labels.append(ekey[:40])
    dpsi_df.index = row_labels

    dpsi_filled = dpsi_df.fillna(0)
    cmap = LinearSegmentedColormap.from_list(
        "dpsi_diverging", [COLOR_DOWN, "white", COLOR_UP], N=256
    )
    vmax = max(abs(dpsi_filled.values.min()), abs(dpsi_filled.values.max()))
    vmax = max(vmax, dpsi_cutoff)
    fig_height = max(8, len(dpsi_df) * 0.25 + 2)
    fig_width = max(6, len(labels) * 1.5 + 4)

    try:
        g = sns.clustermap(
            dpsi_filled,
            cmap=cmap,
            center=0,
            vmin=-vmax,
            vmax=vmax,
            figsize=(fig_width, fig_height),
            linewidths=0.5,
            linecolor="white",
            yticklabels=True,
            xticklabels=True,
            row_cluster=len(dpsi_filled) > 1,
            col_cluster=len(dpsi_filled.columns) > 1,
        )
        g.fig.suptitle(
            f"dPSI Heatmap \u2014 {event_type} Events (top {len(dpsi_df)} by |dPSI|)",
            fontsize=13,
            fontweight="bold",
            y=1.02,
        )
        outpath = outdir / f"heatmap_dpsi_{event_type}.{fig_format}"
        g.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
        plt.close("all")
        print(f"  Saved: {outpath}")
    except Exception as e:
        print(f"  {event_type} heatmap failed: {e}")
        plt.close("all")


# ---------------------------------------------------------------------------
# rMATS event pie chart
# ---------------------------------------------------------------------------


def rmats_event_pie_chart(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Pie chart per condition showing splicing event type distribution."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    names = list(condition_results.keys())
    labels = [condition_labels[n] for n in names]
    n_conds = len(names)

    if n_conds == 0:
        return

    nrows, ncols_grid = grid_dims(n_conds)
    fig, axes = plt.subplots(nrows, ncols_grid, figsize=(5 * ncols_grid, 5 * nrows))
    if n_conds == 1:
        axes_flat = [axes]
    else:
        axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for idx, (name, lbl) in enumerate(zip(names, labels)):
        ax = axes_flat[idx]
        counts = []
        et_labels = []
        colors = []
        for et in RMATS_EVENT_TYPES:
            filt_df = condition_results[name]["rmats_filtered"].get(et, pd.DataFrame())
            n_events = len(filt_df)
            if n_events > 0:
                counts.append(n_events)
                et_labels.append(et)
                colors.append(EVENT_COLORS.get(et, "#888888"))

        if not counts:
            ax.text(0.5, 0.5, "No significant\nevents", ha="center", va="center", fontsize=11, transform=ax.transAxes)
            ax.set_title(lbl, fontsize=12, fontweight="bold")
            ax.axis("off")
            continue

        total = sum(counts)
        ax.pie(
            counts,
            labels=et_labels,
            colors=colors,
            autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct / 100 * total))})",
            startangle=90,
            pctdistance=0.65,
            textprops={"fontsize": 9},
        )
        ax.set_title(f"{lbl}\n(n={total:,} events)", fontsize=12, fontweight="bold")

    for i in range(n_conds, len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.suptitle(
        "rMATS \u2014 Significant Event Type Distribution",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout()
    outpath = outdir / f"rmats_event_type_pie.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Combined DESeq2 + rMATS analyses
# ---------------------------------------------------------------------------


def deseq2_vs_rmats_venn(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    rmats_cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Per-condition Venn: DE genes vs genes with significant splicing events."""
    if not _VENN_AVAILABLE:
        print("  Skipping DE vs AS Venn (matplotlib-venn not installed)")
        return
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS
    r_cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_col_de = cols["gene_name"]
    gene_col_as = r_cols["gene_name"]

    for name, data in condition_results.items():
        lbl = condition_labels[name]
        filt_de = data["deseq2_filtered"].get("all_genes", pd.DataFrame())
        rmats_f = data.get("rmats_filtered", {})

        if len(filt_de) == 0 or not any(len(v) > 0 for v in rmats_f.values()):
            continue

        de_genes = set(filt_de[gene_col_de].dropna().unique())
        as_genes: set = set()
        for et_df in rmats_f.values():
            if len(et_df) > 0 and gene_col_as in et_df.columns:
                as_genes.update(et_df[gene_col_as].dropna().unique())

        fig, ax = plt.subplots(figsize=(7, 7))
        venn2([de_genes, as_genes], set_labels=["DE genes", "AS genes"], ax=ax)
        ax.set_title(f"{lbl} \u2014 DE vs Alternatively Spliced Genes", fontsize=12, fontweight="bold")

        outpath = outdir / f"deseq2_vs_rmats_venn_{name}.{fig_format}"
        fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
        plt.close(fig)
        print(f"  Saved: {outpath}")


def log2fc_vs_dpsi_scatter(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    rmats_cols: dict[str, str] | None = None,
    dpsi_cutoff: float = 0.1,
    log2fc_cutoff: float = 1.0,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Scatter of log2FC vs dPSI for genes that are both DE and AS."""
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS
    r_cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_col_de = cols["gene_name"]
    fc_col = cols["log2fc"]
    gene_col_as = r_cols["gene_name"]
    dpsi_col = r_cols["inclevel_diff"]

    valid_conds = []
    for name, data in condition_results.items():
        filt_de = data["deseq2_filtered"].get("all_genes", pd.DataFrame())
        rmats_f = data.get("rmats_filtered", {})
        if len(filt_de) > 0 and any(len(v) > 0 for v in rmats_f.values()):
            valid_conds.append(name)

    if not valid_conds:
        print("  Skipping log2FC vs dPSI: no conditions with both DE and AS data")
        return

    nrows, ncols_grid = grid_dims(len(valid_conds))
    fig, axes = plt.subplots(
        nrows, ncols_grid, figsize=(6 * ncols_grid, 6 * nrows), squeeze=False
    )
    axes_flat = [axes[r][c] for r in range(nrows) for c in range(ncols_grid)]

    for ax_idx, name in enumerate(valid_conds):
        ax = axes_flat[ax_idx]
        lbl = condition_labels[name]
        data = condition_results[name]
        filt_de = data["deseq2_filtered"]["all_genes"]
        rmats_f = data["rmats_filtered"]

        de_lookup = filt_de.dropna(subset=[gene_col_de, fc_col]).set_index(gene_col_de)[fc_col]
        de_lookup = de_lookup[~de_lookup.index.duplicated(keep="first")]

        plotted = False
        for et in RMATS_EVENT_TYPES:
            et_df = rmats_f.get(et, pd.DataFrame())
            if len(et_df) == 0 or gene_col_as not in et_df.columns:
                continue
            et_df = et_df.dropna(subset=[gene_col_as, dpsi_col]).copy()
            et_df = et_df[et_df[gene_col_as].isin(de_lookup.index)]
            if len(et_df) == 0:
                continue

            x_vals = de_lookup.loc[et_df[gene_col_as].values].values
            y_vals = et_df[dpsi_col].values
            ax.scatter(
                x_vals,
                y_vals,
                c=EVENT_COLORS.get(et, "#888888"),
                s=8,
                alpha=0.65,
                edgecolors="none",
                rasterized=True,
                label=f"{et} ({len(et_df):,})",
            )
            plotted = True

        if not plotted:
            ax.set_visible(False)
            continue

        ax.axhline(0, color="black", lw=0.7)
        ax.axvline(0, color="black", lw=0.7)
        ax.axhline(dpsi_cutoff, color="grey", ls="--", lw=0.6)
        ax.axhline(-dpsi_cutoff, color="grey", ls="--", lw=0.6)
        ax.axvline(log2fc_cutoff, color="grey", ls="--", lw=0.6)
        ax.axvline(-log2fc_cutoff, color="grey", ls="--", lw=0.6)

        ax.set_xlabel("log$_2$ Fold Change (DESeq2)", fontsize=9)
        ax.set_ylabel("$\\Delta$PSI (rMATS)", fontsize=9)
        ax.set_title(f"{lbl}", fontsize=10)
        ax.legend(fontsize=7, markerscale=2)

    for i in range(len(valid_conds), nrows * ncols_grid):
        axes_flat[i].set_visible(False)

    fig.suptitle(
        "Combined: log$_2$FC vs $\\Delta$PSI\n(Genes Both DE and Alternatively Spliced)",
        fontsize=12,
        fontweight="bold",
    )
    plt.tight_layout()
    outpath = outdir / f"log2fc_vs_dpsi_scatter.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ---------------------------------------------------------------------------
# Gene overlap summary
# ---------------------------------------------------------------------------


def gene_overlap_summary(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    rmats_cols: dict[str, str] | None = None,
) -> None:
    """Create a table showing which genes are DE, AS, or both per condition."""
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS
    r_cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_name_col_de = cols.get("gene_name", "gene_name")
    gene_name_col_as = r_cols.get("gene_name", "geneSymbol")

    de_genes_per_cond: dict[str, set] = {}
    as_genes_per_cond: dict[str, set] = {}
    cond_names = list(condition_results.keys())

    for cond_name in cond_names:
        res = condition_results[cond_name]
        de_set: set = set()
        filt = res.get("deseq2_filtered", {}).get("all_genes", pd.DataFrame())
        if not filt.empty and gene_name_col_de in filt.columns:
            de_set = set(filt[gene_name_col_de].dropna().unique())
        de_genes_per_cond[cond_name] = de_set

        as_set: set = set()
        rmats_filt = res.get("rmats_filtered", {})
        for et_df in rmats_filt.values():
            if not et_df.empty and gene_name_col_as in et_df.columns:
                as_set.update(et_df[gene_name_col_as].dropna().unique())
        as_genes_per_cond[cond_name] = as_set

    all_genes: set = set()
    for s in list(de_genes_per_cond.values()) + list(as_genes_per_cond.values()):
        all_genes.update(s)

    if not all_genes:
        print("  Gene overlap summary: no genes found across conditions.")
        return

    rows = []
    for gene in sorted(all_genes):
        row: dict = {"Gene": gene}
        n_conditions = 0
        for cond_name in cond_names:
            lbl = condition_labels.get(cond_name, cond_name)
            is_de = gene in de_genes_per_cond[cond_name]
            is_as = gene in as_genes_per_cond[cond_name]
            if is_de and is_as:
                row[lbl] = "DE+AS"
                n_conditions += 1
            elif is_de:
                row[lbl] = "DE"
                n_conditions += 1
            elif is_as:
                row[lbl] = "AS"
                n_conditions += 1
            else:
                row[lbl] = ""
        row["Conditions"] = n_conditions
        rows.append(row)

    summary_df = pd.DataFrame(rows)
    summary_df = summary_df.sort_values(
        by=["Conditions", "Gene"], ascending=[False, True]
    ).reset_index(drop=True)

    outpath = outdir / "gene_overlap_summary.xlsx"
    summary_df.to_excel(outpath, index=False, freeze_panes=(1, 1))
    print(f"  Saved gene overlap summary: {len(summary_df)} genes across {len(cond_names)} conditions")
    print(f"    -> {outpath}")


# ---------------------------------------------------------------------------
# Summary dashboard
# ---------------------------------------------------------------------------


def summary_dashboard(
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    go_results: dict | None,
    gsea_results: dict | None,
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    rmats_cols: dict[str, str] | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Create a 2x3 summary dashboard with key statistics."""
    if cols is None:
        cols = DEFAULT_DESEQ2_COLS
    r_cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cond_names = list(condition_results.keys())
    cond_labels_list = [condition_labels.get(cn, cn) for cn in cond_names]
    n_conds = len(cond_names)

    gene_name_col_de = cols.get("gene_name", "gene_name")
    gene_name_col_as = r_cols.get("gene_name", "geneSymbol")

    # Panel 1: DEG counts
    deg_up = []
    deg_down = []
    for cn in cond_names:
        filt = condition_results[cn].get("deseq2_filtered", {}).get("all_genes", pd.DataFrame())
        if not filt.empty and "direction" in filt.columns:
            deg_up.append(int((filt["direction"].str.lower() == "up").sum()))
            deg_down.append(int((filt["direction"].str.lower() == "down").sum()))
        else:
            deg_up.append(0)
            deg_down.append(0)

    # Panel 2: Splicing events
    splice_counts = {et: [] for et in RMATS_EVENT_TYPES}
    for cn in cond_names:
        rmats_filt = condition_results[cn].get("rmats_filtered", {})
        for et in RMATS_EVENT_TYPES:
            df = rmats_filt.get(et, pd.DataFrame())
            splice_counts[et].append(len(df))

    # Panel 3: DE vs AS overlap
    de_only_counts = []
    as_only_counts = []
    both_counts = []
    all_de_genes: set = set()
    all_as_genes: set = set()
    all_both_genes: set = set()
    for cn in cond_names:
        filt = condition_results[cn].get("deseq2_filtered", {}).get("all_genes", pd.DataFrame())
        de_set: set = set()
        if not filt.empty and gene_name_col_de in filt.columns:
            de_set = set(filt[gene_name_col_de].dropna().unique())
        as_set: set = set()
        rmats_filt = condition_results[cn].get("rmats_filtered", {})
        for et_df in rmats_filt.values():
            if not et_df.empty and gene_name_col_as in et_df.columns:
                as_set.update(et_df[gene_name_col_as].dropna().unique())
        both = de_set & as_set
        de_only_counts.append(len(de_set - both))
        as_only_counts.append(len(as_set - both))
        both_counts.append(len(both))
        all_de_genes.update(de_set)
        all_as_genes.update(as_set)
        all_both_genes.update(both)

    # Panel 4: GO terms
    go_available = go_results is not None and len(go_results) > 0
    go_up_counts = []
    go_down_counts = []
    if go_available:
        for cn in cond_names:
            cond_go = go_results.get(cn, {})
            if cond_go is None:
                cond_go = {}
            up_df = cond_go.get("up", None)
            dn_df = cond_go.get("down", None)
            go_up_counts.append(
                len(up_df) if up_df is not None and not (isinstance(up_df, pd.DataFrame) and up_df.empty) else 0
            )
            go_down_counts.append(
                len(dn_df) if dn_df is not None and not (isinstance(dn_df, pd.DataFrame) and dn_df.empty) else 0
            )

    # Build figure
    fig = plt.figure(figsize=(18, 10))
    gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.30)

    x = np.arange(n_conds)
    bar_width = 0.35

    # Panel 1: DEG Counts
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.bar(x - bar_width / 2, deg_up, bar_width, label="Up", color=COLOR_UP)
    ax1.bar(x + bar_width / 2, deg_down, bar_width, label="Down", color=COLOR_DOWN)
    ax1.set_xticks(x)
    ax1.set_xticklabels(cond_labels_list, fontsize=9, rotation=15, ha="right")
    ax1.set_ylabel("Number of DEGs")
    ax1.set_title("DEG Counts", fontsize=12, fontweight="bold")
    ax1.legend(fontsize=9)

    # Panel 2: Splicing Events
    ax2 = fig.add_subplot(gs[0, 1])
    bottoms = np.zeros(n_conds)
    for et in RMATS_EVENT_TYPES:
        vals = np.array(splice_counts[et])
        ax2.bar(x, vals, bar_width * 1.5, bottom=bottoms, label=et, color=EVENT_COLORS.get(et, "#888888"))
        bottoms += vals
    ax2.set_xticks(x)
    ax2.set_xticklabels(cond_labels_list, fontsize=9, rotation=15, ha="right")
    ax2.set_ylabel("Number of Events")
    ax2.set_title("Splicing Events by Type", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=8, ncol=2, loc="upper right")

    # Panel 3: DE vs AS Overlap
    ax3 = fig.add_subplot(gs[0, 2])
    w = 0.25
    ax3.bar(x - w, de_only_counts, w, label="DE only", color=COLOR_UP)
    ax3.bar(x, as_only_counts, w, label="AS only", color="#009E73")
    ax3.bar(x + w, both_counts, w, label="DE + AS", color="#CC79A7")
    ax3.set_xticks(x)
    ax3.set_xticklabels(cond_labels_list, fontsize=9, rotation=15, ha="right")
    ax3.set_ylabel("Number of Genes")
    ax3.set_title("DE vs AS Overlap", fontsize=12, fontweight="bold")
    ax3.legend(fontsize=9)

    # Panel 4: GO Terms
    ax4 = fig.add_subplot(gs[1, 0])
    if go_available and (sum(go_up_counts) + sum(go_down_counts)) > 0:
        ax4.bar(x - bar_width / 2, go_up_counts, bar_width, label="Up-regulated", color=COLOR_UP)
        ax4.bar(x + bar_width / 2, go_down_counts, bar_width, label="Down-regulated", color=COLOR_DOWN)
        ax4.set_xticks(x)
        ax4.set_xticklabels(cond_labels_list, fontsize=9, rotation=15, ha="right")
        ax4.set_ylabel("Enriched GO Terms")
        ax4.legend(fontsize=9)
    else:
        ax4.text(0.5, 0.5, "GO Enrichment\nN/A", transform=ax4.transAxes, fontsize=16, ha="center", va="center", color="grey")
        ax4.set_xticks([])
        ax4.set_yticks([])
    ax4.set_title("GO Terms Enriched", fontsize=12, fontweight="bold")

    # Panel 5: placeholder
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.text(0.5, 0.5, "GSEA Pathways\nN/A", transform=ax5.transAxes, fontsize=16, ha="center", va="center", color="grey")
    ax5.set_xticks([])
    ax5.set_yticks([])
    ax5.set_title("GSEA Pathways", fontsize=12, fontweight="bold")

    # Panel 6: Key Numbers
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.set_xticks([])
    ax6.set_yticks([])
    for spine in ax6.spines.values():
        spine.set_visible(False)
    ax6.set_facecolor("#F7F7F7")

    summary_lines = [
        ("Total Unique DEGs", f"{len(all_de_genes):,}"),
        ("Total Unique AS Genes", f"{len(all_as_genes):,}"),
        ("Genes in Both DE + AS", f"{len(all_both_genes):,}"),
        ("Conditions Analyzed", str(n_conds)),
    ]
    ax6.set_title("Key Numbers", fontsize=12, fontweight="bold")
    y_start = 0.88
    y_step = 0.18
    for i, (key, val) in enumerate(summary_lines):
        y_pos = y_start - i * y_step
        ax6.text(0.08, y_pos, key + ":", transform=ax6.transAxes, fontsize=11, fontweight="bold", va="top", ha="left")
        ax6.text(0.92, y_pos, val, transform=ax6.transAxes, fontsize=13, fontweight="bold", va="top", ha="right", color=COLOR_DOWN)

    fig.suptitle("RNA-seq Analysis Summary Dashboard", fontsize=16, fontweight="bold", y=0.98)

    outpath = outdir / f"analysis_summary_dashboard.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")
