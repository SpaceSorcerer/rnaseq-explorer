"""QC analysis: PCA, sample correlation, and top DEG heatmaps.

These analyses require a normalized counts matrix (genes x samples)
and are used for quality control visualization.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from rnaseq_explorer.viz.theme import CONDITION_COLORS


# ---------------------------------------------------------------------------
# Counts matrix loading
# ---------------------------------------------------------------------------


def load_counts_matrix(
    counts_path: str | Path,
    sample_metadata: dict[str, str] | None = None,
    conditions: list[dict] | None = None,
) -> tuple[Optional[pd.DataFrame], dict[str, str]]:
    """Load a normalized counts matrix (genes x samples) for QC plots.

    Parameters
    ----------
    counts_path : str or Path
        Path to normalized_counts.tsv or .csv.
    sample_metadata : dict or None
        {sample_name: condition_label} mapping. Auto-detected if None.
    conditions : list or None
        The CONDITIONS list for auto-detection.

    Returns
    -------
    tuple
        (counts_df, metadata_dict) or (None, {}) on failure.
    """
    if not counts_path:
        print("[INFO] No counts file provided -- skipping PCA, correlation heatmap, top DEG heatmap")
        return None, {}

    path = Path(counts_path)
    if not path.exists():
        print(f"[WARNING] Counts file not found: {counts_path}")
        return None, {}

    ext = path.suffix.lower()
    sep = "\t" if ext in (".tsv", ".tab") else ","
    try:
        with open(path, "r") as fh:
            header_line = fh.readline()
        first_field = header_line.split(sep)[0].strip().strip('"')
        if first_field == "" or first_field == "X":
            df = pd.read_csv(path, sep=sep, index_col=0, quotechar='"')
            print("  [INFO] Detected R row.names format -- using first column as index")
        else:
            df = pd.read_csv(path, sep=sep, quotechar='"')
            if df.columns[0] in ("gene_id", "GeneID", "Geneid", "X") or \
               df.iloc[:, 0].astype(str).str.upper().str.startswith("ENS").mean() > 0.1:
                df = df.set_index(df.columns[0])
                print(f"  [INFO] Set '{df.index.name}' as gene index")
    except Exception as e:
        print(f"[WARNING] Failed to load counts file: {e}")
        return None, {}

    print(f"  Loaded counts matrix: {df.shape[0]:,} genes x {df.shape[1]} samples")

    # Strip Ensembl version suffixes from index
    idx_str = df.index.astype(str)
    ens_frac = idx_str.str.upper().str.startswith("ENS").sum() / max(len(idx_str), 1)
    if ens_frac > 0.1:
        df.index = idx_str.str.replace(r"\.\d+$", "", regex=True)
        print("  [INFO] Stripped Ensembl version suffixes from gene index")

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) < df.shape[1]:
        dropped = [c for c in df.columns if c not in numeric_cols]
        print(f"  [INFO] Dropping non-numeric columns: {dropped}")
        df = df[numeric_cols]

    if df.shape[1] < 2:
        print("[WARNING] Counts matrix has fewer than 2 sample columns, skipping")
        return None, {}

    # Auto-detect sample metadata from CONDITIONS if not provided
    if not sample_metadata and conditions:
        sample_metadata = {}
        for col in df.columns:
            col_lower = str(col).lower()
            for cond in conditions:
                cond_name = cond.get("name", "")
                cond_label = cond.get("label", cond_name)
                if cond_name.lower() in col_lower or cond_label.lower() in col_lower:
                    sample_metadata[col] = cond_label
                    break
        if sample_metadata:
            matched = len(sample_metadata)
            total = len(df.columns)
            print(f"  [INFO] Auto-detected metadata for {matched}/{total} samples from CONDITIONS")
        else:
            print("  [INFO] Could not auto-detect sample metadata")

    if sample_metadata is None:
        sample_metadata = {}

    return df, sample_metadata


# ---------------------------------------------------------------------------
# PCA plot
# ---------------------------------------------------------------------------


def compute_pca(
    counts_df: pd.DataFrame,
    metadata: dict[str, str] | None = None,
    outdir: str | Path | None = None,
    pca_file: str | Path | None = None,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate a PCA scatter plot from counts data or a pre-computed PCA file.

    Parameters
    ----------
    counts_df : pd.DataFrame or None
        Normalized counts matrix (genes x samples).
    metadata : dict or None
        {sample_name: condition_label} for coloring.
    outdir : Path or None
        Directory to save the plot.
    pca_file : str or Path or None
        Path to pre-computed pca_data.csv.
    fig_format : str
        Figure format.
    fig_dpi : int
        Figure DPI.
    """
    if outdir is None:
        return
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    oi_palette = CONDITION_COLORS

    # Mode 1: Pre-computed PCA file
    if pca_file is not None:
        pca_path = Path(pca_file)
        if not pca_path.exists():
            print(f"  [WARNING] PCA file not found: {pca_file}")
            return
        pca_df = pd.read_csv(pca_path)
        required = {"sample", "PC1", "PC2"}
        if not required.issubset(set(pca_df.columns)):
            print(f"  [WARNING] PCA file missing required columns {required}")
            return

        conditions_list = pca_df["condition"].unique() if "condition" in pca_df.columns else ["Unknown"]
        color_map = {c: oi_palette[i % len(oi_palette)] for i, c in enumerate(conditions_list)}

        pc1_var = pca_df["PC1_variance_pct"].iloc[0] if "PC1_variance_pct" in pca_df.columns else None
        pc2_var = pca_df["PC2_variance_pct"].iloc[0] if "PC2_variance_pct" in pca_df.columns else None

        fig, ax = plt.subplots(figsize=(8, 6))
        for cond in conditions_list:
            mask = pca_df["condition"] == cond if "condition" in pca_df.columns else [True] * len(pca_df)
            subset = pca_df[mask]
            ax.scatter(subset["PC1"], subset["PC2"], c=color_map.get(cond, "#999999"),
                       label=cond, s=80, edgecolor="white", linewidth=0.5, zorder=3)

        xlabel = f"PC1 ({pc1_var:.1f}% variance)" if pc1_var is not None else "PC1"
        ylabel = f"PC2 ({pc2_var:.1f}% variance)" if pc2_var is not None else "PC2"
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title("PCA -- Sample Clustering")
        ax.axhline(0, color="#CCCCCC", linewidth=0.8, zorder=1)
        ax.axvline(0, color="#CCCCCC", linewidth=0.8, zorder=1)
        ax.legend(title="Condition", frameon=True)
        plt.tight_layout()
        fname = outdir / f"pca_plot.{fig_format}"
        fig.savefig(fname, dpi=fig_dpi)
        plt.close(fig)
        print(f"  Saved: {fname.name}")
        return

    # Mode 2: Compute PCA from counts
    if counts_df is None or counts_df.shape[1] < 2:
        print("  [INFO] Insufficient data for PCA plot")
        return

    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        print("  [WARNING] scikit-learn not installed, skipping PCA plot")
        return

    log_counts = np.log2(counts_df + 1).T
    scaler = StandardScaler()
    scaled = scaler.fit_transform(log_counts)

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled)
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100

    sample_names = counts_df.columns.tolist()
    if metadata:
        conditions_list = list(dict.fromkeys(metadata.get(s, "Unknown") for s in sample_names))
    else:
        conditions_list = ["Unknown"]
    color_map = {c: oi_palette[i % len(oi_palette)] for i, c in enumerate(conditions_list)}

    fig, ax = plt.subplots(figsize=(8, 6))
    for i, sample in enumerate(sample_names):
        cond = metadata.get(sample, "Unknown") if metadata else "Unknown"
        ax.scatter(pcs[i, 0], pcs[i, 1], c=color_map.get(cond, "#999999"),
                   s=80, edgecolor="white", linewidth=0.5, zorder=3,
                   label=cond if cond not in [ax.get_legend_handles_labels()[1]] else "")

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), title="Condition", frameon=True)

    ax.set_xlabel(f"PC1 ({pc1_var:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2_var:.1f}% variance)")
    ax.set_title("PCA -- Sample Clustering")
    ax.axhline(0, color="#CCCCCC", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#CCCCCC", linewidth=0.8, zorder=1)

    ax.text(0.02, 0.02,
            "PCA computed from log2(counts+1); for publication, use DESeq2 VST (PMID: 25516281)",
            transform=ax.transAxes, fontsize=7, color="#666666", style="italic",
            verticalalignment="bottom")

    plt.tight_layout()
    fname = outdir / f"pca_plot.{fig_format}"
    fig.savefig(fname, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {fname.name}")


# ---------------------------------------------------------------------------
# Sample correlation heatmap
# ---------------------------------------------------------------------------


def compute_sample_correlation(
    counts_df: pd.DataFrame,
    metadata: dict[str, str],
    outdir: str | Path,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate a sample-sample correlation heatmap with hierarchical clustering.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Normalized counts matrix (genes x samples).
    metadata : dict
        {sample_name: condition_label}.
    outdir : str or Path
        Directory to save the plot.
    fig_format : str
        Figure format.
    fig_dpi : int
        Figure DPI.
    """
    if counts_df is None or counts_df.shape[1] < 2:
        print("  [INFO] Insufficient data for correlation heatmap")
        return

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log_counts = np.log2(counts_df + 1)
    corr = log_counts.corr(method="pearson")

    oi_palette = CONDITION_COLORS
    if metadata:
        conditions_list = list(dict.fromkeys(metadata.get(s, "Unknown") for s in corr.columns))
        color_map = {c: oi_palette[i % len(oi_palette)] for i, c in enumerate(conditions_list)}
        col_colors = pd.Series(
            {s: color_map.get(metadata.get(s, "Unknown"), "#999999") for s in corr.columns},
            name="Condition",
        )
    else:
        col_colors = None

    try:
        g = sns.clustermap(
            corr, method="average", metric="euclidean", cmap="viridis",
            vmin=corr.values[np.triu_indices_from(corr.values, k=1)].min() if len(corr) > 1 else 0,
            vmax=1.0, col_colors=col_colors, row_colors=col_colors,
            linewidths=0.5,
            figsize=(max(8, len(corr) * 0.6), max(7, len(corr) * 0.55)),
            dendrogram_ratio=0.12,
            cbar_pos=(0.02, 0.8, 0.03, 0.15),
        )
        g.ax_heatmap.set_title("Sample Correlation (Pearson, log2 counts+1)", pad=20)

        if metadata:
            legend_patches = [mpatches.Patch(color=color_map[c], label=c) for c in conditions_list]
            g.ax_heatmap.legend(
                handles=legend_patches, title="Condition",
                bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9,
            )

        fname = outdir / f"sample_correlation_heatmap.{fig_format}"
        g.savefig(fname, dpi=fig_dpi, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {fname.name}")
    except Exception as e:
        print(f"  [WARNING] Failed to generate correlation heatmap: {e}")


# ---------------------------------------------------------------------------
# Top DEG heatmap
# ---------------------------------------------------------------------------


def compute_top_deg_heatmap(
    counts_df: pd.DataFrame,
    condition_results: dict[str, dict],
    condition_labels: dict[str, str],
    metadata: dict[str, str],
    outdir: str | Path,
    cols: dict[str, str] | None = None,
    top_n: int = 50,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Generate a heatmap of top DEGs (by padj) across all conditions.

    Parameters
    ----------
    counts_df : pd.DataFrame
        Normalized counts matrix (genes x samples).
    condition_results : dict
        Pipeline condition_results structure.
    condition_labels : dict
        Maps condition name -> label.
    metadata : dict
        {sample_name: condition_label}.
    outdir : str or Path
        Directory to save the plot.
    cols : dict or None
        DESeq2 column name mapping.
    top_n : int
        Number of top genes to display.
    fig_format : str
        Figure format.
    fig_dpi : int
        Figure DPI.
    """
    from rnaseq_explorer.engine.deseq2 import DEFAULT_DESEQ2_COLS

    if cols is None:
        cols = DEFAULT_DESEQ2_COLS
    if counts_df is None or counts_df.shape[1] < 2:
        print("  [INFO] Insufficient counts data for top DEG heatmap")
        return

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    padj_col = cols.get("padj", "padj")
    gene_id_col = cols.get("gene_id", "gene_id")
    gene_name_col = cols.get("gene_name", "gene_name")

    all_top_genes: list[str] = []
    for cond_name, data in condition_results.items():
        deg_df = data.get("deseq2_filtered", {}).get("all_genes", pd.DataFrame())
        if len(deg_df) == 0:
            continue
        if padj_col in deg_df.columns and gene_id_col in deg_df.columns:
            sorted_df = deg_df.sort_values(padj_col)
            all_top_genes.extend(sorted_df[gene_id_col].dropna().astype(str).tolist())

    if not all_top_genes:
        print("  [INFO] No DEGs found for top DEG heatmap")
        return

    seen: set[str] = set()
    unique_genes: list[str] = []
    for g in all_top_genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)
    top_genes = unique_genes[:top_n]

    counts_idx = counts_df.index.astype(str)
    mask = counts_idx.isin(top_genes)
    if mask.sum() == 0:
        top_genes_stripped = [g.split(".")[0] for g in top_genes]
        mask = counts_idx.isin(top_genes_stripped)
    if mask.sum() == 0:
        print("  [INFO] None of the top DEGs found in counts matrix -- skipping heatmap")
        return

    subset = counts_df.loc[mask].copy()
    print(f"  Top DEG heatmap: {len(subset)} genes matched in counts matrix")

    from scipy import stats

    log_counts = np.log2(subset + 1)
    z_scored = log_counts.apply(lambda row: stats.zscore(row, nan_policy="omit"), axis=1)
    z_scored = z_scored.clip(-3, 3)

    # Try gene names as row labels
    id_to_name: dict[str, str] = {}
    for cond_name, data in condition_results.items():
        raw_df = data.get("deseq2_raw", pd.DataFrame())
        if gene_id_col in raw_df.columns and gene_name_col in raw_df.columns:
            for _, row in raw_df[[gene_id_col, gene_name_col]].dropna().iterrows():
                gid = str(row[gene_id_col])
                gname = str(row[gene_name_col])
                if gid and gname and gname != "nan":
                    id_to_name[gid] = gname
    if id_to_name:
        z_scored.index = [id_to_name.get(str(g), str(g)) for g in z_scored.index]

    oi_palette = CONDITION_COLORS
    if metadata:
        conditions_unique = list(dict.fromkeys(metadata.get(s, "Unknown") for s in z_scored.columns))
        color_map = {c: oi_palette[i % len(oi_palette)] for i, c in enumerate(conditions_unique)}
        col_colors = pd.Series(
            {s: color_map.get(metadata.get(s, "Unknown"), "#999999") for s in z_scored.columns},
            name="Condition",
        )
    else:
        col_colors = None

    try:
        n_genes = len(z_scored)
        fig_height = max(8, n_genes * 0.25)
        g = sns.clustermap(
            z_scored, cmap="RdBu_r", center=0, vmin=-3, vmax=3,
            col_colors=col_colors, method="ward", metric="euclidean",
            linewidths=0.3,
            figsize=(max(8, len(z_scored.columns) * 0.8), fig_height),
            dendrogram_ratio=(0.1, 0.08),
            cbar_pos=(0.02, 0.8, 0.03, 0.15),
            yticklabels=True,
        )
        g.ax_heatmap.set_title(f"Top {top_n} DEGs -- z-scored log2(counts+1)", pad=20)
        g.ax_heatmap.set_ylabel("")
        g.ax_heatmap.tick_params(axis="y", labelsize=7)

        if metadata:
            legend_patches = [mpatches.Patch(color=color_map[c], label=c) for c in conditions_unique]
            g.ax_heatmap.legend(
                handles=legend_patches, title="Condition",
                bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9,
            )

        fname = outdir / f"top_deg_heatmap.{fig_format}"
        g.savefig(fname, dpi=fig_dpi, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {fname.name}")
    except Exception as e:
        print(f"  [WARNING] Failed to generate top DEG heatmap: {e}")
