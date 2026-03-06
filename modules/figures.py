#!/usr/bin/env python3
"""
Figure generation functions for RNA-seq analysis pipeline.

Contains all publication-quality plotting functions:
- Volcano plots (DEG and splicing)
- MA plots
- Venn diagrams (regular and directional 3-panel)
- Violin plots (log2FC and dPSI distributions)
- Biotype pie charts and bar charts
- Concordance scatter plots
- Heatmaps
- Summary dashboards

All figures use:
- Okabe-Ito color palette for color-blind accessibility
- 300+ DPI for publication quality
- Consistent styling and formatting
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

# Try to import optional dependencies
try:
    from matplotlib_venn import venn2, venn3
    HAS_VENN = True
except ImportError:
    HAS_VENN = False

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False


# ─── Color Palettes ──────────────────────────────────────────────────────────

# Okabe-Ito color palette (color-blind friendly)
COLOR_UP = "#D72638"
COLOR_DOWN = "#1B98E0"
COLOR_NS = "#BFBFBF"

EVENT_COLORS = {
    "SE": "#E64B35", "A3SS": "#4DBBD5", "A5SS": "#00A087",
    "RI": "#3C5488", "MXE": "#F39B7F"
}

BIOTYPE_COLORS = {
    "Protein Coding": "#4C72B0", "lncRNA": "#DD8452",
    "Pseudogene": "#55A868", "Small ncRNA": "#C44E52", "Other": "#8172B2"
}


# ─── Volcano Plots ───────────────────────────────────────────────────────────

def plot_volcano(df_all, df_sig, label, output_path, log2fc_cutoff=0.4, padj_cutoff=0.01):
    """
    Generate publication-quality volcano plot.

    Args:
        df_all: All genes (filtered for NaN in padj, log2fc, basemean)
        df_sig: Significant DEGs with "direction" column (up/down)
        label: Condition label for title
        output_path: Output file path (Path object)
        log2fc_cutoff: log2 fold change threshold for cutoff lines
        padj_cutoff: Adjusted p-value threshold for cutoff line
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # All points (not significant)
    ns = df_all[~df_all.index.isin(df_sig.index)]
    ax.scatter(ns["log2fc"], -np.log10(ns["padj"].clip(lower=1e-300)),
               c=COLOR_NS, s=4, alpha=0.3, rasterized=True, label="Not significant")

    # Significant up
    up = df_sig[df_sig["direction"] == "up"]
    ax.scatter(up["log2fc"], -np.log10(up["padj"].clip(lower=1e-300)),
               c=COLOR_UP, s=8, alpha=0.6, label=f"Up ({len(up):,})")

    # Significant down
    down = df_sig[df_sig["direction"] == "down"]
    ax.scatter(down["log2fc"], -np.log10(down["padj"].clip(lower=1e-300)),
               c=COLOR_DOWN, s=8, alpha=0.6, label=f"Down ({len(down):,})")

    # Cutoff lines
    ax.axhline(-np.log10(padj_cutoff), ls="--", lw=0.8, c="grey", alpha=0.5)
    ax.axvline(log2fc_cutoff, ls="--", lw=0.8, c="grey", alpha=0.5)
    ax.axvline(-log2fc_cutoff, ls="--", lw=0.8, c="grey", alpha=0.5)

    # Label top genes
    if len(df_sig) > 0 and "gene_name" in df_sig.columns:
        top_genes = df_sig.nlargest(10, "log2fc", keep="first")
        bottom_genes = df_sig.nsmallest(10, "log2fc", keep="first")
        to_label = pd.concat([top_genes, bottom_genes]).drop_duplicates()

        for _, row in to_label.iterrows():
            name = row.get("gene_name", row.get("gene_id", ""))
            if pd.notna(name) and not str(name).startswith("ENS"):
                ax.annotate(str(name),
                            (row["log2fc"], -np.log10(max(row["padj"], 1e-300))),
                            fontsize=6, alpha=0.8, ha="center",
                            textcoords="offset points", xytext=(0, 5))

    ax.set_xlabel("log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-log₁₀(padj)", fontsize=12)
    ax.set_title(f"Differential Expression: {label}", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_volcano_by_biotype(df_all, df_sig, label, output_path, biotype_filter,
                             log2fc_cutoff=0.4, padj_cutoff=0.01):
    """
    Generate volcano plot for specific biotype subset.

    Args:
        df_all: All genes (filtered for NaN)
        df_sig: Significant DEGs
        label: Condition label
        output_path: Base output path (will be modified with biotype suffix)
        biotype_filter: Single biotype string or list of biotypes to include
        log2fc_cutoff: log2 fold change threshold for cutoff lines
        padj_cutoff: Adjusted p-value threshold for cutoff line
    """
    if isinstance(biotype_filter, str):
        biotype_filter = [biotype_filter]

    if "biotype_group" not in df_all.columns or "biotype_group" not in df_sig.columns:
        print(f"  Warning: No biotype_group column, skipping biotype-specific volcano")
        return

    # Filter data by biotype
    df_all_filtered = df_all[df_all["biotype_group"].isin(biotype_filter)]
    df_sig_filtered = df_sig[df_sig["biotype_group"].isin(biotype_filter)]

    if len(df_sig_filtered) == 0:
        print(f"  No significant {biotype_filter} genes to plot")
        return

    # Create filename with biotype suffix
    biotype_str = "_".join([b.replace(" ", "").replace("ncRNA", "ncRNA") for b in biotype_filter])
    base_name = output_path.stem.replace("_volcano", "")
    new_filename = f"{base_name}_volcano_{biotype_str}.png"
    full_path = output_path.parent / new_filename

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # All points (not significant)
    ns = df_all_filtered[~df_all_filtered.index.isin(df_sig_filtered.index)]
    if len(ns) > 0:
        ax.scatter(ns["log2fc"], -np.log10(ns["padj"].clip(lower=1e-300)),
                   c=COLOR_NS, s=4, alpha=0.3, rasterized=True, label="Not significant")

    # Significant up
    up = df_sig_filtered[df_sig_filtered["direction"] == "up"]
    if len(up) > 0:
        ax.scatter(up["log2fc"], -np.log10(up["padj"].clip(lower=1e-300)),
                   c=COLOR_UP, s=8, alpha=0.6, label=f"Up ({len(up):,})")

    # Significant down
    down = df_sig_filtered[df_sig_filtered["direction"] == "down"]
    if len(down) > 0:
        ax.scatter(down["log2fc"], -np.log10(down["padj"].clip(lower=1e-300)),
                   c=COLOR_DOWN, s=8, alpha=0.6, label=f"Down ({len(down):,})")

    # Cutoff lines
    ax.axhline(-np.log10(padj_cutoff), ls="--", lw=0.8, c="grey", alpha=0.5)
    ax.axvline(log2fc_cutoff, ls="--", lw=0.8, c="grey", alpha=0.5)
    ax.axvline(-log2fc_cutoff, ls="--", lw=0.8, c="grey", alpha=0.5)

    # Label top genes
    if len(df_sig_filtered) > 0 and "gene_name" in df_sig_filtered.columns:
        n_to_label = min(10, len(df_sig_filtered))
        top_genes = df_sig_filtered.nlargest(n_to_label, "log2fc", keep="first")
        bottom_genes = df_sig_filtered.nsmallest(n_to_label, "log2fc", keep="first")
        to_label = pd.concat([top_genes, bottom_genes]).drop_duplicates()

        for _, row in to_label.iterrows():
            name = row.get("gene_name", row.get("gene_id", ""))
            if pd.notna(name) and not str(name).startswith("ENS"):
                ax.annotate(str(name),
                            (row["log2fc"], -np.log10(max(row["padj"], 1e-300))),
                            fontsize=6, alpha=0.8, ha="center",
                            textcoords="offset points", xytext=(0, 5))

    biotype_label = " + ".join(biotype_filter)
    ax.set_xlabel("log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-log₁₀(padj)", fontsize=12)
    ax.set_title(f"Differential Expression: {label}\n({biotype_label})",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(full_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {full_path.name}")


# ─── MA Plot ─────────────────────────────────────────────────────────────────

def plot_ma(df_all, df_sig, label, output_path, log2fc_cutoff=0.4):
    """
    Generate MA plot (log2FC vs baseMean).

    Args:
        df_all: All genes with log2fc and basemean columns
        df_sig: Significant DEGs with direction column
        label: Condition label for title
        output_path: Output file path
        log2fc_cutoff: log2 fold change threshold for cutoff lines
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    ns = df_all[~df_all.index.isin(df_sig.index)]
    ax.scatter(np.log10(ns["basemean"].clip(lower=0.01)),
               ns["log2fc"], c=COLOR_NS, s=3, alpha=0.2, rasterized=True)

    up = df_sig[df_sig["direction"] == "up"]
    down = df_sig[df_sig["direction"] == "down"]
    ax.scatter(np.log10(up["basemean"].clip(lower=0.01)),
               up["log2fc"], c=COLOR_UP, s=6, alpha=0.5, label=f"Up ({len(up):,})")
    ax.scatter(np.log10(down["basemean"].clip(lower=0.01)),
               down["log2fc"], c=COLOR_DOWN, s=6, alpha=0.5, label=f"Down ({len(down):,})")

    ax.axhline(0, ls="-", lw=0.5, c="black", alpha=0.3)
    ax.axhline(log2fc_cutoff, ls="--", lw=0.5, c="grey", alpha=0.5)
    ax.axhline(-log2fc_cutoff, ls="--", lw=0.5, c="grey", alpha=0.5)

    ax.set_xlabel("log₁₀(baseMean)", fontsize=12)
    ax.set_ylabel("log₂ Fold Change", fontsize=12)
    ax.set_title(f"MA Plot: {label}", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# ─── Biotype Distribution ───────────────────────────────────────────────────

def plot_biotype_distribution(df_sig, label, output_path):
    """
    Bar chart of DEGs by biotype group.

    Args:
        df_sig: Significant DEGs with biotype_group column
        label: Condition label for title
        output_path: Output file path
    """
    if "biotype_group" not in df_sig.columns or len(df_sig) == 0:
        return

    counts = df_sig["biotype_group"].value_counts()
    # Reorder
    order = ["Protein Coding", "lncRNA", "Pseudogene", "Small ncRNA", "Other"]
    counts = counts.reindex([c for c in order if c in counts.index])

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = [BIOTYPE_COLORS.get(b, "#999999") for b in counts.index]
    bars = ax.bar(counts.index, counts.values, color=colors, edgecolor="white", linewidth=0.5)

    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(val), ha="center", va="bottom", fontsize=10)

    ax.set_ylabel("Number of DEGs", fontsize=12)
    ax.set_title(f"DEGs by Biotype: {label}", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_biotype_direction(df_sig, label, output_path):
    """
    Stacked bar chart: biotype x direction.

    Args:
        df_sig: Significant DEGs with biotype_group and direction columns
        label: Condition label for title
        output_path: Output file path
    """
    if len(df_sig) == 0:
        return

    order = ["Protein Coding", "lncRNA", "Pseudogene", "Small ncRNA", "Other"]
    ct = pd.crosstab(df_sig["biotype_group"], df_sig["direction"])
    ct = ct.reindex([c for c in order if c in ct.index])

    fig, ax = plt.subplots(figsize=(8, 5))
    ct_plot = ct.reindex(columns=["up", "down"], fill_value=0)
    ct_plot.plot(kind="bar", stacked=True, color=[COLOR_UP, COLOR_DOWN],
                 ax=ax, edgecolor="white", linewidth=0.5)

    ax.set_ylabel("Number of DEGs", fontsize=12)
    ax.set_title(f"DEGs by Biotype & Direction: {label}", fontsize=13, fontweight="bold")
    ax.legend(["Upregulated", "Downregulated"], fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xticks(rotation=30, ha="right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# ─── rMATS Splicing Plots ───────────────────────────────────────────────────

def plot_splicing_summary(rmats_df, label, output_path, event_types=None):
    """
    Bar chart of significant splicing events by type.

    Args:
        rmats_df: rMATS results with event_type column
        label: Condition label for title
        output_path: Output file path
        event_types: List of event types in order (default: SE, A3SS, A5SS, RI, MXE)
    """
    if len(rmats_df) == 0:
        return

    if event_types is None:
        event_types = ["SE", "A3SS", "A5SS", "RI", "MXE"]

    counts = rmats_df["event_type"].value_counts()
    counts = counts.reindex([e for e in event_types if e in counts.index])

    fig, ax = plt.subplots(figsize=(7, 5))
    colors = [EVENT_COLORS.get(e, "#999") for e in counts.index]
    bars = ax.bar(counts.index, counts.values, color=colors, edgecolor="white")

    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(val), ha="center", va="bottom", fontsize=10)

    ax.set_ylabel("Significant Events", fontsize=12)
    ax.set_title(f"Alternative Splicing Events: {label}", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_splicing_volcano(rmats_df, label, output_path, use_pval=True,
                          pvalue_cutoff=0.01, fdr_cutoff=0.05, inclevel_cutoff=0.1,
                          event_types=None):
    """
    Volcano plot for splicing: IncLevelDifference vs -log10(FDR or PValue).

    Args:
        rmats_df: rMATS results with IncLevelDifference, PValue, FDR, event_type columns
        label: Condition label for title
        output_path: Output file path
        use_pval: Use p-value instead of FDR for y-axis (default: True)
        pvalue_cutoff: P-value threshold for cutoff line
        fdr_cutoff: FDR threshold for cutoff line
        inclevel_cutoff: IncLevelDifference threshold for cutoff lines
        event_types: List of event types (default: SE, A3SS, A5SS, RI, MXE)
    """
    if len(rmats_df) == 0:
        return

    if event_types is None:
        event_types = ["SE", "A3SS", "A5SS", "RI", "MXE"]

    # Choose metric based on configuration
    y_col = "PValue" if use_pval else "FDR"
    y_label = "-log₁₀(p-value)" if use_pval else "-log₁₀(FDR)"
    cutoff_val = pvalue_cutoff if use_pval else fdr_cutoff

    fig, ax = plt.subplots(figsize=(8, 6))

    for event_type in event_types:
        sub = rmats_df[rmats_df["event_type"] == event_type]
        if len(sub) == 0:
            continue
        ax.scatter(sub["IncLevelDifference"],
                   -np.log10(sub[y_col].clip(lower=1e-300)),
                   c=EVENT_COLORS.get(event_type, "#999"),
                   s=15, alpha=0.6, label=f"{event_type} ({len(sub)})")

    ax.axhline(-np.log10(cutoff_val), ls="--", lw=0.8, c="grey", alpha=0.5)
    ax.axvline(inclevel_cutoff, ls="--", lw=0.8, c="grey", alpha=0.5)
    ax.axvline(-inclevel_cutoff, ls="--", lw=0.8, c="grey", alpha=0.5)

    ax.set_xlabel("ΔΨ (IncLevelDifference)", fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(f"Splicing Events: {label}", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# ─── Cross-Condition Comparisons ────────────────────────────────────────────

def plot_overlap_comparison(all_filtered, output_path, conditions):
    """
    Compare DEG overlap across conditions.

    Args:
        all_filtered: Dictionary of {condition_name: filtered_df}
        output_path: Output file path
        conditions: List of condition dictionaries with "name" and "label" keys
    """
    if len(all_filtered) < 2:
        return

    names = list(all_filtered.keys())
    gene_sets = {}
    for name, df in all_filtered.items():
        if "gene_id_base" in df.columns:
            gene_sets[name] = set(df["gene_id_base"])
        elif "gene_id" in df.columns:
            gene_sets[name] = set(df["gene_id"])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: bar chart of DEG counts
    ax = axes[0]
    counts = {n: len(s) for n, s in gene_sets.items()}
    labels = [next(c["label"] for c in conditions if c["name"] == n) for n in names]
    bars = ax.bar(labels, list(counts.values()),
                  color=[COLOR_UP, COLOR_DOWN, "#55A868"][:len(names)],
                  edgecolor="white")
    for bar, val in zip(bars, counts.values()):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                str(val), ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("Number of DEGs", fontsize=12)
    ax.set_title("DEG Counts per Condition", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.sca(ax)
    plt.xticks(rotation=20, ha="right")

    # Right: pairwise overlaps as heatmap
    ax2 = axes[1]
    n = len(names)
    overlap_matrix = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            overlap_matrix[i, j] = len(gene_sets[names[i]] & gene_sets[names[j]])

    short_labels = [next(c["label"] for c in conditions if c["name"] == names[i]).split(" vs ")[0]
                    for i in range(n)]
    im = ax2.imshow(overlap_matrix, cmap="YlOrRd", aspect="auto")
    ax2.set_xticks(range(n))
    ax2.set_yticks(range(n))
    ax2.set_xticklabels(short_labels, rotation=30, ha="right")
    ax2.set_yticklabels(short_labels)

    for i in range(n):
        for j in range(n):
            ax2.text(j, i, str(overlap_matrix[i, j]),
                     ha="center", va="center", fontsize=11, fontweight="bold")

    ax2.set_title("DEG Overlap Matrix", fontsize=13, fontweight="bold")
    plt.colorbar(im, ax=ax2, shrink=0.8)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_concordance_scatter(filtered_dict, cond_a, cond_b, output_path, conditions):
    """
    Scatter plot of log2FC for shared DEGs between two conditions.

    Args:
        filtered_dict: Dictionary of {condition_name: filtered_df}
        cond_a: First condition name
        cond_b: Second condition name
        output_path: Output file path
        conditions: List of condition dictionaries with "name" and "label" keys
    """
    df_a = filtered_dict.get(cond_a)
    df_b = filtered_dict.get(cond_b)
    if df_a is None or df_b is None:
        return

    # Merge on gene_id_base
    col_a = "gene_id_base" if "gene_id_base" in df_a.columns else "gene_id"
    col_b = "gene_id_base" if "gene_id_base" in df_b.columns else "gene_id"

    merged = df_a[[col_a, "log2fc", "gene_name"]].merge(
        df_b[[col_b, "log2fc"]],
        left_on=col_a, right_on=col_b, suffixes=("_a", "_b")
    )

    if len(merged) == 0:
        print(f"  No shared DEGs between {cond_a} and {cond_b}")
        return

    label_a = next(c["label"] for c in conditions if c["name"] == cond_a)
    label_b = next(c["label"] for c in conditions if c["name"] == cond_b)

    fig, ax = plt.subplots(figsize=(7, 7))

    # Color by concordance
    concordant = (merged["log2fc_a"] * merged["log2fc_b"]) > 0
    ax.scatter(merged.loc[concordant, "log2fc_a"],
               merged.loc[concordant, "log2fc_b"],
               c="#55A868", s=20, alpha=0.6, label=f"Concordant ({concordant.sum()})")
    ax.scatter(merged.loc[~concordant, "log2fc_a"],
               merged.loc[~concordant, "log2fc_b"],
               c="#C44E52", s=20, alpha=0.6, label=f"Discordant ({(~concordant).sum()})")

    # Identity line
    lim = max(abs(merged["log2fc_a"]).max(), abs(merged["log2fc_b"]).max()) + 0.5
    ax.plot([-lim, lim], [-lim, lim], ls="--", c="grey", alpha=0.5, lw=0.8)
    ax.axhline(0, ls="-", c="grey", alpha=0.2, lw=0.5)
    ax.axvline(0, ls="-", c="grey", alpha=0.2, lw=0.5)

    # Label top discordant and concordant genes
    if "gene_name" in merged.columns:
        # Label genes with largest absolute fold changes
        merged["max_abs_fc"] = merged[["log2fc_a", "log2fc_b"]].abs().max(axis=1)
        top = merged.nlargest(15, "max_abs_fc")
        for _, row in top.iterrows():
            name = row.get("gene_name", "")
            if pd.notna(name) and not str(name).startswith("ENS"):
                ax.annotate(str(name), (row["log2fc_a"], row["log2fc_b"]),
                            fontsize=6, alpha=0.7, ha="center",
                            textcoords="offset points", xytext=(0, 5))

    # Stats
    r = np.corrcoef(merged["log2fc_a"], merged["log2fc_b"])[0, 1]
    n_conc = concordant.sum()
    n_disc = (~concordant).sum()
    pct_conc = 100 * n_conc / len(merged) if len(merged) > 0 else 0

    ax.set_xlabel(f"log₂FC ({label_a})", fontsize=12)
    ax.set_ylabel(f"log₂FC ({label_b})", fontsize=12)
    ax.set_title(f"Fold Change Concordance\n{label_a} vs {label_b}\n"
                 f"r = {r:.3f}, {n_conc} concordant ({pct_conc:.0f}%), "
                 f"{n_disc} discordant",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_concordance_all_pairs(all_data, filtered_dict, output_path, conditions):
    """
    Scatter of log2FC for ALL shared genes (not just DEGs) between each pair.

    Args:
        all_data: Dictionary of {condition_name: all_genes_df}
        filtered_dict: Dictionary of {condition_name: filtered_df}
        output_path: Output file path
        conditions: List of condition dictionaries with "name" and "label" keys
    """
    names = list(all_data.keys())
    n_pairs = len(names) * (len(names) - 1) // 2
    if n_pairs == 0:
        return

    fig, axes = plt.subplots(1, n_pairs, figsize=(7 * n_pairs, 6))
    if n_pairs == 1:
        axes = [axes]

    pair_idx = 0
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            ax = axes[pair_idx]
            df_a = all_data[names[i]]
            df_b = all_data[names[j]]

            col_a = "gene_id_base" if "gene_id_base" in df_a.columns else "gene_id"
            col_b = "gene_id_base" if "gene_id_base" in df_b.columns else "gene_id"

            merged = df_a[[col_a, "log2fc"]].merge(
                df_b[[col_b, "log2fc"]],
                left_on=col_a, right_on=col_b, suffixes=("_a", "_b")
            )

            label_a = next(c["label"] for c in conditions if c["name"] == names[i])
            label_b = next(c["label"] for c in conditions if c["name"] == names[j])

            ax.scatter(merged["log2fc_a"], merged["log2fc_b"],
                       c=COLOR_NS, s=2, alpha=0.1, rasterized=True)

            # Highlight DEGs
            sig_a = filtered_dict.get(names[i])
            sig_b = filtered_dict.get(names[j])
            if sig_a is not None and sig_b is not None:
                sig_col_a = "gene_id_base" if "gene_id_base" in sig_a.columns else "gene_id"
                sig_col_b = "gene_id_base" if "gene_id_base" in sig_b.columns else "gene_id"
                shared_degs = set(sig_a[sig_col_a]) & set(sig_b[sig_col_b])
                highlight = merged[merged[col_a].isin(shared_degs)]
                if len(highlight) > 0:
                    conc = (highlight["log2fc_a"] * highlight["log2fc_b"]) > 0
                    ax.scatter(highlight.loc[conc, "log2fc_a"],
                               highlight.loc[conc, "log2fc_b"],
                               c="#55A868", s=10, alpha=0.6)
                    ax.scatter(highlight.loc[~conc, "log2fc_a"],
                               highlight.loc[~conc, "log2fc_b"],
                               c="#C44E52", s=10, alpha=0.6)

            r = np.corrcoef(merged["log2fc_a"], merged["log2fc_b"])[0, 1]
            ax.plot([-5, 5], [-5, 5], ls="--", c="grey", alpha=0.4, lw=0.8)
            ax.axhline(0, ls="-", c="grey", alpha=0.2)
            ax.axvline(0, ls="-", c="grey", alpha=0.2)
            ax.set_xlabel(f"log₂FC ({label_a.split(' vs ')[0]})", fontsize=10)
            ax.set_ylabel(f"log₂FC ({label_b.split(' vs ')[0]})", fontsize=10)
            ax.set_title(f"r = {r:.3f}", fontsize=11, fontweight="bold")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            pair_idx += 1

    plt.suptitle("Global Fold Change Correlation (All Genes)", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# ─── Directional Venn Diagrams ──────────────────────────────────────────────

def plot_directional_venn_3panel(df_a, df_b, label_a, label_b, output_path):
    """
    Three-panel Venn diagram showing directional overlap.

    Panel A: UP in both conditions
    Panel B: DOWN in both conditions
    Panel C: All significant (any direction)

    Args:
        df_a: First condition DEGs with direction column
        df_b: Second condition DEGs with direction column
        label_a: First condition label
        label_b: Second condition label
        output_path: Output file path
    """
    if not HAS_VENN:
        print("  Warning: matplotlib-venn not installed, skipping Venn diagrams")
        return

    col_a = "gene_id_base" if "gene_id_base" in df_a.columns else "gene_id"
    col_b = "gene_id_base" if "gene_id_base" in df_b.columns else "gene_id"

    # All significant genes
    genes_a_all = set(df_a[col_a])
    genes_b_all = set(df_b[col_b])

    # UP in each condition
    genes_a_up = set(df_a[df_a["direction"] == "up"][col_a])
    genes_b_up = set(df_b[df_b["direction"] == "up"][col_b])

    # DOWN in each condition
    genes_a_down = set(df_a[df_a["direction"] == "down"][col_a])
    genes_b_down = set(df_b[df_b["direction"] == "down"][col_b])

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel A: Concordant UP
    ax = axes[0]
    v = venn2([genes_a_up, genes_b_up], set_labels=(label_a, label_b), ax=ax)
    ax.set_title(f"Upregulated in Both\n({len(genes_a_up & genes_b_up)} shared)",
                 fontsize=13, fontweight="bold")
    if v:
        for text in v.set_labels:
            if text:
                text.set_fontsize(10)
        for text in v.subset_labels:
            if text:
                text.set_fontsize(11)

    # Panel B: Concordant DOWN
    ax = axes[1]
    v = venn2([genes_a_down, genes_b_down], set_labels=(label_a, label_b), ax=ax)
    ax.set_title(f"Downregulated in Both\n({len(genes_a_down & genes_b_down)} shared)",
                 fontsize=13, fontweight="bold")
    if v:
        for text in v.set_labels:
            if text:
                text.set_fontsize(10)
        for text in v.subset_labels:
            if text:
                text.set_fontsize(11)

    # Panel C: All significant
    ax = axes[2]
    v = venn2([genes_a_all, genes_b_all], set_labels=(label_a, label_b), ax=ax)
    ax.set_title(f"All Significant DEGs\n({len(genes_a_all & genes_b_all)} shared)",
                 fontsize=13, fontweight="bold")
    if v:
        for text in v.set_labels:
            if text:
                text.set_fontsize(10)
        for text in v.subset_labels:
            if text:
                text.set_fontsize(11)

    plt.suptitle(f"Directional Overlap: {label_a} vs {label_b}",
                 fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_splicing_venn_3panel(rmats_a, rmats_b, label_a, label_b, event_type, output_path):
    """
    Three-panel Venn for alternative splicing events.

    Panel A: ALL overlapping events (any change in both)
    Panel B: Both DOWN/SKIPPED (dPSI < 0 in both)
    Panel C: Both UP/INCLUDED (dPSI > 0 in both)

    Args:
        rmats_a: First condition rMATS results
        rmats_b: Second condition rMATS results
        label_a: First condition label
        label_b: Second condition label
        event_type: Splicing event type (SE, A3SS, etc.)
        output_path: Output file path
    """
    if not HAS_VENN:
        print("  Warning: matplotlib-venn not installed, skipping Venn diagrams")
        return

    # Filter by event type
    events_a = rmats_a[rmats_a["event_type"] == event_type]
    events_b = rmats_b[rmats_b["event_type"] == event_type]

    if len(events_a) == 0 or len(events_b) == 0:
        print(f"  No {event_type} events to compare")
        return

    # Use gene symbol as event identifier (or ID if available)
    id_col = "geneSymbol" if "geneSymbol" in events_a.columns else "ID"

    # All events
    all_a = set(events_a[id_col].dropna())
    all_b = set(events_b[id_col].dropna())

    # DOWN/SKIPPED (dPSI < 0)
    down_a = set(events_a[events_a["direction"] == "excluded"][id_col].dropna())
    down_b = set(events_b[events_b["direction"] == "excluded"][id_col].dropna())

    # UP/INCLUDED (dPSI > 0)
    up_a = set(events_a[events_a["direction"] == "included"][id_col].dropna())
    up_b = set(events_b[events_b["direction"] == "included"][id_col].dropna())

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel A: All events
    ax = axes[0]
    v = venn2([all_a, all_b], set_labels=(label_a, label_b), ax=ax)
    ax.set_title(f"All {event_type} Events\n({len(all_a & all_b)} shared)",
                 fontsize=13, fontweight="bold")
    if v:
        for text in v.set_labels:
            if text:
                text.set_fontsize(10)
        for text in v.subset_labels:
            if text:
                text.set_fontsize(11)

    # Panel B: DOWN/SKIPPED in both
    ax = axes[1]
    v = venn2([down_a, down_b], set_labels=(label_a, label_b), ax=ax)
    ax.set_title(f"Skipped/Excluded in Both\n({len(down_a & down_b)} shared)",
                 fontsize=13, fontweight="bold")
    if v:
        for text in v.set_labels:
            if text:
                text.set_fontsize(10)
        for text in v.subset_labels:
            if text:
                text.set_fontsize(11)

    # Panel C: UP/INCLUDED in both
    ax = axes[2]
    v = venn2([up_a, up_b], set_labels=(label_a, label_b), ax=ax)
    ax.set_title(f"Included in Both\n({len(up_a & up_b)} shared)",
                 fontsize=13, fontweight="bold")
    if v:
        for text in v.set_labels:
            if text:
                text.set_fontsize(10)
        for text in v.subset_labels:
            if text:
                text.set_fontsize(11)

    plt.suptitle(f"Splicing Overlap ({event_type}): {label_a} vs {label_b}",
                 fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# ─── Summary Dashboard ──────────────────────────────────────────────────────

def plot_summary_dashboard(summary_df, output_path):
    """
    Multi-panel summary dashboard.

    Args:
        summary_df: Summary DataFrame with columns: Condition, Upregulated, Downregulated,
                    Protein Coding DEGs, lncRNA DEGs, Splicing Events
        output_path: Output file path
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    labels = summary_df["Condition"].tolist()
    x = np.arange(len(labels))

    # Panel 1: DEG counts
    ax = axes[0, 0]
    ax.bar(x - 0.15, summary_df["Upregulated"], 0.3, color=COLOR_UP, label="Upregulated")
    ax.bar(x + 0.15, summary_df["Downregulated"], 0.3, color=COLOR_DOWN, label="Downregulated")
    ax.set_xticks(x)
    ax.set_xticklabels([l.split(" vs ")[0] for l in labels], rotation=20, ha="right")
    ax.set_ylabel("Number of DEGs")
    ax.set_title("Differential Expression", fontweight="bold")
    ax.legend()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel 2: Biotype composition
    ax = axes[0, 1]
    ax.bar(x - 0.15, summary_df["Protein Coding DEGs"], 0.3,
           color=BIOTYPE_COLORS["Protein Coding"], label="Protein Coding")
    ax.bar(x + 0.15, summary_df["lncRNA DEGs"], 0.3,
           color=BIOTYPE_COLORS["lncRNA"], label="lncRNA")
    ax.set_xticks(x)
    ax.set_xticklabels([l.split(" vs ")[0] for l in labels], rotation=20, ha="right")
    ax.set_ylabel("Number of DEGs")
    ax.set_title("DEGs by Biotype", fontweight="bold")
    ax.legend()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel 3: Splicing events
    ax = axes[1, 0]
    ax.bar(x, summary_df["Splicing Events"], color="#3C5488", edgecolor="white")
    for i, val in enumerate(summary_df["Splicing Events"]):
        ax.text(i, val + 1, str(val), ha="center", va="bottom", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels([l.split(" vs ")[0] for l in labels], rotation=20, ha="right")
    ax.set_ylabel("Significant Events")
    ax.set_title("Alternative Splicing (rMATS)", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Panel 4: Summary text table
    ax = axes[1, 1]
    ax.axis("off")
    cell_text = []
    for _, row in summary_df.iterrows():
        cell_text.append([
            row["Condition"].split(" vs ")[0],
            f"{row['DEGs']:,}",
            f"{row['Upregulated']:,}↑ / {row['Downregulated']:,}↓",
            f"{row['Splicing Events']:,}",
        ])
    table = ax.table(cellText=cell_text,
                     colLabels=["Condition", "DEGs", "Up/Down", "Splicing"],
                     cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    ax.set_title("Summary", fontweight="bold", pad=20)

    plt.suptitle("RNA-seq Analysis: MIAT OE vs QKI-KO vs polyQKI-KO",
                 fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# ─── Top DEGs Heatmap ───────────────────────────────────────────────────────

def plot_top_genes_heatmap(all_data, filtered_dict, output_path, conditions, n_top=50):
    """
    Heatmap of top DEGs across all conditions.

    Args:
        all_data: Dictionary of {condition_name: all_genes_df}
        filtered_dict: Dictionary of {condition_name: filtered_df}
        output_path: Output file path
        conditions: List of condition dictionaries with "name" and "label" keys
        n_top: Number of top genes to include (default: 50)
    """
    # Collect all significant gene IDs
    all_sig_genes = set()
    for name, df in filtered_dict.items():
        col = "gene_id_base" if "gene_id_base" in df.columns else "gene_id"
        all_sig_genes.update(df[col].tolist())

    if len(all_sig_genes) == 0:
        return

    # Build matrix of log2FC values
    names = list(all_data.keys())
    rows = []
    for gene in all_sig_genes:
        row = {"gene": gene}
        for name in names:
            df = all_data[name]
            col = "gene_id_base" if "gene_id_base" in df.columns else "gene_id"
            match = df[df[col] == gene]
            if len(match) > 0:
                row[name] = match.iloc[0]["log2fc"]
                if "gene_name" in match.columns:
                    gn = match.iloc[0]["gene_name"]
                    if pd.notna(gn) and not str(gn).startswith("ENS"):
                        row["gene_name"] = gn
            else:
                row[name] = np.nan
        rows.append(row)

    matrix = pd.DataFrame(rows)
    if "gene_name" not in matrix.columns:
        matrix["gene_name"] = matrix["gene"]

    # Select top genes by maximum absolute fold change
    matrix["max_abs"] = matrix[names].abs().max(axis=1)
    matrix = matrix.dropna(subset=names, how="all")
    top = matrix.nlargest(min(n_top, len(matrix)), "max_abs")

    # Plot
    plot_data = top.set_index("gene_name")[names]
    plot_data.columns = [next(c["label"] for c in conditions if c["name"] == n).split(" vs ")[0]
                         for n in names]

    fig, ax = plt.subplots(figsize=(6, max(8, len(plot_data) * 0.25)))

    if HAS_SEABORN:
        sns.heatmap(plot_data, cmap="RdBu_r", center=0, ax=ax,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={"label": "log₂FC"}, linewidths=0.2)
    else:
        im = ax.imshow(plot_data.values, cmap="RdBu_r", aspect="auto",
                       vmin=-plot_data.values[np.isfinite(plot_data.values)].max(),
                       vmax=plot_data.values[np.isfinite(plot_data.values)].max())
        ax.set_xticks(range(len(plot_data.columns)))
        ax.set_xticklabels(plot_data.columns, rotation=30, ha="right")
        ax.set_yticks(range(len(plot_data.index)))
        ax.set_yticklabels(plot_data.index, fontsize=7)
        plt.colorbar(im, ax=ax, label="log₂FC")

    ax.set_title(f"Top {len(plot_data)} DEGs Across Conditions",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# ─── Violin Plots ───────────────────────────────────────────────────────────

def plot_log2fc_violin(filtered_dict, output_path, conditions):
    """
    Violin plot showing log2FC distribution across all conditions.

    Args:
        filtered_dict: Dictionary of {condition_name: filtered_df}
        output_path: Output file path
        conditions: List of condition dictionaries with "name" and "label" keys
    """
    if not HAS_SEABORN:
        print("  Warning: seaborn not installed, skipping violin plots")
        return

    # Collect data
    data = []
    for cond in conditions:
        name = cond["name"]
        df = filtered_dict.get(name)
        if df is not None and len(df) > 0:
            for _, row in df.iterrows():
                data.append({
                    "Condition": cond["label"].split(" vs ")[0],
                    "log2FC": row["log2fc"],
                    "Direction": row["direction"]
                })

    if len(data) == 0:
        print("  No data for violin plot")
        return

    plot_df = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(10, 6))

    sns.violinplot(data=plot_df, x="Condition", y="log2FC",
                   hue="Direction", split=True,
                   palette={"up": COLOR_UP, "down": COLOR_DOWN},
                   ax=ax, inner="quartile")

    ax.axhline(0, ls="--", lw=0.8, c="grey", alpha=0.5)
    ax.set_ylabel("log₂ Fold Change", fontsize=12)
    ax.set_xlabel("Condition", fontsize=12)
    ax.set_title("Distribution of Fold Changes Across Conditions",
                 fontsize=13, fontweight="bold")
    ax.legend(title="Direction", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xticks(rotation=20, ha="right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_dpsi_violin(rmats_dict, output_path, conditions):
    """
    Violin plot showing dPSI distribution for splicing events.

    Args:
        rmats_dict: Dictionary of {condition_name: rmats_df}
        output_path: Output file path
        conditions: List of condition dictionaries with "name" and "label" keys
    """
    if not HAS_SEABORN:
        print("  Warning: seaborn not installed, skipping violin plots")
        return

    # Collect data
    data = []
    for cond in conditions:
        name = cond["name"]
        df = rmats_dict.get(name)
        if df is not None and len(df) > 0:
            for _, row in df.iterrows():
                data.append({
                    "Condition": cond["label"].split(" vs ")[0],
                    "dPSI": row["IncLevelDifference"],
                    "Event Type": row["event_type"],
                    "Direction": row["direction"]
                })

    if len(data) == 0:
        print("  No splicing data for violin plot")
        return

    plot_df = pd.DataFrame(data)

    fig, ax = plt.subplots(figsize=(12, 6))

    sns.violinplot(data=plot_df, x="Event Type", y="dPSI",
                   hue="Condition",
                   palette=[COLOR_UP, COLOR_DOWN, "#55A868"][:len(conditions)],
                   ax=ax, inner="box")

    ax.axhline(0, ls="--", lw=0.8, c="grey", alpha=0.5)
    ax.set_ylabel("ΔΨ (IncLevelDifference)", fontsize=12)
    ax.set_xlabel("Event Type", fontsize=12)
    ax.set_title("Distribution of Splicing Changes Across Conditions",
                 fontsize=13, fontweight="bold")
    ax.legend(title="Condition", fontsize=9, bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")
