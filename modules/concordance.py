"""
Concordance and Cross-Condition Analysis Module
================================================

Functions for analyzing DEG overlap and concordance across multiple conditions:
- Cross-condition overlap visualization
- Concordance calculations (directional agreement)
- Directional splitting (up/down DEGs)
- Set operations (intersections, unions)
- Excel export for cross-condition comparisons
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

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

# ─── Colors ──────────────────────────────────────────────────────────────────

COLOR_UP = "#D72638"
COLOR_DOWN = "#1B98E0"
COLOR_NS = "#BFBFBF"
COLOR_CONCORDANT = "#55A868"
COLOR_DISCORDANT = "#C44E52"


# ─── Cross-Condition Overlap ─────────────────────────────────────────────────

def plot_overlap_comparison(all_filtered, output_path, conditions_config):
    """
    Compare DEG overlap across conditions.

    Args:
        all_filtered: Dict of condition_name -> filtered DEG DataFrame
        output_path: Path to save figure
        conditions_config: List of condition config dicts (with 'label' key)
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
    labels = [conditions_config[i]["label"] for i, n in enumerate(names)]
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

    short_labels = [conditions_config[i]["label"].split(" vs ")[0] for i in range(n)]
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


# ─── Concordance Scatter Plots ───────────────────────────────────────────────

def plot_concordance_scatter(filtered_dict, cond_a, cond_b, output_path, conditions_config):
    """
    Scatter plot of log2FC for shared DEGs between two conditions.

    Args:
        filtered_dict: Dict of condition_name -> filtered DEG DataFrame
        cond_a: Name of first condition
        cond_b: Name of second condition
        output_path: Path to save figure
        conditions_config: List of condition config dicts (with 'label', 'name' keys)
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

    label_a = next(c["label"] for c in conditions_config if c["name"] == cond_a)
    label_b = next(c["label"] for c in conditions_config if c["name"] == cond_b)

    fig, ax = plt.subplots(figsize=(7, 7))

    # Color by concordance
    concordant = (merged["log2fc_a"] * merged["log2fc_b"]) > 0
    ax.scatter(merged.loc[concordant, "log2fc_a"],
               merged.loc[concordant, "log2fc_b"],
               c=COLOR_CONCORDANT, s=20, alpha=0.6, label=f"Concordant ({concordant.sum()})")
    ax.scatter(merged.loc[~concordant, "log2fc_a"],
               merged.loc[~concordant, "log2fc_b"],
               c=COLOR_DISCORDANT, s=20, alpha=0.6, label=f"Discordant ({(~concordant).sum()})")

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


def plot_concordance_all_pairs(all_data, filtered_dict, output_path, conditions_config):
    """
    Scatter of log2FC for ALL shared genes (not just DEGs) between each pair.

    Args:
        all_data: Dict of condition_name -> full DataFrame (all genes)
        filtered_dict: Dict of condition_name -> filtered DEG DataFrame
        output_path: Path to save figure
        conditions_config: List of condition config dicts
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

            label_a = conditions_config[i]["label"]
            label_b = conditions_config[j]["label"]

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
                               c=COLOR_CONCORDANT, s=10, alpha=0.6)
                    ax.scatter(highlight.loc[~conc, "log2fc_a"],
                               highlight.loc[~conc, "log2fc_b"],
                               c=COLOR_DISCORDANT, s=10, alpha=0.6)

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


# ─── Directional Venn Diagrams ───────────────────────────────────────────────

def plot_directional_venn_3panel(df_a, df_b, label_a, label_b, output_path):
    """
    Three-panel Venn diagram showing directional overlap.

    Panel A: UP in both conditions
    Panel B: DOWN in both conditions
    Panel C: All significant (any direction)

    Args:
        df_a: Filtered DEG DataFrame for condition A
        df_b: Filtered DEG DataFrame for condition B
        label_a: Label for condition A
        label_b: Label for condition B
        output_path: Path to save figure
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
        rmats_a: rMATS DataFrame for condition A
        rmats_b: rMATS DataFrame for condition B
        label_a: Label for condition A
        label_b: Label for condition B
        event_type: Event type to filter by (SE, A3SS, A5SS, RI, MXE)
        output_path: Path to save figure
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


# ─── Excel Export Functions ──────────────────────────────────────────────────

def export_directional_overlap_excel(df_a, df_b, cond_a, cond_b, output_path):
    """
    Create Excel workbook with directional overlap sheets.

    Sheets:
    - Summary: Counts per category
    - Concordant_UP: Genes up in both
    - Concordant_DOWN: Genes down in both
    - Discordant: Genes in opposite directions
    - Only_A: Genes only in condition A
    - Only_B: Genes only in condition B

    Args:
        df_a: Filtered DEG DataFrame for condition A
        df_b: Filtered DEG DataFrame for condition B
        cond_a: Name of condition A
        cond_b: Name of condition B
        output_path: Path to save Excel file
    """
    col_a = "gene_id_base" if "gene_id_base" in df_a.columns else "gene_id"
    col_b = "gene_id_base" if "gene_id_base" in df_b.columns else "gene_id"

    genes_a = set(df_a[col_a])
    genes_b = set(df_b[col_b])
    shared = genes_a & genes_b
    only_a = genes_a - genes_b
    only_b = genes_b - genes_a

    # Categorize shared genes
    concordant_up = []
    concordant_down = []
    discordant = []

    for gene in shared:
        row_a = df_a[df_a[col_a] == gene].iloc[0]
        row_b = df_b[df_b[col_b] == gene].iloc[0]

        dir_a = row_a["direction"]
        dir_b = row_b["direction"]

        entry = {
            "gene_id": gene,
            "gene_name": row_a.get("gene_name", gene),
            f"log2FC_{cond_a}": row_a["log2fc"],
            f"padj_{cond_a}": row_a["padj"],
            f"log2FC_{cond_b}": row_b["log2fc"],
            f"padj_{cond_b}": row_b["padj"],
            "biotype": row_a.get("biotype_group", "unknown")
        }

        if dir_a == "up" and dir_b == "up":
            concordant_up.append(entry)
        elif dir_a == "down" and dir_b == "down":
            concordant_down.append(entry)
        else:
            discordant.append(entry)

    # Only A genes
    only_a_list = []
    for gene in only_a:
        row = df_a[df_a[col_a] == gene].iloc[0]
        only_a_list.append({
            "gene_id": gene,
            "gene_name": row.get("gene_name", gene),
            "log2FC": row["log2fc"],
            "padj": row["padj"],
            "direction": row["direction"],
            "biotype": row.get("biotype_group", "unknown")
        })

    # Only B genes
    only_b_list = []
    for gene in only_b:
        row = df_b[df_b[col_b] == gene].iloc[0]
        only_b_list.append({
            "gene_id": gene,
            "gene_name": row.get("gene_name", gene),
            "log2FC": row["log2fc"],
            "padj": row["padj"],
            "direction": row["direction"],
            "biotype": row.get("biotype_group", "unknown")
        })

    # Create Excel workbook
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # Summary sheet
        summary = pd.DataFrame([
            {"Category": "Concordant UP", "Count": len(concordant_up)},
            {"Category": "Concordant DOWN", "Count": len(concordant_down)},
            {"Category": "Discordant", "Count": len(discordant)},
            {"Category": f"Only {cond_a}", "Count": len(only_a_list)},
            {"Category": f"Only {cond_b}", "Count": len(only_b_list)},
            {"Category": "Total Shared", "Count": len(shared)},
        ])
        summary.to_excel(writer, sheet_name="Summary", index=False)

        # Data sheets
        if concordant_up:
            pd.DataFrame(concordant_up).to_excel(writer, sheet_name="Concordant_UP", index=False)
        if concordant_down:
            pd.DataFrame(concordant_down).to_excel(writer, sheet_name="Concordant_DOWN", index=False)
        if discordant:
            pd.DataFrame(discordant).to_excel(writer, sheet_name="Discordant", index=False)
        if only_a_list:
            pd.DataFrame(only_a_list).to_excel(writer, sheet_name=f"Only_{cond_a}", index=False)
        if only_b_list:
            pd.DataFrame(only_b_list).to_excel(writer, sheet_name=f"Only_{cond_b}", index=False)

    print(f"  Exported directional overlap: {output_path.name}")
    print(f"    Concordant UP: {len(concordant_up)}, DOWN: {len(concordant_down)}, Discordant: {len(discordant)}")


def export_shared_degs_excel(df_a, df_b, cond_a, cond_b, output_path):
    """
    Export shared DEGs with concordance information to Excel.

    Args:
        df_a: Filtered DEG DataFrame for condition A
        df_b: Filtered DEG DataFrame for condition B
        cond_a: Name of condition A
        cond_b: Name of condition B
        output_path: Path to save Excel file

    Returns:
        DataFrame of shared genes with concordance info
    """
    col_a = "gene_id_base" if "gene_id_base" in df_a.columns else "gene_id"
    col_b = "gene_id_base" if "gene_id_base" in df_b.columns else "gene_id"

    shared = set(df_a[col_a]) & set(df_b[col_b])
    if len(shared) == 0:
        print(f"  No shared DEGs between {cond_a} and {cond_b}")
        return pd.DataFrame()

    # Get details for shared genes
    shared_rows = []
    for gene in shared:
        row_a = df_a[df_a[col_a] == gene].iloc[0]
        row_b = df_b[df_b[col_b] == gene].iloc[0]

        gn = row_a.get("gene_name", gene)
        if pd.isna(gn) or str(gn).startswith("ENS"):
            gn = row_b.get("gene_name", gene)

        fc_a = row_a["log2fc"]
        fc_b = row_b["log2fc"]
        concordant = (fc_a > 0 and fc_b > 0) or (fc_a < 0 and fc_b < 0)

        shared_rows.append({
            "gene_id": gene,
            "gene_name": gn,
            f"log2FC_{cond_a}": fc_a,
            f"log2FC_{cond_b}": fc_b,
            "concordance": "Concordant" if concordant else "Discordant",
            "direction_a": "up" if fc_a > 0 else "down",
            "direction_b": "up" if fc_b > 0 else "down",
        })

    overlap_df = pd.DataFrame(shared_rows).sort_values(
        f"log2FC_{cond_a}", key=abs, ascending=False
    )
    overlap_df.to_excel(output_path, index=False)

    return overlap_df


# ─── Concordance Metrics ─────────────────────────────────────────────────────

def calculate_concordance_metrics(df_a, df_b):
    """
    Calculate concordance metrics between two DEG datasets.

    Args:
        df_a: Filtered DEG DataFrame for condition A
        df_b: Filtered DEG DataFrame for condition B

    Returns:
        Dict with concordance metrics
    """
    col_a = "gene_id_base" if "gene_id_base" in df_a.columns else "gene_id"
    col_b = "gene_id_base" if "gene_id_base" in df_b.columns else "gene_id"

    genes_a = set(df_a[col_a])
    genes_b = set(df_b[col_b])
    shared = genes_a & genes_b

    if len(shared) == 0:
        return {
            "n_shared": 0,
            "n_concordant": 0,
            "n_discordant": 0,
            "concordance_rate": 0.0,
            "pearson_r": np.nan
        }

    # Count directional concordance
    n_concordant = 0
    n_discordant = 0
    fc_a_vals = []
    fc_b_vals = []

    for gene in shared:
        row_a = df_a[df_a[col_a] == gene].iloc[0]
        row_b = df_b[df_b[col_b] == gene].iloc[0]

        fc_a = row_a["log2fc"]
        fc_b = row_b["log2fc"]

        fc_a_vals.append(fc_a)
        fc_b_vals.append(fc_b)

        if (fc_a > 0 and fc_b > 0) or (fc_a < 0 and fc_b < 0):
            n_concordant += 1
        else:
            n_discordant += 1

    # Calculate Pearson correlation
    r = np.corrcoef(fc_a_vals, fc_b_vals)[0, 1] if len(fc_a_vals) > 1 else np.nan

    return {
        "n_shared": len(shared),
        "n_concordant": n_concordant,
        "n_discordant": n_discordant,
        "concordance_rate": n_concordant / len(shared) if len(shared) > 0 else 0.0,
        "pearson_r": r
    }
