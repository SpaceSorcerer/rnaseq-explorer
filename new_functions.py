"""
=============================================================================
New / Updated Functions for DESeq2 & rMATS Pipeline
=============================================================================
Event-level matching via coordinate keys, opposite-direction Venn panels,
clustered heatmap, event pie chart, and pairwise Excel export.

These functions replace or supplement the existing gene-level versions in
deseq2_rmats_filter_pipeline.py. They use the same globals, imports, and
coding conventions.

Integration notes:
  - _make_event_key() is a new helper — place near _style_venn().
  - Updated functions replace their existing counterparts in-place.
  - New functions (rmats_event_heatmap, rmats_event_pie_chart,
    export_pairwise_workbook) are appended and wired into the multi-condition
    comparison section.
=============================================================================
"""

# ---- Imports identical to main pipeline ---------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from itertools import combinations

from matplotlib_venn import venn2, venn3

try:
    from upsetplot import UpSet, from_memberships
    _UPSET_AVAILABLE = True
except ImportError:
    _UPSET_AVAILABLE = False
    UpSet = None
    from_memberships = None

try:
    from scipy.stats import pearsonr
    _SCIPY_AVAILABLE = True
except ImportError:
    _SCIPY_AVAILABLE = False
    pearsonr = None


# ---- Globals (imported from main pipeline at integration time) ----------
# Duplicated here so the file is self-contained for review.
RMATS_EVENT_TYPES = ["SE", "A3SS", "A5SS", "RI", "MXE"]
RMATS_COLS = {
    "event_id":      "ID",
    "gene_id":       "GeneID",
    "gene_name":     "geneSymbol",
    "pvalue":        "PValue",
    "fdr":           "FDR",
    "inclevel_diff": "IncLevelDifference",
}
DESEQ2_COLS = {
    "gene_id":   "gene_id",
    "gene_name": "gene_name",
    "log2fc":    "log2FoldChange",
    "basemean":  "baseMean",
    "padj":      "padj",
    "pvalue":    "pvalue",
    "biotype":   "biotype",
    "stat":      "stat",
    "lfcSE":     "lfcSE",
}
INCLEVEL_DIFF_CUTOFF = 0.1
FIG_DPI = 300
FIG_FORMAT = "png"
FONT_SIZE = 12
COLOR_UP = "#E69F00"
COLOR_DOWN = "#0072B2"
COLOR_NS = "#BFBFBF"
EVENT_COLORS = {
    "SE":   "#E64B35",
    "A3SS": "#4DBBD5",
    "A5SS": "#00A087",
    "RI":   "#3C5488",
    "MXE":  "#F39B7F",
}

# Coordinate columns that uniquely identify each splicing event
_COORD_COLS = {
    "SE": ["chr", "strand", "exonStart_0base", "exonEnd",
            "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "A3SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",
             "shortES", "shortEE", "flankingES", "flankingEE"],
    "A5SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",
             "shortES", "shortEE", "flankingES", "flankingEE"],
    "RI": ["chr", "strand", "riExonStart_0base", "riExonEnd",
            "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "MXE": ["chr", "strand", "1stExonStart_0base", "1stExonEnd",
            "2ndExonStart_0base", "2ndExonEnd",
            "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
}


# ---------------------------------------------------------------------------
# Helper stubs referenced by updated functions (exist in main pipeline)
# ---------------------------------------------------------------------------
def _style_venn(v, n_sets):
    """Apply Okabe-Ito colors to a matplotlib_venn Venn diagram object."""
    if v is None:
        return
    if n_sets == 2:
        for rid, col in [('10', '#0072B2'), ('01', '#E69F00'), ('11', '#009E73')]:
            patch = v.get_patch_by_id(rid)
            if patch:
                patch.set_color(col)
                patch.set_alpha(0.6)
    elif n_sets == 3:
        for rid, col in [('100', '#0072B2'), ('010', '#E69F00'), ('001', '#56B4E9'),
                         ('110', '#009E73'), ('101', '#CC79A7'), ('011', '#F0E442'),
                         ('111', '#D55E00')]:
            patch = v.get_patch_by_id(rid)
            if patch:
                patch.set_color(col)
                patch.set_alpha(0.6)


def add_count_box(ax, n_up, n_down, total, position="lower left",
                  up_label="Up", down_label="Down"):
    """Add a compact count box in a data-sparse corner of the plot."""
    text = (f"{up_label}: {n_up:,}\n"
            f"{down_label}: {n_down:,}\n"
            f"Total: {total:,}")
    loc = {"upper left": (0.02, 0.98), "upper right": (0.98, 0.98),
           "lower left": (0.02, 0.02), "lower right": (0.98, 0.02)}
    x, y = loc.get(position, (0.02, 0.02))
    ha = "left" if "left" in position else "right"
    va = "top" if "upper" in position else "bottom"
    ax.text(x, y, text, transform=ax.transAxes, fontsize=9,
            va=va, ha=ha, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="grey", alpha=0.9))


def _grid_dims(n):
    """Return (nrows, ncols) for a tight subplot grid of n panels."""
    if n <= 1:
        return (1, 1)
    if n == 2:
        return (1, 2)
    if n == 3:
        return (1, 3)
    if n == 4:
        return (2, 2)
    if n == 5:
        return (2, 3)
    ncols = int(np.ceil(np.sqrt(n)))
    nrows = int(np.ceil(n / ncols))
    return (nrows, ncols)


def _best_gene_key(df):
    """Return the column name that gives the most reliable unique gene identifier."""
    id_col = DESEQ2_COLS.get("gene_id", "")
    name_col = DESEQ2_COLS.get("gene_name", "")
    if id_col and id_col in df.columns:
        sample = df[id_col].dropna().astype(str).head(500)
        ens_frac = sample.str.upper().str.startswith("ENS").sum() / max(len(sample), 1)
        if ens_frac > 0.1:
            return id_col, "Ensembl ID"
    return name_col, "gene name"


# ===================================================================
# 1. _make_event_key() — unique event identifier from coordinates
# ===================================================================

def _make_event_key(df, event_type):
    """Create a unique event key from genomic coordinate columns.

    Returns a pandas Series of strings in the format
    ``chr:strand:col1:col2:...`` for each row in *df*.  If any required
    coordinate column is missing, an empty string is returned for every
    row so callers can filter gracefully.

    Parameters
    ----------
    df : pd.DataFrame
        rMATS filtered or raw DataFrame for a single event type.
    event_type : str
        One of ``RMATS_EVENT_TYPES`` (SE, A3SS, A5SS, RI, MXE).
    """
    coord_cols = _COORD_COLS.get(event_type, [])
    if not coord_cols:
        return pd.Series([""] * len(df), index=df.index)

    # Check that all required columns are present
    missing = [c for c in coord_cols if c not in df.columns]
    if missing:
        return pd.Series([""] * len(df), index=df.index)

    # Build key by concatenating coordinate values with ':' separator
    key = df[coord_cols[0]].astype(str)
    for col in coord_cols[1:]:
        key = key + ":" + df[col].astype(str)
    return key


# ===================================================================
# 2. pairwise_splicing_venns() — EVENT-LEVEL with opposite-direction
# ===================================================================

def pairwise_splicing_venns(condition_results, condition_labels, outdir):
    """Pairwise 5-panel Venn diagrams for splicing events (event-level).

    For each pair of conditions x each event type:
        1. All Significant events
        2. Included in both  (dPSI >= INCLEVEL_DIFF_CUTOFF)
        3. Excluded in both  (dPSI <= -INCLEVEL_DIFF_CUTOFF)
        4. Included in {A}, Excluded in {B}
        5. Excluded in {A}, Included in {B}

    Events are matched by genomic coordinates via _make_event_key().
    """
    outdir = Path(outdir)
    names = list(condition_results.keys())
    dpsi_col = RMATS_COLS["inclevel_diff"]
    gene_col = RMATS_COLS["gene_name"]

    for name_a, name_b in combinations(names, 2):
        label_a = condition_labels[name_a]
        label_b = condition_labels[name_b]

        for et in RMATS_EVENT_TYPES:
            df_a = condition_results[name_a]["rmats_filtered"].get(et)
            df_b = condition_results[name_b]["rmats_filtered"].get(et)

            if df_a is None or len(df_a) == 0 or df_b is None or len(df_b) == 0:
                continue

            # Build event keys
            key_a = _make_event_key(df_a, et)
            key_b = _make_event_key(df_b, et)

            if key_a.eq("").all() or key_b.eq("").all():
                continue

            # Add keys as temporary column for set operations
            df_a = df_a.copy()
            df_b = df_b.copy()
            df_a["_ekey"] = key_a.values
            df_b["_ekey"] = key_b.values

            events_all_a = set(df_a["_ekey"].dropna().unique())
            events_all_b = set(df_b["_ekey"].dropna().unique())

            # Direction subsets
            events_inc_a = set(
                df_a.loc[df_a[dpsi_col] >= INCLEVEL_DIFF_CUTOFF, "_ekey"]
                .dropna().unique())
            events_inc_b = set(
                df_b.loc[df_b[dpsi_col] >= INCLEVEL_DIFF_CUTOFF, "_ekey"]
                .dropna().unique())
            events_exc_a = set(
                df_a.loc[df_a[dpsi_col] <= -INCLEVEL_DIFF_CUTOFF, "_ekey"]
                .dropna().unique())
            events_exc_b = set(
                df_b.loc[df_b[dpsi_col] <= -INCLEVEL_DIFF_CUTOFF, "_ekey"]
                .dropna().unique())

            # Opposite-direction intersections
            inc_a_exc_b = events_inc_a & events_exc_b
            exc_a_inc_b = events_exc_a & events_inc_b

            panels = [
                ("All Significant", events_all_a, events_all_b),
                (f"Included in Both (dPSI \u2265 {INCLEVEL_DIFF_CUTOFF})",
                 events_inc_a, events_inc_b),
                (f"Excluded in Both (dPSI \u2264 \u2212{INCLEVEL_DIFF_CUTOFF})",
                 events_exc_a, events_exc_b),
                (f"Inc {label_a} / Exc {label_b}",
                 events_inc_a & events_all_b, events_exc_b & events_all_a),
                (f"Exc {label_a} / Inc {label_b}",
                 events_exc_a & events_all_b, events_inc_b & events_all_a),
            ]

            fig, axes = plt.subplots(2, 3, figsize=(20, 12))
            axes_flat = axes.flatten()

            for idx, (panel_title, set_a, set_b) in enumerate(panels):
                ax = axes_flat[idx]
                v = venn2([set_a, set_b], set_labels=(label_a, label_b), ax=ax)
                _style_venn(v, 2)
                ax.set_title(f"{panel_title}\n(n={len(set_a | set_b)})",
                             fontsize=11, fontweight="bold")

            # Hide unused 6th panel
            axes_flat[5].set_visible(False)

            fig.suptitle(f"{et} Splicing Events (event-level) \u2014 "
                         f"{label_a} vs {label_b}",
                         fontsize=14, fontweight="bold")
            plt.tight_layout()
            outpath = outdir / f"venn_splicing_{et}_{name_a}_vs_{name_b}.{FIG_FORMAT}"
            fig.savefig(outpath, format=FIG_FORMAT, dpi=FIG_DPI, bbox_inches="tight")
            plt.close(fig)

            print(f"  Saved: {outpath.name} "
                  f"(All: {len(events_all_a - events_all_b)}|"
                  f"{len(events_all_a & events_all_b)}|"
                  f"{len(events_all_b - events_all_a)}, "
                  f"Inc/Exc: {len(inc_a_exc_b)}, "
                  f"Exc/Inc: {len(exc_a_inc_b)})")


# ===================================================================
# 3. pairwise_deg_venns() — 5-panel with opposite-direction panels
# ===================================================================

def pairwise_deg_venns(condition_results, condition_labels, outdir):
    """Pairwise 5-panel Venn diagrams for differentially expressed genes.

    Panels:
        1. All Significant DEGs
        2. Upregulated in both
        3. Downregulated in both
        4. Up in {A}, Down in {B}
        5. Down in {A}, Up in {B}
    """
    outdir = Path(outdir)
    names = list(condition_results.keys())

    for name_a, name_b in combinations(names, 2):
        label_a = condition_labels[name_a]
        label_b = condition_labels[name_b]

        filt_a = condition_results[name_a]["deseq2_filtered"]["all_genes"]
        filt_b = condition_results[name_b]["deseq2_filtered"]["all_genes"]

        key_col_a, _ = _best_gene_key(filt_a)
        key_col_b, _ = _best_gene_key(filt_b)
        if not key_col_a or key_col_a not in filt_a.columns:
            key_col_a = DESEQ2_COLS["gene_name"]
        if not key_col_b or key_col_b not in filt_b.columns:
            key_col_b = DESEQ2_COLS["gene_name"]

        all_a = set(filt_a[key_col_a].dropna().unique())
        all_b = set(filt_b[key_col_b].dropna().unique())
        up_a = set(filt_a.loc[filt_a["direction"] == "up", key_col_a].dropna().unique())
        up_b = set(filt_b.loc[filt_b["direction"] == "up", key_col_b].dropna().unique())
        down_a = set(filt_a.loc[filt_a["direction"] == "down", key_col_a].dropna().unique())
        down_b = set(filt_b.loc[filt_b["direction"] == "down", key_col_b].dropna().unique())

        # Opposite-direction intersections
        up_a_down_b = up_a & down_b
        down_a_up_b = down_a & up_b

        panels = [
            ("All Significant DEGs", all_a, all_b),
            ("Upregulated", up_a, up_b),
            ("Downregulated", down_a, down_b),
            (f"Up in {label_a} / Down in {label_b}", up_a & all_b, down_b & all_a),
            (f"Down in {label_a} / Up in {label_b}", down_a & all_b, up_b & all_a),
        ]

        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        axes_flat = axes.flatten()

        for idx, (panel_title, set_a, set_b) in enumerate(panels):
            ax = axes_flat[idx]
            v = venn2([set_a, set_b], set_labels=(label_a, label_b), ax=ax)
            _style_venn(v, 2)
            ax.set_title(f"{panel_title}\n(n={len(set_a | set_b)})",
                         fontsize=11, fontweight="bold")

        # Hide unused 6th panel
        axes_flat[5].set_visible(False)

        fig.suptitle(f"Differentially Expressed Genes \u2014 "
                     f"{label_a} vs {label_b}",
                     fontsize=14, fontweight="bold")
        plt.tight_layout()
        outpath = outdir / f"venn_deg_{name_a}_vs_{name_b}.{FIG_FORMAT}"
        fig.savefig(outpath, format=FIG_FORMAT, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)

        print(f"  Saved: {outpath.name} "
              f"(All: {len(all_a - all_b)}|{len(all_a & all_b)}|{len(all_b - all_a)}, "
              f"Up: {len(up_a - up_b)}|{len(up_a & up_b)}|{len(up_b - up_a)}, "
              f"Down: {len(down_a - down_b)}|{len(down_a & down_b)}|{len(down_b - down_a)}, "
              f"Up/Down: {len(up_a_down_b)}, Down/Up: {len(down_a_up_b)})")


# ===================================================================
# 4. rmats_cross_condition_venn() — EVENT-LEVEL
# ===================================================================

def rmats_cross_condition_venn(condition_results, condition_labels, outdir):
    """Venn diagrams for significant splicing events across conditions (event-level)."""
    names = list(condition_results.keys())
    labels = [condition_labels[n] for n in names]
    gene_col = RMATS_COLS["gene_name"]

    # --- Event-level Venn (all event types combined) ---
    event_sets = []
    for name in names:
        all_sig_events = set()
        for et, filt_df in condition_results[name]["rmats_filtered"].items():
            keys = _make_event_key(filt_df, et)
            all_sig_events.update(keys[keys != ""].unique())
        event_sets.append(all_sig_events)

    fig, ax = plt.subplots(figsize=(8, 8))
    if len(event_sets) == 2:
        v = venn2(event_sets, set_labels=labels, ax=ax)
        _style_venn(v, 2)
    elif len(event_sets) == 3:
        v = venn3(event_sets, set_labels=labels, ax=ax)
        _style_venn(v, 3)
    ax.set_title("rMATS \u2014 Significant Splicing Events (coordinate-level)",
                 fontsize=13, fontweight="bold")
    outpath = Path(outdir) / f"venn_rmats_events.{FIG_FORMAT}"
    fig.savefig(outpath, format=FIG_FORMAT, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")

    # --- Per-event-type Venns ---
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
                keys = _make_event_key(filt_df, et)
                events = set(keys[keys != ""].unique())
            else:
                events = set()
            sets.append(events)

        if len(sets) == 2:
            v = venn2(sets, set_labels=labels, ax=ax)
            _style_venn(v, 2)
        elif len(sets) == 3:
            v = venn3(sets, set_labels=labels, ax=ax)
            _style_venn(v, 3)
        ax.set_title(f"{et}", fontsize=12, fontweight="bold")

    fig.suptitle("rMATS \u2014 Splicing Event Overlap by Type (coordinate-level)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    outpath = Path(outdir) / f"venn_rmats_events_by_type.{FIG_FORMAT}"
    fig.savefig(outpath, format=FIG_FORMAT, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ===================================================================
# 5. rmats_upset_plot() — EVENT-LEVEL
# ===================================================================

def rmats_upset_plot(rmats_conditions, condition_labels, outdir):
    """UpSet plot for rMATS splicing event sets across 3+ conditions (event-level)."""
    if not _UPSET_AVAILABLE:
        print("  Skipping rMATS UpSet: 'upsetplot' not installed (pip install upsetplot)")
        return
    if len(rmats_conditions) < 3:
        print("  Skipping rMATS UpSet: requires 3+ conditions with rMATS data")
        return

    names = list(rmats_conditions.keys())
    labels = [condition_labels[n] for n in names]

    # Event-level upset (all event types combined)
    event_sets = {}
    for name, lbl in zip(names, labels):
        events = set()
        for et, et_df in rmats_conditions[name]["rmats_filtered"].items():
            if len(et_df) > 0:
                keys = _make_event_key(et_df, et)
                # Prefix with event type to avoid cross-type collisions
                events.update(f"{et}|{k}" for k in keys[keys != ""].unique())
        event_sets[lbl] = events

    all_events = set.union(*event_sets.values()) if event_sets else set()
    if all_events:
        memberships = [tuple(lbl for lbl in labels if ev in event_sets[lbl])
                       for ev in all_events]
        try:
            upset_data = from_memberships(memberships)
            upset_data = upset_data.groupby(
                level=list(range(upset_data.index.nlevels))).sum()
            upset = UpSet(upset_data, show_counts=True, sort_by="cardinality")
            upset.plot()
            plt.suptitle("rMATS UpSet \u2014 Splicing Event Overlap (coordinate-level)",
                         y=1.02, fontsize=12, fontweight="bold")
            outpath = Path(outdir) / f"rmats_upset_events.{FIG_FORMAT}"
            plt.savefig(outpath, format=FIG_FORMAT, bbox_inches="tight")
            plt.close("all")
            print(f"  Saved: {outpath}")
        except Exception as e:
            print(f"  rMATS UpSet (events) failed: {e}")

    # Per-event-type upsets
    for et in RMATS_EVENT_TYPES:
        et_sets = {}
        for name, lbl in zip(names, labels):
            filt_df = rmats_conditions[name]["rmats_filtered"].get(et, pd.DataFrame())
            if len(filt_df) > 0:
                keys = _make_event_key(filt_df, et)
                et_sets[lbl] = set(keys[keys != ""].unique())
        active = {lbl: s for lbl, s in et_sets.items() if s}
        if len(active) < 2:
            continue
        all_et_events = set.union(*active.values())
        all_labels = labels
        memberships = [tuple(lbl for lbl in all_labels
                             if lbl in et_sets and ev in et_sets[lbl])
                       for ev in all_et_events]
        try:
            upset_data = from_memberships(memberships)
            upset_data = upset_data.groupby(
                level=list(range(upset_data.index.nlevels))).sum()
            upset = UpSet(upset_data, show_counts=True, sort_by="cardinality")
            upset.plot()
            plt.suptitle(f"rMATS UpSet \u2014 {et} Event Overlap (coordinate-level)",
                         y=1.02, fontsize=12, fontweight="bold")
            outpath = Path(outdir) / f"rmats_upset_{et}.{FIG_FORMAT}"
            plt.savefig(outpath, format=FIG_FORMAT, bbox_inches="tight")
            plt.close("all")
            print(f"  Saved: {outpath}")
        except Exception as e:
            print(f"  rMATS UpSet ({et}) failed: {e}")


# ===================================================================
# 6. rmats_directional_venn_diagrams() — EVENT-LEVEL
# ===================================================================

def rmats_directional_venn_diagrams(condition_results, condition_labels, outdir):
    """4-panel directional Venn diagrams for splicing events by event type (event-level).

    Panels: A. All Events, B. Concordant Included, C. Concordant Excluded, D. Discordant
    Events are matched by genomic coordinates via _make_event_key().
    """
    names = list(condition_results.keys())
    labels = [condition_labels[n] for n in names]
    dpsi_col = RMATS_COLS["inclevel_diff"]

    if len(names) < 2:
        print("  Directional Venn diagrams require at least 2 conditions")
        return

    for et in RMATS_EVENT_TYPES:
        # Collect dPSI per event key for this event type
        dfs = {}
        for name in names:
            if et in condition_results[name]["rmats_filtered"]:
                df = condition_results[name]["rmats_filtered"][et]
                if len(df) > 0:
                    df = df.copy()
                    df["_ekey"] = _make_event_key(df, et).values
                    df = df[df["_ekey"] != ""]
                    if len(df) > 0:
                        # Use mean dPSI when multiple rows share the same event key
                        event_dpsi = df.groupby("_ekey")[dpsi_col].mean()
                        dfs[name] = event_dpsi

        if len(dfs) < 2:
            print(f"  Skipping {et} directional Venn (insufficient conditions)")
            continue

        # Identify all events per condition
        all_events_per_cond = {n: set(d.index) for n, d in dfs.items()}

        # Find shared events and classify by direction
        shared_events = set.intersection(*all_events_per_cond.values())

        concordant_up = set()
        concordant_down = set()
        discordant = set()

        for ekey in shared_events:
            signs = [np.sign(dfs[name][ekey]) for name in names]
            if all(s > 0 for s in signs):
                concordant_up.add(ekey)
            elif all(s < 0 for s in signs):
                concordant_down.add(ekey)
            else:
                discordant.add(ekey)

        # Validation
        computed_all = len(concordant_up) + len(concordant_down) + len(discordant)
        if computed_all != len(shared_events):
            print(f"  WARNING: Venn math mismatch for {et}: "
                  f"shared={len(shared_events)} but sum={computed_all}")

        # Create 4-panel figure
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        # Panel A: All Events
        ax = axes[0]
        if len(names) == 2:
            v = venn2(list(all_events_per_cond.values()), set_labels=labels, ax=ax)
            _style_venn(v, 2)
        elif len(names) == 3:
            v = venn3(list(all_events_per_cond.values()), set_labels=labels, ax=ax)
            _style_venn(v, 3)
        ax.set_title(
            f"A. All Events (n={sum(len(s) for s in all_events_per_cond.values())})",
            fontsize=12, fontweight="bold")

        # Panel B: Concordant Included
        ax = axes[1]
        concordant_up_per_cond = []
        for name in names:
            events_up = set()
            if name in dfs:
                for ekey in concordant_up:
                    if ekey in dfs[name].index:
                        events_up.add(ekey)
            concordant_up_per_cond.append(events_up)

        if len(names) == 2:
            v = venn2(concordant_up_per_cond, set_labels=labels, ax=ax)
            _style_venn(v, 2)
        elif len(names) == 3:
            v = venn3(concordant_up_per_cond, set_labels=labels, ax=ax)
            _style_venn(v, 3)
        ax.set_title(f"B. Concordant Included (n={len(concordant_up)})",
                     fontsize=12, fontweight="bold")

        # Panel C: Concordant Excluded
        ax = axes[2]
        concordant_down_per_cond = []
        for name in names:
            events_down = set()
            if name in dfs:
                for ekey in concordant_down:
                    if ekey in dfs[name].index:
                        events_down.add(ekey)
            concordant_down_per_cond.append(events_down)

        if len(names) == 2:
            v = venn2(concordant_down_per_cond, set_labels=labels, ax=ax)
            _style_venn(v, 2)
        elif len(names) == 3:
            v = venn3(concordant_down_per_cond, set_labels=labels, ax=ax)
            _style_venn(v, 3)
        ax.set_title(f"C. Concordant Excluded (n={len(concordant_down)})",
                     fontsize=12, fontweight="bold")

        # Panel D: Discordant
        ax = axes[3]
        discordant_per_cond = []
        for name in names:
            events_disc = set()
            if name in dfs:
                for ekey in discordant:
                    if ekey in dfs[name].index:
                        events_disc.add(ekey)
            discordant_per_cond.append(events_disc)

        if len(names) == 2:
            v = venn2(discordant_per_cond, set_labels=labels, ax=ax)
            _style_venn(v, 2)
        elif len(names) == 3:
            v = venn3(discordant_per_cond, set_labels=labels, ax=ax)
            _style_venn(v, 3)
        ax.set_title(f"D. Discordant (n={len(discordant)})",
                     fontsize=12, fontweight="bold")

        fig.suptitle(f"Directional Splicing Overlap \u2014 {et} Events (coordinate-level)",
                     fontsize=14, fontweight="bold", y=0.98)
        plt.tight_layout()

        outpath = Path(outdir) / f"venn_rmats_directional_{et}.{FIG_FORMAT}"
        fig.savefig(outpath, format=FIG_FORMAT, dpi=FIG_DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {outpath} (Inc={len(concordant_up)}, Exc={len(concordant_down)}, "
              f"Disc={len(discordant)}, Total={len(shared_events)})")

    # --- Pairwise directional Venns (venn2) ---
    if len(names) >= 3:
        for name_a, name_b in combinations(names, 2):
            label_a = condition_labels[name_a]
            label_b = condition_labels[name_b]
            pair_labels = [label_a, label_b]

            for et in RMATS_EVENT_TYPES:
                pair_dfs = {}
                for nm in (name_a, name_b):
                    if et in condition_results[nm]["rmats_filtered"]:
                        df = condition_results[nm]["rmats_filtered"][et]
                        if len(df) > 0:
                            df = df.copy()
                            df["_ekey"] = _make_event_key(df, et).values
                            df = df[df["_ekey"] != ""]
                            if len(df) > 0:
                                pair_dfs[nm] = df.groupby("_ekey")[dpsi_col].mean()

                if len(pair_dfs) < 2:
                    continue

                events_a = set(pair_dfs[name_a].index)
                events_b = set(pair_dfs[name_b].index)
                shared = events_a & events_b

                conc_up = set()
                conc_down = set()
                disc = set()
                for ekey in shared:
                    sa = np.sign(pair_dfs[name_a][ekey])
                    sb = np.sign(pair_dfs[name_b][ekey])
                    if sa > 0 and sb > 0:
                        conc_up.add(ekey)
                    elif sa < 0 and sb < 0:
                        conc_down.add(ekey)
                    else:
                        disc.add(ekey)

                # Validation
                computed = len(conc_up) + len(conc_down) + len(disc)
                if computed != len(shared):
                    print(f"  WARNING: Pairwise venn math mismatch for {et} "
                          f"{name_a} vs {name_b}: shared={len(shared)} sum={computed}")

                fig, axes = plt.subplots(2, 2, figsize=(14, 12))
                axes_flat = axes.flatten()

                # Panel A: All Events
                v = venn2([events_a, events_b], set_labels=pair_labels, ax=axes_flat[0])
                _style_venn(v, 2)
                axes_flat[0].set_title(
                    f"A. All Events (n={len(events_a | events_b)})",
                    fontsize=12, fontweight="bold")

                # Panel B: Concordant Included
                cup_a = {e for e in conc_up if e in pair_dfs[name_a].index}
                cup_b = {e for e in conc_up if e in pair_dfs[name_b].index}
                v = venn2([cup_a, cup_b], set_labels=pair_labels, ax=axes_flat[1])
                _style_venn(v, 2)
                axes_flat[1].set_title(
                    f"B. Concordant Included (n={len(conc_up)})",
                    fontsize=12, fontweight="bold")

                # Panel C: Concordant Excluded
                cdn_a = {e for e in conc_down if e in pair_dfs[name_a].index}
                cdn_b = {e for e in conc_down if e in pair_dfs[name_b].index}
                v = venn2([cdn_a, cdn_b], set_labels=pair_labels, ax=axes_flat[2])
                _style_venn(v, 2)
                axes_flat[2].set_title(
                    f"C. Concordant Excluded (n={len(conc_down)})",
                    fontsize=12, fontweight="bold")

                # Panel D: Discordant
                dsc_a = {e for e in disc if e in pair_dfs[name_a].index}
                dsc_b = {e for e in disc if e in pair_dfs[name_b].index}
                v = venn2([dsc_a, dsc_b], set_labels=pair_labels, ax=axes_flat[3])
                _style_venn(v, 2)
                axes_flat[3].set_title(
                    f"D. Discordant (n={len(disc)})",
                    fontsize=12, fontweight="bold")

                fig.suptitle(
                    f"Directional Splicing \u2014 {et} \u2014 "
                    f"{label_a} vs {label_b} (coordinate-level)",
                    fontsize=13, fontweight="bold", y=0.98)
                plt.tight_layout()
                outpath = (Path(outdir) /
                           f"venn_rmats_directional_{et}_{name_a}_vs_{name_b}.{FIG_FORMAT}")
                fig.savefig(outpath, format=FIG_FORMAT, dpi=FIG_DPI, bbox_inches="tight")
                plt.close(fig)
                print(f"  Saved: {outpath.name} (Inc={len(conc_up)}, "
                      f"Exc={len(conc_down)}, Disc={len(disc)})")


# ===================================================================
# 7. pairwise_dpsi_scatter() — EVENT-LEVEL
# ===================================================================

def pairwise_dpsi_scatter(rmats_conditions, condition_labels, outdir):
    """Pairwise scatter of dPSI for shared splicing events (event-level matching)."""
    names = list(rmats_conditions.keys())
    if len(names) < 2:
        return

    dpsi_col = RMATS_COLS["inclevel_diff"]
    gene_col = RMATS_COLS["gene_name"]
    pairs = list(combinations(names, 2))
    nrows, ncols = _grid_dims(len(pairs))

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), squeeze=False)
    axes_flat = [axes[r][c] for r in range(nrows) for c in range(ncols)]

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

            # Build event keys
            keyA = _make_event_key(dfA, et)
            keyB = _make_event_key(dfB, et)
            if keyA.eq("").all() or keyB.eq("").all():
                continue

            dfA = dfA.copy()
            dfB = dfB.copy()
            dfA["_ekey"] = keyA.values
            dfB["_ekey"] = keyB.values

            # Keep gene symbol for tooltip / labeling
            cols_a = ["_ekey", dpsi_col]
            cols_b = ["_ekey", dpsi_col]
            if gene_col in dfA.columns:
                cols_a.append(gene_col)
            if gene_col in dfB.columns:
                cols_b.append(gene_col)

            merged = dfA[cols_a].merge(
                dfB[cols_b], on="_ekey", suffixes=("_A", "_B"))
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
            ax.scatter(sub[f"{dpsi_col}_A"], sub[f"{dpsi_col}_B"],
                       c=EVENT_COLORS.get(et, "#888888"), s=8, alpha=0.6,
                       edgecolors="none", rasterized=True, label=f"{et} ({len(sub):,})")

        lims = [min(combined[f"{dpsi_col}_A"].min(),
                    combined[f"{dpsi_col}_B"].min()) - 0.05,
                max(combined[f"{dpsi_col}_A"].max(),
                    combined[f"{dpsi_col}_B"].max()) + 0.05]
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
            ax.annotate(f"R\u00b2 = {r**2:.3f}", xy=(0, 1), xycoords="axes fraction",
                        xytext=(4, 4), textcoords="offset points",
                        ha="left", va="bottom", fontsize=9,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                  edgecolor="grey", alpha=0.9))

        ax.set_xlabel(f"dPSI  {lblA}", fontsize=9)
        ax.set_ylabel(f"dPSI  {lblB}", fontsize=9)
        ax.set_title(f"{lblA} vs {lblB}", fontsize=10)
        ax.legend(fontsize=7, markerscale=2, loc="lower right")

    for i in range(len(pairs), nrows * ncols):
        axes_flat[i].set_visible(False)

    fig.suptitle("Pairwise dPSI Comparison (Shared Splicing Events, coordinate-level)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    outpath = Path(outdir) / f"pairwise_dpsi_scatter.{FIG_FORMAT}"
    fig.savefig(outpath, format=FIG_FORMAT, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ===================================================================
# 8. rmats_event_heatmap() — clustered heatmap of dPSI across conditions
# ===================================================================

def rmats_event_heatmap(condition_results, condition_labels, event_type, outdir):
    """Clustered heatmap of dPSI values across conditions for a specific event type.

    Parameters
    ----------
    condition_results : dict
        {name: {"rmats_filtered": {et: df, ...}, ...}} structure.
    condition_labels : dict
        {name: label} mapping.
    event_type : str
        One of RMATS_EVENT_TYPES (e.g. "SE").
    outdir : str or Path
        Output directory.
    """
    outdir = Path(outdir)
    names = list(condition_results.keys())
    labels = [condition_labels[n] for n in names]
    dpsi_col = RMATS_COLS["inclevel_diff"]
    gene_col = RMATS_COLS["gene_name"]

    # Collect dPSI per event key for each condition
    event_data = {}
    gene_lookup = {}  # event_key -> geneSymbol

    for name, lbl in zip(names, labels):
        filt_df = condition_results[name]["rmats_filtered"].get(event_type, pd.DataFrame())
        if len(filt_df) == 0:
            continue
        df = filt_df.copy()
        df["_ekey"] = _make_event_key(df, event_type).values
        df = df[df["_ekey"] != ""]
        if len(df) == 0:
            continue

        # Mean dPSI per event key (handles potential duplicates)
        grouped = df.groupby("_ekey").agg(
            dpsi=(dpsi_col, "mean"),
            gene=(gene_col, "first") if gene_col in df.columns else (dpsi_col, "count"),
        )
        event_data[lbl] = grouped["dpsi"]

        # Store gene symbol mapping
        if gene_col in df.columns:
            for ekey, gene in df.groupby("_ekey")[gene_col].first().items():
                if ekey not in gene_lookup and pd.notna(gene):
                    gene_lookup[ekey] = gene

    if len(event_data) < 2:
        print(f"  Skipping {event_type} heatmap (insufficient conditions with data)")
        return

    # Build matrix: rows = events, columns = conditions
    dpsi_df = pd.DataFrame(event_data)

    # Keep only events significant in at least 1 condition (non-NaN in >= 1 col)
    dpsi_df = dpsi_df.dropna(how="all")
    if len(dpsi_df) == 0:
        print(f"  Skipping {event_type} heatmap (no events to plot)")
        return

    # Limit to top 80 events by max |dPSI| across conditions
    max_events = 80
    if len(dpsi_df) > max_events:
        max_abs_dpsi = dpsi_df.abs().max(axis=1)
        top_idx = max_abs_dpsi.nlargest(max_events).index
        dpsi_df = dpsi_df.loc[top_idx]

    # Replace event keys with gene symbols for row labels
    row_labels = []
    for ekey in dpsi_df.index:
        gene = gene_lookup.get(ekey, "")
        if gene:
            # Abbreviate coordinate key for uniqueness
            coord_short = ekey.split(":")[0] + ":" + ekey.split(":")[-1]
            row_labels.append(f"{gene} ({coord_short})")
        else:
            row_labels.append(ekey[:40])
    dpsi_df.index = row_labels

    # Fill NaN with 0 for clustering (event not significant in that condition)
    dpsi_filled = dpsi_df.fillna(0)

    # Create diverging colormap matching pipeline colors (blue-white-orange)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(
        "dpsi_diverging", [COLOR_DOWN, "white", COLOR_UP], N=256)

    # Determine symmetric color limits
    vmax = max(abs(dpsi_filled.values.min()), abs(dpsi_filled.values.max()))
    vmax = max(vmax, INCLEVEL_DIFF_CUTOFF)  # Ensure at least the cutoff is visible

    # Determine figure height based on number of events
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
            dendrogram_ratio=(0.15, 0.05),
            cbar_kws={"label": "dPSI (IncLevelDifference)", "shrink": 0.6},
            linewidths=0.5,
            linecolor="white",
            yticklabels=True,
            xticklabels=True,
            row_cluster=len(dpsi_filled) > 1,
            col_cluster=len(dpsi_filled.columns) > 1,
        )

        g.fig.suptitle(
            f"dPSI Heatmap \u2014 {event_type} Events (coordinate-level, "
            f"top {len(dpsi_df)} by |dPSI|)",
            fontsize=13, fontweight="bold", y=1.02)

        # Adjust row label font size for readability
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_yticklabels(), fontsize=7, rotation=0)
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xticklabels(), fontsize=10, rotation=45, ha="right")

        outpath = outdir / f"heatmap_dpsi_{event_type}.{FIG_FORMAT}"
        g.savefig(outpath, format=FIG_FORMAT, dpi=FIG_DPI, bbox_inches="tight")
        plt.close("all")
        print(f"  Saved: {outpath}")
    except Exception as e:
        print(f"  {event_type} heatmap failed: {e}")
        plt.close("all")


# ===================================================================
# 9. rmats_event_pie_chart() — event type distribution per condition
# ===================================================================

def rmats_event_pie_chart(condition_results, condition_labels, outdir):
    """Pie chart per condition showing splicing event type distribution.

    Parameters
    ----------
    condition_results : dict
        {name: {"rmats_filtered": {et: df, ...}, ...}}.
    condition_labels : dict
        {name: label}.
    outdir : str or Path
        Output directory.
    """
    outdir = Path(outdir)
    names = list(condition_results.keys())
    labels = [condition_labels[n] for n in names]
    n_conds = len(names)

    if n_conds == 0:
        return

    nrows, ncols = _grid_dims(n_conds)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
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
            ax.text(0.5, 0.5, "No significant\nevents", ha="center", va="center",
                    fontsize=11, transform=ax.transAxes)
            ax.set_title(lbl, fontsize=12, fontweight="bold")
            ax.axis("off")
            continue

        total = sum(counts)
        wedges, texts, autotexts = ax.pie(
            counts,
            labels=et_labels,
            colors=colors,
            autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct / 100 * total))})",
            startangle=90,
            pctdistance=0.65,
            textprops={"fontsize": 9},
        )
        for autotext in autotexts:
            autotext.set_fontsize(8)
            autotext.set_fontweight("bold")

        ax.set_title(f"{lbl}\n(n={total:,} events)", fontsize=12, fontweight="bold")

    # Hide unused panels
    for i in range(n_conds, len(axes_flat)):
        axes_flat[i].set_visible(False)

    fig.suptitle("rMATS \u2014 Significant Event Type Distribution",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    outpath = outdir / f"rmats_event_type_pie.{FIG_FORMAT}"
    fig.savefig(outpath, format=FIG_FORMAT, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


# ===================================================================
# 10. export_pairwise_workbook() — Excel export of pairwise comparisons
# ===================================================================

def export_pairwise_workbook(condition_results, condition_labels, outdir):
    """Export pairwise comparison Excel workbooks for DESeq2 and rMATS data.

    For each pair of conditions, creates one .xlsx file with:
      - DESeq2 sheets: shared (all/up/down), opposite-direction, condition-only
      - rMATS sheets (per event type): same breakdown by event coordinates
      - Summary sheet with row counts per category

    Parameters
    ----------
    condition_results : dict
        {name: {"deseq2_filtered": {"all_genes": df}, "rmats_filtered": {et: df}}}.
    condition_labels : dict
        {name: label}.
    outdir : str or Path
        Output directory.
    """
    outdir = Path(outdir)
    names = list(condition_results.keys())
    dpsi_col = RMATS_COLS["inclevel_diff"]
    gene_col = RMATS_COLS["gene_name"]
    lfc_col = DESEQ2_COLS["log2fc"]

    def _short(label):
        """Short label from condition label for column suffixing."""
        parts = label.split(" vs ")
        return parts[0].replace(" ", "_") if parts else label.replace(" ", "_")

    for name_a, name_b in combinations(names, 2):
        label_a = condition_labels[name_a]
        label_b = condition_labels[name_b]
        short_a = _short(label_a)
        short_b = _short(label_b)

        wb_name = f"pairwise_{name_a}_vs_{name_b}.xlsx"
        wb_path = outdir / wb_name
        summary_rows = []

        with pd.ExcelWriter(wb_path, engine="openpyxl") as writer:

            # -------------------------------------------------------
            # DESeq2 sheets
            # -------------------------------------------------------
            filt_a = condition_results[name_a]["deseq2_filtered"]["all_genes"].copy()
            filt_b = condition_results[name_b]["deseq2_filtered"]["all_genes"].copy()

            # Determine best gene key for matching
            key_col_a, _ = _best_gene_key(filt_a)
            key_col_b, _ = _best_gene_key(filt_b)
            if not key_col_a or key_col_a not in filt_a.columns:
                key_col_a = DESEQ2_COLS["gene_name"]
            if not key_col_b or key_col_b not in filt_b.columns:
                key_col_b = DESEQ2_COLS["gene_name"]

            # Ensure consistent key column name for merging
            merge_key = key_col_a
            if key_col_b != key_col_a and key_col_b in filt_b.columns:
                # Rename to match for merge
                filt_b = filt_b.rename(columns={key_col_b: merge_key})

            # Suffix non-key columns
            suffix_a = f"_{short_a}"
            suffix_b = f"_{short_b}"

            # Build gene sets
            all_a = set(filt_a[merge_key].dropna().unique())
            all_b = set(filt_b[merge_key].dropna().unique())

            dir_col = "direction"
            up_a = set(filt_a.loc[filt_a[dir_col] == "up", merge_key].dropna().unique())
            up_b = set(filt_b.loc[filt_b[dir_col] == "up", merge_key].dropna().unique())
            down_a = set(filt_a.loc[filt_a[dir_col] == "down", merge_key].dropna().unique())
            down_b = set(filt_b.loc[filt_b[dir_col] == "down", merge_key].dropna().unique())

            shared_all = all_a & all_b
            shared_up = up_a & up_b
            shared_down = down_a & down_b
            up_a_down_b = up_a & down_b
            down_a_up_b = down_a & up_b
            only_a = all_a - all_b
            only_b = all_b - all_a

            de_categories = {
                "DE_shared_all": shared_all,
                "DE_shared_up": shared_up,
                "DE_shared_down": shared_down,
                f"DE_up_{short_a}_down_{short_b}": up_a_down_b,
                f"DE_down_{short_a}_up_{short_b}": down_a_up_b,
                f"DE_only_{short_a}": only_a,
                f"DE_only_{short_b}": only_b,
            }

            for sheet_name, gene_set in de_categories.items():
                if not gene_set:
                    # Write empty sheet with header
                    empty_df = pd.DataFrame(columns=[merge_key])
                    sheet_label = sheet_name[:31]
                    empty_df.to_excel(writer, sheet_name=sheet_label, index=False)
                    summary_rows.append({"category": sheet_name, "count": 0})
                    continue

                # Subset both conditions' data for these genes
                sub_a = filt_a[filt_a[merge_key].isin(gene_set)].copy()
                sub_b = filt_b[filt_b[merge_key].isin(gene_set)].copy()

                # Rename non-key columns with suffixes
                rename_a = {c: f"{c}{suffix_a}" for c in sub_a.columns if c != merge_key}
                rename_b = {c: f"{c}{suffix_b}" for c in sub_b.columns if c != merge_key}
                sub_a = sub_a.rename(columns=rename_a)
                sub_b = sub_b.rename(columns=rename_b)

                # Merge side by side
                merged = sub_a.merge(sub_b, on=merge_key, how="outer")

                sheet_label = sheet_name[:31]
                merged.to_excel(writer, sheet_name=sheet_label, index=False)
                summary_rows.append({"category": sheet_name, "count": len(merged)})
                print(f"   {sheet_label}: {len(merged):,} genes")

            # -------------------------------------------------------
            # rMATS sheets (per event type)
            # -------------------------------------------------------
            for et in RMATS_EVENT_TYPES:
                df_a = condition_results[name_a]["rmats_filtered"].get(et)
                df_b = condition_results[name_b]["rmats_filtered"].get(et)

                has_a = df_a is not None and len(df_a) > 0
                has_b = df_b is not None and len(df_b) > 0

                if not has_a and not has_b:
                    continue

                # Prepare DataFrames with event keys
                if has_a:
                    df_a = df_a.copy()
                    df_a["_ekey"] = _make_event_key(df_a, et).values
                    df_a = df_a[df_a["_ekey"] != ""]
                else:
                    df_a = pd.DataFrame(columns=["_ekey"])

                if has_b:
                    df_b = df_b.copy()
                    df_b["_ekey"] = _make_event_key(df_b, et).values
                    df_b = df_b[df_b["_ekey"] != ""]
                else:
                    df_b = pd.DataFrame(columns=["_ekey"])

                events_all_a = set(df_a["_ekey"].dropna().unique())
                events_all_b = set(df_b["_ekey"].dropna().unique())

                # Direction subsets
                events_inc_a = set(
                    df_a.loc[df_a[dpsi_col] >= INCLEVEL_DIFF_CUTOFF, "_ekey"]
                    .dropna().unique()) if has_a and dpsi_col in df_a.columns else set()
                events_inc_b = set(
                    df_b.loc[df_b[dpsi_col] >= INCLEVEL_DIFF_CUTOFF, "_ekey"]
                    .dropna().unique()) if has_b and dpsi_col in df_b.columns else set()
                events_exc_a = set(
                    df_a.loc[df_a[dpsi_col] <= -INCLEVEL_DIFF_CUTOFF, "_ekey"]
                    .dropna().unique()) if has_a and dpsi_col in df_a.columns else set()
                events_exc_b = set(
                    df_b.loc[df_b[dpsi_col] <= -INCLEVEL_DIFF_CUTOFF, "_ekey"]
                    .dropna().unique()) if has_b and dpsi_col in df_b.columns else set()

                rmats_categories = {
                    f"{et}_shared_all": events_all_a & events_all_b,
                    f"{et}_shared_included": events_inc_a & events_inc_b,
                    f"{et}_shared_excluded": events_exc_a & events_exc_b,
                    f"{et}_inc_{short_a}_exc_{short_b}": events_inc_a & events_exc_b,
                    f"{et}_exc_{short_a}_inc_{short_b}": events_exc_a & events_inc_b,
                    f"{et}_only_{short_a}": events_all_a - events_all_b,
                    f"{et}_only_{short_b}": events_all_b - events_all_a,
                }

                # Coordinate columns used as merge keys (not suffixed)
                coord_cols = _COORD_COLS.get(et, [])
                # Shared columns that should not be suffixed
                id_cols = [c for c in [RMATS_COLS["gene_id"], RMATS_COLS["gene_name"]]
                           if (has_a and c in df_a.columns) or
                              (has_b and c in df_b.columns)]
                shared_merge_cols = (
                    [c for c in coord_cols
                     if (has_a and c in df_a.columns) or
                        (has_b and c in df_b.columns)]
                    + id_cols
                )

                for sheet_name, event_set in rmats_categories.items():
                    if not event_set:
                        empty_df = pd.DataFrame(columns=["_ekey"])
                        sheet_label = sheet_name[:31]
                        empty_df.to_excel(writer, sheet_name=sheet_label, index=False)
                        summary_rows.append({"category": sheet_name, "count": 0})
                        continue

                    sub_a = df_a[df_a["_ekey"].isin(event_set)].copy() if has_a else pd.DataFrame()
                    sub_b = df_b[df_b["_ekey"].isin(event_set)].copy() if has_b else pd.DataFrame()

                    # Drop the rMATS ID column (run-specific)
                    rmats_id = RMATS_COLS["event_id"]
                    if rmats_id in sub_a.columns:
                        sub_a = sub_a.drop(columns=[rmats_id])
                    if rmats_id in sub_b.columns:
                        sub_b = sub_b.drop(columns=[rmats_id])

                    # Suffix non-shared columns
                    if len(sub_a) > 0:
                        rename_a = {c: f"{short_a}_{c}" for c in sub_a.columns
                                    if c not in shared_merge_cols and c != "_ekey"}
                        sub_a = sub_a.rename(columns=rename_a)

                    if len(sub_b) > 0:
                        rename_b = {c: f"{short_b}_{c}" for c in sub_b.columns
                                    if c not in shared_merge_cols and c != "_ekey"}
                        sub_b = sub_b.rename(columns=rename_b)

                    # Merge on event key
                    if len(sub_a) > 0 and len(sub_b) > 0:
                        merged = sub_a.merge(sub_b.drop(
                            columns=[c for c in shared_merge_cols if c in sub_b.columns],
                            errors="ignore"),
                            on="_ekey", how="outer")
                    elif len(sub_a) > 0:
                        merged = sub_a
                    else:
                        merged = sub_b

                    # Drop the internal key column before export
                    if "_ekey" in merged.columns:
                        merged = merged.drop(columns=["_ekey"])

                    sheet_label = sheet_name[:31]
                    merged.to_excel(writer, sheet_name=sheet_label, index=False)
                    summary_rows.append({"category": sheet_name, "count": len(merged)})
                    print(f"   {sheet_label}: {len(merged):,} events")

            # -------------------------------------------------------
            # Summary sheet
            # -------------------------------------------------------
            if summary_rows:
                summary_df = pd.DataFrame(summary_rows)
                summary_df.to_excel(writer, sheet_name="Summary", index=False)

        print(f"  Saved: {wb_path}")
