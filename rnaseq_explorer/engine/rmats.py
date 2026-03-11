"""rMATS event parsing, filtering, and visualization.

Handles all 5 rMATS event types: SE, A3SS, A5SS, RI, MXE.
Provides loading, column validation, filtering (single/dual mode),
and scatter/volcano/distribution plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from rnaseq_explorer.engine.deseq2 import (
    _resolve_column,
    _strip_ensembl_version,
    load_file,
    validate_columns,
)
from rnaseq_explorer.viz.theme import (
    COLOR_NS,
    COLOR_UP,
    EVENT_COLORS,
    add_count_box,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

RMATS_EVENT_TYPES: list[str] = ["SE", "A3SS", "A5SS", "RI", "MXE"]
RMATS_FILE_SUFFIX: str = ".MATS.JCEC.txt"

DEFAULT_RMATS_COLS: dict[str, str] = {
    "event_id": "ID",
    "gene_id": "GeneID",
    "gene_name": "geneSymbol",
    "pvalue": "PValue",
    "fdr": "FDR",
    "inclevel_diff": "IncLevelDifference",
}

# Column-name aliases for auto-detection
_RMATS_ALIASES: dict[str, list[str]] = {
    "event_id": ["id", "event_id", "eventid"],
    "gene_id": ["geneid", "gene_id", "ensembl_gene_id"],
    "gene_name": ["genesymbol", "gene_name", "gene_symbol", "symbol"],
    "pvalue": ["pvalue", "pval", "p.value", "p_value"],
    "fdr": ["fdr", "padj", "adj.p.val", "q.value", "qvalue"],
    "inclevel_diff": [
        "incleveldifference", "inclevel_diff", "deltapsi", "dpsi", "delta_psi",
    ],
}

# Coordinate columns that uniquely identify each splicing event
COORD_COLS: dict[str, list[str]] = {
    "SE": [
        "chr", "strand", "exonStart_0base", "exonEnd",
        "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
    ],
    "A3SS": [
        "chr", "strand", "longExonStart_0base", "longExonEnd",
        "shortES", "shortEE", "flankingES", "flankingEE",
    ],
    "A5SS": [
        "chr", "strand", "longExonStart_0base", "longExonEnd",
        "shortES", "shortEE", "flankingES", "flankingEE",
    ],
    "RI": [
        "chr", "strand", "riExonStart_0base", "riExonEnd",
        "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
    ],
    "MXE": [
        "chr", "strand", "1stExonStart_0base", "1stExonEnd",
        "2ndExonStart_0base", "2ndExonEnd",
        "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
    ],
}


# ---------------------------------------------------------------------------
# Column validation
# ---------------------------------------------------------------------------


def _validate_rmats_columns(
    df: pd.DataFrame, col_mapping: dict[str, str]
) -> pd.DataFrame:
    """Validate that required rMATS columns exist, trying common aliases.

    Parameters
    ----------
    df : pd.DataFrame
        Raw rMATS data.
    col_mapping : dict
        Expected column mapping (like DEFAULT_RMATS_COLS).

    Returns
    -------
    pd.DataFrame
        Possibly renamed DataFrame.

    Raises
    ------
    ValueError
        If a required column cannot be found.
    """
    alias_map = {
        "FDR": ["FDR", "fdr", "adj.P.Val", "padj", "q-value"],
        "PValue": ["PValue", "pvalue", "P.Value", "p_value", "p-value"],
        "IncLevelDifference": [
            "IncLevelDifference", "IncLevel_Difference", "dPSI", "inc_level_diff",
        ],
        "geneSymbol": ["geneSymbol", "GeneSymbol", "gene_symbol", "geneName", "gene_name"],
        "GeneID": ["GeneID", "gene_id", "Ensembl_ID"],
        "ID": ["ID", "id", "event_id"],
    }
    rename = {}
    for _key, expected in col_mapping.items():
        if not expected or expected in df.columns:
            continue
        aliases = alias_map.get(expected, [])
        found = False
        for alias in aliases:
            if alias in df.columns:
                rename[alias] = expected
                print(f"  [WARN] rMATS column '{alias}' mapped to '{expected}'")
                found = True
                break
        if not found:
            raise ValueError(
                f"Required rMATS column '{expected}' not found. "
                f"Available columns: {list(df.columns)}. "
                f"Check your rMATS output format."
            )
    if rename:
        df = df.rename(columns=rename)
    return df


def normalize_rmats_columns(
    df: pd.DataFrame,
    cols: dict[str, str],
    file_label: str = "rMATS file",
) -> pd.DataFrame:
    """Rename df columns to match the names in *cols*, strip Ensembl versions.

    Parameters
    ----------
    df : pd.DataFrame
        Raw rMATS data.
    cols : dict
        Column name mapping.
    file_label : str
        Label for log messages.

    Returns
    -------
    pd.DataFrame
        DataFrame with normalized column names.
    """
    rename_map = {}
    for key, configured in cols.items():
        if not configured:
            continue
        if configured in df.columns:
            continue
        actual = _resolve_column(df, key, configured, _RMATS_ALIASES, file_label)
        if actual and actual != configured:
            rename_map[actual] = configured

    if rename_map:
        df = df.rename(columns=rename_map)

    gene_id_col = cols.get("gene_id", "")
    if gene_id_col and gene_id_col in df.columns:
        sample = df[gene_id_col].dropna().astype(str).head(500)
        ens_frac = sample.str.upper().str.startswith("ENS").sum() / max(len(sample), 1)
        if ens_frac > 0.1:
            df = df.copy()
            df[gene_id_col] = _strip_ensembl_version(df[gene_id_col])

    return df


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_all_rmats(
    rmats_dir: str | Path,
    cols: dict[str, str] | None = None,
    event_types: list[str] | None = None,
    file_suffix: str = RMATS_FILE_SUFFIX,
) -> dict[str, pd.DataFrame]:
    """Load all rMATS event type files from a directory.

    Parameters
    ----------
    rmats_dir : str or Path
        Directory containing rMATS output files.
    cols : dict or None
        Column name mapping. Defaults to DEFAULT_RMATS_COLS.
    event_types : list[str] or None
        Event types to load. Defaults to all 5.
    file_suffix : str
        File suffix pattern for rMATS files.

    Returns
    -------
    dict
        Mapping of event_type -> DataFrame.
    """
    if cols is None:
        cols = DEFAULT_RMATS_COLS
    if event_types is None:
        event_types = RMATS_EVENT_TYPES

    rmats_dir = Path(rmats_dir)
    all_data: dict[str, pd.DataFrame] = {}

    for event_type in event_types:
        filepath = rmats_dir / f"{event_type}{file_suffix}"
        if not filepath.exists():
            print(f"  WARNING: {filepath.name} not found, skipping {event_type}")
            continue

        df = load_file(filepath, name=f"rMATS {event_type}")
        df = normalize_rmats_columns(df, cols, f"rMATS {event_type}")
        required = [v for v in cols.values() if v is not None]
        validate_columns(df, required, name=f"rMATS {event_type}")

        critical = ["FDR", "PValue", "IncLevelDifference", "geneSymbol"]
        status_parts = [
            f"{c} \u2713" if c in df.columns else f"{c} \u2717" for c in critical
        ]
        print(f"  [INFO] {event_type} columns: {', '.join(status_parts)}")

        gene_col = cols["gene_name"]
        id_col = cols["gene_id"]
        if gene_col in df.columns and id_col in df.columns:
            mask = df[gene_col].isna() | (df[gene_col].str.strip() == "")
            if mask.any():
                df.loc[mask, gene_col] = df.loc[mask, id_col]
                print(
                    f"  [INFO] Filled {mask.sum()} missing geneSymbol values with GeneID"
                )

        df["event_type"] = event_type
        all_data[event_type] = df

    return all_data


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------


def filter_rmats(
    df: pd.DataFrame,
    cols: dict[str, str] | None = None,
    fdr_cutoff: float = 0.05,
    pval_cutoff: float = 0.05,
    dpsi_cutoff: float = 0.1,
    use_fdr: bool = True,
    dual_filter: bool = False,
    event_type: str = "",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply rMATS cutoffs and return (raw, filtered) DataFrames.

    Parameters
    ----------
    df : pd.DataFrame
        Raw rMATS data for one event type.
    cols : dict or None
        Column name mapping.
    fdr_cutoff : float
        FDR threshold.
    pval_cutoff : float
        P-value threshold.
    dpsi_cutoff : float
        Absolute IncLevelDifference threshold.
    use_fdr : bool
        If True, filter by FDR; otherwise by PValue.
    dual_filter : bool
        If True, require both FDR AND PValue thresholds.
    event_type : str
        Event type label for log messages.

    Returns
    -------
    tuple
        (cleaned_raw_df, filtered_df)
    """
    if cols is None:
        cols = DEFAULT_RMATS_COLS

    df = _validate_rmats_columns(df, cols)

    if dual_filter:
        drop_cols = [
            c for c in [cols["fdr"], cols["pvalue"], cols["inclevel_diff"]]
            if c in df.columns
        ]
        df = df.dropna(subset=drop_cols)
        mask = (
            (df[cols["fdr"]] < fdr_cutoff)
            & (df[cols["pvalue"]] < pval_cutoff)
            & (df[cols["inclevel_diff"]].abs() >= dpsi_cutoff)
        )
        filtered = df[mask].copy()
        print(
            f"  {event_type}: {len(df):,} total -> {len(filtered):,} significant "
            f"(FDR < {fdr_cutoff} AND PValue < {pval_cutoff}, |dPSI| >= {dpsi_cutoff})"
        )
    else:
        pval_col = cols["fdr"] if use_fdr else cols["pvalue"]
        pval_threshold = fdr_cutoff if use_fdr else pval_cutoff
        pval_label = "FDR" if use_fdr else "PValue"
        df = df.dropna(subset=[pval_col, cols["inclevel_diff"]])
        mask = (
            (df[pval_col] < pval_threshold)
            & (df[cols["inclevel_diff"]].abs() >= dpsi_cutoff)
        )
        filtered = df[mask].copy()
        print(
            f"  {event_type}: {len(df):,} total -> {len(filtered):,} significant "
            f"({pval_label} < {pval_threshold}, |dPSI| >= {dpsi_cutoff})"
        )

    return df, filtered


# ---------------------------------------------------------------------------
# Event key helper
# ---------------------------------------------------------------------------


def make_event_key(df: pd.DataFrame, event_type: str) -> pd.Series:
    """Create a unique event key from genomic coordinate columns.

    Returns a pandas Series of strings in the format
    ``chr:strand:col1:col2:...`` for each row.

    Parameters
    ----------
    df : pd.DataFrame
        rMATS data for a single event type.
    event_type : str
        One of RMATS_EVENT_TYPES.

    Returns
    -------
    pd.Series
        String keys for each event.
    """
    coord_cols = COORD_COLS.get(event_type, [])
    if not coord_cols:
        return pd.Series([""] * len(df), index=df.index)

    missing = [c for c in coord_cols if c not in df.columns]
    if missing:
        return pd.Series([""] * len(df), index=df.index)

    key = df[coord_cols[0]].astype(str)
    for col in coord_cols[1:]:
        key = key + ":" + df[col].astype(str)
    return key


def parse_inclevel_mean(series: pd.Series) -> pd.Series:
    """Parse IncLevel1/IncLevel2 comma-separated PSI values to mean.

    Parameters
    ----------
    series : pd.Series
        rMATS IncLevel column values (comma-separated PSI per replicate).

    Returns
    -------
    pd.Series
        Mean PSI values.
    """
    def _parse(val):
        if pd.isna(val):
            return np.nan
        parts = str(val).split(",")
        nums = []
        for p in parts:
            p = p.strip()
            if p and p != "NA":
                try:
                    nums.append(float(p))
                except ValueError:
                    pass
        return np.mean(nums) if nums else np.nan

    return series.apply(_parse)


# ---------------------------------------------------------------------------
# Per-condition rMATS visualization
# ---------------------------------------------------------------------------


def rmats_scatter(
    df: pd.DataFrame,
    event_type: str,
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    fdr_cutoff: float = 0.05,
    pval_cutoff: float = 0.05,
    dpsi_cutoff: float = 0.1,
    use_fdr: bool = True,
    dual_filter: bool = False,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Scatter plot: IncLevelDifference vs -log10(pvalue) for one event type.

    Parameters
    ----------
    df : pd.DataFrame
        Raw rMATS data for one event type.
    event_type : str
        Splicing event type (SE, A3SS, etc.).
    outdir : str or Path
        Output directory for the figure.
    rmats_cols : dict or None
        Column name mapping.
    fdr_cutoff, pval_cutoff, dpsi_cutoff : float
        Significance thresholds.
    use_fdr : bool
        If True, use FDR on y-axis; else PValue.
    dual_filter : bool
        If True, require both FDR AND PValue.
    fig_format : str
        Output format (png, svg, pdf).
    fig_dpi : int
        Resolution.
    """
    cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)

    if dual_filter:
        pval_col = cols["fdr"]
        pval_threshold = fdr_cutoff
        pval_label = "FDR"
    else:
        pval_col = cols["fdr"] if use_fdr else cols["pvalue"]
        pval_threshold = fdr_cutoff if use_fdr else pval_cutoff
        pval_label = "FDR" if use_fdr else "PValue"

    data = df.dropna(subset=[pval_col, cols["inclevel_diff"]]).copy()
    data["-log10p"] = -np.log10(data[pval_col].clip(lower=1e-300))

    if dual_filter:
        sig = (
            (data[cols["fdr"]] < fdr_cutoff)
            & (data[cols["pvalue"]] < pval_cutoff)
            & (data[cols["inclevel_diff"]].abs() >= dpsi_cutoff)
        )
    else:
        sig = (
            (data[pval_col] < pval_threshold)
            & (data[cols["inclevel_diff"]].abs() >= dpsi_cutoff)
        )
    data["significant"] = np.where(sig, "Significant", "NS")

    fig, ax = plt.subplots(figsize=(8, 6))

    for status, color in [
        ("NS", COLOR_NS),
        ("Significant", EVENT_COLORS.get(event_type, COLOR_UP)),
    ]:
        subset = data[data["significant"] == status]
        lbl = "NS" if status == "NS" else f"{status} ({len(subset):,})"
        ax.scatter(
            subset[cols["inclevel_diff"]],
            subset["-log10p"],
            c=color,
            s=10,
            alpha=0.5,
            edgecolors="none",
            label=lbl,
            rasterized=True,
        )

    ax.axhline(-np.log10(pval_threshold), color="grey", ls="--", lw=0.8)
    ax.axvline(dpsi_cutoff, color="grey", ls="--", lw=0.8)
    ax.axvline(-dpsi_cutoff, color="grey", ls="--", lw=0.8)

    sig_data = data[data["significant"] == "Significant"]
    n_inc = (sig_data[cols["inclevel_diff"]] >= dpsi_cutoff).sum()
    n_exc = (sig_data[cols["inclevel_diff"]] <= -dpsi_cutoff).sum()
    add_count_box(
        ax,
        n_inc,
        n_exc,
        n_inc + n_exc,
        position="lower left",
        up_label=f"Included (dPSI\u22650.1)",
        down_label=f"Excluded (dPSI\u2264\u22120.1)",
    )

    ax.set_xlabel("$\\Delta$PSI (IncLevelDifference)")
    ax.set_ylabel(f"-log$_{{10}}$ ({pval_label})")
    ax.set_title(
        f"rMATS - {event_type} (Skipped Exon)"
        if event_type == "SE"
        else f"rMATS - {event_type}"
    )
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        fontsize=10,
        markerscale=2,
    )

    outpath = outdir / f"rmats_{event_type}_scatter.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def rmats_combined_volcano(
    all_data: dict[str, pd.DataFrame],
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    fdr_cutoff: float = 0.05,
    pval_cutoff: float = 0.05,
    dpsi_cutoff: float = 0.1,
    use_fdr: bool = True,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Combined scatter: all event types overlaid on one plot.

    Parameters
    ----------
    all_data : dict
        Mapping of event_type -> raw DataFrame.
    outdir : str or Path
        Output directory.
    rmats_cols : dict or None
        Column name mapping.
    fdr_cutoff, pval_cutoff, dpsi_cutoff : float
        Significance thresholds.
    use_fdr : bool
        Whether to use FDR or PValue.
    fig_format : str
        Output format.
    fig_dpi : int
        Resolution.
    """
    cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)

    pval_col = cols["fdr"] if use_fdr else cols["pvalue"]
    pval_threshold = fdr_cutoff if use_fdr else pval_cutoff
    pval_label = "FDR" if use_fdr else "PValue"

    fig, ax = plt.subplots(figsize=(10, 7))

    for event_type, df in all_data.items():
        data = df.dropna(subset=[pval_col, cols["inclevel_diff"]]).copy()
        data["-log10p"] = -np.log10(data[pval_col].clip(lower=1e-300))

        sig = (
            (data[pval_col] < pval_threshold)
            & (data[cols["inclevel_diff"]].abs() >= dpsi_cutoff)
        )
        sig_data = data[sig]
        ns_data = data[~sig]

        ax.scatter(
            ns_data[cols["inclevel_diff"]],
            ns_data["-log10p"],
            c=COLOR_NS,
            s=6,
            alpha=0.15,
            edgecolors="none",
            rasterized=True,
        )
        ax.scatter(
            sig_data[cols["inclevel_diff"]],
            sig_data["-log10p"],
            c=EVENT_COLORS.get(event_type, COLOR_UP),
            s=12,
            alpha=0.7,
            edgecolors="none",
            label=f"{event_type} ({len(sig_data):,})",
            rasterized=True,
        )

    ax.axhline(-np.log10(pval_threshold), color="grey", ls="--", lw=0.8)
    ax.axvline(dpsi_cutoff, color="grey", ls="--", lw=0.8)
    ax.axvline(-dpsi_cutoff, color="grey", ls="--", lw=0.8)

    # Per-event-type breakdown box
    lines = []
    grand_inc = 0
    grand_exc = 0
    for et, df in all_data.items():
        d = df.dropna(subset=[pval_col, cols["inclevel_diff"]])
        s = (d[pval_col] < pval_threshold) & (
            d[cols["inclevel_diff"]].abs() >= dpsi_cutoff
        )
        sig_d = d[s]
        n_inc = int((sig_d[cols["inclevel_diff"]] >= dpsi_cutoff).sum())
        n_exc = int((sig_d[cols["inclevel_diff"]] <= -dpsi_cutoff).sum())
        grand_inc += n_inc
        grand_exc += n_exc
        lines.append(f"{et}: {n_inc:,} inc / {n_exc:,} exc")
    lines.append(
        f"Total: {grand_inc + grand_exc:,} ({grand_inc:,} inc / {grand_exc:,} exc)"
    )
    box_text = "\n".join(lines)
    ax.text(
        0.02,
        0.02,
        box_text,
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
        ha="left",
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="grey",
            alpha=0.9,
        ),
    )

    ax.set_xlabel("$\\Delta$PSI (IncLevelDifference)")
    ax.set_ylabel(f"-log$_{{10}}$ ({pval_label})")
    ax.set_title("rMATS - All Splicing Event Types")
    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        fontsize=10,
        markerscale=2,
        title="Significant Events",
    )

    outpath = outdir / f"rmats_all_events_scatter.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def rmats_event_summary_chart(
    filtered_counts: dict[str, int],
    outdir: str | Path,
    use_fdr: bool = True,
    fdr_cutoff: float = 0.05,
    pval_cutoff: float = 0.05,
    dpsi_cutoff: float = 0.1,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Bar chart comparing significant event counts across all event types.

    Parameters
    ----------
    filtered_counts : dict
        Mapping of event_type -> count of significant events.
    outdir : str or Path
        Output directory.
    use_fdr : bool
        Whether FDR or PValue was used for filtering.
    fdr_cutoff, pval_cutoff, dpsi_cutoff : float
        Thresholds (for labeling the title).
    fig_format : str
        Output format.
    fig_dpi : int
        Resolution.
    """
    outdir = Path(outdir)
    event_types = list(filtered_counts.keys())
    counts = list(filtered_counts.values())
    colors = [EVENT_COLORS.get(et, "#888888") for et in event_types]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(event_types, counts, color=colors, edgecolor="black", linewidth=0.5)

    for bar, val in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(counts) * 0.02,
            str(val),
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

    pval_label = "FDR" if use_fdr else "PValue"
    pval_threshold = fdr_cutoff if use_fdr else pval_cutoff
    ax.set_xlabel("Splicing Event Type")
    ax.set_ylabel("Number of Significant Events")
    ax.set_title(
        f"rMATS - Significant Events by Type\n"
        f"({pval_label} < {pval_threshold}, "
        f"|$\\Delta$PSI| >= {dpsi_cutoff})"
    )

    outpath = outdir / f"rmats_event_type_summary.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def rmats_dpsi_distribution(
    all_filtered: dict[str, pd.DataFrame],
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    dpsi_cutoff: float = 0.1,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Violin/box plot of delta-PSI distributions across event types.

    Parameters
    ----------
    all_filtered : dict
        Mapping of event_type -> filtered DataFrame.
    outdir : str or Path
        Output directory.
    rmats_cols : dict or None
        Column name mapping.
    dpsi_cutoff : float
        IncLevelDifference cutoff for annotation.
    fig_format : str
        Output format.
    fig_dpi : int
        Resolution.
    """
    cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)

    plot_data = []
    for event_type, df in all_filtered.items():
        if len(df) == 0:
            continue
        tmp = df[[cols["inclevel_diff"]]].copy()
        tmp["Event Type"] = event_type
        plot_data.append(tmp)

    if not plot_data:
        print("  No significant events to plot dPSI distribution")
        return

    combined = pd.concat(plot_data, ignore_index=True)
    combined.rename(columns={cols["inclevel_diff"]: "dPSI"}, inplace=True)

    fig, ax = plt.subplots(figsize=(9, 6))
    event_order = [
        et
        for et in RMATS_EVENT_TYPES
        if et in all_filtered and len(all_filtered[et]) > 0
    ]
    palette = [EVENT_COLORS.get(et, "#888888") for et in event_order]

    sns.violinplot(
        data=combined,
        x="Event Type",
        y="dPSI",
        order=event_order,
        palette=palette,
        inner=None,
        alpha=0.3,
        ax=ax,
    )
    sns.stripplot(
        data=combined,
        x="Event Type",
        y="dPSI",
        order=event_order,
        palette=palette,
        size=3,
        alpha=0.6,
        jitter=True,
        ax=ax,
    )

    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(dpsi_cutoff, color="grey", ls="--", lw=0.6, alpha=0.5)
    ax.axhline(-dpsi_cutoff, color="grey", ls="--", lw=0.6, alpha=0.5)

    # Per-event inc/exc counts
    y_max = combined["dPSI"].max()
    y_min = combined["dPSI"].min()
    y_pad = (y_max - y_min) * 0.08
    for i, et in enumerate(event_order):
        et_data = combined[combined["Event Type"] == et]["dPSI"]
        n_inc = int((et_data >= dpsi_cutoff).sum())
        n_exc = int((et_data <= -dpsi_cutoff).sum())
        ax.text(
            i,
            y_max + y_pad,
            f"{n_inc} inc",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color=EVENT_COLORS.get(et, "#333333"),
        )
        ax.text(
            i,
            y_min - y_pad,
            f"{n_exc} exc",
            ha="center",
            va="top",
            fontsize=9,
            fontweight="bold",
            color=EVENT_COLORS.get(et, "#333333"),
        )
    ax.set_ylim(y_min - y_pad * 3, y_max + y_pad * 3)

    ax.set_ylabel("$\\Delta$PSI (IncLevelDifference)")
    ax.set_title("Distribution of $\\Delta$PSI by Event Type (Significant Events)")

    outpath = outdir / f"rmats_dpsi_distribution.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")


def rmats_psi_scatter(
    rmats_raw: dict[str, pd.DataFrame],
    rmats_filtered: dict[str, pd.DataFrame],
    event_type: str,
    outdir: str | Path,
    rmats_cols: dict[str, str] | None = None,
    dpsi_cutoff: float = 0.1,
    fig_format: str = "png",
    fig_dpi: int = 300,
) -> None:
    """Mean PSI (IncLevel1 vs IncLevel2) scatter showing direction of splicing shift.

    Parameters
    ----------
    rmats_raw : dict
        Mapping of event_type -> raw DataFrame.
    rmats_filtered : dict
        Mapping of event_type -> filtered DataFrame.
    event_type : str
        Event type to plot.
    outdir : str or Path
        Output directory.
    rmats_cols : dict or None
        Column name mapping.
    dpsi_cutoff : float
        IncLevelDifference cutoff for offset lines.
    fig_format : str
        Output format.
    fig_dpi : int
        Resolution.
    """
    cols = rmats_cols or DEFAULT_RMATS_COLS
    outdir = Path(outdir)

    df = rmats_raw.get(event_type, pd.DataFrame())
    if len(df) == 0:
        return

    if "IncLevel1" not in df.columns or "IncLevel2" not in df.columns:
        print(
            f"  Skipping PSI scatter ({event_type}): "
            "IncLevel1/IncLevel2 columns not found"
        )
        return

    df = df.copy()
    df["_psi1"] = parse_inclevel_mean(df["IncLevel1"])
    df["_psi2"] = parse_inclevel_mean(df["IncLevel2"])
    df = df.dropna(subset=["_psi1", "_psi2"])
    if len(df) == 0:
        return

    id_col = cols["event_id"]
    filt_df = rmats_filtered.get(event_type, pd.DataFrame())
    sig_ids = (
        set(filt_df[id_col].unique())
        if len(filt_df) > 0 and id_col in filt_df.columns
        else set()
    )
    df["_sig"] = df[id_col].isin(sig_ids) if id_col in df.columns else False

    color = EVENT_COLORS.get(event_type, "#333333")
    fig, ax = plt.subplots(figsize=(6, 6))

    ns_mask = ~df["_sig"]
    sig_mask = df["_sig"]
    ax.scatter(
        df.loc[ns_mask, "_psi1"],
        df.loc[ns_mask, "_psi2"],
        c=COLOR_NS,
        s=5,
        alpha=0.2,
        edgecolors="none",
        rasterized=True,
        label="NS",
    )
    ax.scatter(
        df.loc[sig_mask, "_psi1"],
        df.loc[sig_mask, "_psi2"],
        c=color,
        s=8,
        alpha=0.7,
        edgecolors="none",
        rasterized=True,
        label=f"Significant ({sig_mask.sum():,})",
    )

    # Diagonal and cutoff offset lines
    ax.plot([0, 1], [0, 1], "k--", lw=0.7, alpha=0.5)
    ax.plot(
        [0, 1 - dpsi_cutoff],
        [dpsi_cutoff, 1],
        color="grey",
        ls=":",
        lw=0.6,
        alpha=0.6,
    )
    ax.plot(
        [dpsi_cutoff, 1],
        [0, 1 - dpsi_cutoff],
        color="grey",
        ls=":",
        lw=0.6,
        alpha=0.6,
    )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Mean PSI \u2014 Sample Group 1")
    ax.set_ylabel("Mean PSI \u2014 Sample Group 2")
    ax.set_title(f"PSI Shift \u2014 {event_type}")
    ax.legend(fontsize=9, markerscale=2)

    outpath = outdir / f"rmats_{event_type}_psi_scatter.{fig_format}"
    fig.savefig(outpath, format=fig_format, dpi=fig_dpi)
    plt.close(fig)
    print(f"  Saved: {outpath}")
