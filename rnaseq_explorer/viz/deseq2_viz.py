"""Interactive DEG visualizations using Plotly.

Provides volcano plots, MA plots, p-value distributions, log2FC distributions,
top gene bar charts, and biotype breakdowns for DESeq2 results.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from rnaseq_explorer.viz.theme import (
    PALETTE,
    BIOTYPE_COLORS,
    FONT_SIZE_ANNOTATION,
    setup_plotly_theme,
)


def _classify_genes(
    df: pd.DataFrame,
    log2fc_col: str,
    padj_col: str,
    log2fc_cutoff: float,
    padj_cutoff: float,
) -> pd.Series:
    """Classify genes as Up, Down, or NS based on thresholds.

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results.
    log2fc_col : str
        Column name for log2 fold change.
    padj_col : str
        Column name for adjusted p-value.
    log2fc_cutoff : float
        Absolute log2FC threshold.
    padj_cutoff : float
        Adjusted p-value threshold.

    Returns
    -------
    pd.Series
        Classification labels: 'Up', 'Down', or 'NS'.
    """
    conditions = [
        (df[padj_col] < padj_cutoff) & (df[log2fc_col] >= log2fc_cutoff),
        (df[padj_col] < padj_cutoff) & (df[log2fc_col] <= -log2fc_cutoff),
    ]
    choices = ["Up", "Down"]
    return pd.Series(np.select(conditions, choices, default="NS"), index=df.index)


def volcano_plot(
    df: pd.DataFrame,
    log2fc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    gene_col: str = "gene_name",
    log2fc_cutoff: float = 1.0,
    padj_cutoff: float = 0.05,
    genes_of_interest: Optional[Sequence[str]] = None,
) -> go.Figure:
    """Create an interactive volcano plot of DESeq2 results.

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results with log2FC, padj, and gene name columns.
    log2fc_col : str
        Column name for log2 fold change.
    padj_col : str
        Column name for adjusted p-value.
    gene_col : str
        Column name for gene identifiers.
    log2fc_cutoff : float
        Absolute log2FC threshold for significance.
    padj_cutoff : float
        Adjusted p-value threshold for significance.
    genes_of_interest : Sequence[str], optional
        Gene names to label on the plot.

    Returns
    -------
    go.Figure
        Interactive Plotly volcano plot.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Volcano Plot — No Data")
        return fig

    plot_df = df.dropna(subset=[log2fc_col, padj_col]).copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df[padj_col].clip(lower=1e-300))
    plot_df["direction"] = _classify_genes(plot_df, log2fc_col, padj_col, log2fc_cutoff, padj_cutoff)

    color_map = {"Up": PALETTE["up"], "Down": PALETTE["down"], "NS": PALETTE["neutral"]}
    fig = go.Figure()

    for direction in ["NS", "Down", "Up"]:
        subset = plot_df[plot_df["direction"] == direction]
        if subset.empty:
            continue
        n = len(subset)
        fig.add_trace(
            go.Scattergl(
                x=subset[log2fc_col],
                y=subset["neg_log10_padj"],
                mode="markers",
                marker=dict(color=color_map[direction], size=4, opacity=0.6),
                name=f"{direction} ({n:,})",
                text=subset[gene_col] if gene_col in subset.columns else None,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "log2FC: %{x:.3f}<br>"
                    "padj: %{customdata:.2e}<br>"
                    "<extra></extra>"
                ),
                customdata=subset[padj_col],
            )
        )

    # Threshold lines
    max_y = plot_df["neg_log10_padj"].replace([np.inf], np.nan).max()
    max_y = max_y if pd.notna(max_y) else 10
    fig.add_hline(
        y=-np.log10(padj_cutoff), line_dash="dash", line_color="#666666", line_width=1,
        annotation_text=f"padj={padj_cutoff}", annotation_position="top right",
    )
    fig.add_vline(
        x=log2fc_cutoff, line_dash="dash", line_color="#666666", line_width=1,
    )
    fig.add_vline(
        x=-log2fc_cutoff, line_dash="dash", line_color="#666666", line_width=1,
    )

    # Label genes of interest
    if genes_of_interest and gene_col in plot_df.columns:
        goi_df = plot_df[plot_df[gene_col].isin(genes_of_interest)]
        if not goi_df.empty:
            fig.add_trace(
                go.Scatter(
                    x=goi_df[log2fc_col],
                    y=goi_df["neg_log10_padj"],
                    mode="markers+text",
                    marker=dict(color=PALETTE["highlight"], size=8, symbol="diamond",
                                line=dict(color="black", width=1)),
                    text=goi_df[gene_col],
                    textposition="top center",
                    textfont=dict(size=FONT_SIZE_ANNOTATION),
                    name="Genes of Interest",
                    hovertemplate=(
                        "<b>%{text}</b><br>"
                        "log2FC: %{x:.3f}<br>"
                        "padj: %{customdata:.2e}<br>"
                        "<extra></extra>"
                    ),
                    customdata=goi_df[padj_col],
                )
            )

    n_up = int((plot_df["direction"] == "Up").sum())
    n_down = int((plot_df["direction"] == "Down").sum())
    fig.update_layout(
        title=f"Volcano Plot — {n_up:,} Up, {n_down:,} Down",
        xaxis_title="log₂ Fold Change",
        yaxis_title="-log₁₀ adjusted p-value",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def ma_plot(
    df: pd.DataFrame,
    basemean_col: str = "baseMean",
    log2fc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    gene_col: str = "gene_name",
    log2fc_cutoff: float = 1.0,
    padj_cutoff: float = 0.05,
) -> go.Figure:
    """Create an interactive MA plot (log10 baseMean vs log2FC).

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results.
    basemean_col : str
        Column name for baseMean.
    log2fc_col : str
        Column name for log2 fold change.
    padj_col : str
        Column name for adjusted p-value.
    gene_col : str
        Column name for gene identifiers.
    log2fc_cutoff : float
        Absolute log2FC threshold.
    padj_cutoff : float
        Adjusted p-value threshold.

    Returns
    -------
    go.Figure
        Interactive Plotly MA plot.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="MA Plot — No Data")
        return fig

    plot_df = df.dropna(subset=[basemean_col, log2fc_col, padj_col]).copy()
    plot_df = plot_df[plot_df[basemean_col] > 0].copy()
    plot_df["log10_basemean"] = np.log10(plot_df[basemean_col])
    plot_df["direction"] = _classify_genes(plot_df, log2fc_col, padj_col, log2fc_cutoff, padj_cutoff)

    color_map = {"Up": PALETTE["up"], "Down": PALETTE["down"], "NS": PALETTE["neutral"]}
    fig = go.Figure()

    for direction in ["NS", "Down", "Up"]:
        subset = plot_df[plot_df["direction"] == direction]
        if subset.empty:
            continue
        fig.add_trace(
            go.Scattergl(
                x=subset["log10_basemean"],
                y=subset[log2fc_col],
                mode="markers",
                marker=dict(color=color_map[direction], size=4, opacity=0.6),
                name=f"{direction} ({len(subset):,})",
                text=subset[gene_col] if gene_col in subset.columns else None,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "baseMean: %{customdata:.1f}<br>"
                    "log2FC: %{y:.3f}<br>"
                    "<extra></extra>"
                ),
                customdata=subset[basemean_col],
            )
        )

    fig.add_hline(y=0, line_color="#333333", line_width=1)
    fig.add_hline(y=log2fc_cutoff, line_dash="dash", line_color="#666666", line_width=1)
    fig.add_hline(y=-log2fc_cutoff, line_dash="dash", line_color="#666666", line_width=1)

    fig.update_layout(
        title="MA Plot",
        xaxis_title="log₁₀ Mean Expression",
        yaxis_title="log₂ Fold Change",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def pvalue_distribution(
    df: pd.DataFrame,
    pvalue_col: str = "pvalue",
    padj_col: str = "padj",
) -> go.Figure:
    """Plot overlaid histograms of raw and adjusted p-values.

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results.
    pvalue_col : str
        Column name for raw p-values.
    padj_col : str
        Column name for adjusted p-values.

    Returns
    -------
    go.Figure
        Plotly figure with overlaid histograms.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="P-value Distribution — No Data")
        return fig

    fig = go.Figure()

    if pvalue_col in df.columns:
        pvals = df[pvalue_col].dropna()
        fig.add_trace(
            go.Histogram(
                x=pvals,
                nbinsx=50,
                name="Raw p-value",
                marker_color=PALETTE["accent1"],
                opacity=0.7,
            )
        )

    if padj_col in df.columns:
        padj_vals = df[padj_col].dropna()
        fig.add_trace(
            go.Histogram(
                x=padj_vals,
                nbinsx=50,
                name="Adjusted p-value",
                marker_color=PALETTE["accent4"],
                opacity=0.7,
            )
        )

    fig.update_layout(
        title="P-value Distribution",
        xaxis_title="P-value",
        yaxis_title="Count",
        barmode="overlay",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def log2fc_distribution(
    df: pd.DataFrame,
    log2fc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    padj_cutoff: float = 0.05,
) -> go.Figure:
    """Plot histogram of log2FC split by significance.

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results.
    log2fc_col : str
        Column name for log2 fold change.
    padj_col : str
        Column name for adjusted p-value.
    padj_cutoff : float
        Adjusted p-value threshold to classify significance.

    Returns
    -------
    go.Figure
        Plotly figure with stacked histograms.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="log2FC Distribution — No Data")
        return fig

    plot_df = df.dropna(subset=[log2fc_col, padj_col]).copy()
    sig = plot_df[plot_df[padj_col] < padj_cutoff]
    ns = plot_df[plot_df[padj_col] >= padj_cutoff]

    fig = go.Figure()

    if not ns.empty:
        fig.add_trace(
            go.Histogram(
                x=ns[log2fc_col],
                nbinsx=60,
                name=f"Not significant ({len(ns):,})",
                marker_color=PALETTE["neutral"],
                opacity=0.6,
            )
        )

    if not sig.empty:
        fig.add_trace(
            go.Histogram(
                x=sig[log2fc_col],
                nbinsx=60,
                name=f"Significant ({len(sig):,})",
                marker_color=PALETTE["up"],
                opacity=0.7,
            )
        )

    fig.update_layout(
        title="log₂ Fold Change Distribution",
        xaxis_title="log₂ Fold Change",
        yaxis_title="Count",
        barmode="overlay",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def top_genes_bar(
    df: pd.DataFrame,
    log2fc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    gene_col: str = "gene_name",
    n: int = 20,
) -> go.Figure:
    """Bar chart of top N upregulated and downregulated genes by |log2FC|.

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results (typically already significance-filtered).
    log2fc_col : str
        Column name for log2 fold change.
    padj_col : str
        Column name for adjusted p-value.
    gene_col : str
        Column name for gene identifiers.
    n : int
        Number of top genes per direction.

    Returns
    -------
    go.Figure
        Horizontal bar chart.
    """
    setup_plotly_theme()

    if df.empty or gene_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Top Genes — No Data")
        return fig

    plot_df = df.dropna(subset=[log2fc_col, gene_col]).copy()
    plot_df["abs_log2fc"] = plot_df[log2fc_col].abs()

    up = plot_df[plot_df[log2fc_col] > 0].nlargest(n, "abs_log2fc")
    down = plot_df[plot_df[log2fc_col] < 0].nlargest(n, "abs_log2fc")

    combined = pd.concat([up, down]).sort_values(log2fc_col)

    colors = [PALETTE["up"] if v > 0 else PALETTE["down"] for v in combined[log2fc_col]]

    fig = go.Figure(
        go.Bar(
            y=combined[gene_col],
            x=combined[log2fc_col],
            orientation="h",
            marker_color=colors,
            hovertemplate=(
                "<b>%{y}</b><br>"
                "log2FC: %{x:.3f}<br>"
                "<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title=f"Top {n} Upregulated & Downregulated Genes",
        xaxis_title="log₂ Fold Change",
        yaxis_title="",
        height=max(400, len(combined) * 22),
    )

    return fig


def biotype_breakdown(
    df: pd.DataFrame,
    biotype_col: str = "biotype_group",
    direction_col: str = "direction",
) -> go.Figure:
    """Stacked bar chart of biotype distribution per direction (Up/Down).

    Parameters
    ----------
    df : pd.DataFrame
        DESeq2 results with biotype and direction columns.
    biotype_col : str
        Column name for biotype group.
    direction_col : str
        Column name for direction (Up/Down).

    Returns
    -------
    go.Figure
        Stacked bar chart.
    """
    setup_plotly_theme()

    if df.empty or biotype_col not in df.columns or direction_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No biotype data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Biotype Breakdown — No Data")
        return fig

    ct = pd.crosstab(df[direction_col], df[biotype_col])

    fig = go.Figure()
    for biotype in ct.columns:
        color = BIOTYPE_COLORS.get(biotype, PALETTE["neutral"])
        fig.add_trace(
            go.Bar(
                x=ct.index,
                y=ct[biotype],
                name=biotype,
                marker_color=color,
                hovertemplate=f"<b>{biotype}</b><br>%{{x}}: %{{y:,}}<extra></extra>",
            )
        )

    fig.update_layout(
        title="Biotype Distribution by Direction",
        xaxis_title="Direction",
        yaxis_title="Gene Count",
        barmode="stack",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig
