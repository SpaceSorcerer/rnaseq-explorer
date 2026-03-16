"""Interactive splicing visualizations using Plotly.

Provides dPSI volcano plots, event type pie charts, dPSI distributions,
top splicing event bar charts, and gene-by-event-count summaries for rMATS results.
"""

from __future__ import annotations


import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from rnaseq_explorer.viz.theme import (
    PALETTE,
    EVENT_COLORS,
    setup_plotly_theme,
)


def dpsi_volcano(
    df: pd.DataFrame,
    dpsi_col: str = "IncLevelDifference",
    fdr_col: str = "FDR",
    gene_col: str = "GeneID",
    event_type_col: str = "event_type",
    dpsi_cutoff: float = 0.1,
    fdr_cutoff: float = 0.05,
) -> go.Figure:
    """Create a dPSI volcano plot colored by event type.

    Parameters
    ----------
    df : pd.DataFrame
        rMATS results with dPSI, FDR, gene, and event type columns.
    dpsi_col : str
        Column name for delta-PSI (IncLevelDifference).
    fdr_col : str
        Column name for FDR.
    gene_col : str
        Column name for gene identifier.
    event_type_col : str
        Column name for event type (SE, MXE, A3SS, A5SS, RI).
    dpsi_cutoff : float
        Absolute dPSI threshold.
    fdr_cutoff : float
        FDR threshold.

    Returns
    -------
    go.Figure
        Interactive dPSI volcano plot.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="dPSI Volcano — No Data")
        return fig

    plot_df = df.dropna(subset=[dpsi_col, fdr_col]).copy()
    plot_df["neg_log10_fdr"] = -np.log10(plot_df[fdr_col].clip(lower=1e-300))

    fig = go.Figure()

    event_types = plot_df[event_type_col].unique() if event_type_col in plot_df.columns else ["Unknown"]

    for etype in sorted(event_types):
        if event_type_col in plot_df.columns:
            subset = plot_df[plot_df[event_type_col] == etype]
        else:
            subset = plot_df
        if subset.empty:
            continue

        color = EVENT_COLORS.get(etype, PALETTE["neutral"])
        gene_text = subset[gene_col] if gene_col in subset.columns else None

        fig.add_trace(
            go.Scattergl(
                x=subset[dpsi_col],
                y=subset["neg_log10_fdr"],
                mode="markers",
                marker=dict(color=color, size=5, opacity=0.7),
                name=f"{etype} ({len(subset):,})",
                text=gene_text,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "dPSI: %{x:.3f}<br>"
                    "FDR: %{customdata:.2e}<br>"
                    "<extra></extra>"
                ),
                customdata=subset[fdr_col],
            )
        )

    # Threshold lines
    fig.add_hline(
        y=-np.log10(fdr_cutoff), line_dash="dash", line_color="#666666", line_width=1,
        annotation_text=f"FDR={fdr_cutoff}", annotation_position="top right",
    )
    fig.add_vline(x=dpsi_cutoff, line_dash="dash", line_color="#666666", line_width=1)
    fig.add_vline(x=-dpsi_cutoff, line_dash="dash", line_color="#666666", line_width=1)

    fig.update_layout(
        title="Splicing dPSI Volcano Plot",
        xaxis_title="ΔPSI (IncLevelDifference)",
        yaxis_title="-log₁₀ FDR",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def event_type_pie(
    df: pd.DataFrame,
    event_type_col: str = "event_type",
) -> go.Figure:
    """Pie chart showing distribution of splicing event types.

    Parameters
    ----------
    df : pd.DataFrame
        rMATS results.
    event_type_col : str
        Column name for event type.

    Returns
    -------
    go.Figure
        Pie chart of event type counts.
    """
    setup_plotly_theme()

    if df.empty or event_type_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No event type data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Event Types — No Data")
        return fig

    counts = df[event_type_col].value_counts()
    colors = [EVENT_COLORS.get(et, PALETTE["neutral"]) for et in counts.index]

    fig = go.Figure(
        go.Pie(
            labels=counts.index,
            values=counts.values,
            marker=dict(colors=colors),
            textinfo="label+percent+value",
            hovertemplate="<b>%{label}</b><br>Count: %{value:,}<br>%{percent}<extra></extra>",
        )
    )

    fig.update_layout(title="Splicing Event Type Distribution")

    return fig


def dpsi_distribution(
    df: pd.DataFrame,
    dpsi_col: str = "IncLevelDifference",
    event_type_col: str = "event_type",
) -> go.Figure:
    """Histogram of dPSI values faceted by event type.

    Parameters
    ----------
    df : pd.DataFrame
        rMATS results.
    dpsi_col : str
        Column name for dPSI.
    event_type_col : str
        Column name for event type.

    Returns
    -------
    go.Figure
        Faceted histogram of dPSI distributions.
    """
    setup_plotly_theme()

    if df.empty or dpsi_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No dPSI data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="dPSI Distribution — No Data")
        return fig

    event_types = sorted(df[event_type_col].unique()) if event_type_col in df.columns else ["All"]
    n_types = len(event_types)

    fig = make_subplots(
        rows=1, cols=n_types,
        subplot_titles=event_types,
        shared_yaxes=True,
    )

    for i, etype in enumerate(event_types, 1):
        if event_type_col in df.columns:
            subset = df[df[event_type_col] == etype]
        else:
            subset = df

        color = EVENT_COLORS.get(etype, PALETTE["accent1"])
        fig.add_trace(
            go.Histogram(
                x=subset[dpsi_col],
                nbinsx=40,
                marker_color=color,
                opacity=0.8,
                name=etype,
                showlegend=True,
            ),
            row=1, col=i,
        )

    fig.update_layout(
        title="ΔPSI Distribution by Event Type",
        height=350,
        legend=dict(orientation="h", yanchor="bottom", y=1.08, xanchor="right", x=1),
    )
    fig.update_xaxes(title_text="ΔPSI", row=1, col=1)
    fig.update_yaxes(title_text="Count", row=1, col=1)

    return fig


def top_splicing_events(
    df: pd.DataFrame,
    dpsi_col: str = "IncLevelDifference",
    fdr_col: str = "FDR",
    gene_col: str = "GeneID",
    n: int = 20,
) -> go.Figure:
    """Bar chart of top N splicing events ranked by |dPSI|.

    Parameters
    ----------
    df : pd.DataFrame
        rMATS results.
    dpsi_col : str
        Column name for dPSI.
    fdr_col : str
        Column name for FDR.
    gene_col : str
        Column name for gene identifier.
    n : int
        Number of top events to show.

    Returns
    -------
    go.Figure
        Horizontal bar chart of top events.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Top Splicing Events — No Data")
        return fig

    plot_df = df.dropna(subset=[dpsi_col]).copy()
    plot_df["abs_dpsi"] = plot_df[dpsi_col].abs()
    top = plot_df.nlargest(n, "abs_dpsi").sort_values(dpsi_col)

    # Build labels
    if gene_col in top.columns:
        labels = top[gene_col].astype(str)
    else:
        labels = top.index.astype(str)

    # Make labels unique by appending index if duplicated
    if labels.duplicated().any():
        labels = [f"{g} [{i}]" for i, g in enumerate(labels)]

    colors = [PALETTE["up"] if v > 0 else PALETTE["down"] for v in top[dpsi_col]]

    fig = go.Figure(
        go.Bar(
            y=labels,
            x=top[dpsi_col].values,
            orientation="h",
            marker_color=colors,
            hovertemplate="<b>%{y}</b><br>dPSI: %{x:.3f}<extra></extra>",
        )
    )

    fig.update_layout(
        title=f"Top {n} Splicing Events by |ΔPSI|",
        xaxis_title="ΔPSI",
        yaxis_title="",
        height=max(400, n * 22),
    )

    return fig


def genes_by_event_count(
    df: pd.DataFrame,
    gene_col: str = "GeneID",
    event_type_col: str = "event_type",
    n: int = 30,
) -> go.Figure:
    """Stacked bar chart: genes ranked by number of splicing events.

    Parameters
    ----------
    df : pd.DataFrame
        rMATS results.
    gene_col : str
        Column name for gene identifier.
    event_type_col : str
        Column name for event type.
    n : int
        Number of top genes to show.

    Returns
    -------
    go.Figure
        Stacked horizontal bar chart.
    """
    setup_plotly_theme()

    if df.empty or gene_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Genes by Event Count — No Data")
        return fig

    if event_type_col in df.columns:
        ct = pd.crosstab(df[gene_col], df[event_type_col])
    else:
        ct = df.groupby(gene_col).size().to_frame("Events")

    ct["total"] = ct.sum(axis=1)
    ct = ct.nlargest(n, "total").drop(columns=["total"]).sort_values(
        ct.columns[0], ascending=True
    )

    fig = go.Figure()

    for col in ct.columns:
        color = EVENT_COLORS.get(col, PALETTE["neutral"])
        fig.add_trace(
            go.Bar(
                y=ct.index,
                x=ct[col],
                name=col,
                orientation="h",
                marker_color=color,
                hovertemplate=f"<b>%{{y}}</b><br>{col}: %{{x:,}}<extra></extra>",
            )
        )

    fig.update_layout(
        title=f"Top {n} Genes by Splicing Event Count",
        xaxis_title="Number of Events",
        yaxis_title="",
        barmode="stack",
        height=max(400, n * 22),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig
