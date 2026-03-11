"""Interactive QC visualizations using Plotly.

Provides PCA plots, sample correlation heatmaps, and top DEG expression
heatmaps for quality control assessment of RNA-seq experiments.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from rnaseq_explorer.viz.theme import (
    PALETTE,
    CONDITION_COLORS,
    FONT_SIZE_TITLE,
    condition_color_map,
    setup_plotly_theme,
)


def pca_plot(
    pca_df: pd.DataFrame,
    pc1_col: str = "PC1",
    pc2_col: str = "PC2",
    var1: float = 0.0,
    var2: float = 0.0,
    sample_col: str = "sample",
    condition_col: str = "condition",
) -> go.Figure:
    """Interactive PCA scatter plot colored by condition.

    Parameters
    ----------
    pca_df : pd.DataFrame
        DataFrame with PC coordinates and sample metadata.
    pc1_col : str
        Column name for PC1 values.
    pc2_col : str
        Column name for PC2 values.
    var1 : float
        Variance explained by PC1 (0-100 scale).
    var2 : float
        Variance explained by PC2 (0-100 scale).
    sample_col : str
        Column name for sample identifiers.
    condition_col : str
        Column name for condition labels.

    Returns
    -------
    go.Figure
        Interactive PCA scatter plot.
    """
    setup_plotly_theme()

    if pca_df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No PCA data available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="PCA Plot — No Data")
        return fig

    conditions = pca_df[condition_col].unique() if condition_col in pca_df.columns else ["Unknown"]
    cmap = condition_color_map(conditions)

    fig = go.Figure()

    for cond in conditions:
        if condition_col in pca_df.columns:
            subset = pca_df[pca_df[condition_col] == cond]
        else:
            subset = pca_df
        color = cmap.get(cond, PALETTE["neutral"])

        sample_text = subset[sample_col] if sample_col in subset.columns else None

        fig.add_trace(
            go.Scatter(
                x=subset[pc1_col],
                y=subset[pc2_col],
                mode="markers+text",
                marker=dict(color=color, size=10, line=dict(color="white", width=1)),
                text=sample_text,
                textposition="top center",
                textfont=dict(size=9),
                name=str(cond),
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    f"PC1: %{{x:.2f}}<br>"
                    f"PC2: %{{y:.2f}}<br>"
                    "<extra></extra>"
                ),
            )
        )

    x_title = f"PC1 ({var1:.1f}% variance)" if var1 > 0 else "PC1"
    y_title = f"PC2 ({var2:.1f}% variance)" if var2 > 0 else "PC2"

    fig.update_layout(
        title="Principal Component Analysis",
        xaxis_title=x_title,
        yaxis_title=y_title,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def correlation_heatmap(
    corr_matrix: pd.DataFrame,
    sample_labels: Optional[Sequence[str]] = None,
) -> go.Figure:
    """Annotated correlation heatmap for sample-sample comparisons.

    Parameters
    ----------
    corr_matrix : pd.DataFrame
        Square correlation matrix (samples x samples).
    sample_labels : Sequence[str], optional
        Labels for samples. Uses matrix index/columns if not provided.

    Returns
    -------
    go.Figure
        Annotated heatmap.
    """
    setup_plotly_theme()

    if corr_matrix.empty:
        fig = go.Figure()
        fig.add_annotation(text="No correlation data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Correlation Heatmap — No Data")
        return fig

    labels = sample_labels if sample_labels is not None else corr_matrix.columns.tolist()

    # Round for display
    z_text = corr_matrix.round(3).values.tolist()

    fig = go.Figure(
        go.Heatmap(
            z=corr_matrix.values,
            x=labels,
            y=labels,
            colorscale="RdBu_r",
            zmin=corr_matrix.values.min(),
            zmax=1.0,
            text=z_text,
            texttemplate="%{text:.3f}",
            textfont=dict(size=8),
            colorbar=dict(title="Pearson r"),
            hovertemplate="%{x} vs %{y}<br>r = %{z:.4f}<extra></extra>",
        )
    )

    n = len(labels)
    fig.update_layout(
        title="Sample-Sample Correlation",
        xaxis=dict(tickangle=45, tickfont=dict(size=9)),
        yaxis=dict(tickfont=dict(size=9), autorange="reversed"),
        height=max(400, n * 35 + 150),
        width=max(500, n * 35 + 200),
    )

    return fig


def top_deg_heatmap(
    expr_matrix: pd.DataFrame,
    gene_list: Sequence[str],
    sample_labels: Optional[Sequence[str]] = None,
    condition_labels: Optional[Sequence[str]] = None,
) -> go.Figure:
    """Z-scored expression heatmap for selected DEGs.

    Parameters
    ----------
    expr_matrix : pd.DataFrame
        Normalized expression matrix (genes x samples).
    gene_list : Sequence[str]
        Genes to include in the heatmap.
    sample_labels : Sequence[str], optional
        Sample display labels.
    condition_labels : Sequence[str], optional
        Condition annotations per sample (for color bar).

    Returns
    -------
    go.Figure
        Z-scored heatmap.
    """
    setup_plotly_theme()

    if expr_matrix.empty or not gene_list:
        fig = go.Figure()
        fig.add_annotation(text="No expression data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Top DEG Heatmap — No Data")
        return fig

    # Filter to available genes
    available = [g for g in gene_list if g in expr_matrix.index]
    if not available:
        fig = go.Figure()
        fig.add_annotation(text="No matching genes found", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Top DEG Heatmap — No Matching Genes")
        return fig

    sub = expr_matrix.loc[available].copy()

    # Z-score per gene (row)
    row_means = sub.mean(axis=1)
    row_stds = sub.std(axis=1).replace(0, 1)
    z_scored = sub.sub(row_means, axis=0).div(row_stds, axis=0)

    x_labels = sample_labels if sample_labels is not None else z_scored.columns.tolist()

    fig = go.Figure(
        go.Heatmap(
            z=z_scored.values,
            x=x_labels,
            y=z_scored.index.tolist(),
            colorscale=[
                [0, "#0072B2"],
                [0.5, "#FFFFFF"],
                [1, "#E69F00"],
            ],
            zmid=0,
            colorbar=dict(title="Z-score"),
            hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Z-score: %{z:.2f}<extra></extra>",
        )
    )

    n_genes = len(available)
    fig.update_layout(
        title=f"Expression Heatmap ({n_genes} genes)",
        xaxis=dict(tickangle=45, tickfont=dict(size=9)),
        yaxis=dict(tickfont=dict(size=9)),
        height=max(400, n_genes * 16 + 150),
    )

    return fig
