"""Interactive cross-condition comparison visualizations using Plotly.

Provides direction concordance heatmaps, log2FC scatter plots, and
overlap bar charts for comparing DEG results across conditions.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from rnaseq_explorer.viz.theme import (
    PALETTE,
    CONDITION_COLORS,
    setup_plotly_theme,
)


def direction_concordance_heatmap(
    concordance_matrix: pd.DataFrame,
    condition_labels: Optional[Sequence[str]] = None,
) -> go.Figure:
    """Pearson correlation heatmap of log2FC across conditions.

    Parameters
    ----------
    concordance_matrix : pd.DataFrame
        Square matrix of pairwise correlations between conditions.
    condition_labels : Sequence[str], optional
        Display labels for conditions.

    Returns
    -------
    go.Figure
        Annotated concordance heatmap.
    """
    setup_plotly_theme()

    if concordance_matrix.empty:
        fig = go.Figure()
        fig.add_annotation(text="No concordance data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Direction Concordance — No Data")
        return fig

    labels = (
        list(condition_labels)
        if condition_labels is not None
        else concordance_matrix.columns.tolist()
    )

    z_text = concordance_matrix.round(3).values.tolist()

    fig = go.Figure(
        go.Heatmap(
            z=concordance_matrix.values,
            x=labels,
            y=labels,
            colorscale="RdBu_r",
            zmin=-1,
            zmax=1,
            text=z_text,
            texttemplate="%{text:.3f}",
            textfont=dict(size=10),
            colorbar=dict(title="Pearson r"),
            hovertemplate="%{x} vs %{y}<br>r = %{z:.4f}<extra></extra>",
        )
    )

    n = len(labels)
    fig.update_layout(
        title="log₂FC Direction Concordance",
        xaxis=dict(tickangle=45, tickfont=dict(size=10)),
        yaxis=dict(tickfont=dict(size=10), autorange="reversed"),
        height=max(400, n * 60 + 150),
        width=max(500, n * 60 + 200),
    )

    return fig


def log2fc_scatter(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    gene_col: str = "gene_name",
    log2fc_col: str = "log2FoldChange",
    cond1_name: str = "Condition 1",
    cond2_name: str = "Condition 2",
    padj_col: Optional[str] = "padj",
    padj_cutoff: float = 0.05,
) -> go.Figure:
    """Pairwise scatter of log2FC between two conditions with fit line.

    Parameters
    ----------
    df1 : pd.DataFrame
        DESeq2 results for condition 1.
    df2 : pd.DataFrame
        DESeq2 results for condition 2.
    gene_col : str
        Column name for gene identifier (used for merging).
    log2fc_col : str
        Column name for log2 fold change.
    cond1_name : str
        Display name for condition 1.
    cond2_name : str
        Display name for condition 2.
    padj_col : str, optional
        Column for adjusted p-value (for coloring).
    padj_cutoff : float
        Threshold for significance coloring.

    Returns
    -------
    go.Figure
        Scatter plot with fit line.
    """
    setup_plotly_theme()

    if df1.empty or df2.empty:
        fig = go.Figure()
        fig.add_annotation(text="Insufficient data for comparison", showarrow=False, font=dict(size=16))
        fig.update_layout(title="log₂FC Scatter — No Data")
        return fig

    if gene_col not in df1.columns or gene_col not in df2.columns:
        fig = go.Figure()
        fig.add_annotation(text=f"Column '{gene_col}' not found", showarrow=False, font=dict(size=16))
        fig.update_layout(title="log₂FC Scatter — Missing Column")
        return fig

    merged = pd.merge(
        df1[[gene_col, log2fc_col]].rename(columns={log2fc_col: "x"}),
        df2[[gene_col, log2fc_col]].rename(columns={log2fc_col: "y"}),
        on=gene_col,
        how="inner",
    ).dropna()

    if merged.empty:
        fig = go.Figure()
        fig.add_annotation(text="No overlapping genes", showarrow=False, font=dict(size=16))
        fig.update_layout(title="log₂FC Scatter — No Overlap")
        return fig

    # Classify significance in both conditions
    sig_both = False
    if padj_col and padj_col in df1.columns and padj_col in df2.columns:
        merged_sig = pd.merge(
            df1[[gene_col, padj_col]].rename(columns={padj_col: "padj1"}),
            df2[[gene_col, padj_col]].rename(columns={padj_col: "padj2"}),
            on=gene_col,
            how="inner",
        )
        merged = pd.merge(merged, merged_sig, on=gene_col, how="left")
        merged["sig_category"] = "NS in both"
        merged.loc[
            (merged["padj1"] < padj_cutoff) & (merged["padj2"] >= padj_cutoff), "sig_category"
        ] = f"Sig in {cond1_name} only"
        merged.loc[
            (merged["padj1"] >= padj_cutoff) & (merged["padj2"] < padj_cutoff), "sig_category"
        ] = f"Sig in {cond2_name} only"
        merged.loc[
            (merged["padj1"] < padj_cutoff) & (merged["padj2"] < padj_cutoff), "sig_category"
        ] = "Sig in both"
        sig_both = True

    fig = go.Figure()

    if sig_both:
        color_map = {
            "NS in both": PALETTE["neutral"],
            f"Sig in {cond1_name} only": CONDITION_COLORS[0],
            f"Sig in {cond2_name} only": CONDITION_COLORS[1],
            "Sig in both": PALETTE["highlight"],
        }
        for cat in ["NS in both", f"Sig in {cond1_name} only", f"Sig in {cond2_name} only", "Sig in both"]:
            subset = merged[merged["sig_category"] == cat]
            if subset.empty:
                continue
            fig.add_trace(
                go.Scattergl(
                    x=subset["x"],
                    y=subset["y"],
                    mode="markers",
                    marker=dict(color=color_map.get(cat, PALETTE["neutral"]), size=4, opacity=0.6),
                    name=f"{cat} ({len(subset):,})",
                    text=subset[gene_col],
                    hovertemplate=(
                        "<b>%{text}</b><br>"
                        f"{cond1_name} log2FC: %{{x:.3f}}<br>"
                        f"{cond2_name} log2FC: %{{y:.3f}}<br>"
                        "<extra></extra>"
                    ),
                )
            )
    else:
        fig.add_trace(
            go.Scattergl(
                x=merged["x"],
                y=merged["y"],
                mode="markers",
                marker=dict(color=PALETTE["accent1"], size=4, opacity=0.6),
                name=f"Genes ({len(merged):,})",
                text=merged[gene_col],
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    f"{cond1_name}: %{{x:.3f}}<br>"
                    f"{cond2_name}: %{{y:.3f}}<br>"
                    "<extra></extra>"
                ),
            )
        )

    # Fit line
    valid = merged[["x", "y"]].dropna()
    if len(valid) > 2:
        corr = valid["x"].corr(valid["y"])
        m, b = np.polyfit(valid["x"], valid["y"], 1)
        x_range = np.linspace(valid["x"].min(), valid["x"].max(), 100)
        fig.add_trace(
            go.Scatter(
                x=x_range,
                y=m * x_range + b,
                mode="lines",
                line=dict(color="#333333", width=1.5, dash="dash"),
                name=f"Fit (r={corr:.3f})",
                showlegend=True,
            )
        )

    # Identity line
    lim = max(abs(merged["x"].max()), abs(merged["y"].max()), abs(merged["x"].min()), abs(merged["y"].min()))
    fig.add_trace(
        go.Scatter(
            x=[-lim, lim], y=[-lim, lim],
            mode="lines",
            line=dict(color="#CCCCCC", width=1, dash="dot"),
            name="Identity",
            showlegend=True,
        )
    )

    fig.update_layout(
        title=f"log₂FC: {cond1_name} vs {cond2_name}",
        xaxis_title=f"log₂FC — {cond1_name}",
        yaxis_title=f"log₂FC — {cond2_name}",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def overlap_bar(
    overlap_data: pd.DataFrame,
    condition_labels: Optional[Sequence[str]] = None,
) -> go.Figure:
    """Stacked bar chart showing unique and shared genes per condition.

    Parameters
    ----------
    overlap_data : pd.DataFrame
        DataFrame with columns like 'condition', 'unique', 'shared'.
        Alternatively, accepts a matrix where rows=conditions and
        columns represent overlap categories.
    condition_labels : Sequence[str], optional
        Display labels for conditions.

    Returns
    -------
    go.Figure
        Stacked bar chart.
    """
    setup_plotly_theme()

    if overlap_data.empty:
        fig = go.Figure()
        fig.add_annotation(text="No overlap data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Gene Overlap — No Data")
        return fig

    # Detect format
    if "condition" in overlap_data.columns:
        # Long format with condition, unique, shared columns
        conditions = overlap_data["condition"].tolist()
        if condition_labels is not None:
            conditions = list(condition_labels)

        fig = go.Figure()

        if "unique" in overlap_data.columns:
            fig.add_trace(
                go.Bar(
                    x=conditions,
                    y=overlap_data["unique"],
                    name="Unique",
                    marker_color=PALETTE["accent1"],
                    hovertemplate="%{x}<br>Unique: %{y:,}<extra></extra>",
                )
            )
        if "shared" in overlap_data.columns:
            fig.add_trace(
                go.Bar(
                    x=conditions,
                    y=overlap_data["shared"],
                    name="Shared",
                    marker_color=PALETTE["highlight"],
                    hovertemplate="%{x}<br>Shared: %{y:,}<extra></extra>",
                )
            )
    else:
        # Matrix format: rows=conditions, columns=categories
        conditions = (
            list(condition_labels)
            if condition_labels is not None
            else overlap_data.index.tolist()
        )

        fig = go.Figure()
        for i, col in enumerate(overlap_data.columns):
            color = CONDITION_COLORS[i % len(CONDITION_COLORS)]
            fig.add_trace(
                go.Bar(
                    x=conditions,
                    y=overlap_data[col],
                    name=str(col),
                    marker_color=color,
                    hovertemplate=f"%{{x}}<br>{col}: %{{y:,}}<extra></extra>",
                )
            )

    fig.update_layout(
        title="Gene Overlap by Condition",
        xaxis_title="Condition",
        yaxis_title="Gene Count",
        barmode="stack",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig
