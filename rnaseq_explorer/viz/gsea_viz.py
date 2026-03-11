"""Interactive enrichment visualizations using Plotly.

Provides NES bar charts, enrichment dot plots, ORA dot plots, leading edge
extraction, and enrichment comparison plots for GSEA and ORA results.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from rnaseq_explorer.viz.theme import (
    PALETTE,
    CATEGORY_COLORS,
    FONT_SIZE_TITLE,
    FONT_SIZE_ANNOTATION,
    setup_plotly_theme,
)


def _truncate_name(name: str, max_len: int = 55) -> str:
    """Truncate a pathway name for display.

    Parameters
    ----------
    name : str
        Full pathway name.
    max_len : int
        Maximum display length.

    Returns
    -------
    str
        Truncated name with ellipsis if needed.
    """
    if len(name) <= max_len:
        return name
    return name[: max_len - 3] + "..."


def nes_bar_chart(
    gsea_results: pd.DataFrame,
    n: int = 20,
    fdr_cutoff: float = 0.25,
    nes_col: str = "NES",
    fdr_col: str = "FDR q-val",
    term_col: str = "Term",
) -> go.Figure:
    """Horizontal bar chart of top pathways ranked by NES.

    Parameters
    ----------
    gsea_results : pd.DataFrame
        GSEA results with NES and FDR columns.
    n : int
        Number of top pathways per direction.
    fdr_cutoff : float
        FDR threshold for filtering.
    nes_col : str
        Column name for normalized enrichment score.
    fdr_col : str
        Column name for FDR.
    term_col : str
        Column name for pathway/term name.

    Returns
    -------
    go.Figure
        Horizontal bar chart of NES values.
    """
    setup_plotly_theme()

    if gsea_results.empty:
        fig = go.Figure()
        fig.add_annotation(text="No GSEA results available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GSEA NES Bar Chart — No Data")
        return fig

    # Normalize column names (gseapy uses different conventions)
    df = gsea_results.copy()
    for candidate in ["NES", "nes"]:
        if candidate in df.columns:
            df = df.rename(columns={candidate: nes_col})
            break
    for candidate in ["FDR q-val", "fdr", "FDR", "NOM p-val"]:
        if candidate in df.columns and fdr_col not in df.columns:
            df = df.rename(columns={candidate: fdr_col})
            break
    for candidate in ["Term", "term", "Name", "name", "pathway"]:
        if candidate in df.columns and term_col not in df.columns:
            df = df.rename(columns={candidate: term_col})
            break

    if nes_col not in df.columns or fdr_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="Required columns not found", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GSEA NES — Missing Columns")
        return fig

    sig = df[df[fdr_col] <= fdr_cutoff].copy()
    if sig.empty:
        sig = df.copy()

    pos = sig[sig[nes_col] > 0].nlargest(n, nes_col)
    neg = sig[sig[nes_col] < 0].nsmallest(n, nes_col)
    combined = pd.concat([neg, pos]).sort_values(nes_col)

    if term_col in combined.columns:
        labels = [_truncate_name(t) for t in combined[term_col]]
    else:
        labels = combined.index.astype(str).tolist()

    colors = [PALETTE["up"] if v > 0 else PALETTE["down"] for v in combined[nes_col]]

    fig = go.Figure(
        go.Bar(
            y=labels,
            x=combined[nes_col].values,
            orientation="h",
            marker_color=colors,
            hovertemplate=(
                "<b>%{y}</b><br>"
                "NES: %{x:.3f}<br>"
                "<extra></extra>"
            ),
        )
    )

    fig.add_vline(x=0, line_color="#333333", line_width=1)

    fig.update_layout(
        title=f"GSEA Top Pathways (FDR ≤ {fdr_cutoff})",
        xaxis_title="Normalized Enrichment Score (NES)",
        yaxis_title="",
        height=max(400, len(combined) * 24),
    )

    return fig


def enrichment_dot_plot(
    gsea_results: pd.DataFrame,
    n: int = 30,
    fdr_cutoff: float = 0.25,
    nes_col: str = "NES",
    fdr_col: str = "FDR q-val",
    term_col: str = "Term",
    gene_count_col: str = "Gene %",
    size_col: Optional[str] = None,
) -> go.Figure:
    """Dot plot: x=NES, y=pathway, size=gene overlap, color=-log10(FDR).

    Parameters
    ----------
    gsea_results : pd.DataFrame
        GSEA results.
    n : int
        Number of top pathways to show.
    fdr_cutoff : float
        FDR threshold for filtering.
    nes_col : str
        Column name for NES.
    fdr_col : str
        Column name for FDR.
    term_col : str
        Column name for pathway name.
    gene_count_col : str
        Column name for gene count or percentage.
    size_col : str, optional
        Explicit column for dot size. Falls back to gene_count_col.

    Returns
    -------
    go.Figure
        Enrichment dot plot.
    """
    setup_plotly_theme()

    if gsea_results.empty:
        fig = go.Figure()
        fig.add_annotation(text="No GSEA results available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Enrichment Dot Plot — No Data")
        return fig

    df = gsea_results.copy()
    # Auto-detect columns
    col_map = {
        nes_col: ["NES", "nes"],
        fdr_col: ["FDR q-val", "fdr", "FDR", "NOM p-val"],
        term_col: ["Term", "term", "Name", "name", "pathway"],
    }
    for target, candidates in col_map.items():
        if target not in df.columns:
            for c in candidates:
                if c in df.columns:
                    df = df.rename(columns={c: target})
                    break

    if nes_col not in df.columns or fdr_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="Required columns not found", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Enrichment Dot Plot — Missing Columns")
        return fig

    sig = df[df[fdr_col] <= fdr_cutoff].copy()
    if sig.empty:
        sig = df.copy()

    top = sig.reindex(sig[nes_col].abs().nlargest(n).index)
    top["neg_log10_fdr"] = -np.log10(top[fdr_col].clip(lower=1e-300))

    # Determine size values
    size_source = size_col or gene_count_col
    if size_source in top.columns:
        sizes = top[size_source].fillna(1)
        # Normalize to reasonable marker sizes
        sizes = 8 + 20 * (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-9)
    else:
        sizes = 12

    if term_col in top.columns:
        labels = [_truncate_name(t) for t in top[term_col]]
    else:
        labels = top.index.astype(str).tolist()

    fig = go.Figure(
        go.Scatter(
            x=top[nes_col],
            y=labels,
            mode="markers",
            marker=dict(
                size=sizes,
                color=top["neg_log10_fdr"],
                colorscale="Viridis",
                colorbar=dict(title="-log₁₀ FDR"),
                line=dict(color="#333333", width=0.5),
            ),
            hovertemplate=(
                "<b>%{y}</b><br>"
                "NES: %{x:.3f}<br>"
                "-log10(FDR): %{marker.color:.2f}<br>"
                "<extra></extra>"
            ),
        )
    )

    fig.add_vline(x=0, line_color="#999999", line_width=0.8)

    fig.update_layout(
        title=f"Enrichment Dot Plot (Top {n})",
        xaxis_title="Normalized Enrichment Score (NES)",
        yaxis_title="",
        height=max(400, len(top) * 22),
    )

    return fig


def leading_edge_table(
    gsea_results: pd.DataFrame,
    pathway_name: str,
    term_col: str = "Term",
    lead_edge_col: str = "Lead_genes",
) -> pd.DataFrame:
    """Extract leading edge genes for a specific pathway.

    Parameters
    ----------
    gsea_results : pd.DataFrame
        GSEA results with leading edge gene information.
    pathway_name : str
        Name of the pathway to query.
    term_col : str
        Column name for pathway/term name.
    lead_edge_col : str
        Column name for leading edge gene list.

    Returns
    -------
    pd.DataFrame
        DataFrame with gene names from the leading edge.
    """
    if gsea_results.empty:
        return pd.DataFrame(columns=["gene"])

    # Auto-detect term column
    for candidate in [term_col, "Term", "term", "Name", "name", "pathway"]:
        if candidate in gsea_results.columns:
            term_col = candidate
            break

    # Auto-detect leading edge column
    for candidate in [lead_edge_col, "Lead_genes", "lead_genes", "genes", "Gene"]:
        if candidate in gsea_results.columns:
            lead_edge_col = candidate
            break

    if term_col not in gsea_results.columns or lead_edge_col not in gsea_results.columns:
        return pd.DataFrame(columns=["gene"])

    row = gsea_results[gsea_results[term_col] == pathway_name]
    if row.empty:
        return pd.DataFrame(columns=["gene"])

    genes_str = row.iloc[0][lead_edge_col]
    if pd.isna(genes_str) or not isinstance(genes_str, str):
        return pd.DataFrame(columns=["gene"])

    genes = [g.strip() for g in genes_str.split(";") if g.strip()]
    if not genes:
        # Try comma separation
        genes = [g.strip() for g in str(genes_str).split(",") if g.strip()]

    return pd.DataFrame({"gene": genes})


def ora_dot_plot(
    ora_results: pd.DataFrame,
    n: int = 20,
    fdr_cutoff: float = 0.05,
    overlap_col: str = "Overlap_ratio",
    padj_col: str = "Adjusted P-value",
    term_col: str = "Term",
    score_col: str = "Combined Score",
) -> go.Figure:
    """Dot plot for ORA results: x=overlap, y=term, size=score, color=-log10(padj).

    Parameters
    ----------
    ora_results : pd.DataFrame
        Over-representation analysis results.
    n : int
        Number of top terms to show.
    fdr_cutoff : float
        Adjusted p-value threshold.
    overlap_col : str
        Column name for overlap ratio.
    padj_col : str
        Column name for adjusted p-value.
    term_col : str
        Column name for term.
    score_col : str
        Column name for combined/enrichment score.

    Returns
    -------
    go.Figure
        ORA dot plot.
    """
    setup_plotly_theme()

    if ora_results.empty:
        fig = go.Figure()
        fig.add_annotation(text="No ORA results available", showarrow=False, font=dict(size=16))
        fig.update_layout(title="ORA Dot Plot — No Data")
        return fig

    df = ora_results.copy()

    # Auto-detect columns
    for target, candidates in {
        padj_col: ["Adjusted P-value", "adjusted_pvalue", "padj", "FDR", "q-value"],
        term_col: ["Term", "term", "Name", "name", "native"],
        score_col: ["Combined Score", "combined_score", "Odds Ratio", "odds_ratio"],
    }.items():
        if target not in df.columns:
            for c in candidates:
                if c in df.columns:
                    df = df.rename(columns={c: target})
                    break

    if padj_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="Required columns not found", showarrow=False, font=dict(size=16))
        fig.update_layout(title="ORA Dot Plot — Missing Columns")
        return fig

    sig = df[df[padj_col] <= fdr_cutoff].copy()
    if sig.empty:
        sig = df.head(n).copy()

    top = sig.head(n).copy()
    top["neg_log10_padj"] = -np.log10(top[padj_col].clip(lower=1e-300))

    # x-axis: overlap ratio or combined score
    if overlap_col in top.columns:
        x_vals = top[overlap_col]
        x_title = "Overlap Ratio"
    elif score_col in top.columns:
        x_vals = top[score_col]
        x_title = "Combined Score"
    else:
        x_vals = top["neg_log10_padj"]
        x_title = "-log₁₀ adjusted p-value"

    # Size: combined score
    if score_col in top.columns:
        sizes = top[score_col].fillna(1).abs()
        sizes = 8 + 20 * (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-9)
    else:
        sizes = 12

    if term_col in top.columns:
        labels = [_truncate_name(t) for t in top[term_col]]
    else:
        labels = top.index.astype(str).tolist()

    fig = go.Figure(
        go.Scatter(
            x=x_vals,
            y=labels,
            mode="markers",
            marker=dict(
                size=sizes,
                color=top["neg_log10_padj"],
                colorscale="Viridis",
                colorbar=dict(title="-log₁₀ padj"),
                line=dict(color="#333333", width=0.5),
            ),
            hovertemplate=(
                "<b>%{y}</b><br>"
                f"{x_title}: %{{x:.3f}}<br>"
                "-log10(padj): %{marker.color:.2f}<br>"
                "<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title=f"ORA Enrichment (Top {n})",
        xaxis_title=x_title,
        yaxis_title="",
        height=max(400, len(top) * 22),
    )

    return fig


def enrichment_comparison(
    up_results: pd.DataFrame,
    down_results: pd.DataFrame,
    n: int = 15,
    nes_col: str = "NES",
    fdr_col: str = "FDR q-val",
    term_col: str = "Term",
    fdr_cutoff: float = 0.25,
) -> go.Figure:
    """Side-by-side NES comparison for up- and down-regulated gene enrichment.

    Parameters
    ----------
    up_results : pd.DataFrame
        GSEA results for upregulated genes.
    down_results : pd.DataFrame
        GSEA results for downregulated genes.
    n : int
        Number of top pathways per set.
    nes_col : str
        Column name for NES.
    fdr_col : str
        Column name for FDR.
    term_col : str
        Column name for pathway name.
    fdr_cutoff : float
        FDR threshold.

    Returns
    -------
    go.Figure
        Side-by-side bar chart.
    """
    setup_plotly_theme()

    def _prep(df: pd.DataFrame, label: str) -> pd.DataFrame:
        if df.empty:
            return pd.DataFrame()
        out = df.copy()
        # Auto-detect columns
        for target, candidates in {
            nes_col: ["NES", "nes"],
            fdr_col: ["FDR q-val", "fdr", "FDR"],
            term_col: ["Term", "term", "Name", "name"],
        }.items():
            if target not in out.columns:
                for c in candidates:
                    if c in out.columns:
                        out = out.rename(columns={c: target})
                        break
        if nes_col in out.columns and fdr_col in out.columns:
            sig = out[out[fdr_col] <= fdr_cutoff]
            if sig.empty:
                sig = out
            top = sig.reindex(sig[nes_col].abs().nlargest(n).index)
            top["source"] = label
            return top
        return pd.DataFrame()

    up_prep = _prep(up_results, "Upregulated")
    down_prep = _prep(down_results, "Downregulated")

    if up_prep.empty and down_prep.empty:
        fig = go.Figure()
        fig.add_annotation(text="No enrichment results", showarrow=False, font=dict(size=16))
        fig.update_layout(title="Enrichment Comparison — No Data")
        return fig

    fig = go.Figure()

    for source_df, name, color in [
        (up_prep, "Upregulated", PALETTE["up"]),
        (down_prep, "Downregulated", PALETTE["down"]),
    ]:
        if source_df.empty:
            continue
        labels = (
            [_truncate_name(t) for t in source_df[term_col]]
            if term_col in source_df.columns
            else source_df.index.astype(str).tolist()
        )
        fig.add_trace(
            go.Bar(
                y=labels,
                x=source_df[nes_col].values,
                orientation="h",
                name=name,
                marker_color=color,
                opacity=0.8,
                hovertemplate="<b>%{y}</b><br>NES: %{x:.3f}<extra></extra>",
            )
        )

    fig.add_vline(x=0, line_color="#333333", line_width=1)

    fig.update_layout(
        title="Enrichment Comparison: Up vs Down Genes",
        xaxis_title="NES",
        yaxis_title="",
        barmode="group",
        height=max(450, n * 28),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig
