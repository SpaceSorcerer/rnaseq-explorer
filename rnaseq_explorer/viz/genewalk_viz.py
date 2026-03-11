"""GeneWalk result exploration visualizations using Plotly.

Provides volcano plots, per-gene GO bar charts, bipartite network graphs,
gene-by-GO heatmaps, GO domain pie charts, and gene summary rankings.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from rnaseq_explorer.viz.theme import (
    PALETTE,
    CATEGORY_COLORS,
    CONDITION_COLORS,
    FONT_SIZE_TITLE,
    FONT_SIZE_ANNOTATION,
    setup_plotly_theme,
)


def gw_volcano(
    df: pd.DataFrame,
    sim_col: str = "sim",
    padj_col: str = "gene_padj",
    gene_col: str = "hgnc_symbol",
    padj_cutoff: float = 0.1,
) -> go.Figure:
    """GeneWalk volcano plot: similarity vs -log10(padj).

    Parameters
    ----------
    df : pd.DataFrame
        GeneWalk results.
    sim_col : str
        Column name for similarity score.
    padj_col : str
        Column name for adjusted p-value.
    gene_col : str
        Column name for gene symbol.
    padj_cutoff : float
        Adjusted p-value threshold.

    Returns
    -------
    go.Figure
        Interactive scatter plot.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No GeneWalk data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GeneWalk Volcano — No Data")
        return fig

    plot_df = df.dropna(subset=[sim_col, padj_col]).copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df[padj_col].clip(lower=1e-300))
    plot_df["significant"] = plot_df[padj_col] < padj_cutoff

    fig = go.Figure()

    for is_sig, label, color in [
        (False, "NS", PALETTE["neutral"]),
        (True, f"padj < {padj_cutoff}", PALETTE["accent2"]),
    ]:
        subset = plot_df[plot_df["significant"] == is_sig]
        if subset.empty:
            continue
        gene_text = subset[gene_col] if gene_col in subset.columns else None
        fig.add_trace(
            go.Scattergl(
                x=subset[sim_col],
                y=subset["neg_log10_padj"],
                mode="markers",
                marker=dict(color=color, size=5, opacity=0.7),
                name=f"{label} ({len(subset):,})",
                text=gene_text,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "Similarity: %{x:.3f}<br>"
                    "padj: %{customdata:.2e}<br>"
                    "<extra></extra>"
                ),
                customdata=subset[padj_col],
            )
        )

    fig.add_hline(
        y=-np.log10(padj_cutoff), line_dash="dash", line_color="#666666", line_width=1,
        annotation_text=f"padj={padj_cutoff}", annotation_position="top right",
    )

    fig.update_layout(
        title="GeneWalk: Similarity vs Significance",
        xaxis_title="Similarity Score",
        yaxis_title="-log₁₀ adjusted p-value",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def gw_gene_bar(
    df: pd.DataFrame,
    gene: str,
    padj_col: str = "gene_padj",
    go_name_col: str = "go_name",
    sim_col: str = "sim",
    n: int = 15,
) -> go.Figure:
    """Bar chart of top GO terms for a single gene.

    Parameters
    ----------
    df : pd.DataFrame
        GeneWalk results.
    gene : str
        Gene symbol to query.
    padj_col : str
        Column name for adjusted p-value.
    go_name_col : str
        Column name for GO term name.
    sim_col : str
        Column name for similarity score.
    n : int
        Number of top GO terms to show.

    Returns
    -------
    go.Figure
        Horizontal bar chart of GO terms.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No GeneWalk data", showarrow=False, font=dict(size=16))
        fig.update_layout(title=f"GeneWalk: {gene} — No Data")
        return fig

    # Find gene in any likely column
    gene_col_candidates = ["hgnc_symbol", "gene", "gene_name", "GeneID"]
    gene_data = pd.DataFrame()
    for col in gene_col_candidates:
        if col in df.columns:
            gene_data = df[df[col].astype(str).str.upper() == gene.upper()]
            if not gene_data.empty:
                break

    if gene_data.empty:
        fig = go.Figure()
        fig.add_annotation(text=f"Gene '{gene}' not found", showarrow=False, font=dict(size=16))
        fig.update_layout(title=f"GeneWalk: {gene} — Not Found")
        return fig

    if go_name_col not in gene_data.columns or sim_col not in gene_data.columns:
        fig = go.Figure()
        fig.add_annotation(text="Required columns missing", showarrow=False, font=dict(size=16))
        fig.update_layout(title=f"GeneWalk: {gene} — Missing Columns")
        return fig

    top = gene_data.nlargest(n, sim_col)
    labels = top[go_name_col].astype(str).tolist()
    # Truncate long names
    labels = [l[:50] + "..." if len(l) > 50 else l for l in labels]

    colors = []
    if padj_col in top.columns:
        for p in top[padj_col]:
            colors.append(PALETTE["accent2"] if p < 0.1 else PALETTE["neutral"])
    else:
        colors = [PALETTE["accent2"]] * len(top)

    fig = go.Figure(
        go.Bar(
            y=labels,
            x=top[sim_col].values,
            orientation="h",
            marker_color=colors,
            hovertemplate="<b>%{y}</b><br>Similarity: %{x:.3f}<extra></extra>",
        )
    )

    fig.update_layout(
        title=f"GeneWalk: Top GO Terms for {gene}",
        xaxis_title="Similarity Score",
        yaxis_title="",
        height=max(350, n * 24),
    )

    return fig


def gw_network(
    df: pd.DataFrame,
    padj_cutoff: float = 0.1,
    min_sim: float = 0.1,
    padj_col: str = "gene_padj",
    gene_col: str = "hgnc_symbol",
    go_name_col: str = "go_name",
    sim_col: str = "sim",
    max_nodes: int = 100,
) -> go.Figure:
    """Bipartite network: genes connected to GO terms by similarity.

    Uses a spring layout to position gene and GO term nodes.

    Parameters
    ----------
    df : pd.DataFrame
        GeneWalk results.
    padj_cutoff : float
        Filter to significant associations.
    min_sim : float
        Minimum similarity to include an edge.
    padj_col : str
        Column for adjusted p-value.
    gene_col : str
        Column for gene symbol.
    go_name_col : str
        Column for GO term name.
    sim_col : str
        Column for similarity score.
    max_nodes : int
        Maximum number of nodes to display (for performance).

    Returns
    -------
    go.Figure
        Interactive network graph.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No GeneWalk data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GeneWalk Network — No Data")
        return fig

    required = [gene_col, go_name_col, sim_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        fig = go.Figure()
        fig.add_annotation(
            text=f"Missing columns: {', '.join(missing)}", showarrow=False, font=dict(size=16)
        )
        fig.update_layout(title="GeneWalk Network — Missing Columns")
        return fig

    # Filter
    filt = df.copy()
    if padj_col in filt.columns:
        filt = filt[filt[padj_col] < padj_cutoff]
    filt = filt[filt[sim_col] >= min_sim]

    if filt.empty:
        fig = go.Figure()
        fig.add_annotation(text="No significant associations", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GeneWalk Network — No Significant Results")
        return fig

    # Limit for performance
    filt = filt.nlargest(max_nodes, sim_col)

    # Build graph positions using simple spring layout
    genes = filt[gene_col].unique().tolist()
    go_terms = filt[go_name_col].unique().tolist()

    # Simple circular layout: genes on left, GO terms on right
    node_positions = {}
    for i, g in enumerate(genes):
        y = i / max(len(genes) - 1, 1) if len(genes) > 1 else 0.5
        node_positions[g] = (0.0, y)
    for i, t in enumerate(go_terms):
        y = i / max(len(go_terms) - 1, 1) if len(go_terms) > 1 else 0.5
        node_positions[t] = (1.0, y)

    # Edge traces
    edge_x, edge_y = [], []
    for _, row in filt.iterrows():
        g_pos = node_positions.get(row[gene_col])
        t_pos = node_positions.get(row[go_name_col])
        if g_pos and t_pos:
            edge_x.extend([g_pos[0], t_pos[0], None])
            edge_y.extend([g_pos[1], t_pos[1], None])

    fig = go.Figure()

    # Edges
    fig.add_trace(
        go.Scatter(
            x=edge_x, y=edge_y,
            mode="lines",
            line=dict(width=0.5, color="#999999"),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    # Gene nodes
    gene_x = [node_positions[g][0] for g in genes]
    gene_y = [node_positions[g][1] for g in genes]
    fig.add_trace(
        go.Scatter(
            x=gene_x, y=gene_y,
            mode="markers+text",
            marker=dict(color=PALETTE["up"], size=12, line=dict(color="white", width=1)),
            text=genes,
            textposition="middle left",
            textfont=dict(size=9),
            name="Genes",
            hovertemplate="<b>%{text}</b><extra>Gene</extra>",
        )
    )

    # GO term nodes
    go_x = [node_positions[t][0] for t in go_terms]
    go_y = [node_positions[t][1] for t in go_terms]
    go_labels = [t[:35] + "..." if len(t) > 35 else t for t in go_terms]
    fig.add_trace(
        go.Scatter(
            x=go_x, y=go_y,
            mode="markers+text",
            marker=dict(color=PALETTE["accent1"], size=10, symbol="square",
                        line=dict(color="white", width=1)),
            text=go_labels,
            textposition="middle right",
            textfont=dict(size=8),
            name="GO Terms",
            hovertemplate="<b>%{text}</b><extra>GO Term</extra>",
        )
    )

    fig.update_layout(
        title=f"GeneWalk Network (padj < {padj_cutoff}, sim ≥ {min_sim})",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
        height=max(500, (len(genes) + len(go_terms)) * 12),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def gw_heatmap(
    df: pd.DataFrame,
    padj_cutoff: float = 0.1,
    padj_col: str = "gene_padj",
    gene_col: str = "hgnc_symbol",
    go_name_col: str = "go_name",
    sim_col: str = "sim",
    max_genes: int = 30,
    max_terms: int = 30,
) -> go.Figure:
    """Gene x GO term similarity heatmap.

    Parameters
    ----------
    df : pd.DataFrame
        GeneWalk results.
    padj_cutoff : float
        Filter to significant associations.
    padj_col : str
        Column for adjusted p-value.
    gene_col : str
        Column for gene symbol.
    go_name_col : str
        Column for GO term name.
    sim_col : str
        Column for similarity score.
    max_genes : int
        Maximum genes to display.
    max_terms : int
        Maximum GO terms to display.

    Returns
    -------
    go.Figure
        Interactive heatmap.
    """
    setup_plotly_theme()

    if df.empty:
        fig = go.Figure()
        fig.add_annotation(text="No GeneWalk data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GeneWalk Heatmap — No Data")
        return fig

    required = [gene_col, go_name_col, sim_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        fig = go.Figure()
        fig.add_annotation(text=f"Missing: {', '.join(missing)}", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GeneWalk Heatmap — Missing Columns")
        return fig

    filt = df.copy()
    if padj_col in filt.columns:
        filt = filt[filt[padj_col] < padj_cutoff]

    if filt.empty:
        fig = go.Figure()
        fig.add_annotation(text="No significant associations", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GeneWalk Heatmap — No Significant Results")
        return fig

    # Pivot to matrix
    pivot = filt.pivot_table(
        index=gene_col, columns=go_name_col, values=sim_col, aggfunc="mean"
    ).fillna(0)

    # Limit dimensions
    if pivot.shape[0] > max_genes:
        row_means = pivot.mean(axis=1).nlargest(max_genes).index
        pivot = pivot.loc[row_means]
    if pivot.shape[1] > max_terms:
        col_means = pivot.mean(axis=0).nlargest(max_terms).index
        pivot = pivot[col_means]

    # Truncate long GO term names
    new_cols = [c[:40] + "..." if len(c) > 40 else c for c in pivot.columns]

    fig = go.Figure(
        go.Heatmap(
            z=pivot.values,
            x=new_cols,
            y=pivot.index.tolist(),
            colorscale="Viridis",
            colorbar=dict(title="Similarity"),
            hovertemplate="Gene: %{y}<br>GO: %{x}<br>Similarity: %{z:.3f}<extra></extra>",
        )
    )

    fig.update_layout(
        title=f"Gene × GO Term Similarity (padj < {padj_cutoff})",
        xaxis=dict(tickangle=45, tickfont=dict(size=8)),
        yaxis=dict(tickfont=dict(size=9)),
        height=max(400, pivot.shape[0] * 18 + 150),
    )

    return fig


def gw_domain_pie(
    df: pd.DataFrame,
    padj_cutoff: float = 0.1,
    padj_col: str = "gene_padj",
    domain_col: str = "go_domain",
) -> go.Figure:
    """Pie chart of GO domain distribution (BP, CC, MF).

    Parameters
    ----------
    df : pd.DataFrame
        GeneWalk results.
    padj_cutoff : float
        Filter to significant associations.
    padj_col : str
        Column for adjusted p-value.
    domain_col : str
        Column for GO domain.

    Returns
    -------
    go.Figure
        Pie chart.
    """
    setup_plotly_theme()

    if df.empty or domain_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No GO domain data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GO Domain Distribution — No Data")
        return fig

    filt = df.copy()
    if padj_col in filt.columns:
        filt = filt[filt[padj_col] < padj_cutoff]

    if filt.empty:
        filt = df.copy()

    counts = filt[domain_col].value_counts()

    # Map to friendly names
    domain_names = {
        "biological_process": "Biological Process",
        "cellular_component": "Cellular Component",
        "molecular_function": "Molecular Function",
        "BP": "Biological Process",
        "CC": "Cellular Component",
        "MF": "Molecular Function",
    }
    labels = [domain_names.get(d, d) for d in counts.index]

    # Use category colors where possible
    domain_to_cat = {
        "biological_process": "BP", "Biological Process": "BP", "BP": "BP",
        "cellular_component": "CC", "Cellular Component": "CC", "CC": "CC",
        "molecular_function": "MF", "Molecular Function": "MF", "MF": "MF",
    }
    colors = [
        CATEGORY_COLORS.get(domain_to_cat.get(d, ""), PALETTE["neutral"])
        for d in counts.index
    ]

    fig = go.Figure(
        go.Pie(
            labels=labels,
            values=counts.values,
            marker=dict(colors=colors),
            textinfo="label+percent+value",
            hovertemplate="<b>%{label}</b><br>Count: %{value:,}<br>%{percent}<extra></extra>",
        )
    )

    fig.update_layout(title=f"GO Domain Distribution (padj < {padj_cutoff})")

    return fig


def gw_gene_summary(
    df: pd.DataFrame,
    padj_cutoff: float = 0.1,
    metric: str = "count",
    padj_col: str = "gene_padj",
    gene_col: str = "hgnc_symbol",
    sim_col: str = "sim",
    n: int = 25,
) -> go.Figure:
    """Gene ranking by number of significant GO terms or mean similarity.

    Parameters
    ----------
    df : pd.DataFrame
        GeneWalk results.
    padj_cutoff : float
        Filter to significant associations.
    metric : str
        "count" for number of GO terms, "mean_sim" for average similarity.
    padj_col : str
        Column for adjusted p-value.
    gene_col : str
        Column for gene symbol.
    sim_col : str
        Column for similarity score.
    n : int
        Number of top genes to show.

    Returns
    -------
    go.Figure
        Horizontal bar chart.
    """
    setup_plotly_theme()

    if df.empty or gene_col not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="No GeneWalk data", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GeneWalk Gene Summary — No Data")
        return fig

    filt = df.copy()
    if padj_col in filt.columns:
        filt = filt[filt[padj_col] < padj_cutoff]

    if filt.empty:
        fig = go.Figure()
        fig.add_annotation(text="No significant results", showarrow=False, font=dict(size=16))
        fig.update_layout(title="GeneWalk Gene Summary — No Significant Results")
        return fig

    if metric == "mean_sim" and sim_col in filt.columns:
        agg = filt.groupby(gene_col)[sim_col].mean().nlargest(n).sort_values()
        x_title = "Mean Similarity"
        title = f"Top {n} Genes by Mean Similarity"
    else:
        agg = filt.groupby(gene_col).size().nlargest(n).sort_values()
        x_title = "Number of Significant GO Terms"
        title = f"Top {n} Genes by GO Term Count"

    fig = go.Figure(
        go.Bar(
            y=agg.index.tolist(),
            x=agg.values,
            orientation="h",
            marker_color=PALETTE["accent2"],
            hovertemplate="<b>%{y}</b><br>%{x:.2f}<extra></extra>",
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title=x_title,
        yaxis_title="",
        height=max(400, n * 22),
    )

    return fig
