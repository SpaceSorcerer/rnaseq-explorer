"""GeneWalk visualization page for RNA-seq Explorer.

Displays GeneWalk functional annotation results with network, heatmap,
and summary visualizations.
"""

from __future__ import annotations

import streamlit as st
import pandas as pd

from rnaseq_explorer.viz.genewalk_viz import (
    gw_volcano,
    gw_gene_bar,
    gw_network,
    gw_heatmap,
    gw_domain_pie,
    gw_gene_summary,
)


def render(settings: dict) -> None:
    """Render the GeneWalk page.

    Parameters
    ----------
    settings : dict
        Settings dict from sidebar.render_sidebar().
    """
    st.title("GeneWalk Functional Annotation")

    gw_df: pd.DataFrame | None = st.session_state.get("genewalk_data")

    if gw_df is None or gw_df.empty:
        st.info("Upload GeneWalk results in the sidebar to view functional annotations.")
        return

    # Auto-detect columns
    sim_col = _detect_col(gw_df, ["sim", "similarity"]) or "sim"
    padj_col = _detect_col(gw_df, ["gene_padj", "padj"]) or "gene_padj"
    gene_col = _detect_col(gw_df, ["hgnc_symbol", "gene", "gene_name"]) or "hgnc_symbol"
    go_name_col = _detect_col(gw_df, ["go_name", "GO_name", "go_term"]) or "go_name"

    gw_padj_cutoff = st.slider(
        "GeneWalk padj cutoff",
        min_value=0.01,
        max_value=0.25,
        value=0.1,
        step=0.01,
        key="gw_padj",
    )

    st.markdown("---")

    # ---- Volcano + Domain Pie ----
    col1, col2 = st.columns([2, 1])

    with col1:
        st.subheader("Similarity vs Significance")
        st.plotly_chart(
            gw_volcano(gw_df, sim_col=sim_col, padj_col=padj_col, gene_col=gene_col, padj_cutoff=gw_padj_cutoff),
            use_container_width=True,
        )

    with col2:
        st.subheader("GO Domain Distribution")
        st.plotly_chart(
            gw_domain_pie(gw_df, padj_cutoff=gw_padj_cutoff, padj_col=padj_col),
            use_container_width=True,
        )

    st.markdown("---")

    # ---- Gene Summary ----
    col3, col4 = st.columns(2)

    with col3:
        st.subheader("Gene Summary (by GO Term Count)")
        st.plotly_chart(
            gw_gene_summary(
                gw_df, padj_cutoff=gw_padj_cutoff, metric="count",
                padj_col=padj_col, gene_col=gene_col, sim_col=sim_col,
            ),
            use_container_width=True,
        )

    with col4:
        st.subheader("Gene Summary (by Mean Similarity)")
        st.plotly_chart(
            gw_gene_summary(
                gw_df, padj_cutoff=gw_padj_cutoff, metric="mean_sim",
                padj_col=padj_col, gene_col=gene_col, sim_col=sim_col,
            ),
            use_container_width=True,
        )

    st.markdown("---")

    # ---- Per-Gene Exploration ----
    st.subheader("Per-Gene GO Terms")

    if gene_col in gw_df.columns:
        gene_options = sorted(gw_df[gene_col].dropna().unique())
        selected_gene = st.selectbox("Select Gene", gene_options, key="gw_gene_select")

        if selected_gene:
            st.plotly_chart(
                gw_gene_bar(
                    gw_df, gene=selected_gene,
                    padj_col=padj_col, go_name_col=go_name_col, sim_col=sim_col,
                ),
                use_container_width=True,
            )

    st.markdown("---")

    # ---- Network ----
    st.subheader("Gene-GO Term Network")
    st.plotly_chart(
        gw_network(
            gw_df, padj_cutoff=gw_padj_cutoff,
            padj_col=padj_col, gene_col=gene_col, go_name_col=go_name_col, sim_col=sim_col,
        ),
        use_container_width=True,
    )

    st.markdown("---")

    # ---- Heatmap ----
    st.subheader("Gene x GO Term Heatmap")
    st.plotly_chart(
        gw_heatmap(
            gw_df, padj_cutoff=gw_padj_cutoff,
            padj_col=padj_col, gene_col=gene_col, go_name_col=go_name_col, sim_col=sim_col,
        ),
        use_container_width=True,
    )

    # ---- Data Table ----
    st.markdown("---")
    st.subheader("Data Table")
    st.dataframe(gw_df, use_container_width=True, height=400)

    csv = gw_df.to_csv(index=False)
    st.download_button(
        "Download as CSV",
        csv,
        file_name="genewalk_results.csv",
        mime="text/csv",
        key="gw_download",
    )


def _detect_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name found in df from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
