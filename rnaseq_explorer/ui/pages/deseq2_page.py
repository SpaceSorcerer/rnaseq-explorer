"""DESeq2 visualization page for RNA-seq Explorer.

Displays all DEG visualizations with interactive controls and data table.
"""

from __future__ import annotations

import streamlit as st
import pandas as pd

from rnaseq_explorer.viz.deseq2_viz import (
    volcano_plot,
    ma_plot,
    pvalue_distribution,
    log2fc_distribution,
    top_genes_bar,
    biotype_breakdown,
)


def render(settings: dict) -> None:
    """Render the DESeq2 analysis page.

    Parameters
    ----------
    settings : dict
        Settings dict from sidebar.render_sidebar().
    """
    st.title("Differential Expression (DESeq2)")

    deseq2_df: pd.DataFrame | None = st.session_state.get("deseq2_data")

    if deseq2_df is None or deseq2_df.empty:
        st.info("Upload DESeq2 results in the sidebar to view visualizations.")
        return

    # Auto-detect columns
    log2fc_col = _detect_col(deseq2_df, ["log2FoldChange", "log2fc", "logFC"]) or "log2FoldChange"
    padj_col = _detect_col(deseq2_df, ["padj", "pvalue", "p_value"]) or "padj"
    gene_col = _detect_col(deseq2_df, ["gene_name", "Gene", "gene", "hgnc_symbol"]) or "gene_name"
    basemean_col = _detect_col(deseq2_df, ["baseMean", "basemean", "AveExpr"]) or "baseMean"
    pvalue_col = _detect_col(deseq2_df, ["pvalue", "PValue", "p_value"]) or "pvalue"

    # Genes of interest input
    goi_text = st.text_input(
        "Genes of Interest (comma-separated)",
        value="",
        help="Highlight specific genes on the volcano plot.",
    )
    genes_of_interest = (
        [g.strip() for g in goi_text.split(",") if g.strip()] if goi_text else None
    )

    st.markdown("---")

    # ---- Volcano + MA ----
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Volcano Plot")
        st.plotly_chart(
            volcano_plot(
                deseq2_df,
                log2fc_col=log2fc_col,
                padj_col=padj_col,
                gene_col=gene_col,
                log2fc_cutoff=settings["log2fc_cutoff"],
                padj_cutoff=settings["padj_cutoff"],
                genes_of_interest=genes_of_interest,
            ),
            use_container_width=True,
        )

    with col2:
        st.subheader("MA Plot")
        st.plotly_chart(
            ma_plot(
                deseq2_df,
                basemean_col=basemean_col,
                log2fc_col=log2fc_col,
                padj_col=padj_col,
                gene_col=gene_col,
                log2fc_cutoff=settings["log2fc_cutoff"],
                padj_cutoff=settings["padj_cutoff"],
            ),
            use_container_width=True,
        )

    # ---- Distributions ----
    col3, col4 = st.columns(2)

    with col3:
        st.subheader("P-value Distribution")
        st.plotly_chart(
            pvalue_distribution(deseq2_df, pvalue_col=pvalue_col, padj_col=padj_col),
            use_container_width=True,
        )

    with col4:
        st.subheader("log₂FC Distribution")
        st.plotly_chart(
            log2fc_distribution(
                deseq2_df,
                log2fc_col=log2fc_col,
                padj_col=padj_col,
                padj_cutoff=settings["padj_cutoff"],
            ),
            use_container_width=True,
        )

    # ---- Top Genes + Biotype ----
    col5, col6 = st.columns(2)

    with col5:
        st.subheader("Top Genes")
        st.plotly_chart(
            top_genes_bar(
                deseq2_df,
                log2fc_col=log2fc_col,
                padj_col=padj_col,
                gene_col=gene_col,
                n=settings["n_top_genes"],
            ),
            use_container_width=True,
        )

    with col6:
        biotype_col = _detect_col(deseq2_df, ["biotype_group", "biotype", "gene_biotype"])
        direction_col = _detect_col(deseq2_df, ["direction", "Direction"])
        if biotype_col and direction_col:
            st.subheader("Biotype Breakdown")
            st.plotly_chart(
                biotype_breakdown(deseq2_df, biotype_col=biotype_col, direction_col=direction_col),
                use_container_width=True,
            )

    # ---- Data Table ----
    st.markdown("---")
    st.subheader("Data Table")

    # Filter options
    show_sig_only = st.checkbox("Show significant genes only", value=False)

    display_df = deseq2_df.copy()
    if show_sig_only and padj_col in display_df.columns and log2fc_col in display_df.columns:
        display_df = display_df[
            (display_df[padj_col] < settings["padj_cutoff"])
            & (display_df[log2fc_col].abs() >= settings["log2fc_cutoff"])
        ]

    st.dataframe(display_df, use_container_width=True, height=400)

    # Export
    csv = display_df.to_csv(index=False)
    st.download_button(
        "Download as CSV",
        csv,
        file_name="deseq2_filtered.csv",
        mime="text/csv",
    )


def _detect_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name found in df from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
