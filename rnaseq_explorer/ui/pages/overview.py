"""Overview page for RNA-seq Explorer.

Displays summary metrics and quick-look charts for uploaded data.
"""

from __future__ import annotations

import streamlit as st
import pandas as pd

from rnaseq_explorer.viz.deseq2_viz import volcano_plot, top_genes_bar
from rnaseq_explorer.viz.gsea_viz import nes_bar_chart


def render(settings: dict) -> None:
    """Render the overview page.

    Parameters
    ----------
    settings : dict
        Settings dict from sidebar.render_sidebar().
    """
    st.title("Overview")

    deseq2_df: pd.DataFrame | None = st.session_state.get("deseq2_data")
    rmats_df: pd.DataFrame | None = st.session_state.get("rmats_data")
    gsea_df: pd.DataFrame | None = st.session_state.get("gsea_data")

    if deseq2_df is None and rmats_df is None:
        st.info("Upload data using the sidebar to get started.")
        return

    # ---- Summary Metrics ----
    st.subheader("Summary Metrics")

    cols = st.columns(4)

    if deseq2_df is not None and not deseq2_df.empty:
        padj_col = _detect_col(deseq2_df, ["padj", "pvalue", "p_value"])
        log2fc_col = _detect_col(deseq2_df, ["log2FoldChange", "log2fc", "logFC"])

        total_genes = len(deseq2_df)

        if padj_col and log2fc_col:
            sig = deseq2_df[
                (deseq2_df[padj_col] < settings["padj_cutoff"])
                & (deseq2_df[log2fc_col].abs() >= settings["log2fc_cutoff"])
            ]
            n_up = int((sig[log2fc_col] > 0).sum())
            n_down = int((sig[log2fc_col] < 0).sum())
        else:
            n_up = n_down = 0

        cols[0].metric("Total Genes", f"{total_genes:,}")
        cols[1].metric("Up-regulated", f"{n_up:,}")
        cols[2].metric("Down-regulated", f"{n_down:,}")
    else:
        cols[0].metric("Total Genes", "—")
        cols[1].metric("Up-regulated", "—")
        cols[2].metric("Down-regulated", "—")

    if rmats_df is not None and not rmats_df.empty:
        fdr_col = _detect_col(rmats_df, ["FDR", "fdr"])
        dpsi_col = _detect_col(rmats_df, ["IncLevelDifference", "dPSI"])
        if fdr_col and dpsi_col:
            n_splice = int(
                ((rmats_df[fdr_col] < settings["fdr_cutoff"])
                 & (rmats_df[dpsi_col].abs() >= settings["dpsi_cutoff"])).sum()
            )
        else:
            n_splice = len(rmats_df)
        cols[3].metric("Splicing Events", f"{n_splice:,}")
    else:
        cols[3].metric("Splicing Events", "—")

    st.markdown("---")

    # ---- Quick-Look Charts ----
    st.subheader("Quick-Look Charts")

    if deseq2_df is not None and not deseq2_df.empty:
        log2fc_col = _detect_col(deseq2_df, ["log2FoldChange", "log2fc", "logFC"])
        padj_col = _detect_col(deseq2_df, ["padj", "pvalue"])
        gene_col = _detect_col(deseq2_df, ["gene_name", "Gene", "gene", "hgnc_symbol"])

        if log2fc_col and padj_col:
            col_a, col_b = st.columns(2)

            with col_a:
                st.plotly_chart(
                    volcano_plot(
                        deseq2_df,
                        log2fc_col=log2fc_col,
                        padj_col=padj_col,
                        gene_col=gene_col or "gene_name",
                        log2fc_cutoff=settings["log2fc_cutoff"],
                        padj_cutoff=settings["padj_cutoff"],
                    ),
                    use_container_width=True,
                )

            with col_b:
                st.plotly_chart(
                    top_genes_bar(
                        deseq2_df,
                        log2fc_col=log2fc_col,
                        padj_col=padj_col,
                        gene_col=gene_col or "gene_name",
                        n=10,
                    ),
                    use_container_width=True,
                )

    if gsea_df is not None and not gsea_df.empty:
        st.plotly_chart(
            nes_bar_chart(gsea_df, n=10, fdr_cutoff=settings["fdr_cutoff"]),
            use_container_width=True,
        )


def _detect_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name found in df from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
