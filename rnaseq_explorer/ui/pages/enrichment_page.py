"""Enrichment visualization page for RNA-seq Explorer.

Displays GSEA and ORA visualizations with database selection and leading edge tables.
"""

from __future__ import annotations

import streamlit as st
import pandas as pd

from rnaseq_explorer.viz.gsea_viz import (
    nes_bar_chart,
    enrichment_dot_plot,
    leading_edge_table,
    ora_dot_plot,
    enrichment_comparison,
)


def render(settings: dict) -> None:
    """Render the enrichment analysis page.

    Parameters
    ----------
    settings : dict
        Settings dict from sidebar.render_sidebar().
    """
    st.title("Pathway Enrichment")

    gsea_df: pd.DataFrame | None = st.session_state.get("gsea_data")
    ora_df: pd.DataFrame | None = st.session_state.get("ora_data")

    if gsea_df is None and ora_df is None:
        st.info("Upload GSEA or ORA results in the sidebar to view enrichment visualizations.")
        return

    tab_gsea, tab_ora, tab_compare = st.tabs(["GSEA", "ORA", "Comparison"])

    # ---- GSEA Tab ----
    with tab_gsea:
        if gsea_df is not None and not gsea_df.empty:
            # Database filter
            db_col = _detect_col(gsea_df, ["database", "db", "Gene_set", "gene_set"])
            if db_col and db_col in gsea_df.columns:
                databases = sorted(gsea_df[db_col].unique())
                selected_db = st.selectbox("Database", databases, key="gsea_db")
                filtered_gsea = gsea_df[gsea_df[db_col] == selected_db]
            else:
                filtered_gsea = gsea_df

            col1, col2 = st.columns(2)

            with col1:
                st.subheader("NES Bar Chart")
                st.plotly_chart(
                    nes_bar_chart(
                        filtered_gsea,
                        n=settings["n_top_genes"],
                        fdr_cutoff=settings["fdr_cutoff"],
                    ),
                    use_container_width=True,
                )

            with col2:
                st.subheader("Enrichment Dot Plot")
                st.plotly_chart(
                    enrichment_dot_plot(
                        filtered_gsea,
                        n=30,
                        fdr_cutoff=settings["fdr_cutoff"],
                    ),
                    use_container_width=True,
                )

            # Leading edge table
            st.markdown("---")
            st.subheader("Leading Edge Genes")

            term_col = _detect_col(filtered_gsea, ["Term", "term", "Name", "pathway"])
            if term_col and term_col in filtered_gsea.columns:
                pathway_options = filtered_gsea[term_col].tolist()
                if pathway_options:
                    selected_pathway = st.selectbox(
                        "Select Pathway",
                        pathway_options,
                        key="le_pathway",
                    )
                    le_df = leading_edge_table(filtered_gsea, selected_pathway, term_col=term_col)
                    if not le_df.empty:
                        st.dataframe(le_df, use_container_width=True)
                    else:
                        st.info("No leading edge gene data available for this pathway.")

            # Full data table
            st.markdown("---")
            st.subheader("GSEA Results Table")
            st.dataframe(filtered_gsea, use_container_width=True, height=300)
        else:
            st.info("Upload GSEA results to view this tab.")

    # ---- ORA Tab ----
    with tab_ora:
        if ora_df is not None and not ora_df.empty:
            # Database filter
            db_col = _detect_col(ora_df, ["database", "source", "Gene_set", "db"])
            if db_col and db_col in ora_df.columns:
                databases = sorted(ora_df[db_col].unique())
                selected_db = st.selectbox("Database", databases, key="ora_db")
                filtered_ora = ora_df[ora_df[db_col] == selected_db]
            else:
                filtered_ora = ora_df

            st.subheader("ORA Dot Plot")
            st.plotly_chart(
                ora_dot_plot(
                    filtered_ora,
                    n=settings["n_top_genes"],
                    fdr_cutoff=settings["padj_cutoff"],
                ),
                use_container_width=True,
            )

            st.markdown("---")
            st.subheader("ORA Results Table")
            st.dataframe(filtered_ora, use_container_width=True, height=300)
        else:
            st.info("Upload ORA results to view this tab.")

    # ---- Comparison Tab ----
    with tab_compare:
        if gsea_df is not None and not gsea_df.empty:
            nes_col = _detect_col(gsea_df, ["NES", "nes"])
            if nes_col:
                up_df = gsea_df[gsea_df[nes_col] > 0]
                down_df = gsea_df[gsea_df[nes_col] < 0]

                st.subheader("Up vs Down Enrichment Comparison")
                st.plotly_chart(
                    enrichment_comparison(up_df, down_df, n=15, fdr_cutoff=settings["fdr_cutoff"]),
                    use_container_width=True,
                )
            else:
                st.info("NES column not found in GSEA results.")
        else:
            st.info("Upload GSEA results to view enrichment comparison.")


def _detect_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name found in df from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
