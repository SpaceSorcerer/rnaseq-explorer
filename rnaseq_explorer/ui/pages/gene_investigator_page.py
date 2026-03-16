"""Gene Investigator page for RNA-seq Explorer.

Provides a gene search interface that aggregates evidence across all
uploaded data types and displays a visual evidence card.
"""

from __future__ import annotations

import streamlit as st
import pandas as pd

from rnaseq_explorer.viz.gene_investigator import (
    investigate_gene,
    gene_evidence_card,
)


def render(settings: dict) -> None:
    """Render the Gene Investigator page.

    Parameters
    ----------
    settings : dict
        Settings dict from sidebar.render_sidebar().
    """
    st.title("Gene Investigator")
    st.markdown(
        "Search for a gene to aggregate all available evidence: "
        "differential expression, splicing, enrichment, and functional annotations."
    )

    # ---- Gene Search ----
    # Build autocomplete list from available data
    gene_options: list[str] = []

    deseq2_df = st.session_state.get("deseq2_data")
    if deseq2_df is not None and not deseq2_df.empty:
        gene_col = _detect_col(deseq2_df, ["gene_name", "Gene", "gene", "hgnc_symbol"])
        if gene_col:
            gene_options.extend(deseq2_df[gene_col].dropna().unique().tolist())

    rmats_df = st.session_state.get("rmats_data")
    if rmats_df is not None and not rmats_df.empty:
        gene_col = _detect_col(rmats_df, ["GeneID", "geneSymbol", "gene_name"])
        if gene_col:
            gene_options.extend(rmats_df[gene_col].dropna().unique().tolist())

    gene_options = sorted(set(str(g) for g in gene_options if pd.notna(g)))

    col_search, col_btn = st.columns([3, 1])
    with col_search:
        if gene_options:
            gene_name = st.selectbox(
                "Gene name",
                options=[""] + gene_options,
                index=0,
                help="Select or type a gene name.",
                key="gene_search",
            )
        else:
            gene_name = st.text_input(
                "Gene name",
                value="",
                help="Type a gene name to search.",
                key="gene_search_text",
            )

    with col_btn:
        st.markdown("<br>", unsafe_allow_html=True)
        search_clicked = st.button("Investigate", type="primary", key="investigate_btn")

    if not gene_name or not search_clicked:
        st.info("Enter a gene name and click **Investigate** to see all available evidence.")
        return

    st.markdown("---")

    # ---- Gather Evidence ----
    with st.spinner(f"Gathering evidence for {gene_name}..."):
        try:
            evidence = investigate_gene(
                gene_name=gene_name,
                deseq2_results=st.session_state.get("deseq2_data"),
                gsea_results=st.session_state.get("gsea_data"),
                ora_results=st.session_state.get("ora_data"),
                rmats_results=st.session_state.get("rmats_data"),
                genewalk_results=st.session_state.get("genewalk_data"),
            )
        except Exception as e:
            st.warning(f"Could not gather evidence for {gene_name}: {e}. Check that your data has the expected columns.")
            return

    # ---- Evidence Card ----
    try:
        figures, summary_text = gene_evidence_card(evidence)
    except Exception as e:
        st.warning(f"Could not build evidence card: {e}. Check that your data has the expected columns.")
        return

    # Summary header
    deg_info = evidence.get("deg", {})
    if deg_info:
        cols = st.columns(4)
        log2fc = deg_info.get("log2fc")
        padj = deg_info.get("padj")
        direction = deg_info.get("direction", "—")
        biotype = deg_info.get("biotype", "—")

        cols[0].metric("log₂FC", f"{log2fc:.3f}" if log2fc is not None else "—")
        cols[1].metric("padj", f"{padj:.2e}" if padj is not None else "—")
        cols[2].metric("Direction", str(direction))
        cols[3].metric("Biotype", str(biotype))
    else:
        st.warning(f"No differential expression data found for {gene_name}.")

    # Evidence counts
    n_gsea = len(evidence.get("gsea", []))
    n_ora = len(evidence.get("ora", []))
    n_splice = len(evidence.get("splicing", []))
    n_gw = len(evidence.get("genewalk", []))

    ecols = st.columns(4)
    ecols[0].metric("GSEA Pathways", n_gsea)
    ecols[1].metric("ORA Terms", n_ora)
    ecols[2].metric("Splicing Events", n_splice)
    ecols[3].metric("GeneWalk GO Terms", n_gw)

    st.markdown("---")

    # Render figures
    for fig in figures:
        try:
            st.plotly_chart(fig, use_container_width=True)
        except Exception as e:
            st.warning(f"Could not render evidence figure: {e}.")

    # Text summary
    with st.expander("Full Evidence Summary (text)", expanded=False):
        st.text(summary_text)

    # Detail tables
    if evidence.get("gsea"):
        with st.expander(f"GSEA Pathways ({n_gsea})", expanded=False):
            st.dataframe(pd.DataFrame(evidence["gsea"]), use_container_width=True)

    if evidence.get("ora"):
        with st.expander(f"ORA Terms ({n_ora})", expanded=False):
            st.dataframe(pd.DataFrame(evidence["ora"]), use_container_width=True)

    if evidence.get("splicing"):
        with st.expander(f"Splicing Events ({n_splice})", expanded=False):
            st.dataframe(pd.DataFrame(evidence["splicing"]), use_container_width=True)

    if evidence.get("genewalk"):
        with st.expander(f"GeneWalk GO Terms ({n_gw})", expanded=False):
            st.dataframe(pd.DataFrame(evidence["genewalk"]), use_container_width=True)


def _detect_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name found in df from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
