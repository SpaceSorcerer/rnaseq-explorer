"""Cross-condition comparison page for RNA-seq Explorer.

Displays concordance heatmaps, log2FC scatter plots, and overlap bar charts
when multiple DESeq2 result sets are available.
"""

from __future__ import annotations

import streamlit as st
import pandas as pd

from rnaseq_explorer.viz.cross_condition_viz import (
    direction_concordance_heatmap,
    log2fc_scatter,
    overlap_bar,
)


def render(settings: dict) -> None:
    """Render the cross-condition comparison page.

    Parameters
    ----------
    settings : dict
        Settings dict from sidebar.render_sidebar().
    """
    st.title("Cross-Condition Comparison")

    # Check for multiple condition datasets in session_state
    condition_data: dict[str, pd.DataFrame] = st.session_state.get("condition_datasets", {})

    if len(condition_data) < 2:
        st.info(
            "Cross-condition comparison requires at least 2 DESeq2 result sets. "
            "Upload multiple result files, or load a multi-condition dataset."
        )

        # Offer manual upload of additional conditions
        st.subheader("Upload Additional Conditions")
        new_name = st.text_input("Condition name", value="", key="new_cond_name")
        new_file = st.file_uploader(
            "DESeq2 results for this condition",
            type=["csv", "tsv", "xlsx"],
            key="new_cond_file",
        )

        if new_name and new_file:
            try:
                new_df = _load_file(new_file)
                if new_df is not None:
                    if "condition_datasets" not in st.session_state:
                        st.session_state["condition_datasets"] = {}
                    st.session_state["condition_datasets"][new_name] = new_df
                    st.success(f"Added condition '{new_name}' ({len(new_df):,} genes).")
                    st.rerun()
            except Exception as e:
                st.error(f"Error loading file: {e}")

        # Also add main deseq2 data if available
        deseq2_df = st.session_state.get("deseq2_data")
        if deseq2_df is not None and not deseq2_df.empty:
            if "condition_datasets" not in st.session_state:
                st.session_state["condition_datasets"] = {}
            if "Primary" not in st.session_state["condition_datasets"]:
                st.session_state["condition_datasets"]["Primary"] = deseq2_df
                st.info("Primary DESeq2 data added. Upload at least one more condition above.")
        return

    condition_names = list(condition_data.keys())
    st.success(f"Loaded {len(condition_names)} conditions: {', '.join(condition_names)}")

    # Auto-detect shared columns
    first_df = list(condition_data.values())[0]
    gene_col = _detect_col(first_df, ["gene_name", "Gene", "gene", "hgnc_symbol"]) or "gene_name"
    log2fc_col = _detect_col(first_df, ["log2FoldChange", "log2fc", "logFC"]) or "log2FoldChange"

    st.markdown("---")

    # ---- Direction Concordance Heatmap ----
    st.subheader("Direction Concordance")

    try:
        # Build log2FC matrix: genes x conditions
        fc_frames = {}
        for cond_name, cond_df in condition_data.items():
            lfc = _detect_col(cond_df, ["log2FoldChange", "log2fc", "logFC"])
            gc = _detect_col(cond_df, ["gene_name", "Gene", "gene", "hgnc_symbol"])
            if lfc and gc:
                fc_frames[cond_name] = cond_df.set_index(gc)[lfc]

        if len(fc_frames) >= 2:
            fc_matrix = pd.DataFrame(fc_frames).dropna()
            concordance = fc_matrix.corr(method="pearson")

            st.plotly_chart(
                direction_concordance_heatmap(concordance, condition_labels=condition_names),
                use_container_width=True,
            )
        else:
            st.warning("Not enough shared genes for concordance analysis.")
    except Exception as e:
        st.error(f"Concordance computation failed: {e}")

    st.markdown("---")

    # ---- Pairwise log2FC Scatter ----
    st.subheader("Pairwise log₂FC Scatter")

    col1, col2 = st.columns(2)
    with col1:
        cond1 = st.selectbox("Condition 1", condition_names, index=0, key="scatter_c1")
    with col2:
        cond2_options = [c for c in condition_names if c != cond1]
        cond2 = st.selectbox("Condition 2", cond2_options, index=0, key="scatter_c2")

    if cond1 and cond2:
        df1 = condition_data[cond1]
        df2 = condition_data[cond2]
        gc1 = _detect_col(df1, ["gene_name", "Gene", "gene"]) or gene_col
        gc2 = _detect_col(df2, ["gene_name", "Gene", "gene"]) or gene_col
        lfc1 = _detect_col(df1, ["log2FoldChange", "log2fc"]) or log2fc_col
        lfc2 = _detect_col(df2, ["log2FoldChange", "log2fc"]) or log2fc_col

        # Standardize column names for merge
        df1_std = df1.rename(columns={gc1: "gene_name", lfc1: "log2FoldChange"})
        df2_std = df2.rename(columns={gc2: "gene_name", lfc2: "log2FoldChange"})

        try:
            st.plotly_chart(
                log2fc_scatter(
                    df1_std, df2_std,
                    gene_col="gene_name",
                    log2fc_col="log2FoldChange",
                    cond1_name=cond1,
                    cond2_name=cond2,
                ),
                use_container_width=True,
            )
        except Exception as e:
            st.warning(f"Could not render log2FC Scatter: {e}. Check that your data has the expected columns.")

    st.markdown("---")

    # ---- Overlap Bar ----
    st.subheader("Gene Overlap")

    try:
        sig_genes = {}
        for cond_name, cond_df in condition_data.items():
            padj = _detect_col(cond_df, ["padj", "pvalue"])
            lfc = _detect_col(cond_df, ["log2FoldChange", "log2fc"])
            gc = _detect_col(cond_df, ["gene_name", "Gene", "gene"])
            if padj and lfc and gc:
                sig = cond_df[
                    (cond_df[padj] < settings["padj_cutoff"])
                    & (cond_df[lfc].abs() >= settings["log2fc_cutoff"])
                ]
                sig_genes[cond_name] = set(sig[gc].tolist())

        if len(sig_genes) >= 2:
            all_sig = set()
            for s in sig_genes.values():
                all_sig |= s

            overlap_rows = []
            for cond_name, genes in sig_genes.items():
                other_genes = set()
                for other_name, other_set in sig_genes.items():
                    if other_name != cond_name:
                        other_genes |= other_set
                shared = len(genes & other_genes)
                unique = len(genes) - shared
                overlap_rows.append({
                    "condition": cond_name,
                    "unique": unique,
                    "shared": shared,
                })

            overlap_df = pd.DataFrame(overlap_rows)
            st.plotly_chart(
                overlap_bar(overlap_df),
                use_container_width=True,
            )
        else:
            st.info("Need significance data from at least 2 conditions for overlap analysis.")
    except Exception as e:
        st.error(f"Overlap computation failed: {e}")


def _load_file(uploaded_file) -> pd.DataFrame | None:
    """Load an uploaded file into a DataFrame."""
    name = uploaded_file.name.lower()
    if name.endswith(".csv"):
        return pd.read_csv(uploaded_file)
    elif name.endswith(".tsv") or name.endswith(".txt"):
        return pd.read_csv(uploaded_file, sep="\t")
    elif name.endswith(".xlsx"):
        return pd.read_excel(uploaded_file)
    return None


def _detect_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name found in df from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
