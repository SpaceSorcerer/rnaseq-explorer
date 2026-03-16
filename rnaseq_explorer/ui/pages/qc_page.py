"""QC visualization page for RNA-seq Explorer.

Displays PCA plots, sample correlation heatmaps, and top DEG heatmaps.
Requires normalized count matrix upload.
"""

from __future__ import annotations

import streamlit as st
import pandas as pd
import numpy as np

from rnaseq_explorer.viz.qc_viz import (
    pca_plot,
    correlation_heatmap,
    top_deg_heatmap,
)


def render(settings: dict) -> None:
    """Render the QC page.

    Parameters
    ----------
    settings : dict
        Settings dict from sidebar.render_sidebar().
    """
    st.title("Quality Control")

    counts_df: pd.DataFrame | None = st.session_state.get("counts_data")
    deseq2_df: pd.DataFrame | None = st.session_state.get("deseq2_data")

    if counts_df is None or counts_df.empty:
        st.info(
            "Upload a normalized counts matrix (genes x samples) in the sidebar "
            "to view QC visualizations."
        )
        return

    # Ensure numeric matrix
    numeric_df = counts_df.select_dtypes(include=[np.number])
    if numeric_df.empty:
        st.error("Counts matrix has no numeric columns.")
        return

    # Condition mapping
    st.sidebar.markdown("---")
    st.sidebar.subheader("QC: Sample Conditions")
    sample_names = numeric_df.columns.tolist()

    # Auto-detect or ask for condition labels
    condition_input = st.text_area(
        "Condition labels (one per sample, same order as columns)",
        value="\n".join(sample_names),
        height=150,
        help="Enter one condition label per line, matching sample column order.",
    )
    conditions = [c.strip() for c in condition_input.strip().split("\n") if c.strip()]
    if len(conditions) != len(sample_names):
        st.warning(
            f"Number of conditions ({len(conditions)}) does not match "
            f"number of samples ({len(sample_names)}). Using sample names as conditions."
        )
        conditions = sample_names

    # ---- PCA ----
    st.subheader("Principal Component Analysis")

    try:
        from sklearn.decomposition import PCA as SkPCA

        # Transpose: rows=samples, cols=genes
        X = numeric_df.T.fillna(0).values
        pca_model = SkPCA(n_components=min(2, X.shape[1]))
        coords = pca_model.fit_transform(X)
        var_explained = pca_model.explained_variance_ratio_ * 100

        pca_df = pd.DataFrame({
            "PC1": coords[:, 0],
            "PC2": coords[:, 1] if coords.shape[1] > 1 else 0,
            "sample": sample_names,
            "condition": conditions,
        })

        st.plotly_chart(
            pca_plot(
                pca_df,
                var1=var_explained[0],
                var2=var_explained[1] if len(var_explained) > 1 else 0,
                sample_col="sample",
                condition_col="condition",
            ),
            use_container_width=True,
        )
    except ImportError:
        st.warning("Install scikit-learn (`pip install scikit-learn`) for PCA.")
    except Exception as e:
        st.error(f"PCA computation failed: {e}")

    st.markdown("---")

    # ---- Correlation Heatmap ----
    st.subheader("Sample Correlation")

    try:
        corr = numeric_df.corr(method="pearson")
        st.plotly_chart(
            correlation_heatmap(corr, sample_labels=sample_names),
            use_container_width=True,
        )
    except Exception as e:
        st.warning(f"Could not render Correlation Heatmap: {e}. Check that your counts matrix has numeric columns.")

    st.markdown("---")

    # ---- Top DEG Heatmap ----
    st.subheader("Top DEG Expression Heatmap")

    if deseq2_df is not None and not deseq2_df.empty:
        padj_col = _detect_col(deseq2_df, ["padj", "pvalue"])
        log2fc_col = _detect_col(deseq2_df, ["log2FoldChange", "log2fc"])
        gene_col = _detect_col(deseq2_df, ["gene_name", "Gene", "gene"])

        if padj_col and log2fc_col and gene_col:
            sig = deseq2_df[
                (deseq2_df[padj_col] < settings["padj_cutoff"])
            ].nlargest(50, log2fc_col, "all")
            gene_list = sig[gene_col].tolist()

            # Set index if needed
            expr = counts_df.copy()
            if gene_col in expr.columns:
                expr = expr.set_index(gene_col)

            try:
                st.plotly_chart(
                    top_deg_heatmap(
                        expr.select_dtypes(include=[np.number]),
                        gene_list=gene_list,
                        sample_labels=sample_names,
                        condition_labels=conditions,
                    ),
                    use_container_width=True,
                )
            except Exception as e:
                st.warning(f"Could not render Top DEG Heatmap: {e}. Check that your counts matrix and DESeq2 results have matching gene names.")
        else:
            st.info("DESeq2 results required for DEG heatmap.")
    else:
        st.info("Upload DESeq2 results to see top DEG expression heatmap.")


def _detect_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name found in df from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
