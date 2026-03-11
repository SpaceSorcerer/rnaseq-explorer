"""Shared sidebar configuration for RNA-seq Explorer.

Renders file upload widgets, threshold sliders, condition selectors,
and export controls in the Streamlit sidebar.
"""

from __future__ import annotations

from typing import Any

import streamlit as st


def render_sidebar() -> dict[str, Any]:
    """Render the sidebar and return a dict of user settings.

    Returns
    -------
    dict[str, Any]
        Dictionary with keys:
        - deseq2_file: UploadedFile or None
        - rmats_file: UploadedFile or None
        - genewalk_file: UploadedFile or None
        - counts_file: UploadedFile or None
        - gsea_file: UploadedFile or None
        - ora_file: UploadedFile or None
        - log2fc_cutoff: float
        - padj_cutoff: float
        - dpsi_cutoff: float
        - fdr_cutoff: float
        - dark_mode: bool
        - n_top_genes: int
    """
    with st.sidebar:
        st.title("RNA-seq Explorer")
        st.markdown("---")

        # ---- File Uploads ----
        st.header("Data Upload")

        deseq2_file = st.file_uploader(
            "DESeq2 Results",
            type=["csv", "tsv", "xlsx", "txt"],
            help="Upload DESeq2 differential expression results.",
            key="deseq2_upload",
        )

        rmats_file = st.file_uploader(
            "rMATS Results",
            type=["csv", "tsv", "xlsx", "txt"],
            help="Upload rMATS alternative splicing results (combined or single event type).",
            key="rmats_upload",
        )

        genewalk_file = st.file_uploader(
            "GeneWalk Results",
            type=["csv", "tsv", "xlsx", "txt"],
            help="Upload GeneWalk functional annotation results.",
            key="genewalk_upload",
        )

        counts_file = st.file_uploader(
            "Normalized Counts",
            type=["csv", "tsv", "xlsx", "txt"],
            help="Normalized expression matrix (genes x samples) for QC plots.",
            key="counts_upload",
        )

        with st.expander("Enrichment Results", expanded=False):
            gsea_file = st.file_uploader(
                "GSEA Results",
                type=["csv", "tsv", "xlsx", "txt"],
                help="Upload GSEA prerank results.",
                key="gsea_upload",
            )

            ora_file = st.file_uploader(
                "ORA Results",
                type=["csv", "tsv", "xlsx", "txt"],
                help="Upload ORA (g:Profiler / Enrichr) results.",
                key="ora_upload",
            )

        st.markdown("---")

        # ---- Thresholds ----
        st.header("Thresholds")

        log2fc_cutoff = st.slider(
            "log₂FC Cutoff",
            min_value=0.0,
            max_value=5.0,
            value=1.0,
            step=0.1,
            help="Absolute log2 fold change threshold for significance.",
            key="log2fc_cutoff",
        )

        padj_cutoff = st.slider(
            "Adjusted p-value Cutoff",
            min_value=0.001,
            max_value=0.10,
            value=0.05,
            step=0.005,
            format="%.3f",
            help="Adjusted p-value threshold for DEG significance.",
            key="padj_cutoff",
        )

        dpsi_cutoff = st.slider(
            "ΔPSI Cutoff",
            min_value=0.0,
            max_value=0.5,
            value=0.1,
            step=0.01,
            help="Absolute delta-PSI threshold for splicing significance.",
            key="dpsi_cutoff",
        )

        fdr_cutoff = st.slider(
            "FDR Cutoff (Splicing / GSEA)",
            min_value=0.01,
            max_value=0.50,
            value=0.05,
            step=0.01,
            help="FDR threshold for rMATS splicing and GSEA pathway significance.",
            key="fdr_cutoff",
        )

        st.markdown("---")

        # ---- Display Settings ----
        st.header("Settings")

        n_top_genes = st.slider(
            "Top Genes to Display",
            min_value=5,
            max_value=50,
            value=20,
            step=5,
            key="n_top_genes",
        )

        dark_mode = st.toggle(
            "Dark Mode",
            value=False,
            key="dark_mode",
        )

        st.markdown("---")

        # ---- Export ----
        st.header("Export")
        st.info("Use the camera icon on each chart to export as PNG.")

    return {
        "deseq2_file": deseq2_file,
        "rmats_file": rmats_file,
        "genewalk_file": genewalk_file,
        "counts_file": counts_file,
        "gsea_file": gsea_file,
        "ora_file": ora_file,
        "log2fc_cutoff": log2fc_cutoff,
        "padj_cutoff": padj_cutoff,
        "dpsi_cutoff": dpsi_cutoff,
        "fdr_cutoff": fdr_cutoff,
        "dark_mode": dark_mode,
        "n_top_genes": n_top_genes,
    }
