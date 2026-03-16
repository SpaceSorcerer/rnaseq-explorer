"""Splicing visualization page for RNA-seq Explorer.

Displays rMATS alternative splicing visualizations with event type filtering.
"""

from __future__ import annotations

import io

import streamlit as st
import pandas as pd

from rnaseq_explorer.viz.rmats_viz import (
    dpsi_volcano,
    event_type_pie,
    dpsi_distribution,
    top_splicing_events,
    genes_by_event_count,
)


def render(settings: dict) -> None:
    """Render the splicing analysis page.

    Parameters
    ----------
    settings : dict
        Settings dict from sidebar.render_sidebar().
    """
    st.title("Alternative Splicing (rMATS)")

    rmats_df: pd.DataFrame | None = st.session_state.get("rmats_data")

    if rmats_df is None or rmats_df.empty:
        st.info("Upload rMATS results in the sidebar to view visualizations.")
        return

    # Auto-detect columns
    dpsi_col = _detect_col(rmats_df, ["IncLevelDifference", "dPSI"]) or "IncLevelDifference"
    fdr_col = _detect_col(rmats_df, ["FDR", "fdr"]) or "FDR"
    gene_col = _detect_col(rmats_df, ["GeneID", "geneSymbol", "gene_name"]) or "GeneID"
    event_type_col = _detect_col(rmats_df, ["event_type", "EventType", "type"]) or "event_type"

    # Event type filter
    if event_type_col in rmats_df.columns:
        all_types = sorted(rmats_df[event_type_col].unique())
        selected_types = st.multiselect(
            "Event Types",
            options=all_types,
            default=all_types,
            help="Filter to specific splicing event types.",
        )
        filtered_df = rmats_df[rmats_df[event_type_col].isin(selected_types)]
    else:
        filtered_df = rmats_df

    st.markdown("---")

    # ---- Volcano + Pie ----
    col1, col2 = st.columns([2, 1])

    with col1:
        st.subheader("dPSI Volcano Plot")
        try:
            st.plotly_chart(
                dpsi_volcano(
                    filtered_df,
                    dpsi_col=dpsi_col,
                    fdr_col=fdr_col,
                    gene_col=gene_col,
                    event_type_col=event_type_col,
                    dpsi_cutoff=settings["dpsi_cutoff"],
                    fdr_cutoff=settings["fdr_cutoff"],
                ),
                use_container_width=True,
            )
        except Exception as e:
            st.warning(f"Could not render dPSI Volcano Plot: {e}. Check that your data has the expected columns.")

    with col2:
        st.subheader("Event Types")
        try:
            st.plotly_chart(
                event_type_pie(filtered_df, event_type_col=event_type_col),
                use_container_width=True,
            )
        except Exception as e:
            st.warning(f"Could not render Event Type Pie: {e}. Check that your data has the expected columns.")

    # ---- Distribution ----
    st.subheader("ΔPSI Distribution")
    try:
        st.plotly_chart(
            dpsi_distribution(filtered_df, dpsi_col=dpsi_col, event_type_col=event_type_col),
            use_container_width=True,
        )
    except Exception as e:
        st.warning(f"Could not render dPSI Distribution: {e}. Check that your data has the expected columns.")

    # ---- Top Events + Gene Counts ----
    col3, col4 = st.columns(2)

    with col3:
        st.subheader("Top Splicing Events")
        try:
            st.plotly_chart(
                top_splicing_events(
                    filtered_df,
                    dpsi_col=dpsi_col,
                    fdr_col=fdr_col,
                    gene_col=gene_col,
                    n=settings["n_top_genes"],
                ),
                use_container_width=True,
            )
        except Exception as e:
            st.warning(f"Could not render Top Splicing Events: {e}. Check that your data has the expected columns.")

    with col4:
        st.subheader("Genes by Event Count")
        try:
            st.plotly_chart(
                genes_by_event_count(
                    filtered_df,
                    gene_col=gene_col,
                    event_type_col=event_type_col,
                ),
                use_container_width=True,
            )
        except Exception as e:
            st.warning(f"Could not render Genes by Event Count: {e}. Check that your data has the expected columns.")

    # ---- Data Table ----
    st.markdown("---")
    st.subheader("Data Table")

    show_sig = st.checkbox("Show significant events only", value=False, key="splice_sig_only")
    display_df = filtered_df.copy()
    if show_sig and fdr_col in display_df.columns and dpsi_col in display_df.columns:
        display_df = display_df[
            (display_df[fdr_col] < settings["fdr_cutoff"])
            & (display_df[dpsi_col].abs() >= settings["dpsi_cutoff"])
        ]

    st.dataframe(display_df, use_container_width=True, height=400)

    with st.expander("Export", expanded=False):
        csv = display_df.to_csv(index=False)
        st.download_button(
            "Download as CSV",
            csv,
            file_name="rmats_filtered.csv",
            mime="text/csv",
            key="splice_csv",
        )

        try:
            excel_buf = io.BytesIO()
            display_df.to_excel(excel_buf, index=False, engine="openpyxl")
            excel_buf.seek(0)
            st.download_button(
                "Download as Excel",
                excel_buf,
                file_name="rmats_filtered.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                key="splice_xlsx",
            )
        except Exception as e:
            st.error(f"Excel export failed: {e}")


def _detect_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first column name found in df from candidates."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
