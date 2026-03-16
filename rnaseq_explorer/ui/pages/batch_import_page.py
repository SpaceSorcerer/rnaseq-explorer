"""Batch Import page -- load Phase 3 pipeline output directory."""

import streamlit as st
from pathlib import Path


def render(settings: dict):
    """Render the Batch Import page."""
    st.header("Import Phase 3 Pipeline Results")
    st.markdown(
        "Load a complete Phase 3 output directory to explore all conditions, "
        "figures, enrichment results, and cross-condition comparisons interactively."
    )

    # Directory input
    col1, col2 = st.columns([3, 1])
    with col1:
        default_path = st.session_state.get("phase3_directory", "")
        directory = st.text_input(
            "Phase 3 output directory path:",
            value=default_path,
            placeholder="F:\\RNA-SEQ-ANALYSIS\\MIAT-KD-RNAseq\\Phase3_4condition_2026_03_15_...",
            help="Enter the full path to a Phase 3 pipeline output directory.",
        )
    with col2:
        st.markdown("<br>", unsafe_allow_html=True)  # vertical alignment
        load_btn = st.button("Load Results", type="primary", use_container_width=True)

    if load_btn and directory:
        directory = directory.strip().strip('"').strip("'")
        st.session_state["phase3_directory"] = directory

        with st.spinner("Scanning directory..."):
            try:
                from rnaseq_explorer.engine.batch_loader import load_phase3_output
                data = load_phase3_output(directory)
                st.session_state["phase3_data"] = data

                # Also populate condition_datasets for cross-condition page compatibility
                condition_datasets = {}
                for cond_name, cd in data["condition_data"].items():
                    if cd["deseq2"] is not None:
                        condition_datasets[cond_name] = cd["deseq2"]
                st.session_state["condition_datasets"] = condition_datasets

                st.success(f"Loaded {data['summary']['n_conditions']} conditions from {Path(directory).name}")
            except FileNotFoundError as e:
                st.error(f"Directory not found: {e}")
                return
            except Exception as e:
                st.error(f"Error loading results: {e}")
                return

    # Show loaded data summary
    if "phase3_data" in st.session_state:
        data = st.session_state["phase3_data"]
        summary = data["summary"]

        st.divider()
        st.subheader("Loaded Results Summary")

        # Metric cards (WSF pattern)
        cols = st.columns(5)
        cols[0].metric("Conditions", summary["n_conditions"])
        cols[1].metric("Total Figures", summary["n_figures"])
        cols[2].metric("Prism Files", summary["n_prism_files"])
        cols[3].metric("GSEA", "Yes" if summary["has_gsea"] else "No")
        cols[4].metric("DE-Splicing Shortlist", "Yes" if summary["has_shortlist"] else "No")

        # Per-condition details
        st.subheader("Conditions")
        for cond_name in data["conditions"]:
            cd = data["condition_data"][cond_name]
            deg_count = summary.get(f"deg_count_{cond_name}", "?")
            n_figs = len(cd["figures"])
            rmats_types = list(cd["rmats"].keys())

            with st.expander(f"{cond_name} — {deg_count} DEGs, {n_figs} figures"):
                if cd["deseq2"] is not None:
                    st.caption(f"DESeq2: {len(cd['deseq2'])} genes")
                    st.dataframe(cd["deseq2"].head(10), use_container_width=True, height=200)
                if rmats_types:
                    st.caption(f"rMATS event types: {', '.join(rmats_types)}")

        # Navigation hint
        st.info(
            "Results loaded. Navigate to **Cross-Condition**, **Enrichment**, "
            "or **Gene Investigator** tabs to explore interactively."
        )
