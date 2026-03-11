"""Main Streamlit entry point for RNA-seq Explorer.

Multi-page app with sidebar configuration and tabbed navigation.
Run with: streamlit run rnaseq_explorer/ui/app.py
"""

from __future__ import annotations

import io

import pandas as pd
import streamlit as st

from rnaseq_explorer.ui.styles import inject_css
from rnaseq_explorer.ui.sidebar import render_sidebar
from rnaseq_explorer.viz.theme import setup_plotly_theme

# ---------------------------------------------------------------------------
# Page configuration (must be first Streamlit call)
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="RNA-seq Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner="Loading data...")
def _load_uploaded_file(file_bytes: bytes, file_name: str) -> pd.DataFrame | None:
    """Parse an uploaded file into a DataFrame.

    Parameters
    ----------
    file_bytes : bytes
        Raw file content.
    file_name : str
        Original filename (used to detect format).

    Returns
    -------
    pd.DataFrame or None
        Parsed DataFrame, or None on failure.
    """
    try:
        name_lower = file_name.lower()
        if name_lower.endswith(".csv"):
            return pd.read_csv(io.BytesIO(file_bytes))
        elif name_lower.endswith((".tsv", ".txt")):
            return pd.read_csv(io.BytesIO(file_bytes), sep="\t")
        elif name_lower.endswith(".xlsx"):
            return pd.read_excel(io.BytesIO(file_bytes))
        else:
            # Try CSV first, then TSV
            try:
                return pd.read_csv(io.BytesIO(file_bytes))
            except Exception:
                return pd.read_csv(io.BytesIO(file_bytes), sep="\t")
    except Exception as e:
        st.error(f"Failed to load {file_name}: {e}")
        return None


def _store_uploaded_data(settings: dict) -> None:
    """Parse uploaded files and store in session_state.

    Parameters
    ----------
    settings : dict
        Settings dict from render_sidebar().
    """
    file_keys = {
        "deseq2_file": "deseq2_data",
        "rmats_file": "rmats_data",
        "genewalk_file": "genewalk_data",
        "counts_file": "counts_data",
        "gsea_file": "gsea_data",
        "ora_file": "ora_data",
    }

    for file_key, state_key in file_keys.items():
        uploaded = settings.get(file_key)
        if uploaded is not None:
            file_bytes = uploaded.getvalue()
            df = _load_uploaded_file(file_bytes, uploaded.name)
            if df is not None:
                st.session_state[state_key] = df
        elif state_key not in st.session_state:
            st.session_state[state_key] = None


# ---------------------------------------------------------------------------
# Main application
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the RNA-seq Explorer Streamlit app."""
    # Render sidebar and get settings
    settings = render_sidebar()

    # Inject custom CSS
    st.markdown(inject_css(dark_mode=settings.get("dark_mode", False)), unsafe_allow_html=True)

    # Setup Plotly theme
    setup_plotly_theme(dark_mode=settings.get("dark_mode", False))

    # Load uploaded data into session_state
    _store_uploaded_data(settings)

    # Import pages
    from rnaseq_explorer.ui.pages import overview
    from rnaseq_explorer.ui.pages import deseq2_page
    from rnaseq_explorer.ui.pages import splicing_page
    from rnaseq_explorer.ui.pages import enrichment_page
    from rnaseq_explorer.ui.pages import qc_page
    from rnaseq_explorer.ui.pages import cross_condition_page
    from rnaseq_explorer.ui.pages import genewalk_page
    from rnaseq_explorer.ui.pages import gene_investigator_page

    # Build tab list based on available data
    tab_names = ["Overview", "DESeq2", "Splicing", "Enrichment"]

    has_counts = st.session_state.get("counts_data") is not None
    has_genewalk = st.session_state.get("genewalk_data") is not None
    has_multi = len(st.session_state.get("condition_datasets", {})) >= 2

    if has_counts:
        tab_names.append("QC")
    tab_names.append("Cross-Condition")
    if has_genewalk:
        tab_names.append("GeneWalk")
    tab_names.append("Gene Investigator")

    tabs = st.tabs(tab_names)

    tab_idx = 0

    with tabs[tab_idx]:
        overview.render(settings)
    tab_idx += 1

    with tabs[tab_idx]:
        deseq2_page.render(settings)
    tab_idx += 1

    with tabs[tab_idx]:
        splicing_page.render(settings)
    tab_idx += 1

    with tabs[tab_idx]:
        enrichment_page.render(settings)
    tab_idx += 1

    if has_counts:
        with tabs[tab_idx]:
            qc_page.render(settings)
        tab_idx += 1

    with tabs[tab_idx]:
        cross_condition_page.render(settings)
    tab_idx += 1

    if has_genewalk:
        with tabs[tab_idx]:
            genewalk_page.render(settings)
        tab_idx += 1

    with tabs[tab_idx]:
        gene_investigator_page.render(settings)


if __name__ == "__main__":
    main()
