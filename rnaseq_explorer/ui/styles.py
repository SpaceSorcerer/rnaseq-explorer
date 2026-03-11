"""CSS styling for the RNA-seq Explorer Streamlit app.

Provides custom CSS for metric cards, sidebar, tabs, and dark/light mode support.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Custom CSS injected into the Streamlit app
# ---------------------------------------------------------------------------

CUSTOM_CSS = """
<style>
/* ---- Metric Cards ---- */
div[data-testid="stMetric"] {
    background-color: #f8f9fa;
    border-radius: 8px;
    padding: 12px 16px;
    border-left: 4px solid #999999;
    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
}

/* Up-regulated card */
div[data-testid="stMetric"]:has(label:contains("Up")) {
    border-left-color: #D55E00;
}

/* Down-regulated card */
div[data-testid="stMetric"]:has(label:contains("Down")) {
    border-left-color: #0072B2;
}

/* Total card */
div[data-testid="stMetric"]:has(label:contains("Total")) {
    border-left-color: #009E73;
}

/* Splicing card */
div[data-testid="stMetric"]:has(label:contains("Splic")) {
    border-left-color: #E69F00;
}

/* ---- Sidebar ---- */
section[data-testid="stSidebar"] {
    background-color: #fafafa;
    border-right: 1px solid #e0e0e0;
}

section[data-testid="stSidebar"] h1,
section[data-testid="stSidebar"] h2,
section[data-testid="stSidebar"] h3 {
    color: #333333;
    font-weight: 600;
}

section[data-testid="stSidebar"] .stSlider label {
    font-size: 0.9rem;
}

/* ---- Tabs ---- */
button[data-baseweb="tab"] {
    font-size: 0.95rem;
    font-weight: 500;
    padding: 8px 16px;
}

button[data-baseweb="tab"][aria-selected="true"] {
    border-bottom: 3px solid #0072B2;
    color: #0072B2;
}

/* ---- DataFrames ---- */
div[data-testid="stDataFrame"] {
    border-radius: 6px;
    overflow: hidden;
}

/* ---- Dark mode overrides ---- */
@media (prefers-color-scheme: dark) {
    div[data-testid="stMetric"] {
        background-color: #1e1e1e;
        border-left-color: #666666;
    }

    section[data-testid="stSidebar"] {
        background-color: #121212;
        border-right-color: #333333;
    }

    section[data-testid="stSidebar"] h1,
    section[data-testid="stSidebar"] h2,
    section[data-testid="stSidebar"] h3 {
        color: #e0e0e0;
    }
}

/* ---- Utility classes ---- */
.info-box {
    background-color: #e8f4fd;
    border-left: 4px solid #0072B2;
    padding: 12px 16px;
    border-radius: 4px;
    margin: 8px 0;
}

.warning-box {
    background-color: #fff3cd;
    border-left: 4px solid #E69F00;
    padding: 12px 16px;
    border-radius: 4px;
    margin: 8px 0;
}

.success-box {
    background-color: #d4edda;
    border-left: 4px solid #009E73;
    padding: 12px 16px;
    border-radius: 4px;
    margin: 8px 0;
}

/* ---- Hide Streamlit boilerplate ---- */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>
"""

DARK_MODE_CSS = """
<style>
.stApp {
    background-color: #0e1117;
    color: #e0e0e0;
}

div[data-testid="stMetric"] {
    background-color: #1e1e1e;
}

section[data-testid="stSidebar"] {
    background-color: #121212;
    border-right-color: #333333;
}
</style>
"""


def inject_css(dark_mode: bool = False) -> str:
    """Return the CSS string to inject into the Streamlit app.

    Parameters
    ----------
    dark_mode : bool
        If True, include additional dark mode overrides.

    Returns
    -------
    str
        HTML/CSS string for st.markdown(unsafe_allow_html=True).
    """
    css = CUSTOM_CSS
    if dark_mode:
        css += DARK_MODE_CSS
    return css
