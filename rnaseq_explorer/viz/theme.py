"""Unified visualization theme for RNA-seq Explorer.

Provides consistent styling across all charts:
- Colorblind-safe palette (Okabe-Ito)
- Consistent font stack
- Standard axis formatting
- Dark/light mode support for Plotly and Matplotlib
"""

from __future__ import annotations

from typing import Sequence

import matplotlib.pyplot as plt
import seaborn as sns

try:
    import plotly.graph_objects as go
    import plotly.io as pio

    _PLOTLY_AVAILABLE = True
except ImportError:
    _PLOTLY_AVAILABLE = False

# ---------------------------------------------------------------------------
# Okabe-Ito colorblind-safe palette
# ---------------------------------------------------------------------------

PALETTE = {
    "up": "#D55E00",  # Vermillion (upregulated)
    "down": "#0072B2",  # Blue (downregulated)
    "neutral": "#999999",  # Gray (non-significant)
    "highlight": "#E69F00",  # Orange (highlighted/selected)
    "accent1": "#56B4E9",  # Sky blue
    "accent2": "#009E73",  # Bluish green
    "accent3": "#F0E442",  # Yellow
    "accent4": "#CC79A7",  # Reddish purple
}

# Legacy color aliases (matching original pipeline defaults)
COLOR_UP = "#E69F00"  # upregulated / included (Okabe-Ito orange)
COLOR_DOWN = "#0072B2"  # downregulated / excluded (Okabe-Ito blue)
COLOR_NS = "#BFBFBF"  # not significant

CONDITION_COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#E69F00",
    "#56B4E9",
    "#CC79A7",
    "#F0E442",
    "#999999",
]

# rMATS event type colors
EVENT_COLORS = {
    "SE": "#E64B35",
    "A3SS": "#4DBBD5",
    "A5SS": "#00A087",
    "RI": "#3C5488",
    "MXE": "#F39B7F",
}

# Biotype group colors
BIOTYPE_COLORS = {
    "Protein Coding": "#4C72B0",
    "lncRNA": "#DD8452",
    "Pseudogene": "#55A868",
    "Small ncRNA": "#C44E52",
    "Other": "#8172B2",
}

# ORA database category colors (Okabe-Ito)
CATEGORY_COLORS = {
    "BP": "#0072B2",
    "CC": "#E69F00",
    "MF": "#009E73",
    "KEGG": "#CC79A7",
    "Reactome": "#56B4E9",
}

# Venn diagram colors
VENN_COLORS_2 = {"10": "#0072B2", "01": "#E69F00", "11": "#009E73"}
VENN_COLORS_3 = {
    "100": "#0072B2",
    "010": "#E69F00",
    "001": "#56B4E9",
    "110": "#009E73",
    "101": "#CC79A7",
    "011": "#F0E442",
    "111": "#D55E00",
}

# ---------------------------------------------------------------------------
# Font configuration
# ---------------------------------------------------------------------------

FONT_FAMILY = "sans-serif"
FONT_SIZE_DEFAULT = 12
FONT_SIZE_TITLE = 14
FONT_SIZE_LABEL = 12
FONT_SIZE_TICK = 10
FONT_SIZE_ANNOTATION = 9
FONT_SIZE_SMALL = 8

# ---------------------------------------------------------------------------
# Matplotlib theme setup
# ---------------------------------------------------------------------------


def setup_matplotlib_style(
    dpi: int = 300,
    font_size: int = FONT_SIZE_DEFAULT,
    color_up: str = COLOR_UP,
    color_down: str = COLOR_DOWN,
    color_ns: str = COLOR_NS,
) -> None:
    """Set publication-quality matplotlib defaults.

    Parameters
    ----------
    dpi : int
        Resolution for saved figures.
    font_size : int
        Base font size for all text elements.
    color_up : str
        Hex color for upregulated genes.
    color_down : str
        Hex color for downregulated genes.
    color_ns : str
        Hex color for non-significant genes.
    """
    plt.rcParams.update(
        {
            "figure.dpi": dpi,
            "font.size": font_size,
            "font.family": FONT_FAMILY,
            "axes.linewidth": 1.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "figure.facecolor": "white",
            "savefig.bbox": "tight",
            "savefig.dpi": dpi,
        }
    )
    sns.set_palette("deep")


# ---------------------------------------------------------------------------
# Plotly theme setup
# ---------------------------------------------------------------------------


def setup_plotly_theme(dark_mode: bool = False) -> None:
    """Register and activate custom Plotly templates for RNA-seq Explorer.

    Parameters
    ----------
    dark_mode : bool
        If True, use dark background theme. Default is light.
    """
    if not _PLOTLY_AVAILABLE:
        return

    light_template = go.layout.Template(
        layout=go.Layout(
            font=dict(family=FONT_FAMILY, size=FONT_SIZE_DEFAULT),
            title=dict(font=dict(size=FONT_SIZE_TITLE + 2)),
            plot_bgcolor="white",
            paper_bgcolor="white",
            xaxis=dict(
                showgrid=True,
                gridcolor="#E5E5E5",
                gridwidth=0.5,
                linecolor="#333333",
                linewidth=1,
                zeroline=True,
                zerolinecolor="#CCCCCC",
                zerolinewidth=0.8,
            ),
            yaxis=dict(
                showgrid=True,
                gridcolor="#E5E5E5",
                gridwidth=0.5,
                linecolor="#333333",
                linewidth=1,
                zeroline=True,
                zerolinecolor="#CCCCCC",
                zerolinewidth=0.8,
            ),
            colorway=CONDITION_COLORS,
            hovermode="closest",
        )
    )

    dark_template = go.layout.Template(
        layout=go.Layout(
            font=dict(family=FONT_FAMILY, size=FONT_SIZE_DEFAULT, color="#E0E0E0"),
            title=dict(font=dict(size=FONT_SIZE_TITLE + 2, color="#FFFFFF")),
            plot_bgcolor="#1E1E1E",
            paper_bgcolor="#121212",
            xaxis=dict(
                showgrid=True,
                gridcolor="#333333",
                gridwidth=0.5,
                linecolor="#666666",
                linewidth=1,
                zeroline=True,
                zerolinecolor="#444444",
                zerolinewidth=0.8,
            ),
            yaxis=dict(
                showgrid=True,
                gridcolor="#333333",
                gridwidth=0.5,
                linecolor="#666666",
                linewidth=1,
                zeroline=True,
                zerolinecolor="#444444",
                zerolinewidth=0.8,
            ),
            colorway=CONDITION_COLORS,
            hovermode="closest",
        )
    )

    pio.templates["rnaseq_light"] = light_template
    pio.templates["rnaseq_dark"] = dark_template
    pio.templates.default = "rnaseq_dark" if dark_mode else "rnaseq_light"


# ---------------------------------------------------------------------------
# Helper functions for consistent chart creation
# ---------------------------------------------------------------------------


def get_color_map(
    color_up: str = COLOR_UP,
    color_down: str = COLOR_DOWN,
    color_ns: str = COLOR_NS,
) -> dict[str, str]:
    """Return the standard Up/Down/NS color map.

    Parameters
    ----------
    color_up : str
        Hex color for upregulated genes.
    color_down : str
        Hex color for downregulated genes.
    color_ns : str
        Hex color for non-significant genes.

    Returns
    -------
    dict
        Mapping of status labels to hex colors.
    """
    return {"Up": color_up, "Down": color_down, "NS": color_ns}


def condition_color_map(labels: Sequence[str]) -> dict[str, str]:
    """Build a color map from condition labels to Okabe-Ito colors.

    Parameters
    ----------
    labels : Sequence[str]
        Condition labels.

    Returns
    -------
    dict
        Mapping of label -> hex color.
    """
    return {
        label: CONDITION_COLORS[i % len(CONDITION_COLORS)]
        for i, label in enumerate(labels)
    }


def add_count_box(
    ax,
    n_up: int,
    n_down: int,
    total: int,
    position: str = "lower left",
    up_label: str = "Up",
    down_label: str = "Down",
) -> None:
    """Add a compact count box in a data-sparse corner of a matplotlib axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to annotate.
    n_up : int
        Count of upregulated/included items.
    n_down : int
        Count of downregulated/excluded items.
    total : int
        Total significant count.
    position : str
        Corner position: 'upper left', 'upper right', 'lower left', 'lower right'.
    up_label : str
        Label for the up count.
    down_label : str
        Label for the down count.
    """
    text = f"{up_label}: {n_up:,}\n{down_label}: {n_down:,}\nTotal: {total:,}"
    loc = {
        "upper left": (0.02, 0.98),
        "upper right": (0.98, 0.98),
        "lower left": (0.02, 0.02),
        "lower right": (0.98, 0.02),
    }
    x, y = loc.get(position, (0.02, 0.02))
    ha = "left" if "left" in position else "right"
    va = "top" if "upper" in position else "bottom"
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=FONT_SIZE_ANNOTATION,
        va=va,
        ha=ha,
        fontweight="bold",
        bbox=dict(
            boxstyle="round,pad=0.3",
            facecolor="white",
            edgecolor="grey",
            alpha=0.9,
        ),
    )


def style_venn(v, n_sets: int) -> None:
    """Apply Okabe-Ito colors to a matplotlib_venn Venn diagram object.

    Parameters
    ----------
    v : VennDiagram
        The venn2 or venn3 return object.
    n_sets : int
        Number of sets (2 or 3).
    """
    if v is None:
        return
    colors = VENN_COLORS_2 if n_sets == 2 else VENN_COLORS_3
    for rid, col in colors.items():
        patch = v.get_patch_by_id(rid)
        if patch:
            patch.set_color(col)
            patch.set_alpha(0.6)


def grid_dims(n: int) -> tuple[int, int]:
    """Return (nrows, ncols) for a tight subplot grid of *n* panels.

    Parameters
    ----------
    n : int
        Number of panels.

    Returns
    -------
    tuple
        (nrows, ncols)
    """
    import math

    if n <= 1:
        return (1, 1)
    if n == 2:
        return (1, 2)
    if n == 3:
        return (1, 3)
    if n == 4:
        return (2, 2)
    if n == 5:
        return (2, 3)
    ncols = int(math.ceil(math.sqrt(n)))
    nrows = int(math.ceil(n / ncols))
    return (nrows, ncols)


def diverging_cmap():
    """Return a blue-white-orange diverging colormap (colorblind-safe).

    Returns
    -------
    LinearSegmentedColormap
        Matplotlib colormap from Okabe-Ito blue through white to orange.
    """
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "blue_white_orange", [COLOR_DOWN, "#FFFFFF", COLOR_UP]
    )
