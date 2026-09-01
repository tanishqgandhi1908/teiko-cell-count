"""Colours and matplotlib defaults shared by every figure.

Blue/orange are colourblind-safe and clear 3:1 contrast on the chart surface.
"""

from __future__ import annotations

RESPONSE_COLORS = {
    "Responder": "#2a78d6",      # blue
    "Non-responder": "#eb6834",  # orange
}
RESPONSE_ORDER = ["Responder", "Non-responder"]

SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
TEXT_MUTED = "#7a7a75"
GRID = "#e6e5e1"

FONT_STACK = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]


def apply_matplotlib_theme() -> None:
    """Light background, faint gridlines, no top/right spines."""
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.family": "sans-serif",
            "font.sans-serif": FONT_STACK,
            "text.color": TEXT_PRIMARY,
            "axes.labelcolor": TEXT_SECONDARY,
            "axes.edgecolor": GRID,
            "axes.linewidth": 1.0,
            "axes.grid": True,
            "axes.grid.axis": "y",
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.color": TEXT_SECONDARY,
            "ytick.color": TEXT_SECONDARY,
            "xtick.bottom": False,
            "ytick.left": False,
            "legend.frameon": False,
            "figure.dpi": 110,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
        }
    )
