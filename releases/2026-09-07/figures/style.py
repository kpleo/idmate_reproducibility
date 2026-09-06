"""Shared print style for the focused IDMate manuscript figures."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

MM = 1 / 25.4
WIDTH_MM = 183
COLORS = {
    "ink": "#111111", "blue": "#0072B2", "orange": "#D55E00",
    "green": "#009E73", "reference": "#626262", "boundary": "#A4272C",
    "pale_blue": "#EAF3F9", "pale_green": "#EAF5F0", "pale_orange": "#FCF0E8",
}


def use_style():
    for family in ("Arial", "Times New Roman"):
        try:
            font_manager.findfont(family, fallback_to_default=False)
        except ValueError as error:
            raise RuntimeError(f"Required figure font unavailable: {family}. Install it locally before rendering.") from error
    sans, serif = "Arial", "Times New Roman"
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [sans],
        "font.size": 7.6, "font.weight": "bold", "text.color": COLORS["ink"],
        "axes.labelsize": 8.0, "axes.labelweight": "bold",
        "axes.labelcolor": COLORS["ink"], "axes.edgecolor": COLORS["ink"],
        "axes.linewidth": 0.75, "axes.spines.top": True, "axes.spines.right": True,
        "axes.grid": False, "axes.unicode_minus": True,
        "xtick.labelsize": 7.0, "ytick.labelsize": 7.0,
        "xtick.color": COLORS["ink"], "ytick.color": COLORS["ink"],
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "xtick.major.size": 2.6, "ytick.major.size": 2.6,
        "legend.fontsize": 7.2, "legend.frameon": False,
        "lines.linewidth": 1.1, "lines.markersize": 4.0,
        "mathtext.fontset": "custom", "mathtext.rm": serif,
        "mathtext.it": serif + ":italic", "mathtext.bf": serif + ":bold",
        "mathtext.cal": "STIXGeneral", "mathtext.default": "it",
        "text.usetex": False, "svg.fonttype": "none", "pdf.fonttype": 42,
        "ps.fonttype": 42, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def panel_label(ax, letter, x=-0.12, y=1.035):
    return ax.text(x, y, letter, transform=ax.transAxes, fontsize=10,
                   fontweight="bold", ha="left", va="bottom", color=COLORS["ink"])


def save_figure(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    fig.savefig(stem.with_suffix(".pdf"), metadata={"Creator": "IDMate scientific plotting"})
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".png"), dpi=400)
    return stem.with_suffix(".pdf")
