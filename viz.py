"""Shared plot style: one palette, thin marks, recessive axes.

Roles, not raw colours, are used in the notebook:
  OBS    - observed data          (categorical slot 1, blue)
  MODEL  - Markov model / simulated (slot 2, orange)
  THEORY - closed-form theory     (slot 3, aqua)
  STATE  - Dry / Light / Heavy, an ordered light->dark blue ramp
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

OBS = "#2a78d6"
MODEL = "#eb6834"
THEORY = "#1baf7a"
EXTRA = "#4a3aa7"
STATE = ["#86b6ef", "#2a78d6", "#104281"]   # Dry, Light, Heavy

INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"

SEQ = LinearSegmentedColormap.from_list(
    "seq_blue", ["#f4f8fe", "#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"])


def use():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
        "font.size": 10,
        "text.color": INK,
        "axes.labelcolor": INK_2,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2,
        "lines.linewidth": 2,
        "lines.markersize": 6,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "patch.linewidth": 0,
    })


def save(fig, name):
    fig.savefig(f"figures/{name}.png")
    plt.show()
