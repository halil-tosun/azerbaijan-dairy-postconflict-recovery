"""
08a_make_figure1_framework.py
===============================
Generates Figure 1 (the integrated analytical-framework schematic) as a
high-resolution (600 DPI) diagram. Unlike Figures 2-6, this figure does
not depend on any output/*.csv file -- it is a conceptual diagram of the
four-stage empirical framework described in Section 2.3 of the
manuscript, rendered programmatically for full reproducibility (it was
previously a hand-drawn, uncoded image).

Run standalone, or via run_all.py alongside 08_make_figures.py.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

FIG_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

DPI = 600

NAVY = "#1b3a4b"
TEAL = "#2d6e7e"
SLATE = "#5c7a89"
SAND = "#c9a15e"
LIGHT = "#f4f6f7"
INK = "#1f2933"


def _box(ax, x, y, w, h, text, facecolor, edgecolor=INK, fontsize=9.3,
         fontweight="normal", textcolor=INK, lw=1.3, radius=0.12):
    b = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.02,rounding_size={radius}",
                        linewidth=lw, edgecolor=edgecolor, facecolor=facecolor, zorder=2)
    ax.add_patch(b)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
            fontweight=fontweight, color=textcolor, zorder=3, linespacing=1.4,
            family="DejaVu Sans")


def _arrow(ax, x, y1, y2, color=SLATE, style="-|>"):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle=style, mutation_scale=16,
                                  linewidth=1.6, color=color, zorder=1))


def _stage_label(ax, x, y, num, text, color):
    ax.text(x, y, f"STAGE {num}  \u2013  {text}", fontsize=7.6, fontweight="bold", color=color,
            ha="left", va="center", family="DejaVu Sans")


def figure1_framework():
    fig, ax = plt.subplots(figsize=(6.4, 8.8), dpi=DPI)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 15.9)
    ax.axis("off")

    ax.text(5, 15.55, "Analytical Framework of the Empirical Analysis", ha="center", va="center",
            fontsize=12.4, fontweight="bold", color=NAVY, family="DejaVu Sans")
    ax.plot([0.9, 9.1], [15.15, 15.15], color=NAVY, linewidth=1.4)

    _box(ax, 1.0, 13.45, 8.0, 1.15,
         "District-level Panel Data\nAzerbaijan dairy sector, 2000\u20132024\n"
         "67 districts \u00b7 1,523 district-year observations",
         facecolor=LIGHT, edgecolor=NAVY, fontsize=9.0, fontweight="bold", textcolor=NAVY, lw=1.6)
    _arrow(ax, 5, 13.45, 12.75)

    _stage_label(ax, 1.0, 12.55, 1, "Technical Efficiency Estimation", NAVY)
    _box(ax, 1.0, 10.75, 3.85, 1.55,
         "Bootstrap DEA\n(VRS / CRS)\nSimar\u2013Wilson\nbias correction\n200 replications",
         facecolor="#dce8ea", edgecolor=TEAL, fontsize=8.5)
    _box(ax, 5.15, 10.75, 3.85, 1.55,
         "Stochastic Frontier\nAnalysis (SFA)\nTranslog, ML\nJondrow / Battese\u2013Coelli TE",
         facecolor="#dce8ea", edgecolor=TEAL, fontsize=8.5)

    ax.add_patch(FancyArrowPatch((2.925, 10.75), (5, 9.85), arrowstyle="-|>", mutation_scale=15,
                                  linewidth=1.5, color=SLATE, connectionstyle="arc3,rad=-0.18"))
    ax.add_patch(FancyArrowPatch((7.075, 10.75), (5, 9.85), arrowstyle="-|>", mutation_scale=15,
                                  linewidth=1.5, color=SLATE, connectionstyle="arc3,rad=0.18"))

    _stage_label(ax, 1.0, 9.75, 2, "Productivity Change Analysis", TEAL)
    _box(ax, 1.0, 8.15, 8.0, 1.35,
         "Bootstrap Malmquist Productivity Index\n"
         "TFP change = Technical Efficiency Change (TEC) \u00d7 Technological Change (TC)\n"
         "13 economic regions \u00b7 300 bootstrap replications",
         facecolor="#eef2df", edgecolor=SAND, fontsize=8.7, fontweight="bold", textcolor=INK)
    _arrow(ax, 5, 8.15, 7.35)

    _stage_label(ax, 1.0, 7.25, 3, "Determinants of Technical Efficiency", SAND)
    _box(ax, 1.0, 5.75, 8.0, 1.15,
         "Simar\u2013Wilson (2007) Bootstrap Truncated Regression\n"
         "Cattle stock \u00b7 fodder area \u00b7 milk cost \u00b7 milk profitability",
         facecolor="#f6ecd9", edgecolor=SAND, fontsize=8.7)
    _arrow(ax, 5, 5.75, 4.95)

    _stage_label(ax, 1.0, 4.85, 4, "Post-Conflict Policy Evaluation", NAVY)
    _box(ax, 1.0, 3.35, 8.0, 1.15,
         "Two-Way Fixed-Effects Event Study & Difference-in-Differences\n"
         "Treated districts (n = 9) vs. comparison districts (n = 58)",
         facecolor="#dde3ea", edgecolor=NAVY, fontsize=8.7)
    _arrow(ax, 5, 3.35, 2.45)

    _box(ax, 0.8, 0.75, 8.4, 1.5,
         "Integrated Assessment of Post-Conflict\nAgricultural Systems Recovery",
         facecolor=NAVY, edgecolor=NAVY, fontsize=10.8, fontweight="bold", textcolor="white", lw=0)

    plt.tight_layout(pad=0.3)
    plt.savefig(os.path.join(FIG_DIR, 'Figure1_Analytical_Framework.png'),
                dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print('Figure 1 (analytical framework) saved at', DPI, 'DPI')


def run():
    figure1_framework()


if __name__ == '__main__':
    run()
