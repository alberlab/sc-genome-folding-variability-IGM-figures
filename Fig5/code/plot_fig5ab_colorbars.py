#!/usr/bin/env python3
"""Colour bars for the ChimeraX renders of Figure 5a / 5b (same palettes and ranges as code/chimerax/fig5ab_scene.py):
  5a  blue-white-red over 0.2-1.5, SD(dist to CT) (IQR bar dropped, 2026-11)
  5b  blue-yellow over 0.8-2.3, nuclear locale entropy
Horizontal bars with the end values, transparent background.
Outputs: output/figures/Fig5a/colorbar_SD.{png,pdf}, output/figures/Fig5b/colorbar_entropy.{png,pdf}
"""
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from matplotlib.cm import ScalarMappable
from fig5_common import FIG_DIR, DISTCT_NORM, DEFAULT_DISTCT
if DISTCT_NORM != DEFAULT_DISTCT:      # colour bars of the e2c ChimeraX renders: not redrawn for r / Req
    print(f"skipped (DISTCT_NORM={DISTCT_NORM}: Fig5a/5b renders use the e2c values)"); raise SystemExit(0)

BARS = [("Fig5a", "colorbar_SD", ["blue", "white", "red"], 0.2, 1.5, "SD(dist to CT)"),
        ("Fig5b", "colorbar_entropy", ["blue", "yellow"], 0.8, 2.3,
         "Nuclear locale entropy (bits)")]
for panel, name, colors, lo, hi, label in BARS:
    cmap = mcol.LinearSegmentedColormap.from_list(name, colors)
    fig, ax = plt.subplots(figsize=(2.6, 0.9))
    cb = fig.colorbar(ScalarMappable(norm=mcol.Normalize(lo, hi), cmap=cmap), cax=ax, orientation="horizontal")
    cb.set_ticks([lo, hi]); cb.ax.tick_params(labelsize=11, length=3); cb.outline.set_linewidth(0.6)
    cb.set_label(label, fontsize=11, labelpad=2)
    fig.subplots_adjust(left=0.08, right=0.92, bottom=0.55, top=0.85)
    (FIG_DIR / panel).mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / panel / f"{name}.{ext}", dpi=600, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote output/figures/{panel}/{name}.{{png,pdf}}")
