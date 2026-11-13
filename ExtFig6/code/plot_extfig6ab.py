#!/usr/bin/env python3
"""Extended Data Fig. 6a (top row and bottom row of the manuscript panel): per 200 kb bin, SD(dist to CT) (top; IQR version
dropped, 2026-11) and the nuclear locale entropy (bottom) against the radial position RAD, four stage models (1 x 4), points coloured by
the bin's dominant locale in that model; Spearman correlation per panel. 13,315 autosomal non-centromeric bins.
Locale = ENV v2.
Outputs: output/figures/ExtFig6a_top_SDdistCT_vs_RAD.*, ExtFig6a_bottom_entropy_vs_RAD.*, output/source_data/ExtFig6a_stats.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import spearmanr
from extfig6_common import CELLS, CELL_TITLE, LOC_COLORS, LOC_DISP, V_LABEL, load_stage, save, save_table

D = {c: load_stage(c) for c in CELLS}
DRAW = ["rep", "act", "env", "nuc", "spe"]                          # spe (smallest group) drawn on top
LEG = ["spe", "act", "nuc", "rep", "env"]
FS = 16; rows = []
allx = np.concatenate([D[c].norm_rad.to_numpy() for c in CELLS]); x0, x1 = allx.min(), allx.max(); px = 0.03 * (x1 - x0)
for col, ylab, name, panel in [("sd", V_LABEL["sd"], "ExtFig6a_top_SDdistCT_vs_RAD", "6a"),
                               ("entropy_bit", "Nuclear locale entropy (bits)", "ExtFig6a_bottom_entropy_vs_RAD", "6b")]:
    ally = np.concatenate([D[c][col].to_numpy() for c in CELLS]); y0, y1 = np.percentile(ally, [0.1, 99.9]); py = 0.03 * (y1 - y0)
    fig, axs = plt.subplots(1, len(CELLS), figsize=(5.0 * len(CELLS) + 1.6, 5.0), sharex=True, sharey=True)
    for ax, c in zip(axs, CELLS):
        d = D[c]
        for g in DRAW:
            k = d.dom == g
            if k.any():
                ax.scatter(d.loc[k, "norm_rad"], d.loc[k, col], s=6, c=LOC_COLORS[g], alpha=0.55, linewidths=0,
                           rasterized=True, label=LOC_DISP[g])
        r, p = spearmanr(d.norm_rad, d[col])
        rows.append(dict(panel=panel, metric=ylab, cell=c, n=len(d), spearman=r, p=p))
        ax.text(0.04, 0.96, f"Spearman ρ = {r:.3f}", transform=ax.transAxes, ha="left", va="top", fontsize=FS - 2,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
        ax.set_title(CELL_TITLE[c], fontsize=FS + 2); ax.set_xlabel("Radial dist. (RAD)", fontsize=FS)
        ax.set_xlim(x0 - px, x1 + px); ax.set_ylim(y0 - py, y1 + py); ax.tick_params(labelsize=FS - 3)
    axs[0].set_ylabel(ylab, fontsize=FS)
    h, l = axs[0].get_legend_handles_labels(); order = [l.index(LOC_DISP[g]) for g in LEG if LOC_DISP[g] in l]
    fig.legend([h[i] for i in order], [l[i] for i in order], title="dominant locale", loc="center left",
               bbox_to_anchor=(0.995, 0.5), frameon=False, fontsize=FS - 2, title_fontsize=FS - 2, markerscale=3)
    fig.tight_layout(rect=(0, 0, 0.995, 1)); save(fig, name)
S = pd.DataFrame(rows); save_table(S, "ExtFig6a_stats.tsv", float_format="%.4g")
print(S.pivot(index="metric", columns="cell", values="spearman")[CELLS].round(3).to_string())
