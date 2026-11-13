#!/usr/bin/env python3
"""Extended Data Fig. 6e, solid-colour variant (2026-11-11): locale entropy vs SD of the distance to the CT core (IQR dropped),
split by dominant locale (rows spe, env, NL-act, NL-rep; columns = stages). Every bin light grey; the row's dominant-locale
bins in the locale colour, opaque, without locale-strength shading and without colour bars. Panel box = n and Spearman r.
Dominant locale taken in each stage; 13,315 autosomal non-centromeric bins; locale = ENV v2; stages ES, DE, HB, HLC.
Outputs: output/figures/ExtFig6e_entropy_vs_SD_by_dominant_locale_solid.*, output/source_data/ExtFig6e_solid_stats.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import spearmanr
from extfig6_common import CELLS, CELL_TITLE, CELL_COLOR, LOC_COLORS, LOC_DISP, V_LABEL, load_stage, save, save_table

D = {c: load_stage(c) for c in CELLS}
ROWS = ["spe", "env", "act", "rep"]
GREY, CAP, AMIN = "#d9d9d9", 4.0, 0.08
FS_T, FS_L, FS_K, FS_B = 22, 17, 13, 13
PW, PH, CBW = 4.6, 4.3, 1.2


def alpha_cmap(rgb, n=256):
    a = np.clip(np.linspace(0, 1, n), AMIN, 1)
    return mcolors.ListedColormap(np.clip((1 - a)[:, None] + a[:, None] * np.asarray(rgb), 0, 1))


rows = []
for v in ("sd",):                                                       # SD only (2026-11)
    xa = np.concatenate([D[c][v] for c in CELLS]); ya = np.concatenate([D[c].entropy_bit for c in CELLS])
    (xlo, xhi), (ylo, yhi) = np.quantile(xa, [0.001, 0.999]), np.quantile(ya, [0.001, 0.999])
    XL = (xlo - 0.03 * (xhi - xlo), xhi + 0.03 * (xhi - xlo)); YL = (ylo - 0.03 * (yhi - ylo), yhi + 0.03 * (yhi - ylo))
    fig, axs = plt.subplots(len(ROWS), len(CELLS), figsize=(PW * len(CELLS), PH * len(ROWS)), sharex=True, sharey=True,
                            constrained_layout=True, squeeze=False)
    for ri, key in enumerate(ROWS):
        rgb = np.array(mcolors.to_rgb(LOC_COLORS[key])); name = LOC_DISP[key]
        for ci, c in enumerate(CELLS):
            ax = axs[ri, ci]; d = D[c]; hit = (d.dom == key).to_numpy()
            x, y, st = d[v].to_numpy(), d.entropy_bit.to_numpy(), d.strength.to_numpy()
            ax.scatter(x[~hit], y[~hit], s=5, c=GREY, alpha=0.25, edgecolors="none", rasterized=True)
            ax.scatter(x[hit], y[hit], s=8, c=[rgb], edgecolors="none", rasterized=True)   # opaque, no strength shading
            rho = spearmanr(x[hit], y[hit])[0]
            ax.text(0.03, 0.97, f"{name}: {int(hit.sum())}\nSpearman r = {rho:+.2f}", transform=ax.transAxes, ha="left", va="top",
                    fontsize=FS_B, bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
            ax.set_xlim(XL); ax.set_ylim(YL); ax.tick_params(labelsize=FS_K)
            if ri == 0:
                ax.set_title(CELL_TITLE[c], fontsize=FS_T, fontweight="bold", color=CELL_COLOR[c])
            if ci == 0:
                ax.set_ylabel(f"{name}\nlocale entropy (bit)", fontsize=FS_L, color=LOC_COLORS[key], fontweight="bold")
            if ri == len(ROWS) - 1:
                ax.set_xlabel(V_LABEL[v], fontsize=FS_L)
            rows.append(dict(metric=V_LABEL[v], locale=name, cell=c, n=int(hit.sum()), spearman=rho,
                             median_strength=float(np.median(st[hit])), median_entropy=float(np.median(y[hit])),
                             median_metric=float(np.median(x[hit]))))
    from matplotlib.lines import Line2D                                # legend: colour = dominant locale (2026-11-11)
    hd = [Line2D([], [], ls="", marker="o", ms=11, mfc=LOC_COLORS[k], mec="none", label=f"{LOC_DISP[k]}-dominant") for k in ROWS]
    hd.append(Line2D([], [], ls="", marker="o", ms=11, mfc=GREY, mec="none", label="other bins"))
    fig.legend(handles=hd, loc="outside lower center", ncol=len(hd), frameon=False, fontsize=FS_L,
               title="Dominant locale of the bin in that stage", title_fontsize=FS_L)
    save(fig, f"ExtFig6e_entropy_vs_{v.upper()}_by_dominant_locale_solid", dpi=200)
S = pd.DataFrame(rows); save_table(S, "ExtFig6e_solid_stats.tsv", float_format="%.4g")
print(S.pivot_table(index=["metric", "locale"], columns="cell", values="spearman")[CELLS].round(2).to_string())
