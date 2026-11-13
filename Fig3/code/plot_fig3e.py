#!/usr/bin/env python3
"""Figure 3e (IGM row): chr2, Mean(dist to CT) vs the stage's own bulk Eig1 at 1 Mb, one panel per stage model (1 x 4),
coloured by SD(dist to CT) (IQR version dropped, 2026-11). Colour range shared by the four stages
(2nd-98th percentile of the pooled chr2 values). No regression line. Spearman correlation over all plotted chr2 bins.
Outputs: output/figures/Fig3e_IGM_Mean_vs_Eig1_chr2_colorbySD.*, output/source_data/Fig3e_chr2_1Mb.tsv, Fig3e_stats.tsv
"""
import numpy as np, pandas as pd, matplotlib, matplotlib.ticker, matplotlib.pyplot as plt
from fig3_common import CT_SCALE, CELLS, CELL_TITLE, BWR, V_LABEL, load_200kb, load_e1_1mb, to_1mb, spear_text, prange, save, save_table

FS, TS = 20, 18
tabs = []
for c in CELLS:
    m = to_1mb(load_200kb(c), ["mean", "sd", "iqr"]).merge(load_e1_1mb(c), on=["chr", "start"], how="left", validate="1:1")
    m = m[(m.chr == "chr2") & np.isfinite(m.E1_1mb)]
    tabs.append(m.assign(cell=c))
T = pd.concat(tabs, ignore_index=True)
save_table(T[["cell", "chr", "start", "end", "n_200kb", "E1_1mb", "mean", "sd", "iqr"]], "Fig3e_chr2_1Mb.tsv")

stats = []
for v in ("sd",):                                                       # SD only (2026-11)
    crange = prange(T[v])
    fig, axs = plt.subplots(1, 4, figsize=(24, 5.6), sharex=True, sharey=True, constrained_layout=True)
    for ax, c in zip(axs, CELLS):
        t = T[T.cell == c]
        sc = ax.scatter(t.E1_1mb, t["mean"], c=t[v], cmap=BWR, vmin=crange[0], vmax=crange[1], s=40,
                        edgecolor="white", linewidths=0.3)
        r, p = spear_text(ax, t.E1_1mb, t["mean"], loc="ll", fontsize=FS - 2)
        ax.set_title(CELL_TITLE[c], fontsize=FS + 2, fontweight="bold")
        ax.set_xlim(-1.2, 1.4); ax.set_ylim(-1.5 * CT_SCALE, 0.5 * CT_SCALE)          # shows every chr2 bin of all four stages
        ax.set_xlabel("Eig1", fontsize=FS); ax.tick_params(labelsize=TS)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        stats.append(dict(color_by=v.upper(), cell=c, n=len(t), spearman=r, p=p))
    axs[0].set_ylabel(V_LABEL["mean"], fontsize=FS)
    cb = fig.colorbar(sc, ax=axs, fraction=0.025, pad=0.02); cb.set_label(V_LABEL[v], fontsize=FS); cb.ax.tick_params(labelsize=TS)
    save(fig, f"Fig3e_IGM_Mean_vs_Eig1_chr2_colorby{v.upper()}")
S = pd.DataFrame(stats); save_table(S, "Fig3e_stats.tsv", float_format="%.4g")
print(S.to_string(index=False))
