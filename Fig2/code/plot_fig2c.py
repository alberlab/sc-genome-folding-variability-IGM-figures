#!/usr/bin/env python3
"""Figure 2c: dist to CT vs bulk compartment Eig1 at 1 Mb; variability metric V = SD(dist to CT) (IQR version dropped).
  left    Mean(dist to CT) vs Eig1, chr2 only, colour = V
  middle  V vs Eig1, genome-wide (x axis limited to [-1.5, 1.2]: hides 3 bins with E1 < -1.5; statistics use all bins)
  right   V vs Mean(dist to CT), genome-wide, colour = Eig1
Outputs: output/figures/Fig2c_{left,middle,right}_*, output/source_data/Fig2c_h1_1Mb.tsv, Fig2c_stats.tsv
Model values are the 1 Mb mean of the 200 kb bin statistics. Spearman correlation on all plotted bins.
"""
import numpy as np, pandas as pd, matplotlib, matplotlib.ticker, matplotlib.pyplot as plt
from fig2_common import CT_SCALE, BWR, V_LABEL, load_200kb, load_e1_1mb, to_1mb, spear_text, prange, lim, save, save_table

d = load_200kb()
m = to_1mb(d, ["mean", "sd", "iqr"]).merge(load_e1_1mb(), on=["chr", "start"], how="left", validate="1:1")
m = m[np.isfinite(m.E1_1mb)].reset_index(drop=True)
save_table(m[["chr", "start", "end", "n_200kb", "E1_1mb", "mean", "iqr", "sd"]], "Fig2c_h1_1Mb.tsv")
print(f"1 Mb bins with model + E1: {len(m)} (chr2: {(m.chr == 'chr2').sum()})")

E1_C = (-1, 1)                    # colour range of Eig1
FS, TS = 20, 18                   # label / tick font sizes
stats = []


def panel(ax, x, y, c=None, crange=None, clabel=None, fig=None, xlab="", ylab="", xl=None, yl=None, loc="ul"):
    if c is None:
        ax.scatter(x, y, c="lightgray", s=40, edgecolor="white", linewidths=0.3)
    else:
        sc = ax.scatter(x, y, c=c, cmap=BWR, vmin=crange[0], vmax=crange[1], s=40, edgecolor="white", linewidths=0.3)
        cb = fig.colorbar(sc, ax=ax); cb.set_label(clabel, fontsize=FS); cb.ax.tick_params(labelsize=TS)
    r, p = spear_text(ax, x, y, loc=loc, fontsize=FS)
    ax.set_xlabel(xlab, fontsize=FS); ax.set_ylabel(ylab, fontsize=FS); ax.tick_params(labelsize=TS)
    if xl: ax.set_xlim(xl)
    if yl: ax.set_ylim(yl)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    return r, p


NEW_NAME = {"p1_Mean_vs_Eig1_chr2_colorbySD": "Fig2c_left_Mean_vs_Eig1_chr2_colorbySD",   # manuscript 2026-11 numbering
            "p2_SD_vs_Eig1": "Fig2c_middle_SD_vs_Eig1", "p3_SD_vs_Mean_colorbyEig1": "Fig2c_right_SD_vs_Mean_colorbyEig1"}
for v in ("sd",):                                                       # SD only (2026-11)
    V = V_LABEL[v]; vr = prange(m[v]); g = m[m.chr == "chr2"]; yl2 = lim(m[v])
    specs = [
        ("p1_Mean_vs_Eig1_chr2_colorby" + v.upper(), (6, 5),
         dict(x=g.E1_1mb, y=g["mean"], c=g[v], crange=vr, clabel=V, xlab="Eig1", ylab=V_LABEL["mean"],
              xl=(-1, 1), yl=(-1.5 * CT_SCALE, 0.5 * CT_SCALE), loc="ml")),
        (f"p2_{v.upper()}_vs_Eig1", (5, 5),
         dict(x=m.E1_1mb, y=m[v], xlab="Eig1", ylab=V, xl=(-1.5, 1.2), yl=(yl2[0], yl2[1] + 0.25 * (yl2[1] - yl2[0])))),
        (f"p3_{v.upper()}_vs_Mean_colorbyEig1", (6, 5),
         dict(x=m["mean"], y=m[v], c=m.E1_1mb, crange=E1_C, clabel="Eig1", xlab=V_LABEL["mean"], ylab=V,
              xl=lim(m["mean"]), yl=lim(m[v]), loc="ur")),
    ]
    for name, fs, kw in specs:
        fig, ax = plt.subplots(figsize=fs); r, p = panel(ax, fig=fig, **kw); fig.tight_layout()
        save(fig, NEW_NAME[name])
        stats.append(dict(metric=v.upper(), panel=name, n=len(kw["x"]), spearman=r, p=p))

S = pd.DataFrame(stats); save_table(S, "Fig2c_stats.tsv", float_format="%.4g")
print(S.to_string(index=False))
