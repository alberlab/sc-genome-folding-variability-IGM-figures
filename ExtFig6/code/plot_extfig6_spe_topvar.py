#!/usr/bin/env python3
"""Extended Data Fig. 6b (manuscript panel = (B), SD version): speckle-dominant regions and the most variable regions: for each stage model, the 10 % of 200 kb bins with the highest
SD (IQR version dropped, 2026-11) of the distance to the CT core ("top 10 %": metric >= its 90th percentile over
the 13,315 autosomal non-centromeric bins of that model) vs the bins whose dominant nuclear locale is speckle
(spe-dominant: argmax of the five exclusive locale frequencies in that model). Re-implements cell 55 of the original
notebook (nuclear_body_annotation_entropy_allab_absolute_dist-NLrenamed.ipynb, "bar_speDominant_pct_in_topCV10"),
with cisratio CV replaced by SD / IQR(dist to CT) and the occ0.6 models.
  (A) % of spe-dominant regions that are in the top 10 %  = n_spe_top / n_spe; random expectation 10 % (grey dashed line)
  (B) % of top-10 % regions that are spe-dominant         = n_spe_top / n_top; random expectation = the stage's
      fraction of spe-dominant bins (grey dashed segment on each bar)
Locale = ENV v2.
Outputs: output/figures/ExtFig6b_top10_SDdistCT_fraction_spe_dominant.* (B, manuscript), ExtFig6b_alt_spe_in_top10_SDdistCT.* (A),
         output/source_data/ExtFig6b_spe_top10.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from extfig6_common import CELLS, CELL_TITLE, LOC_COLORS, V_LABEL, load_stage, save, save_table

TOP = 0.10
D = {c: load_stage(c) for c in CELLS}
rows = []
for v in ("sd",):                                                       # SD only (2026-11)
    for c in CELLS:
        d = D[c]; thr = d[v].quantile(1 - TOP); top = (d[v] >= thr).to_numpy(); spe = (d.dom == "spe").to_numpy()
        n_spe, n_top, n_st = int(spe.sum()), int(top.sum()), int((spe & top).sum())
        rows.append(dict(metric=V_LABEL[v], cell=c, stage=CELL_TITLE[c], n_bins=len(d), threshold=thr, n_top=n_top,
                         n_spe=n_spe, n_spe_top=n_st, pct_spe_in_top=100 * n_st / n_spe, pct_top_is_spe=100 * n_st / n_top,
                         pct_spe_all=100 * n_spe / len(d)))
T = pd.DataFrame(rows); save_table(T, "ExtFig6b_spe_top10.tsv", float_format="%.4f")
print(T[["metric", "stage", "n_spe_top", "n_spe", "pct_spe_in_top", "n_top", "pct_top_is_spe", "pct_spe_all"]].round(1).to_string(index=False))


def bar(t, col, ylabel, title, name, expected):
    fig, ax = plt.subplots(figsize=(6.2, 4.6)); x = np.arange(len(t))
    ax.bar(x, t[col], width=0.6, color=LOC_COLORS["spe"], edgecolor="black", linewidth=0.6)
    top = t[col].max()
    for i, val in enumerate(t[col]):
        ax.text(i, val + 0.012 * top, f"{val:.1f}%", ha="center", va="bottom", fontsize=11)
    if np.isscalar(expected):
        h = ax.axhline(expected, color="#bdbdbd", ls="--", lw=1.4, zorder=3)
        lab = f"expected {expected:.0f}% (random)"
    else:                                             # stage-specific expectation: one dashed segment per bar
        for i, e in enumerate(expected):
            h, = ax.plot([i - 0.36, i + 0.36], [e, e], color="#bdbdbd", ls="--", lw=1.6, zorder=3)
        lab = "expected (random) = % of all bins\nthat are spe-dominant"
    ax.legend([h], [lab], loc="upper left", frameon=False, fontsize=10, handlelength=2.2)
    ax.set_xticks(x, t.stage, fontsize=12, fontweight="bold"); ax.tick_params(axis="y", labelsize=11)
    ax.set_ylim(0, top * 1.25); ax.set_ylabel(ylabel, fontsize=12); ax.set_title(title, fontsize=12, fontweight="bold")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout(); save(fig, name)


for v in ("sd",):                                                       # SD only (2026-11)
    t = T[T.metric == V_LABEL[v]]
    bar(t, "pct_spe_in_top", f"% of spe-dominant regions\nin top 10% {V_LABEL[v]}",
        f"spe-dominant regions: fraction in top 10% {V_LABEL[v]}", f"ExtFig6b_alt_spe_in_top10_{v.upper()}distCT", 100 * TOP)
    bar(t, "pct_top_is_spe", f"% of top 10% {V_LABEL[v]} regions\nthat are spe-dominant",
        f"top 10% {V_LABEL[v]} regions: fraction spe-dominant", f"ExtFig6b_top10_{v.upper()}distCT_fraction_spe_dominant",
        t.pct_spe_all.to_numpy())
