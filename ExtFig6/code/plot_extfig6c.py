#!/usr/bin/env python3
"""Extended Data Figure 6c: nuclear locale entropy vs SD(dist to CT), 2 x 4 — top row all autosomes merged, bottom row
chr20; columns ES DE HB HLC.
Points = 200 kb autosomal non-centromeric bins coloured by dominant locale (argmax of the five exclusive frequencies),
Spearman rho per panel (in the panel title). Panels in one row share the x range (pooled 0.1-99.9 % of that row's SD
values + 3 %); y = 0.45-2.35 bits everywhere. Locale = ENV v2 (env seeds = outermost 20 % of all beads + 2.5 r).
Split from the earlier combined script plot_new_panels.py (Fig5d / old ExtFig2d / ExtFig6c); drawing code unchanged.
Data: ../../data (shared; read through ../../Fig5/code/fig5_common.py).
Outputs: ../output/figures/ExtFig6c_entropy_vs_SDdistCT_merged_chr20.*, ../output/source_data/ExtFig6c_spearman.tsv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "Fig5" / "code"))       # ../../Fig5/code/fig5_common.py (reads ../../data)
import fig5_common as fc  # noqa: E402

OUT = HERE.parent / fc.OUT.name; FIG, SRC = OUT / "figures", OUT / "source_data"   # this folder's output/ or output_<LOCALE_VERSION>/
FIG.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)
D = {c: fc.load_locale(c).merge(fc.read_ct(fc.CT_DIR / f"{c}_occ06_distCT_top30_2p5r_stats.tsv"), on=["chr", "start"],
                                validate="1:1") for c in fc.CELLS}
for c in fc.CELLS:
    assert len(D[c]) == 13315
DRAW = ["rep", "act", "env", "nuc", "spe"]; LEG = ["spe", "act", "nuc", "rep", "env"]
YL = (0.45, 2.35); FS = 16; rows = []


def panel(ax, d, dot, tag):
    for g in DRAW:
        k = d.dom5 == g
        if k.any():
            ax.scatter(d.loc[k, "sd"], d.loc[k, "entropy_bit"], s=dot, c=fc.LOC_COLORS[g], alpha=0.55, linewidths=0,
                       rasterized=True, label=fc.LOC_DISP[g])
    r, p = spearmanr(d["sd"], d.entropy_bit)
    ax.set_ylim(*YL); ax.tick_params(labelsize=FS - 3)          # rho goes into the panel title (no box over points)
    rows.append(dict(figure=tag[0], cell=tag[1], subset=tag[2], n=len(d), spearman=r, p=p))
    return r


def xlim(frames):
    x = np.concatenate([f["sd"].to_numpy() for f in frames]); x0, x1 = np.percentile(x, [0.1, 99.9]); px = 0.03 * (x1 - x0)
    return x0 - px, x1 + px


def legend(fig, ax, x=0.995):
    h, l = ax.get_legend_handles_labels(); order = [l.index(fc.LOC_DISP[g]) for g in LEG if fc.LOC_DISP[g] in l]
    fig.legend([h[i] for i in order], [l[i] for i in order], title="dominant locale", loc="center left",
               bbox_to_anchor=(x, 0.5), frameon=False, fontsize=FS - 2, title_fontsize=FS - 2, markerscale=2)


def save(fig, name):
    for e in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig); print("wrote", name)


def grid(name, layout, dot_of):
    """layout = list of rows; each row = (row label, subset key, [cells]); subset key 'all' or a chromosome."""
    sub = lambda c, s: D[c] if s == "all" else D[c][D[c].chr == s]
    nr, nc = len(layout), max(len(cs) for _, _, cs in layout)
    fig, axs = plt.subplots(nr, nc, figsize=(5.0 * nc + 1.6, 5.0 * nr), sharey=True, squeeze=False)
    for i, (rlab, s, cs) in enumerate(layout):
        XL = xlim([sub(c, s) for c in cs])                              # x range shared within a row
        for j, c in enumerate(cs):
            ax = axs[i, j]; r = panel(ax, sub(c, s), dot_of(s), (name, fc.CELL_TITLE[c], s)); ax.set_xlim(*XL)
            ax.set_title(f"{fc.CELL_TITLE[c]}, {rlab}  (Spearman ρ = {r:.3f})", fontsize=FS + 1)
            if i == nr - 1:
                ax.set_xlabel(f"SD(dist to CT{fc.CT_UNIT})", fontsize=FS)
        axs[i, 0].set_ylabel("Nuclear locale entropy (bits)", fontsize=FS)
    legend(fig, axs[0, 0]); fig.tight_layout(rect=(0, 0, 0.995, 1)); save(fig, name)


dot = lambda s: 6 if s == "all" else 22
grid("ExtFig6c_entropy_vs_SDdistCT_merged_chr20", [("all autosomes", "all", fc.CELLS), ("chr20", "chr20", fc.CELLS)], dot)
S = pd.DataFrame(rows); S.to_csv(SRC / "ExtFig6c_spearman.tsv", sep="\t", index=False, float_format="%.4g")
print(S.to_string(index=False))
