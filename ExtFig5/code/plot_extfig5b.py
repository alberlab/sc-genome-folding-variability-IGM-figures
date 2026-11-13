#!/usr/bin/env python3
"""Extended Data Fig. 5b (IGM model part): chr2 territory volume and compactness in the four stage models (each structure x
each chr2 copy: n = 2000 per stage). Same definitions and style as Fig. 4f (chr11), chr2 drawn in dark red.
  CT volume ratio  = convex-hull volume of the chr2 CT-core beads (top 30 % strict-cis beads at 2.5 r, the core used for
                     dist to CT) / nuclear volume (4/3 pi a b c of the model's ellipsoid)
  Compactness      = bead radius / Rg, Rg = radius of gyration of all beads of that chr2 copy
Outputs: output/figures/ExtFig5b_chr2_{CT_volume_ratio,compactness}.*, output/source_data/ExtFig5b_summary.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from extfig5_common import DATA, CELLS, CELL_TITLE, save, save_table

COLOR = "darkred"
V = pd.concat([pd.read_csv(DATA / f"{c}_occ06_chr2_ct_volfrac_compactness.tsv", sep="\t").assign(cell=c) for c in CELLS],
              ignore_index=True)
for c in CELLS:
    assert (V.cell == c).sum() == 2000, c
rows = []
for col, lab in (("ct_volume_fraction", "CT volume ratio"), ("compactness_r_over_rg", "Compactness")):
    for c in CELLS:
        v = V.loc[V.cell == c, col].dropna().to_numpy()
        q1, q3 = np.percentile(v, [25, 75])
        rows.append(dict(metric=lab, cell=c, n=len(v), median=np.median(v), q25=q1, q75=q3, mean=v.mean()))
save_table(pd.DataFrame(rows), "ExtFig5b_summary.tsv", float_format="%.5f")
print(pd.DataFrame(rows).pivot(index="metric", columns="cell", values="median")[CELLS].round(4).to_string())


def violin(col, title, ylabel, name):
    data = [V.loc[V.cell == c, col].dropna().to_numpy() for c in CELLS]
    fig, ax = plt.subplots(figsize=(3.2, 4.8))
    parts = ax.violinplot(data, positions=range(1, len(CELLS) + 1), widths=0.72, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(COLOR); b.set_edgecolor("0.35"); b.set_linewidth(1.0); b.set_alpha(0.95)
    for k, v in enumerate(data, start=1):
        ax.plot([k - 0.16, k + 0.16], [np.median(v)] * 2, color="0.75", lw=1.2)
    lo, hi = ax.get_ylim()
    for k, v in enumerate(data, start=1):
        ax.text(k, lo + 1.04 * (hi - lo), f"n={len(v)}", rotation=90, ha="center", va="bottom", fontsize=16)
    ax.set_ylim(lo, lo + 1.42 * (hi - lo))                        # headroom for the n labels above the violins
    ax.set_xticks(range(1, len(CELLS) + 1), [CELL_TITLE[c] for c in CELLS], rotation=90, fontsize=22)
    ax.set_ylabel(ylabel, fontsize=22); ax.set_title(title, fontsize=22); ax.tick_params(axis="y", labelsize=18)
    save(fig, name)


violin("ct_volume_fraction", "CT volume ratio", "CT vol. / nuclei vol", "ExtFig5b_chr2_CT_volume_ratio")
violin("compactness_r_over_rg", "Compactness score", "chr Rg$^{-1}$ (normalized)", "ExtFig5b_chr2_compactness")
