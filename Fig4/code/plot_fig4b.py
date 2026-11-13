#!/usr/bin/env python3
"""Figure 4b: variability of the distance to the chromosome-territory core across radial shells, four stage models.
For each model the 13,315 autosomal non-centromeric 200 kb bins are split into 10 equal-count shells by their mean
radial position (ellipsoid-normalised; 1 = interior, 10 = periphery); the line is the mean IQR (or SD) of dist to CT in
each shell, the band its 95 % confidence interval (mean +- 1.96 SE over bins). SD only (IQR dropped, 2026-11).
Outputs: output/figures/Fig4b_SDdistCT_vs_radial_decile.*, output/source_data/Fig4b_shells.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from fig4_common import (DATA, SHARED, CT_DIR, CELLS, CELL_TITLE, CELL_COLOR, V_LABEL, domain_bins, read_ct, read_radial, radial_deciles,
                         shell_summary, save, save_table)

rows = []
D = {}
for c in CELLS:
    d = domain_bins(c).merge(read_ct(CT_DIR / f"{c}_occ06_distCT_top30_2p5r_stats.tsv"), on=["chr", "start"], validate="1:1") \
        .merge(read_radial(SHARED / "radial" / f"{c}_occ06_radial_200kb.bed"), on=["chr", "start"], validate="1:1")
    d["shell"] = radial_deciles(d); D[c] = d
    for v in ("sd",):                                 # SD only (2026-11)
        s = shell_summary(d, v); s.insert(0, "metric", V_LABEL[v]); s.insert(0, "cell", c)
        s["norm_rad_mean"] = d.groupby("shell").norm_rad.mean().values
        rows.append(s)
S = pd.concat(rows, ignore_index=True); save_table(S, "Fig4b_shells.tsv")

for v in ("sd",):                                 # SD only (2026-11)
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    for c in CELLS:
        s = S[(S.cell == c) & (S.metric == V_LABEL[v])]
        ax.plot(s.shell, s["mean"], lw=2, color=CELL_COLOR[c], label=CELL_TITLE[c])
        ax.fill_between(s.shell, s["mean"] - s.ci95, s["mean"] + s.ci95, color=CELL_COLOR[c], alpha=0.18, lw=0)
    ax.set_xlim(1, 10); ax.set_xticks(range(1, 11))
    ax.set_xlabel("Radial position decile", fontsize=14); ax.set_ylabel(f"Mean {V_LABEL[v]} in shell", fontsize=14)
    ax.text(0.0, -0.13, "Interior", transform=ax.transAxes, ha="left", va="top", color="tab:blue", fontsize=13)
    ax.text(1.0, -0.13, "Periphery", transform=ax.transAxes, ha="right", va="top", color="tab:blue", fontsize=13)
    ax.set_title(f"{V_LABEL[v]} radial distribution", fontsize=15)
    ax.tick_params(labelsize=12); ax.legend(frameon=False, ncol=2, fontsize=12)
    fig.tight_layout(); save(fig, f"Fig4b_{v.upper()}distCT_vs_radial_decile")
print(S.pivot_table(index=["metric", "shell"], columns="cell", values="mean")[CELLS].round(3).to_string())
