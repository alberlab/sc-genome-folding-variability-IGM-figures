#!/usr/bin/env python3
"""Figure 4j: chromosome-intermingling Shannon index across radial shells, four stage models.
Shannon index of a bead in one structure: neighbours = beads within 2.5 bead radii (centre to centre), counted by
chromosome identity (homologous copies pooled), H = -sum p ln p (0 if no neighbour). Per 200 kb bin: mean over the two
copies x 1000 structures. Bins split into 10 equal-count radial shells per model as in Fig. 4b; line = shell mean,
band = 95 % CI.
Outputs: output/figures/Fig4j_Shannon_vs_radial_decile.*, output/source_data/Fig4j_shells.tsv
"""
import pandas as pd, matplotlib.pyplot as plt
from fig4_common import (DATA, SHARED, CELLS, CELL_TITLE, CELL_COLOR, domain_bins, read_radial, radial_deciles, shell_summary, save,
                         save_table)

rows = []
for c in CELLS:
    sh = pd.read_csv(DATA / f"{c}_occ06_Shannon_haploid_200kb.tsv", sep="\t")[["chr", "start", "shannon_mean"]]
    d = domain_bins(c).merge(sh, on=["chr", "start"], validate="1:1") \
        .merge(read_radial(SHARED / "radial" / f"{c}_occ06_radial_200kb.bed"), on=["chr", "start"], validate="1:1")
    d["shell"] = radial_deciles(d)
    s = shell_summary(d, "shannon_mean"); s.insert(0, "cell", c); s["norm_rad_mean"] = d.groupby("shell").norm_rad.mean().values
    rows.append(s)
S = pd.concat(rows, ignore_index=True); save_table(S, "Fig4j_shells.tsv")

fig, ax = plt.subplots(figsize=(7.6, 5.2))
for c in CELLS:
    s = S[S.cell == c]
    ax.plot(s.shell, s["mean"], lw=2, color=CELL_COLOR[c], label=CELL_TITLE[c])
    ax.fill_between(s.shell, s["mean"] - s.ci95, s["mean"] + s.ci95, color=CELL_COLOR[c], alpha=0.18, lw=0)
ax.set_xlim(1, 10); ax.set_xticks(range(1, 11))
ax.set_xlabel("Radial position decile", fontsize=14); ax.set_ylabel("Mean Shannon index in shell", fontsize=14)
ax.text(0.0, -0.13, "Interior", transform=ax.transAxes, ha="left", va="top", color="tab:blue", fontsize=13)
ax.text(1.0, -0.13, "Periphery", transform=ax.transAxes, ha="right", va="top", color="tab:blue", fontsize=13)
ax.set_title("Shannon index radial distribution", fontsize=15)
ax.tick_params(labelsize=12); ax.legend(frameon=False, ncol=2, fontsize=12)
fig.tight_layout(); save(fig, "Fig4j_Shannon_vs_radial_decile")
print(S.pivot(index="shell", columns="cell", values="mean")[CELLS].round(3).to_string())
