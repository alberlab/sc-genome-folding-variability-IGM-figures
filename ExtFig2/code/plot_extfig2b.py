#!/usr/bin/env python3
"""Extended Data Fig. 2b: H1 model variability of the distance to the CT core vs radial position (RAD), 200 kb bins,
density scatter, SD(dist to CT) only (left panel; the IQR panel of the earlier version is dropped). RAD = mean normalised radial position of the bin (0 = nuclear centre,
1 = envelope; ellipsoid-normalised). 13,315 autosomal non-centromeric bins. Colour = per-point KDE density (panel maximum
= 1); axes span the 0.5-99.5 % range of the plotted values (points outside counted in the source data). Spearman r.
Right panel (SD(RAD) vs SD(dist to CT)): plot_extfig2b_right.py.
Outputs: output/figures/ExtFig2b_left_SDdistCT_vs_RAD.*, output/source_data/ExtFig2b_stats.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from extfig2_common import V_LABEL, load_h1_200kb, density_panel, density_colorbar, spear_box, axis_range, save, save_table

d = load_h1_200kb(); x = d.norm_rad.to_numpy(float); xl = axis_range(x)
FS_T, FS_L, FS_K = 22, 20, 16; rows = []
fig, axs = plt.subplots(1, 1, figsize=(6.4, 5.2), squeeze=False); axs = axs[0]
for ax, v in zip(axs, ["sd"]):
    y = d[v].to_numpy(float); yl = axis_range(y)
    density_panel(ax, x, y)
    r, p = spear_box(ax, x, y, pos=(0.96, 0.95, "right", "top"))
    ax.set_xlim(xl); ax.set_ylim(yl)
    rows.append(dict(model_metric=V_LABEL[v], n_200kb_bins=len(d), spearman=r, p=p,
                     n_outside_axes=int(((x < xl[0]) | (x > xl[1]) | (y < yl[0]) | (y > yl[1])).sum())))
    ax.set_xlabel("RAD", fontsize=FS_L); ax.set_ylabel(V_LABEL[v], fontsize=FS_L); ax.tick_params(labelsize=FS_K)
fig.suptitle("ES model", fontsize=FS_T, y=1.02)
density_colorbar(fig, axs, fs=FS_L - 2)
save(fig, "ExtFig2b_left_SDdistCT_vs_RAD")
S = pd.DataFrame(rows); save_table(S, "ExtFig2b_stats.tsv", float_format="%.4g")
print(S.round(3).to_string(index=False))
