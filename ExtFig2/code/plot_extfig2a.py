#!/usr/bin/env python3
"""Extended Data Fig. 2a: H1 model variability of the distance to the CT core vs single-cell Hi-C variability of the cis
contact ratio, 1 Mb bins, density scatter (1 x 2: droplet-based | plate-based scHi-C, model SD only; the earlier 2 x 2
version also had the IQR column).
x  model SD(dist to CT): mean of the autosomal non-centromeric 200 kb bin values inside each 1 Mb bin
y  scHi-C SD(cisratio): column cisratio_S.D. of Supplementary Table 2 (sheet H1-platebased-HiC-batch1) / Table 3
   (sheet "H1-dropletHiC >0.3kb filter")
Bins: 1 Mb bins with >= 1 autosomal non-centromeric 200 kb model bin and a finite scHi-C value; chr2:32-33 Mb (a
droplet trans-contact hotspot) is removed from both assays. Colour = per-point KDE density (panel maximum = 1); axes span
the 0.5-99.5 % range of the plotted values (points outside are counted in the source data). Spearman r per panel.
Outputs: output/figures/ExtFig2a_model_vs_schic_SDcisratio.*, output/source_data/ExtFig2a_{stats,points}.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from extfig2_common import (PLATE_XLSX, PLATE_SHEET, DROPLET_XLSX, DROPLET_SHEET, V_LABEL, load_h1_200kb, to_1mb,
                            read_schic_table, density_panel, density_colorbar, spear_box, axis_range, save, save_table)

HOTSPOT = ("chr2", 32_000_000)
m1 = to_1mb(load_h1_200kb(), ["sd", "iqr"])
ASSAYS = [("droplet", "Droplet-based scHi-C", DROPLET_XLSX, DROPLET_SHEET), ("plate", "Plate-based scHi-C", PLATE_XLSX, PLATE_SHEET)]   # manuscript order
pts, rows = [], []
FS_T, FS_L, FS_K = 22, 20, 16
fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2), squeeze=False)    # SD only: one row, one column per assay
fig.subplots_adjust(wspace=0.42)
for i, (key, title, xlsx, sheet) in enumerate(ASSAYS):
    s = read_schic_table(xlsx, sheet, ["cisratio_S.D."]).rename(columns={"cisratio_S.D.": "schic_sd"})
    d = m1.merge(s[["chr", "start", "schic_sd"]], on=["chr", "start"], how="inner", validate="1:1")
    d = d[np.isfinite(d.schic_sd) & ~((d.chr == HOTSPOT[0]) & (d.start == HOTSPOT[1]))].reset_index(drop=True)
    pts.append(d.assign(assay=key)[["assay", "chr", "start", "n_200kb", "sd", "iqr", "schic_sd"]])
    yl = axis_range(d.schic_sd)
    for j, v in enumerate(["sd"]):
        ax = axs[j, i]; x, y = d[v].to_numpy(float), d.schic_sd.to_numpy(float)
        density_panel(ax, x, y, s=6)
        r, p = spear_box(ax, x, y)
        xl = axis_range(x); ax.set_xlim(xl); ax.set_ylim(yl)
        out = int(((x < xl[0]) | (x > xl[1]) | (y < yl[0]) | (y > yl[1])).sum())
        rows.append(dict(assay=key, sheet=sheet, model_metric=V_LABEL[v], n_1mb_bins=len(d), spearman=r, p=p,
                         n_outside_axes=out))
        ax.set_xlabel(f"Model {V_LABEL[v]}", fontsize=FS_L); ax.tick_params(labelsize=FS_K)
        ax.set_ylabel(f"{title.split()[0]}\nSD(cisratio)", fontsize=FS_L)
        ax.set_title(title, fontsize=FS_T)
density_colorbar(fig, axs, fs=FS_L - 2)
save(fig, "ExtFig2a_model_vs_schic_SDcisratio")
S = pd.DataFrame(rows); save_table(S, "ExtFig2a_stats.tsv", float_format="%.4g")
save_table(pd.concat(pts, ignore_index=True), "ExtFig2a_points.tsv")
print(S[["assay", "model_metric", "n_1mb_bins", "spearman", "n_outside_axes"]].round(3).to_string(index=False))
