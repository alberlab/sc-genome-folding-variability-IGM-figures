#!/usr/bin/env python3
"""Figure 2b: chr11 tracks at 1 Mb — bulk Hi-C, plate-based and droplet scHi-C cis/total ratio, and the H1 model's
distance to the CT core. scHi-C values: Supplementary Tables 2 (plate, sheet H1-platebased-HiC-batch1) and 3 (droplet,
sheet "H1-dropletHiC >0.3kb filter"); bulk cis/total from data/h1_bulk_cis_total_1Mb.tsv. A 1 Mb bin is masked (grey) when
it contains a centromeric 200 kb bin or its Eig1 is NaN; the same mask is applied to every track (as in Fig. 2e).
Outputs: output/figures/Fig2b_chr11_tracks.*, output/source_data/Fig2b_chr11_1Mb.tsv
"""
import numpy as np, pandas as pd
from fig2_common import (CT_UNIT, BULK_TSV, PLATE_XLSX, PLATE_SHEET, DROPLET_XLSX, DROPLET_SHEET, read_schic_table, chr_tracks_1mb,
                         plot_tracks, save_table)

CHR = "chr11"
t, spans = chr_tracks_1mb(CHR)
pl = read_schic_table(PLATE_XLSX, PLATE_SHEET, ["cisratio_mean", "cisratio_S.D."])
dr = read_schic_table(DROPLET_XLSX, DROPLET_SHEET, ["cisratio_mean", "cisratio_S.D."])
bk = pd.read_csv(BULK_TSV, sep="\t").rename(columns={"chrom": "chr"})
for df, pre in ((pl, "plate"), (dr, "droplet")):
    s = df[df.chr == CHR].set_index("start")
    t[f"{pre}_mean"] = s["cisratio_mean"]; t[f"{pre}_sd"] = s["cisratio_S.D."]
t["bulk_cis_total"] = bk[bk.chr == CHR].set_index("start")["cis_total_rmdiag"]
SC = ["plate_mean", "droplet_mean", "plate_sd", "droplet_sd", "bulk_cis_total"]
t.loc[t.masked, SC + ["E1_1mb"]] = np.nan
save_table(t.reset_index().rename(columns={"index": "start"}).assign(chr=CHR)[
    ["chr", "start", "end", "masked", "E1_1mb", "bulk_cis_total", "plate_mean", "droplet_mean", "mean", "plate_sd",
     "droplet_sd", "sd", "iqr"]], "Fig2b_chr11_1Mb.tsv")

MEAN_C, VAR_C = "tab:blue", "tab:red"
TRACKS = [("E1_1mb", "Bulk Eig1", "Eig1", "black", None),
          ("bulk_cis_total", "Bulk cis/total ratio", "cis/total", "dimgray", None),
          ("plate_mean", "Plate-based Mean(cisratio)", "Mean", MEAN_C, (0.6, 0.85)),
          ("droplet_mean", "Droplet Mean(cisratio)", "Mean", MEAN_C, (0.6, 0.85)),       # droplet: 0.3 kb filter
          ("mean", f"Model Mean(dist to CT{CT_UNIT})", "Mean", MEAN_C, None),
          ("plate_sd", "Plate-based SD(cisratio)", "SD", VAR_C, None),
          ("droplet_sd", "Droplet SD(cisratio)", "SD", VAR_C, None),
          ("sd", f"Model SD(dist to CT{CT_UNIT})", "SD", VAR_C, None)]                                 # IQR track removed (SD only, 2026-11)
for c, (lo, hi) in (("plate_mean", (0.6, 0.85)), ("droplet_mean", (0.6, 0.85))):
    out = t[c][(t[c] < lo) | (t[c] > hi)]
    if len(out):
        print(f"note: {c} outside the plotted range {lo}-{hi} at", ", ".join(f"{s / 1e6:g} Mb ({v:.3f})" for s, v in out.items()))
plot_tracks(t, spans, TRACKS, CHR, f"H1 {CHR}: bulk Hi-C, scHi-C and model tracks (1 Mb)", f"Fig2b_{CHR}_tracks")
