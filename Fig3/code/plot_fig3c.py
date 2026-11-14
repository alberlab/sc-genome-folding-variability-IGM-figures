#!/usr/bin/env python3
"""Figure 3c (IGM rows): chr2 + chr11 concatenated tracks at 1 Mb for each stage model — Mean and SD of dist to CT.
Style follows the original compact track panels (6.2 in wide, 1.05 in per row, font size 9): chr2 and chr11 drawn
end to end on one x axis (Mb, restarting at 0 for chr11), solid line at the junction, grey dashed reference line
(Mean: 0 = CT core surface; SD / IQR: median over the four stages). 1 Mb values are the mean of the non-centromeric 200 kb
bins; 1 Mb bins containing a centromeric 200 kb bin are left blank. y ranges are shared by the four stages.
Outputs: output/figures/Fig3c_IGM_{stage}_chr2_chr11_distCT.* (one per stage), Fig3c_IGM_all_stages_chr2_chr11_distCT.*
         (combined, columns ES | DE | HB | iHLC [blank, not modelled] | HLC), output/source_data/Fig3c_chr2_chr11_1Mb.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
from fig3_common import CELLS, CELL_TITLE, V_LABEL, load_200kb, save, save_table

CHRS = ("chr2", "chr11")
ROWS = [("mean", "Mean(dist\nto CT)", "tab:blue"), ("sd", "SD(dist\nto CT)", "tab:red")]          # IQR row removed (SD only, 2026-11)
FIG_W, ROW_H, LW, FS = 6.2, 1.05, 1.0, 9

tabs = []
for c in CELLS:
    d = load_200kb(c, keep_cen=True)
    d = d[d.chr.isin(CHRS)].copy(); d["start_1mb"] = (d.start // 1_000_000) * 1_000_000
    m = d[d.label != "cen"].groupby(["chr", "start_1mb"])[["mean", "sd", "iqr"]].mean()
    cen = d[d.label == "cen"].groupby(["chr", "start_1mb"]).size().index
    grid = d.groupby(["chr", "start_1mb"]).end.max().rename("end").to_frame()
    t = grid.join(m); t["masked"] = t.index.isin(cen)
    t.loc[t.masked, ["mean", "sd", "iqr"]] = np.nan
    tabs.append(t.reset_index().rename(columns={"start_1mb": "start"}).assign(cell=c))
T = pd.concat(tabs, ignore_index=True)
save_table(T[["cell", "chr", "start", "end", "masked", "mean", "sd", "iqr"]], "Fig3c_chr2_chr11_1Mb.tsv")

chr_len = T.groupby("chr").end.max().to_dict()
len_a = chr_len[CHRS[0]]; total = len_a + chr_len[CHRS[1]]
offset = {CHRS[0]: 0, CHRS[1]: len_a}
ylim, yref = {}, {}
for col, _, _ in ROWS:                                     # shared over the four stages
    v = T[col].dropna().to_numpy(); lo, hi = np.percentile(v, [0.5, 99.5]); pad = 0.05 * (hi - lo)
    ylim[col] = (lo - pad, hi + pad); yref[col] = 0.0 if col == "mean" else float(np.median(v))
ticks = [(offset[ch] + p, f"{p / 1e6:.0f}") for ch in CHRS for p in range(0, int(chr_len[ch]), 50_000_000)]
print("chromosome lengths (Mb):", {k: round(v / 1e6, 1) for k, v in chr_len.items()},
      "| masked 1 Mb bins:", T[T.cell == "h1"].groupby("chr").masked.sum().to_dict())


def draw_column(axes, t, ylabels=True):
    """Draw one stage's tracks into a column of axes (Mean / SD)."""
    for ax, (col, ylab, color) in zip(axes, ROWS):
        for ch in CHRS:                                    # chromosomes drawn separately: no line across the junction
            s = t[t.chr == ch]
            ax.plot(s.start + offset[ch], s[col], lw=LW, color=color)
        ax.axhline(yref[col], ls="--", lw=1.0, color="gray")
        ax.axvline(len_a, ls="-", lw=1.0, color="black")   # chr2 | chr11 junction (solid)
        ax.set_ylim(*(ylim[col][::-1] if col == "mean" else ylim[col]))   # Mean flipped: large values at the bottom (2026-11-13)
        if ylabels:
            ax.set_ylabel(ylab, fontsize=FS)
        ax.tick_params(labelsize=FS - 1)
        for sp in ax.spines.values():
            sp.set_linewidth(0.9)
    for ch, xc in ((CHRS[0], len_a / 2), (CHRS[1], len_a + chr_len[CHRS[1]] / 2)):
        axes[0].text(xc, 1.04, ch, transform=axes[0].get_xaxis_transform(), ha="center", va="bottom", fontsize=FS)
    axes[-1].set_xlim(0, total); axes[-1].margins(x=0)
    axes[-1].xaxis.set_major_locator(FixedLocator([p for p, _ in ticks]))
    axes[-1].xaxis.set_major_formatter(FixedFormatter([l for _, l in ticks]))
    axes[-1].set_xlabel("Genomic position (Mb)", fontsize=FS)


# ---- one figure per stage
for c in CELLS:
    fig, axes = plt.subplots(len(ROWS), 1, figsize=(FIG_W, len(ROWS) * ROW_H), sharex=True)
    draw_column(axes, T[T.cell == c])
    fig.suptitle(f"IGM model — {CELL_TITLE[c]}", fontsize=FS + 1, y=1.0)
    fig.tight_layout(h_pad=0.4)
    save(fig, f"Fig3c_IGM_{CELL_TITLE[c]}_chr2_chr11_distCT")

# ---- combined figure in the layout of the original Fig. 3c: ES | DE | HB | iHLC (no model; left blank) | HLC
COLS = ["h1", "DE", "HB", None, "mHEP"]
fig, axes = plt.subplots(len(ROWS), len(COLS), figsize=(FIG_W * len(COLS), len(ROWS) * ROW_H + 0.5), sharex="col",
                         sharey="row", squeeze=False)
for j, c in enumerate(COLS):
    if c is None:                                          # iHLC column kept empty
        for ax in axes[:, j]:
            ax.set_axis_off()
        continue
    draw_column(axes[:, j], T[T.cell == c], ylabels=(j == 0))
    axes[0, j].set_title(CELL_TITLE[c], fontsize=FS + 2, fontweight="bold", pad=16)
    for ax in axes[:, j]:
        ax.tick_params(labelleft=True)                     # keep y tick labels on every column (as the original)
fig.supylabel("IGM", fontsize=FS + 2, fontweight="bold", x=0.002)
fig.tight_layout(h_pad=0.4, w_pad=1.2)
save(fig, "Fig3c_IGM_all_stages_chr2_chr11_distCT")
