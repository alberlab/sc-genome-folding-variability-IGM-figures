#!/usr/bin/env python3
"""Extended Data Fig. 2h (old name 2e): chr2, 3 representative imaging-Variable loci, per-structure dist to CT in the H1 model
(same layout as Fig. 2h).
  left   stacked bars, inside (dist >= 0) / outside fraction per locus
  right  histogram (density, 22 bins over [-3, 3]) of the three loci pooled
Per-structure values (data/h1_occ06_chr2_3loci_distCT.npz): 2 copies x 1000 structures per 200 kb bin; a locus pools
all its 200 kb bins (32-33 Mb 5 bins, 2-3.8 Mb 9, 100.4-102 Mb 8).
Outputs: output/figures/ExtFig2h_chr2_3loci_bar_hist.*, output/source_data/ExtFig2h_chr2_3loci.tsv
"""
import numpy as np, pandas as pd, matplotlib, matplotlib.pyplot as plt, seaborn as sns
from matplotlib.gridspec import GridSpec
from extfig2_common import CT_SCALE, CT_UNIT, CT_TSV, LOCI_NPZ, save, save_table

LOCI = [(0, "chr2:32-33 Mb"), (1, "chr2:2-3.8 Mb"), (2, "chr2:100.4-102 Mb")]
INSIDE_COLOR, OUTSIDE_COLOR, VAR_COLOR = "#4C78A8", "#8B1A1A", "#0D47A1"
XR, NB = (-3, 3), 22
if CT_SCALE != 1.0:
    XR = (-3 * CT_SCALE, 3 * CT_SCALE)

z = np.load(LOCI_NPZ, allow_pickle=True)
ct = pd.read_csv(CT_TSV, sep="\t")                                    # loci bins must exist in the dist-to-CT table
key = set(zip(ct.chr, ct.start))
assert all((c, s) in key for c, s in zip(z["chrom"], z["start"]))
vals, rows = {}, []
for i, name in LOCI:
    cols = np.flatnonzero(z["locus"] == i); v = z["dist"][:, cols].ravel(); v = v[np.isfinite(v)]
    vals[name] = v
    rows.append(dict(locus=name, group="Variable", n_bins=len(cols), n=len(v), inside_pct=100 * (v >= 0).mean(),
                     median=np.median(v), pct_below_m3=100 * (v < XR[0]).mean()))
T = pd.DataFrame(rows)
P = np.concatenate(list(vals.values()))
T.loc[len(T)] = dict(locus="pooled Variable", group="Variable", n_bins=np.nan, n=len(P), inside_pct=100 * (P >= 0).mean(),
                     median=np.median(P), pct_below_m3=100 * (P < XR[0]).mean())
save_table(T, "ExtFig2h_chr2_3loci.tsv", float_format="%.4g")
print(T.to_string(index=False))

sns.set_theme(style="ticks")
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
fig = plt.figure(figsize=(11, 5.2)); gs = GridSpec(1, 2, width_ratios=[0.8, 1.2], wspace=0.35, figure=fig)

# ---- left: stacked inside / outside bars
ax = fig.add_subplot(gs[0]); xs = np.arange(len(LOCI)); ins = T.inside_pct[:len(LOCI)].to_numpy()
ax.bar(xs, ins, color=INSIDE_COLOR, edgecolor="black", lw=1.0, label="Inside")
ax.bar(xs, 100 - ins, bottom=ins, color=OUTSIDE_COLOR, edgecolor="black", lw=1.0, label="Outside")
for xi, (_, r) in zip(xs, T[:len(LOCI)].iterrows()):
    ax.text(xi, ins[xi] / 2, f"{ins[xi]:.0f}%", ha="center", va="center", color="white", fontsize=20)
    ax.text(xi, 101.5, f"n={r.n / 1000:g}k", ha="center", va="bottom", fontsize=15)
ax.plot([-0.35, len(LOCI) - 1 + 0.35], [113, 113], color=VAR_COLOR, lw=2, clip_on=False)   # group bracket
ax.text(np.mean(xs), 114, "Variable", color=VAR_COLOR, ha="center", va="bottom", fontsize=20, fontweight="bold")
ax.axhline(50, ls="--", color="black", alpha=0.6, lw=0.8)
ax.set_xticks(xs, [n for _, n in LOCI], rotation=35, ha="right", rotation_mode="anchor", fontsize=18)
ax.set_xlim(-0.6, len(LOCI) - 0.4)
ax.set_ylim(0, 108); ax.set_yticks(range(0, 101, 20)); ax.set_ylabel("Fraction of loci (%)", fontsize=20); ax.tick_params(axis="y", labelsize=18)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.42), ncol=2, fontsize=18, frameon=False)

# ---- right: pooled histogram of the three loci
ax2 = fig.add_subplot(gs[1]); bins = np.linspace(*XR, NB + 1)
vc = P[(P >= XR[0]) & (P <= XR[1])]
ax2.hist(vc, bins=bins, density=True, alpha=0.30, color=VAR_COLOR, edgecolor=VAR_COLOR, lw=0.6, label="Variable")
ax2.axvline(0, ls=":", color="black", lw=1); ax2.set_ylim(0, ax2.get_ylim()[1] * 1.3)
r = T.iloc[-1]
ax2.text(0.03, 0.97, f"Variable\nin {r.inside_pct:.0f}% / out {100 - r.inside_pct:.0f}%", transform=ax2.transAxes,
         color=VAR_COLOR, fontsize=18, va="top", ha="left")
ax2.set_xlim(*XR); ax2.set_title(f"Combined (values < {XR[0]:.3g} not shown: {r.pct_below_m3:.2f}%)", fontsize=18)
ax2.set_xlabel("Signed dist to CT edge (normed)" if not CT_UNIT else f"Signed dist to CT edge (/ {CT_UNIT[3:]})", fontsize=20); ax2.set_ylabel("Density", fontsize=20); ax2.tick_params(labelsize=18)
ax2.legend(loc="upper right", fontsize=18, frameon=False)
fig.suptitle("H1 model — chr2 representative loci, dist to CT\n(CT core: top 30 % strict-cis beads, 2.5 r; inside = dist ≥ 0)",
             fontsize=22, y=1.12)
save(fig, "ExtFig2h_chr2_3loci_bar_hist")
