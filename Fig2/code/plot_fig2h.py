#!/usr/bin/env python3
"""Figure 2h: chr11, 5 representative loci (imaging-defined Variable / Consistent groups), per-structure dist to CT.
  left   stacked bars, inside (dist >= 0) / outside fraction per locus
  right  combined histogram (density, 22 bins over [-3, 3]) of Variable vs Consistent loci
Per-structure values (data/h1_occ06_chr11_5loci_distCT.npz): 2 copies x 1000 structures per 200 kb bin; a locus pools
all its 200 kb bins (119M 5 bins, 11M 5, 107M 10, 41M 5, 88M 2).
Outputs: output/figures/Fig2h_chr11_5loci_bar_hist.*, output/source_data/Fig2h_chr11_5loci.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
from matplotlib.gridspec import GridSpec
from fig2_common import CT_SCALE, CT_UNIT, CT_TSV, LOCI_NPZ, save, save_table

LOCI = [(0, "chr11:119-120 Mb", "Variable"), (1, "chr11:11-12 Mb", "Variable"), (2, "chr11:107-109 Mb", "Variable"),
        (3, "chr11:41-42 Mb", "Consistent"), (4, "chr11:88-88.4 Mb", "Consistent")]
INSIDE_COLOR, OUTSIDE_COLOR = "#4C78A8", "#8B1A1A"
GC = {"Variable": "#0D47A1", "Consistent": "#1B5E20"}
XR, NB = (-3, 3), 22
if CT_SCALE != 1.0:
    XR = (-3 * CT_SCALE, 3 * CT_SCALE)

z = np.load(LOCI_NPZ, allow_pickle=True)
ct = pd.read_csv(CT_TSV, sep="\t")                                    # loci bins must exist in the dist-to-CT table
key = set(zip(ct.chr, ct.start))
assert all((c, s) in key for c, s in zip(z["chrom"], z["start"]))
vals, rows = {}, []
for i, name, grp in LOCI:
    cols = np.flatnonzero(z["locus"] == i); v = z["dist"][:, cols].ravel(); v = v[np.isfinite(v)]
    vals[name] = v
    rows.append(dict(locus=name, group=grp, n_bins=len(cols), n=len(v), inside_pct=100 * (v >= 0).mean(),
                     median=np.median(v), pct_below_m3=100 * (v < XR[0]).mean()))
T = pd.DataFrame(rows)
G = {g: np.concatenate([vals[n] for _, n, gg in LOCI if gg == g]) for g in GC}
for g, v in G.items():
    T.loc[len(T)] = dict(locus=f"pooled {g}", group=g, n_bins=np.nan, n=len(v), inside_pct=100 * (v >= 0).mean(),
                         median=np.median(v), pct_below_m3=100 * (v < XR[0]).mean())
save_table(T, "Fig2h_chr11_5loci.tsv", float_format="%.4g")
print(T.to_string(index=False))

sns.set_theme(style="ticks")
import matplotlib; matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                                               "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
fig = plt.figure(figsize=(12.5, 5.2)); gs = GridSpec(1, 2, width_ratios=[1.05, 1.2], wspace=0.35, figure=fig)

# ---- left: stacked inside / outside bars
ax = fig.add_subplot(gs[0]); xs = np.arange(len(LOCI)); ins = T.inside_pct[:len(LOCI)].to_numpy()
ax.bar(xs, ins, color=INSIDE_COLOR, edgecolor="black", lw=1.0, label="Inside")
ax.bar(xs, 100 - ins, bottom=ins, color=OUTSIDE_COLOR, edgecolor="black", lw=1.0, label="Outside")
for xi, (_, r) in zip(xs, T[:len(LOCI)].iterrows()):
    ax.text(xi, ins[xi] / 2, f"{ins[xi]:.0f}%", ha="center", va="center", color="white", fontsize=20)
    ax.text(xi, 101.5, f"n={r.n / 1000:g}k", ha="center", va="bottom", fontsize=15)
for g in GC:                                                          # group bracket above the bars
    idx = [k for k, (_, _, gg) in enumerate(LOCI) if gg == g]
    ax.plot([idx[0] - 0.35, idx[-1] + 0.35], [113, 113], color=GC[g], lw=2, clip_on=False)
    ax.text(np.mean(idx), 114, g, color=GC[g], ha="center", va="bottom", fontsize=20, fontweight="bold")
ax.axhline(50, ls="--", color="black", alpha=0.6, lw=0.8)
ax.set_xticks(xs, [n for _, n, _ in LOCI], rotation=35, ha="right", rotation_mode="anchor", fontsize=18)
ax.set_ylim(0, 108); ax.set_yticks(range(0, 101, 20)); ax.set_ylabel("Fraction of loci (%)", fontsize=20); ax.tick_params(axis="y", labelsize=18)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.42), ncol=2, fontsize=18, frameon=False)

# ---- right: combined Variable vs Consistent histogram
ax2 = fig.add_subplot(gs[1]); bins = np.linspace(*XR, NB + 1)
for g in ("Consistent", "Variable"):
    v = G[g]; vc = v[(v >= XR[0]) & (v <= XR[1])]
    ax2.hist(vc, bins=bins, density=True, alpha=0.30, color=GC[g], edgecolor=GC[g], lw=0.6, label=g)
ax2.axvline(0, ls=":", color="black", lw=1); ax2.set_ylim(0, ax2.get_ylim()[1] * 1.3)
for yy, g in ((0.97, "Consistent"), (0.79, "Variable")):
    r = T[T.locus == f"pooled {g}"].iloc[0]
    ax2.text(0.03, yy, f"{g}\nin {r.inside_pct:.0f}% / out {100 - r.inside_pct:.0f}%", transform=ax2.transAxes,
             color=GC[g], fontsize=18, va="top", ha="left")
clip = 100 * np.mean(np.concatenate(list(G.values())) < XR[0])
ax2.set_xlim(*XR); ax2.set_title(f"Combined (values < {XR[0]:.3g} not shown: {clip:.2f}%)", fontsize=18)
ax2.set_xlabel("Signed dist to CT edge (normed)" if not CT_UNIT else f"Signed dist to CT edge (/ {CT_UNIT[3:]})", fontsize=20); ax2.set_ylabel("Density", fontsize=20); ax2.tick_params(labelsize=18)
ax2.legend(loc="upper right", fontsize=18, frameon=False)
fig.suptitle("H1 model — chr11 representative loci, dist to CT\n(CT core: top 30 % strict-cis beads, 2.5 r; inside = dist ≥ 0)",
             fontsize=22, y=1.12)
save(fig, "Fig2h_chr11_5loci_bar_hist")
