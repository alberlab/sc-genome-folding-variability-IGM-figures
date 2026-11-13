#!/usr/bin/env python3
"""Extended Data Fig. 2c (bar plots): nuclear-locale frequencies of 3 representative 200 kb loci in the ES (H1) model,
one panel per locus, horizontal bars of the five locale classes. Loci: high speckle frequency chr12:121.4-121.6 Mb,
high envelope frequency chr21:10.0-10.2 Mb, evenly distributed chr7:57.2-57.4 Mb. Frequency = fraction of the valid
observations (2 copies x 1000 structures) in each class (P_* of the locale table); beads near two or three landmarks
(overlap classes, 4-7 % here) are not drawn, so the five bars need not sum to 1. Box = locale entropy (bits) of the bin.
The 3D structures of the four stage models (beads coloured by locale) are the cmm/pdb files in data/structures/.
Outputs: output/figures/ExtFig2c_example_loci_locale_freq.*, output/source_data/ExtFig2c_example_loci_locale_freq.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from extfig2_common import LOC5, LOC_COLORS, LOC_DISP, locale_table, save, save_table

LOCI = [("High spe probability", "Low locale entropy", "chr12", 121_400_000),
        ("High env probability", "Low locale entropy", "chr21", 10_000_000),
        ("Evenly distributed at different locales", "High locale entropy", "chr7", 57_200_000)]
ORDER = ["spe", "act", "rep", "nuc", "env"]                          # top -> bottom
t = locale_table("h1").set_index(["chr", "start"]); rows = []
for typ, _, ch, s in LOCI:
    r = t.loc[(ch, s)]; assert r.label != "cen"
    f = {k: float(r[f"P_{k}"]) for k in LOC5}
    rows.append(dict(locus_type=typ, region=f"{ch}:{s}-{s + 200_000}", model="ES (H1)", entropy_bit=float(r.entropy_bit),
                     n_valid=int(r.n_valid), **{f"P_{k}": f[k] for k in ORDER}, P_other=1 - sum(f.values())))
T = pd.DataFrame(rows); save_table(T, "ExtFig2c_example_loci_locale_freq.tsv", float_format="%.4f")
print(T.drop(columns=["region", "model", "n_valid"]).round(3).to_string(index=False))

FS = 13
fig, axs = plt.subplots(1, len(LOCI), figsize=(4.0 * len(LOCI), 3.2))
y = np.arange(len(ORDER))
for ax, (typ, ent, ch, s), (_, d) in zip(axs, LOCI, T.iterrows()):
    v = [d[f"P_{k}"] for k in ORDER]
    ax.barh(y, v, color=[LOC_COLORS[k] for k in ORDER], edgecolor="#333333", lw=0.5, height=0.75)
    ax.set_yticks(y, [LOC_DISP[k] for k in ORDER], fontsize=FS - 1); ax.invert_yaxis()
    xmax = max(0.6, np.ceil((max(v) + 0.1) / 0.2) * 0.2)
    ax.set_xlim(0, xmax); ax.set_xticks(np.arange(0, xmax + 1e-9, 0.2)); ax.tick_params(axis="x", labelsize=FS - 2)
    ax.text(0.96, 0.5, f"entropy: {d.entropy_bit:.2f}", transform=ax.transAxes, ha="right", va="center",
            fontsize=FS - 3, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#888888", lw=0.8))
    ax.set_xlabel("Locale probability", fontsize=FS)
    ax.set_title(f"{typ}\n{ent}", fontsize=FS - 1, loc="left")
    ax.text(0, 1.30, f"{ch}:{s / 1e6:.2f}-{(s + 200_000) / 1e6:.2f}Mb", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=FS, fontweight="bold")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
fig.tight_layout(w_pad=1.5); save(fig, "ExtFig2c_example_loci_locale_freq")
