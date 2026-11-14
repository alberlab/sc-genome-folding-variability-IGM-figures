#!/usr/bin/env python3
"""Figure 2d: nuclear locale entropy (5 exclusive classes, bits) of the H1 model.
  left    entropy by dominant locale, box plot (h1 has no nuc-dominant bin -> 4 groups, order spe / NL-act / NL-rep / env)
  middle  entropy vs SD(dist to CT), chr2, 200 kb, colour = dominant locale (alt: genome-wide)
  right   entropy vs Eig1, colour = dominant locale: 200 kb bins vs bulk 200 kb E1 (manuscript panel);
          alt = 1 Mb mean of 200 kb entropy vs the 1 Mb E1 used in Fig 2c (dominant locale at 1 Mb = argmax of the 1 Mb-mean P_*)
Locale = ENV v2 (env seeds = outermost 20 % of all beads + 2.5 r).
Outputs: output/figures/Fig2d_*, output/source_data/Fig2d_*.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from fig2_common import (TEXT_POS, P5, LOC_COLORS, LOC_DISP, V_LABEL, load_200kb, load_e1_1mb, to_1mb, spear_text, lim, save,
                         save_table)

ORDER = ["spe", "act", "env", "rep"]
BOX_ORDER = ["spe", "act", "rep", "env"]            # panel 1 only: NL-act and NL-rep in the middle (2026-11-04)
ENT = "Nuclear locale entropy (bits)"
d = load_200kb()
cnt = d.dom5.value_counts().reindex(["spe", "nuc", "env", "act", "rep"]).fillna(0).astype(int)
print("dominant locale counts (200 kb):", cnt.to_dict())
assert cnt["nuc"] == 0, "h1 now has nuc-dominant bins - revisit the 4-group box plot"
save_table(d, "Fig2d_h1_200kb.tsv")
stats = []

# ---- panel 1: entropy by dominant locale
fig, ax = plt.subplots(figsize=(4.8, 5))
data = [d.loc[d.dom5 == g, "entropy_bit"].to_numpy() for g in BOX_ORDER]
bp = ax.boxplot(data, widths=0.62, patch_artist=True, showfliers=False, whis=1.5,
                medianprops=dict(color="black", lw=1.4))
for patch, g in zip(bp["boxes"], BOX_ORDER):
    patch.set_facecolor(LOC_COLORS[g]); patch.set_alpha(0.75); patch.set_edgecolor("black")
ax.set_xticks(range(1, len(BOX_ORDER) + 1), [f"{LOC_DISP[g]}\n(n={len(x)})" for g, x in zip(BOX_ORDER, data)], rotation=45,
              ha="right", rotation_mode="anchor")
ax.set_ylabel(ENT); ax.set_xlabel("dominant locale")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig, "Fig2d_left_entropy_by_dominant_locale")
for g, x in zip(BOX_ORDER, data):
    stats.append(dict(panel="p1_box", group=LOC_DISP[g], n=len(x), median=np.median(x), mean=x.mean()))


PEARSON = {}


def scatter(ax, df, xcol, xlab, s, loc="lr", xl=None, dense=True, pearson=False):
    """dense = many points (200 kb): translucent, no edge, rasterised; else opaque with a thin white edge."""
    x = df[xcol].to_numpy(float); y = df.entropy_bit.to_numpy(float); ok = np.isfinite(x) & np.isfinite(y)
    for g in ["rep", "act", "env", "spe"]:                   # spe on top (smallest group)
        k = ok & (df.dom5.to_numpy() == g)
        if k.sum():
            ax.scatter(x[k], y[k], s=s, c=LOC_COLORS[g], alpha=0.55 if dense else 0.8, linewidths=0 if dense else 0.4,
                       edgecolor="white", rasterized=dense, label=LOC_DISP[g])
    if pearson:                                              # middle panel: Spearman + Pearson (2026-11-13)
        from scipy.stats import spearmanr, pearsonr
        r, p = spearmanr(x[ok], y[ok]); pr, pp = pearsonr(x[ok], y[ok]); tx, ty, ha, va = TEXT_POS[loc]
        fmt = lambda q: f"{q:.2e}" if q > 0 else "< 1e-300"
        ax.text(tx, ty, f"Spearman r = {r:.2f}, p = {fmt(p)}\nPearson r = {pr:.2f}, p = {fmt(pp)}\nn = {int(ok.sum())}",
                transform=ax.transAxes, fontsize=10, ha=ha, va=va, bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2))
        PEARSON[xlab, int(ok.sum())] = (pr, pp)
    else:
        r, p = spear_text(ax, x[ok], y[ok], loc=loc)
    ax.set_xlim(xl if xl else lim(x[ok])); ax.set_ylim(lim(y[ok])); ax.set_xlabel(xlab); ax.set_ylabel(ENT)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return r, p, int(ok.sum())


def legend(ax):
    h, l = ax.get_legend_handles_labels(); idx = [l.index(LOC_DISP[g]) for g in ORDER if LOC_DISP[g] in l]
    ax.legend([h[i] for i in idx], [l[i] for i in idx], frameon=False, markerscale=max(1.0, np.sqrt(50 / h[0].get_sizes()[0])),
              loc="center left", bbox_to_anchor=(1.0, 0.5), title="dominant locale")


# ---- panel 2: entropy vs IQR / SD (200 kb), genome-wide and chr2 only (dot size 2.5x the 2026-10-16 version)
S_200KB, S_1MB = 40, 75
for v in ("sd",):                                                       # SD only (2026-11)
    fig, ax = plt.subplots(figsize=(6.2, 5)); r, p, n = scatter(ax, d, v, V_LABEL[v], s=S_200KB, pearson=True); legend(ax)
    fig.tight_layout(); save(fig, f"Fig2d_middle_alt_entropy_vs_{v.upper()}_genomewide")
    stats.append(dict(panel=f"p2_entropy_vs_{v.upper()}_200kb", n=n, spearman=r, p=p, pearson=PEARSON[V_LABEL[v], n][0], pearson_p=PEARSON[V_LABEL[v], n][1]))
    c2 = d[d.chr == "chr2"].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.2, 5)); r, p, n = scatter(ax, c2, v, V_LABEL[v], s=S_200KB, dense=False, pearson=True); legend(ax)
    ax.set_title("chr2 (200 kb)"); fig.tight_layout(); save(fig, f"Fig2d_middle_entropy_vs_{v.upper()}_chr2")
    stats.append(dict(panel=f"p2_entropy_vs_{v.upper()}_200kb_chr2", n=n, spearman=r, p=p, pearson=PEARSON[V_LABEL[v], n][0], pearson_p=PEARSON[V_LABEL[v], n][1]))

# ---- panel 3: entropy vs Eig1 (200 kb and 1 Mb)
fig, ax = plt.subplots(figsize=(6.2, 5)); r, p, n = scatter(ax, d, "E1_200kb", "Eig1 (200 kb)", s=6, loc="ll", xl=(-1.5, 1.2)); legend(ax)
fig.tight_layout(); save(fig, "Fig2d_right_entropy_vs_Eig1_200kb")
stats.append(dict(panel="p3_entropy_vs_Eig1_200kb", n=n, spearman=r, p=p))

m = to_1mb(d, ["entropy_bit"] + P5).merge(load_e1_1mb(), on=["chr", "start"], how="left", validate="1:1")
m = m[np.isfinite(m.E1_1mb)].reset_index(drop=True)
m["dom5"] = np.array([c[2:] for c in P5])[np.argmax(m[P5].to_numpy(float), axis=1)]
save_table(m, "Fig2d_h1_1Mb.tsv")
print("dominant locale counts (1 Mb):", m.dom5.value_counts().to_dict())
fig, ax = plt.subplots(figsize=(6.2, 5)); r, p, n = scatter(ax, m, "E1_1mb", "Eig1", s=S_1MB, loc="ll", xl=(-1.5, 1.2), dense=False); legend(ax)
fig.tight_layout(); save(fig, "Fig2d_right_alt_entropy_vs_Eig1_1Mb")
stats.append(dict(panel="p3_entropy_vs_Eig1_1mb", n=n, spearman=r, p=p))

S = pd.DataFrame(stats); save_table(S, "Fig2d_stats.tsv", float_format="%.4g")
print(S.to_string(index=False))
