#!/usr/bin/env python3
"""Extended Data Fig. 6d: SD of the distance to the CT core (top row) and nuclear locale entropy (bottom row) across the
four stage models, by dominant locale. Columns = dominant locale spe, NL-act, env, NL-rep (bins grouped by their dominant
locale in each stage, so a bin can move between groups; nuc-dominant bins, <= 8 per stage, omitted). Each panel: one box
per stage (median, IQR, whiskers 1.5 x IQR, no fliers), dashed line joining the medians, n under each box.
Bottom row: each differentiated stage vs ES, two-sided Mann-Whitney U (unpaired: a bin's group is its dominant locale in
that stage, so the groups differ between stages), brackets with p, for all four locales.
One figure combining the two earlier ones: the SD row of draw() ("ExtFig6d_box_SD_by_locale") and
draw_entropy_tests(groups = all four) ("ExtFig6d_box_entropy_4loc_MWU") of the earlier plot_extfig6d.py; the box / bracket
drawing code is unchanged, only placed on one 2 x 4 grid (shared y axis within each row). Locale = ENV v2.
Outputs: output/figures/ExtFig6d_box_SD_entropy_by_dominant_locale_MWU.*, output/source_data/ExtFig6d_box.tsv,
         ExtFig6d_entropy_MWU_tests.tsv
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from scipy.stats import mannwhitneyu
from extfig6_common import CELLS, CELL_TITLE, LOC_COLORS, LOC_DISP, V_LABEL, load_stage, save, save_table

D = {c: load_stage(c) for c in CELLS}
GROUPS = ["spe", "act", "env", "rep"]
FS = 15; rows = []; test_rows = []


def boxes(ax, vals, g):
    """Box style of the earlier draw() / draw_entropy_tests(); returns the medians."""
    pos = np.arange(len(CELLS)); lum = np.dot(to_rgb(LOC_COLORS[g]), [0.299, 0.587, 0.114])
    ax.boxplot(vals, positions=pos, widths=0.7, patch_artist=True, showfliers=False,
               medianprops=dict(color="white" if lum < 0.5 else "black", lw=2.2),
               boxprops=dict(facecolor=LOC_COLORS[g], edgecolor="black", lw=1.2),
               whiskerprops=dict(color="black", lw=1.2), capprops=dict(color="black", lw=1.2))
    meds = [float(np.median(x)) for x in vals]
    ax.plot(pos, meds, ls="--", lw=1.4, color="0.25", marker="o", ms=4, mfc="0.25", zorder=3)
    ax.set_xticks(pos, [f"{CELL_TITLE[c]}\nn={len(x)}" for c, x in zip(CELLS, vals)], fontsize=FS - 5)
    ax.tick_params(axis="y", labelsize=FS - 4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return meds


fig, axs = plt.subplots(2, len(GROUPS), figsize=(3.6 * len(GROUPS), 4.2 + 4.9), sharey="row",
                        gridspec_kw=dict(height_ratios=[4.2, 4.9]))

# ---- top row: SD(dist to CT)  (= draw("sd", [("sd", V_LABEL["sd"])]) of the earlier script)
for gi, g in enumerate(GROUPS):
    ax = axs[0, gi]; vals = [D[c].loc[D[c].dom == g, "sd"].to_numpy() for c in CELLS]; meds = boxes(ax, vals, g)
    ax.set_title(LOC_DISP[g], fontsize=FS + 2, fontweight="bold", color=LOC_COLORS[g])
    if gi == 0:
        ax.set_ylabel(V_LABEL["sd"], fontsize=FS)
    for c, x, m in zip(CELLS, vals, meds):
        rows.append(dict(version=V_LABEL["sd"], metric=V_LABEL["sd"], locale=LOC_DISP[g], cell=c, n=len(x), median=m,
                         q25=float(np.percentile(x, 25)), q75=float(np.percentile(x, 75))))

# ---- bottom row: locale entropy + Mann-Whitney U vs ES  (= draw_entropy_tests(groups=GROUPS) of the earlier script)
ref = CELLS[0]; others = CELLS[1:]
vals = {g: [D[c].loc[D[c].dom == g, "entropy_bit"].to_numpy() for c in CELLS] for g in GROUPS}
wl = {g: [(np.percentile(x, 25) - 1.5 * (np.percentile(x, 75) - np.percentile(x, 25)),
           np.percentile(x, 75) + 1.5 * (np.percentile(x, 75) - np.percentile(x, 25))) for x in vals[g]] for g in GROUPS}
lo = min(max(w[0], min(x)) for g in GROUPS for w, x in zip(wl[g], vals[g]))
hi = max(min(w[1], max(x)) for g in GROUPS for w, x in zip(wl[g], vals[g])); span = hi - lo
for gi, g in enumerate(GROUPS):
    ax = axs[1, gi]; v = vals[g]; meds = boxes(ax, v, g)
    for k, c in enumerate(others, start=1):                        # bracket ES - <stage>, stacked upwards
        u, p = mannwhitneyu(v[k], v[0], alternative="two-sided")
        yb = hi + (0.06 + 0.13 * (k - 1)) * span
        ax.plot([0, 0, k, k], [yb, yb + 0.03 * span, yb + 0.03 * span, yb], color="black", lw=1.1)
        ax.text(k / 2, yb + 0.04 * span, f"p = {p:.1e}" if p >= 1e-300 else "p < 1e-300", ha="center", va="bottom",
                fontsize=FS - 4)
        test_rows.append(dict(locale=LOC_DISP[g], comparison=f"{CELL_TITLE[c]} vs {CELL_TITLE[ref]}", test="Mann-Whitney U, two-sided",
                              n_ref=len(v[0]), n_stage=len(v[k]), median_ref=meds[0], median_stage=meds[k],
                              median_diff=meds[k] - meds[0], U=float(u), p=float(p)))
    if gi == 0:
        ax.set_ylabel("Locale entropy (bits)", fontsize=FS)
    for c, x, m in zip(CELLS, v, meds):
        rows.append(dict(version=V_LABEL["sd"], metric="Locale entropy (bits)", locale=LOC_DISP[g], cell=c, n=len(x), median=m,
                         q25=float(np.percentile(x, 25)), q75=float(np.percentile(x, 75))))
axs[1, 0].set_ylim(lo - 0.05 * span, hi + (0.06 + 0.13 * len(others) + 0.04) * span)
fig.tight_layout(); save(fig, "ExtFig6d_box_SD_entropy_by_dominant_locale_MWU")

S = pd.DataFrame(rows); save_table(S, "ExtFig6d_box.tsv", float_format="%.4g")
T = pd.DataFrame(test_rows); save_table(T, "ExtFig6d_entropy_MWU_tests.tsv", float_format="%.4g")
pd.set_option("display.width", 200)
print(S.pivot_table(index=["metric", "locale"], columns="cell", values="median")[CELLS].round(3).to_string())
print(T.to_string(index=False))
