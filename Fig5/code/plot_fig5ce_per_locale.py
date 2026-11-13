#!/usr/bin/env python3
"""Fig. 5c and 5e with the five locales shown separately (not merged into active / inactive), updated envelope locale.
Same data, counting and style as the earlier plot_fig5c.py / plot_fig5e.py (fig5_common.load_locale reads the ENV v2
locale table data/locale/<LOCALE_VERSION>/<cell>_occ06_locale_N565.tsv), colours = fig5_common.LOC_COLORS.
  5c  composition of ALL bead instances (2 copies x 1000 structures x 13,315 bins): spe, NL-act, nuc, NL-rep, env, other
      (other = near two or three landmarks, or missing). Every value is written; segments < 3 % are labelled at the right.
  5e  Sankey of the dominant locale (argmax of the five exclusive frequencies) ES -> DE -> HB -> HLC; node label = % of
      13,315 bins; nodes too small for an inside label are labelled at the right in the locale colour.
Outputs: ../output/figures/Fig5c_locale_composition_stacked.*, Fig5e_dominant_locale_sankey.*,
         ../output/source_data/Fig5c_per_locale.tsv, Fig5e_per_locale_{node_fractions,transitions}.tsv
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                                        # fig5_common.py of this folder
import fig5_common as fc  # noqa: E402

OUT = fc.OUT; FIG, SRC = OUT / "figures", OUT / "source_data"   # output/ or output_<LOCALE_VERSION>/
FIG.mkdir(parents=True, exist_ok=True); SRC.mkdir(parents=True, exist_ok=True)
COL, DISP, CELLS, T = fc.LOC_COLORS, fc.LOC_DISP, fc.CELLS, fc.CELL_TITLE
ORDER = ["spe", "act", "nuc", "rep", "env"]                          # top -> bottom


def save(fig, name):
    for e in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig); print("wrote", name)


def txtcol(c):
    return "white" if np.mean(plt.matplotlib.colors.to_rgb(c)) < 0.6 else "black"


# ---------------- 5c: composition of all bead instances, per locale
rows = []
for c in CELLS:
    t = fc.load_locale(c); tot = t.n_obs.sum()
    cnt = {p: float((t[p] * t.n_valid).sum()) for p in fc.P_ALL}
    r = {k: 100 * cnt[f"P_{k}"] / tot for k in ["spe", "act", "nuc", "rep", "env"]}
    r["others"] = 100 - sum(r.values())
    rows.append(dict(cell=c, n_bead_instances=int(tot), **{f"pct_{k}": v for k, v in r.items()}))
S = pd.DataFrame(rows); S.to_csv(SRC / "Fig5c_per_locale.tsv", sep="\t", index=False, float_format="%.3f")
print(S.round(2).to_string(index=False))

STACK = ["others"] + ORDER[::-1]                                     # bottom -> top: other, env, NL-rep, nuc, NL-act, spe
fig, ax = plt.subplots(figsize=(8.2, 6.4))
x = np.arange(len(CELLS)); bottom = np.zeros(len(CELLS)); BW = 0.62
for k in STACK:
    v = S[f"pct_{k}"].to_numpy()
    ax.bar(x, v, bottom=bottom, width=BW, color=COL[k], edgecolor="white", lw=0.6, label=DISP[k])
    for xi, (b, vi) in enumerate(zip(bottom, v)):
        if vi >= 3.0:
            ax.text(xi, b + vi / 2, f"{vi:.1f}", ha="center", va="center", fontsize=12, color=txtcol(COL[k]))
        else:                                                        # small segment: value at the right of the bar
            ax.text(xi + BW / 2 + 0.03, b + vi / 2, f"{vi:.1f}", ha="left", va="center", fontsize=10.5, color=COL[k],
                    fontweight="bold")
    bottom += v
ax.set_xticks(x, [T[c] for c in CELLS], fontsize=14)
ax.set_ylabel("Composition (%)", fontsize=14); ax.set_ylim(0, 105); ax.set_xlim(-0.55, len(CELLS) - 0.3)
ax.tick_params(axis="y", labelsize=12)
h, l = ax.get_legend_handles_labels()
ax.legend(h[::-1], l[::-1], title="Locale", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=12,
          title_fontsize=12)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout(); save(fig, "Fig5c_locale_composition_stacked")

# ---------------- 5e: Sankey of the dominant locale, per locale
L = {}
for c in CELLS:
    t = fc.load_locale(c)[["chr", "start", "dom5"]]
    L[c] = t.rename(columns={"dom5": c})
M = L[CELLS[0]]
for c in CELLS[1:]:
    M = M.merge(L[c], on=["chr", "start"], validate="1:1")
N = len(M); assert N == 13315
frac = pd.DataFrame({c: M[c].value_counts().reindex(ORDER).fillna(0) / N for c in CELLS})
assert np.allclose(frac.sum(), 1.0)
frac.rename(index=DISP).rename_axis("locale").reset_index().to_csv(SRC / "Fig5e_per_locale_node_fractions.tsv", sep="\t",
                                                                    index=False, float_format="%.5f")
trans = []
for a, b in zip(CELLS[:-1], CELLS[1:]):
    ct = pd.crosstab(M[a], M[b]).reindex(index=ORDER, columns=ORDER).fillna(0).astype(int)
    for g in ORDER:
        for h in ORDER:
            if ct.loc[g, h]:
                trans.append(dict(stage_from=a, stage_to=b, locale_from=DISP[g], locale_to=DISP[h], n_bins=int(ct.loc[g, h]),
                                  g=g, h=h))
    print(f"{a} -> {b}: {100 * np.trace(ct.to_numpy()) / N:.1f} % of bins keep their dominant locale")
TR = pd.DataFrame(trans); TR.drop(columns=["g", "h"]).to_csv(SRC / "Fig5e_per_locale_transitions.tsv", sep="\t", index=False)
print((100 * frac).round(2).to_string())

GAP, W = 0.012, 0.2
xs = np.arange(len(CELLS)) * 1.0
top = {}
for c in CELLS:
    y = 1.0
    for g in ORDER:
        top[(c, g)] = y; y -= frac.loc[g, c] + (GAP if frac.loc[g, c] > 0 else 0)


def band(ax, x0, x1, y0a, y0b, y1a, y1b, color):
    t = np.linspace(0, 1, 60); s = (1 - np.cos(np.pi * t)) / 2; xx = x0 + (x1 - x0) * t
    ax.fill_between(xx, y0a + (y1a - y0a) * s, y0b + (y1b - y0b) * s, color=color, alpha=0.35, lw=0)


fig, ax = plt.subplots(figsize=(12, 6.5))
for i, (a, b) in enumerate(zip(CELLS[:-1], CELLS[1:])):
    out_off = {g: top[(a, g)] for g in ORDER}; in_off = {h: top[(b, h)] for h in ORDER}
    sub = TR[(TR.stage_from == a) & (TR.stage_to == b)]
    for g in ORDER:
        for h in ORDER:
            n = sub[(sub.g == g) & (sub.h == h)].n_bins.sum()
            if n == 0:
                continue
            w = n / N
            band(ax, xs[i] + W / 2, xs[i + 1] - W / 2, out_off[g], out_off[g] - w, in_off[h], in_off[h] - w, COL[g])
            out_off[g] -= w; in_off[h] -= w
for j, c in enumerate(CELLS):
    for g in ORDER:
        f = frac.loc[g, c]
        if f <= 0:
            continue
        ax.add_patch(Rectangle((xs[j] - W / 2, top[(c, g)] - f), W, f, facecolor=COL[g], edgecolor="none"))
        lab = f"{100 * f:.2f}%" if f < 0.01 else f"{100 * f:.1f}%"
        if f >= 0.035:
            ax.text(xs[j], top[(c, g)] - f / 2, lab, ha="center", va="center", fontsize=11, color=txtcol(COL[g]))
        else:                                                        # small node: label at the right in the locale colour
            ax.text(xs[j] + W / 2 + 0.02, top[(c, g)] - f / 2, lab, ha="left", va="center", fontsize=10, color=COL[g],
                    fontweight="bold")
    ax.text(xs[j], 1.03, T[c], ha="center", va="bottom", fontsize=15)
ax.set_xlim(xs[0] - 0.3, xs[-1] + 0.35); ax.set_ylim(-0.08, 1.08); ax.axis("off")
handles = [Rectangle((0, 0), 1, 1, color=COL[g]) for g in ORDER]
ax.legend(handles, [DISP[g] for g in ORDER], title="Dominant locale", loc="center left", bbox_to_anchor=(1.0, 0.5),
          frameon=False, fontsize=12, title_fontsize=12)
fig.tight_layout(); save(fig, "Fig5e_dominant_locale_sankey")
