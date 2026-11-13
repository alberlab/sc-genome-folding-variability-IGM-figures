#!/usr/bin/env python3
"""Extended Data Fig. 2d: ES model, nuclear locale entropy vs SD(dist to CT) over all autosomes, split into two locale
classes by the bin's dominant locale (argmax of the five exclusive frequencies): active = spe or NL-act dominant (left),
repressive = env, NL-rep or nuc dominant (right). 1 x 2.

Drawing = the 2-class SD mode of ~/plots/entropy_CV_vs_SDCT_diagnosis/code/plot_entropy_vs_IQR_by_locale_strength.py
(arguments `sd 2class <locale>`; the manuscript panel is its ES column): every bin light grey; the class's bins in the class
colour with per-point opacity = clip((strength - 1) / (CAP - 1), AMIN, 1), strength = top1 / top2 of the renormalised
5-class profile (CAP = 4, AMIN = 0.08), drawn weakest first; n, class composition and the within-class Spearman r at the
top right; one strength colour bar per panel. Axis ranges = pooled 0.1-99.9 % of the four stage models + 3 % (as the
original 2 x 4 figure, so the ES panels keep their axes). Statistics for all four stages go to the source data.
Changes vs that script: ES only, 1 x 2 (class = column instead of row); locale = ENV v2 (the manuscript panel used ENV v1,
env seeds = outermost 20 % of the non-gap beads: ES active r = -0.1362 (n 6,524), repressive r = +0.3733 (n 6,791);
reproduced exactly with the v1 tables in ../data/locale/<v1 name>/ and LOCALE_VERSION=<v1 name>).
CT core = top 30 % strict-cis cisratio, 2.5 r (same table as every other panel).
Data: ../../data/distCT/e2c/<cell>_occ06_distCT_top30_2p5r_stats.tsv and ../../data/locale/<LOCALE_VERSION>/<cell>_occ06_locale_N565.tsv, 4 stages
(default LOCALE_VERSION v2_ab4r = ENV v2; other versions write to output_<LOCALE_VERSION>/).
Outputs: output/figures/ExtFig2d_ES_entropy_vs_SDdistCT_active_repressive.*, output/source_data/ExtFig2d_stats.tsv
"""
import os
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import spearmanr

matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT.parent / "data"                                           # repo-level shared four-stage tables
LOCALE_VERSION = os.environ.get("LOCALE_VERSION", "v2_ab4r")            # locale tables: DATA/locale/<LOCALE_VERSION>/
LOC_DIR = DATA / "locale" / LOCALE_VERSION
if not LOC_DIR.is_dir():
    raise FileNotFoundError(f"LOCALE_VERSION={LOCALE_VERSION}: no locale tables in {LOC_DIR}")
OUT = ROOT / ("output" if LOCALE_VERSION == "v2_ab4r" else f"output_{LOCALE_VERSION}")
# dist-to-CT normalisation is switchable (2026-11-12): DISTCT_NORM = e2c (default, the shipped tables: signed distance /
# core centroid-to-surface distance) | r (signed distance / bead radius) | Req (signed distance / equivalent-sphere radius
# (3 V_hull / 4 pi)^(1/3) of the same core hull). Non-default tables: data/distCT/<norm>/; outputs: <output dir>_distCT-<norm>/.
DEFAULT_DISTCT = "e2c"
DISTCT_NORM = os.environ.get("DISTCT_NORM", DEFAULT_DISTCT)
CT_UNIT = {"e2c": "", "r": " / r", "Req": " / R_eq"}.get(DISTCT_NORM)
if CT_UNIT is None:
    raise ValueError(f"DISTCT_NORM={DISTCT_NORM}: expected e2c, r or Req")
CT_DIR = DATA / "distCT" / DISTCT_NORM                               # default data/distCT/e2c/
if not CT_DIR.is_dir():
    raise FileNotFoundError(f"DISTCT_NORM={DISTCT_NORM}: no dist-to-CT tables in {CT_DIR}")
if DISTCT_NORM != DEFAULT_DISTCT:
    OUT = OUT.parent / f"{OUT.name}_distCT-{DISTCT_NORM}"
FIG_DIR, SRC_DIR = OUT / "figures", OUT / "source_data"
FIG_DIR.mkdir(parents=True, exist_ok=True); SRC_DIR.mkdir(parents=True, exist_ok=True)

CELLS = ["h1", "DE", "HB", "mHEP"]; SHOW = "h1"; CELL_TITLE = {"h1": "ES", "DE": "DE", "HB": "HB", "mHEP": "HLC"}
P5 = ["P_spe", "P_nuc", "P_env", "P_act", "P_rep"]; NAMES5 = np.array(["spe", "nuc", "env", "act", "rep"])
CLASSES = [("active", "active"), ("repressive", "repressive")]
MEMBERS = {"active": ["spe", "act"], "repressive": ["env", "rep", "nuc"]}
DISP5 = {"spe": "spe", "act": "NL-act", "env": "env", "rep": "NL-rep", "nuc": "nuc"}
LC = {"active": "#c0392b", "repressive": "#4a5fa8"}
GREY, CAP, AMIN = "#d9d9d9", 4.0, 0.08
FS_T, FS_L, FS_K, FS_B = 22, 17, 13, 13
PW, PH = 4.6, 4.3


def load(cell):
    lo = pd.read_csv(LOC_DIR / f"{cell}_occ06_locale_N565.tsv", sep="\t", usecols=["chr", "start", "label", "entropy_bit"] + P5)
    ct = pd.read_csv(CT_DIR / f"{cell}_occ06_distCT_top30_2p5r_stats.tsv", sep="\t", usecols=["chr", "start", "iqr_signed_dist",
                                                                                            "sd_signed_dist"])
    d = lo.merge(ct, on=["chr", "start"], validate="1:1")
    d = d[d.chr.isin([f"chr{i}" for i in range(1, 23)]) & (d.label.astype(str) != "cen")].dropna(
        subset=P5 + ["entropy_bit", "iqr_signed_dist"])
    P = d[P5].to_numpy(float); rs = P.sum(axis=1); assert (rs > 0).all()
    Pn = P / rs[:, None]; srt = np.sort(Pn, axis=1)[:, ::-1]
    d["dom"] = NAMES5[np.argmax(Pn, axis=1)]
    d["strength"] = srt[:, 0] / np.clip(srt[:, 1], 1e-12, None)
    assert len(d) == 13315, (cell, len(d))
    return d


def alpha_cmap(rgb, n=256):
    t = np.linspace(0, 1, n); a = np.clip(t, AMIN, 1)
    return mcolors.ListedColormap(np.clip((1 - a)[:, None] * 1.0 + a[:, None] * np.asarray(rgb), 0, 1))


D = {c: load(c) for c in CELLS}
XCOL = "sd_signed_dist"
xa = np.concatenate([D[c][XCOL] for c in CELLS]); ya = np.concatenate([D[c].entropy_bit for c in CELLS])
(xlo, xhi), (ylo, yhi) = np.quantile(xa, [0.001, 0.999]), np.quantile(ya, [0.001, 0.999])
XL = (xlo - 0.03 * (xhi - xlo), xhi + 0.03 * (xhi - xlo)); YL = (ylo - 0.03 * (yhi - ylo), yhi + 0.03 * (yhi - ylo))

rows = []
for cell in CELLS:                                                       # statistics for all four stages
    d = D[cell]; x, y, st = d[XCOL].to_numpy(), d.entropy_bit.to_numpy(), d.strength.to_numpy()
    for key, name in CLASSES:
        hit = d.dom.isin(MEMBERS[key]).to_numpy()
        rows.append(dict(locale_class=name, cell=cell, stage=CELL_TITLE[cell], n=int(hit.sum()),
                         composition=", ".join(f"{DISP5[m]} {int((d.dom == m).sum())}" for m in MEMBERS[key]),
                         spearman_SD_entropy=spearmanr(x[hit], y[hit])[0], median_strength=float(np.median(st[hit])),
                         frac_strength_ge2=float((st[hit] >= 2).mean()), median_entropy=float(np.median(y[hit])),
                         median_SD=float(np.median(x[hit]))))
T = pd.DataFrame(rows)

fig, axs = plt.subplots(1, len(CLASSES), figsize=(PW * len(CLASSES) + 1.2, PH), sharex=True, sharey=True,
                        constrained_layout=True)
d = D[SHOW]; x, y, st = d[XCOL].to_numpy(), d.entropy_bit.to_numpy(), d.strength.to_numpy()
for ci, (key, name) in enumerate(CLASSES):
    ax = axs[ci]; rgb = np.array(mcolors.to_rgb(LC[key])); hit = d.dom.isin(MEMBERS[key]).to_numpy()
    ax.scatter(x[~hit], y[~hit], s=5, c=GREY, alpha=0.25, edgecolors="none", rasterized=True)
    a = np.clip((st[hit] - 1) / (CAP - 1), AMIN, 1); o = np.argsort(a)          # weakest first
    rgba = np.column_stack([np.tile(rgb, (int(hit.sum()), 1)), a])[o]
    ax.scatter(x[hit][o], y[hit][o], s=8, c=rgba, edgecolors="none", rasterized=True)
    r = T[(T.cell == SHOW) & (T.locale_class == name)].iloc[0]
    ax.text(0.97, 0.97, f"{name}: {r.n}\n({r.composition})\nSpearman r = {r.spearman_SD_entropy:+.2f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=FS_B)
    ax.set_xlim(XL); ax.set_ylim(YL); ax.tick_params(labelsize=FS_K, length=3.5, width=0.8)
    ax.set_title(f"{name} ({' + '.join(DISP5[m] for m in MEMBERS[key])} dominant)", fontsize=FS_L, fontweight="bold",
                 color=LC[key])
    ax.set_xlabel(f"SD(dist to CT{CT_UNIT})", fontsize=FS_L)
    if ci == 0:
        ax.set_ylabel("Nuclear locale entropy (bit)", fontsize=FS_L)
    sm = plt.cm.ScalarMappable(norm=mcolors.Normalize(1, CAP), cmap=alpha_cmap(rgb)); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.02, aspect=18)
    cb.set_ticks([1, 2, 3, 4]); cb.ax.set_yticklabels(["1", "2", "3", "≥4"], fontsize=FS_K)
    cb.set_label("strength (top1/top2)", fontsize=FS_K)
fig.suptitle(f"{CELL_TITLE[SHOW]} model, all autosomes (200 kb bins)", fontsize=FS_T - 4, fontweight="bold")
for e in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"ExtFig2d_ES_entropy_vs_SDdistCT_active_repressive.{e}", dpi=200, bbox_inches="tight")
plt.close(fig)
T.to_csv(SRC_DIR / "ExtFig2d_stats.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 200)
print(T[["locale_class", "stage", "n", "spearman_SD_entropy", "median_entropy", "median_SD"]].round(4).to_string(index=False))
print("wrote output/figures/ExtFig2d_ES_entropy_vs_SDdistCT_active_repressive.{png,pdf}, output/source_data/ExtFig2d_stats.tsv")
