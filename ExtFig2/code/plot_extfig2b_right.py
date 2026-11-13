#!/usr/bin/env python3
"""Extended Data Fig. 2b, right panel: SD(RAD) vs SD(dist to CT), genome-wide, per 200 kb bin (13,315 autosomal
non-centromeric bins), ES model (from ~/plots/sdrad_vs_sddistct_4stage/code/plot_sdrad_vs_sddistct_density.py, ES version;
the axis ranges are still taken from the pooled values of the four stage models, so the panel is identical to the earlier
one). Values from data/SDrad_vs_SDdistCT_4stage_perbin.tsv (columns stage chr start sd sd_rad ...; SD(RAD) = pooled SD of
the ellipsoid-normalised radial position over 2 copies x 1000 structures, as SD(dist to CT)). Style as Extended Data Fig. 2a (SD(dist to CT) vs SD(cisratio)): per-point gaussian
KDE density normalised to the panel maximum, jet colormap, dense points on top; axes = 0.5-99.5 % range of the pooled
values (shared by the panels; points outside are not drawn but are in the statistics); Spearman r and n per panel.
Outputs: output/figures/ExtFig2b_right_SDRAD_vs_SDdistCT_ES.{png,pdf}, output/source_data/ExtFig2b_right_stats.tsv
"""
import os
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.stats import gaussian_kde, spearmanr

plt.rcParams.update({"pdf.fonttype": 42, "font.family": "Arial"})
ROOT = Path(__file__).resolve().parents[1]
LOCALE_VERSION = os.environ.get("LOCALE_VERSION", "v2_ab4r")            # no locale input here; output dir follows the folder
OUT = ROOT / ("output" if LOCALE_VERSION == "v2_ab4r" else f"output_{LOCALE_VERSION}")
# dist-to-CT normalisation is switchable (2026-11-12): DISTCT_NORM = e2c (default, the shipped tables: signed distance /
# core centroid-to-surface distance) | r (signed distance / bead radius) | Req (signed distance / equivalent-sphere radius
# (3 V_hull / 4 pi)^(1/3) of the same core hull). Non-default tables: data/distCT/<norm>/; outputs: <output dir>_distCT-<norm>/.
DEFAULT_DISTCT = "e2c"
DISTCT_NORM = os.environ.get("DISTCT_NORM", DEFAULT_DISTCT)
CT_UNIT = {"e2c": "", "r": " / r", "Req": " / R_eq"}.get(DISTCT_NORM)
if CT_UNIT is None:
    raise ValueError(f"DISTCT_NORM={DISTCT_NORM}: expected e2c, r or Req")
CT_DIR = ROOT.parent / "data" / "distCT" / DISTCT_NORM     # shared tables, default data/distCT/e2c/
if not CT_DIR.is_dir():
    raise FileNotFoundError(f"DISTCT_NORM={DISTCT_NORM}: no dist-to-CT tables in {CT_DIR}")
if DISTCT_NORM != DEFAULT_DISTCT:
    OUT = OUT.parent / f"{OUT.name}_distCT-{DISTCT_NORM}"
FIG_DIR, SRC_DIR = OUT / "figures", OUT / "source_data"
FIG_DIR.mkdir(parents=True, exist_ok=True); SRC_DIR.mkdir(parents=True, exist_ok=True)
T = pd.read_csv(CT_DIR / "SDrad_vs_SDdistCT_4stage_perbin.tsv", sep="\t")
STAGES = ["ES", "DE", "HB", "HLC"]
FS_T, FS_L, FS_K, FS_B = 22, 20, 16, 17


def axis_range(v, lo=0.5, hi=99.5, pad=0.04):
    a, b = np.percentile(np.asarray(v, float), [lo, hi]); w = (b - a) * pad
    return a - w, b + w


XL, YL = axis_range(T["sd"]), axis_range(T.sd_rad)
Z = {}
for st in STAGES:
    g = T[T.stage == st]; x, y = g["sd"].to_numpy(), g.sd_rad.to_numpy(); assert len(g) == 13315
    xy = np.vstack([x, y]); z = gaussian_kde(xy)(xy); Z[st] = (x, y, z / z.max())


def panel(ax, st):
    x, y, z = Z[st]; o = np.argsort(z)
    ax.scatter(x[o], y[o], c=z[o], cmap="jet", vmin=0, vmax=1, s=4, edgecolor="none", rasterized=True)
    r = spearmanr(x, y)[0]
    ax.text(0.96, 0.04, f"Spearman r = {r:.2f}\nn = {len(x):,}", transform=ax.transAxes, fontsize=FS_B, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    ax.set_xlim(XL); ax.set_ylim(YL); ax.tick_params(labelsize=FS_K)
    ax.set_title(st, fontsize=FS_T, fontweight="bold"); ax.set_xlabel(f"SD(dist to CT{CT_UNIT})", fontsize=FS_L)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return r


rows = []
for tag, sts, size in (("ES", ["ES"], (6.4, 5.4)),):                     # ES panel only (the 1 x 4 version is not used)
    fig, axs = plt.subplots(1, len(sts), figsize=size, sharex=True, sharey=True, squeeze=False)
    for ax, st in zip(axs[0], sts):
        r = panel(ax, st); print(tag, st, f"Spearman {r:.3f}")
        rows.append(dict(stage=st, n_200kb_bins=len(Z[st][0]), spearman_SDdistCT_SDRAD=r))
    axs[0, 0].set_ylabel("SD(RAD)", fontsize=FS_L)
    cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap="jet"), ax=axs[0].tolist(), fraction=0.025 if len(sts) > 1 else 0.05, pad=0.02)
    cb.set_label("relative point density (KDE)", fontsize=FS_K); cb.ax.tick_params(labelsize=FS_K - 2)
    for e in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"ExtFig2b_right_SDRAD_vs_SDdistCT_{tag}.{e}", dpi=300, bbox_inches="tight")
    plt.close(fig)
pd.DataFrame(rows).to_csv(SRC_DIR / "ExtFig2b_right_stats.tsv", sep="\t", index=False, float_format="%.4g")
print("wrote output/figures/ExtFig2b_right_SDRAD_vs_SDdistCT_ES.{png,pdf}, output/source_data/ExtFig2b_right_stats.tsv")
