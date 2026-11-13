"""Shared inputs, constants and helpers for the Figure 3 model panels (IGM population models of four differentiation
stages: H1 ESC, DE, HB, mHEP; occupancy 0.6, 200 kb, 1000 diploid structures each).

dist to CT   Signed distance of a locus to its chromosome-territory (CT) core, normalised by the core centroid-to-
             surface distance (dimensionless; positive = inside). CT core = convex hull of the top 30 % beads of each
             chromosome copy ranked by strict-cis contact ratio (2.5 bead radii). Per 200 kb bin the 2 copies x 1000
             structures (2000 values) give Mean, SD and IQR.
locale       Nuclear locale entropy (bits) over speckle / nucleolus / envelope / NL-act / NL-rep (5 / 6 / 5 bead radii).
Eig1         Bulk Hi-C compartment eigenvector of the same stage, 1 Mb (column E1 of Supplementary Table 2).
Bins         Autosomal, non-centromeric 200 kb bins (13,315 per stage). 1 Mb values = mean of these 200 kb statistics.
"""
import os
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import spearmanr

matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "sans-serif",
                            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"                       # per-figure input tables (only those unique to this figure)
SHARED = ROOT.parent / "data"              # repo-level shared inputs: locale/, distCT/, radial/, schic/
# Locale input is switchable (2026-11): data/locale/<LOCALE_VERSION>/<cell>_occ06_locale_N565.tsv. Default "v2_ab4r" = ENV v2
# (env seeds = outermost 20 % of all beads + 2.5 r; NL-act / NL-rep by A/B ratio within radius_factor 2 = 4 r). Any other
# version writes to output_<LOCALE_VERSION>/ instead of output/.
DEFAULT_LOCALE = "v2_ab4r"
LOCALE_VERSION = os.environ.get("LOCALE_VERSION", DEFAULT_LOCALE)
LOC_DIR = SHARED / "locale" / LOCALE_VERSION
if not LOC_DIR.is_dir():
    raise FileNotFoundError(f"LOCALE_VERSION={LOCALE_VERSION}: no locale tables in {LOC_DIR}")
OUT = ROOT / ("output" if LOCALE_VERSION == DEFAULT_LOCALE else f"output_{LOCALE_VERSION}")
# dist-to-CT normalisation is switchable (2026-11-12): DISTCT_NORM = e2c (default, the shipped tables: signed distance /
# core centroid-to-surface distance) | r (signed distance / bead radius) | Req (signed distance / equivalent-sphere radius
# (3 V_hull / 4 pi)^(1/3) of the same core hull). Non-default tables: data/distCT/<norm>/; outputs: <output dir>_distCT-<norm>/.
DEFAULT_DISTCT = "e2c"
DISTCT_NORM = os.environ.get("DISTCT_NORM", DEFAULT_DISTCT)
CT_UNIT = {"e2c": "", "r": " / r", "Req": " / R_eq"}.get(DISTCT_NORM)
if CT_UNIT is None:
    raise ValueError(f"DISTCT_NORM={DISTCT_NORM}: expected e2c, r or Req")
CT_DIR = SHARED / "distCT" / DISTCT_NORM     # default e2c tables in data/distCT/e2c/
if not CT_DIR.is_dir():
    raise FileNotFoundError(f"DISTCT_NORM={DISTCT_NORM}: no dist-to-CT tables in {CT_DIR}")
if DISTCT_NORM != DEFAULT_DISTCT:
    OUT = OUT.parent / f"{OUT.name}_distCT-{DISTCT_NORM}"


def _ct_scale():
    """1.0 for e2c; otherwise median over cells of median SD(dist to CT, norm) / median SD(dist to CT, e2c)."""
    if DISTCT_NORM == DEFAULT_DISTCT:
        return 1.0
    import pandas as _pd, numpy as _np
    r = [_pd.read_csv(f, sep="\t").sd_signed_dist.median() / _pd.read_csv(SHARED / "distCT" / DEFAULT_DISTCT / f.name, sep="\t").sd_signed_dist.median()
         for f in sorted(CT_DIR.glob("*_occ06_distCT_top30_2p5r_stats.tsv"))]
    return float(_np.median(r))


CT_SCALE = _ct_scale()          # multiply fixed dist-to-CT axis limits by this (1.0 for the default e2c tables)
FIG_DIR, SRC_DIR = OUT / "figures", OUT / "source_data"
FIG_DIR.mkdir(parents=True, exist_ok=True); SRC_DIR.mkdir(parents=True, exist_ok=True)

CELLS = ["h1", "DE", "HB", "mHEP"]
CELL_TITLE = {"h1": "ES", "DE": "DE", "HB": "HB", "mHEP": "HLC"}               # stage names used in the manuscript
E1_SHEET = {"h1": "H1", "DE": "DE", "HB": "HB", "mHEP": "HLC"}
PLATE_XLSX = SHARED / "schic" / "Table 2-Table 3-Cis_total annotation in plate-based single-cell Hi-C.xlsx"
BWR = mcolors.LinearSegmentedColormap.from_list("blue_white_red", ["blue", "white", "red"])
V_LABEL = {"mean": f"Mean(dist to CT{CT_UNIT})", "iqr": f"IQR(dist to CT{CT_UNIT})", "sd": f"SD(dist to CT{CT_UNIT})"}
AUTOS = [f"chr{i}" for i in range(1, 23)]


def ct_tsv(cell):
    return CT_DIR / f"{cell}_occ06_distCT_top30_2p5r_stats.tsv"


def locale_tsv(cell):
    return LOC_DIR / f"{cell}_occ06_locale_N565.tsv"


def load_200kb(cell, keep_cen=False):
    """Per-200kb-bin table of one stage: chr start end label mean sd iqr entropy_bit (autosomes; cen removed unless keep_cen)."""
    ct = pd.read_csv(ct_tsv(cell), sep="\t").rename(columns={"mean_signed_dist": "mean", "sd_signed_dist": "sd",
                                                            "iqr_signed_dist": "iqr"})
    lc = pd.read_csv(locale_tsv(cell), sep="\t")
    assert len(ct) == len(lc) == 15453, (cell, len(ct), len(lc))
    d = ct[["chr", "start", "end", "mean", "sd", "iqr"]].merge(lc[["chr", "start", "label", "entropy_bit"]],
                                                                on=["chr", "start"], how="inner", validate="1:1")
    assert len(d) == 15453, f"{cell}: dist-to-CT and locale bin grids do not align"
    d = d[d.chr.isin(AUTOS)].copy()
    if not keep_cen:
        d = d[d.label != "cen"].copy()
        assert len(d) == 13315, (cell, len(d))
    return d.reset_index(drop=True)


def load_e1_1mb(cell):
    t = pd.read_excel(PLATE_XLSX, sheet_name=f"{E1_SHEET[cell]}-platebased-HiC-batch1").rename(columns={"chrom": "chr"})
    t["chr"] = t["chr"].astype(str)
    assert (t.start % 1_000_000 == 0).all()
    return t[["chr", "start", "end", "E1"]].rename(columns={"E1": "E1_1mb"})


def to_1mb(d, cols):
    """Mean of the (non-centromeric) 200 kb bins inside each 1 Mb bin (start // 1e6)."""
    g = d.assign(start=(d.start // 1_000_000) * 1_000_000)
    out = g.groupby(["chr", "start"])[cols].mean()
    out["n_200kb"] = g.groupby(["chr", "start"]).size()
    return out.reset_index()


TEXT_POS = {"ul": (0.05, 0.95, "left", "top"), "ur": (0.95, 0.95, "right", "top"), "ll": (0.05, 0.05, "left", "bottom"),
            "ml": (0.05, 0.52, "left", "top"), "lr": (0.95, 0.05, "right", "bottom")}


def spear_text(ax, x, y, fontsize=10, n=True, loc="ul"):
    r, p = spearmanr(x, y); tx, ty, ha, va = TEXT_POS[loc]
    ax.text(tx, ty, f"Spearman r = {r:.2f}\n" + (f"p-value = {p:.2e}" if p > 0 else "p < 1e-300") + (f"\nn = {len(x)}" if n else ""),
            transform=ax.transAxes, fontsize=fontsize, ha=ha, va=va,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2))
    return r, p


def prange(v, lo=2, hi=98):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return tuple(np.percentile(v, [lo, hi]))


def save(fig, name, dpi=300):
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote output/figures/{name}.{{pdf,png}}")


def save_table(df, name, float_format="%.6g"):
    df.to_csv(SRC_DIR / name, sep="\t", index=False, float_format=float_format)
    print(f"wrote output/source_data/{name}")
