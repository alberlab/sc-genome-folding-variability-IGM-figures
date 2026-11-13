"""Shared inputs, constants and helpers for the Extended Data Figure 2 panels (IGM population models, occupancy 0.6,
200 kb, 1000 diploid structures each; panels a, b, e use the H1 ESC model, panel c the four stage models).

Definitions used throughout
---------------------------
dist to CT   Signed distance of a locus to its chromosome-territory (CT) core, normalised by the core centroid-to-
             surface distance (dimensionless; positive = inside the CT core). CT core = convex hull of the top 30 % beads
             of each chromosome copy ranked by strict-cis contact ratio (contacts within 2.5 bead radii, same chromosome
             and same copy). Per 200 kb bin, 2 copies x 1000 structures = 2000 values are pooled; the input table holds
             their mean, SD and IQR (75th - 25th percentile).
RAD          Mean normalised radial position of the bin, sqrt((x/a)^2 + (y/b)^2 + (z/c)^2) with a, b, c = the semi-axes of
             the model's ellipsoidal nucleus (0 = centre, 1 = envelope), 2 copies x 1000 structures.
locale       Nuclear locale of each bead in each structure: speckle (spe), nucleolus (nuc), nuclear envelope (env) within
             5 / 6 / 5 bead radii, otherwise NL-act / NL-rep by A/B ratio > 0.5 or not; beads near two or three landmarks
             form the overlap classes. Per-bin frequencies P_* = fraction of the valid observations (2 copies x 1000).
Bins         Autosomal 200 kb bins that are not centromeric (label != "cen"): 13,315 bins.
1 Mb values  Mean of the non-centromeric 200 kb bin statistics inside each 1 Mb bin (start // 1 Mb).
"""
import os
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.stats import gaussian_kde, spearmanr

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

CT_TSV = CT_DIR / "h1_occ06_distCT_top30_2p5r_stats.tsv"
RAD_BED = SHARED / "radial" / "h1_occ06_radial_200kb.bed"
LOCI_NPZ = CT_DIR / "h1_occ06_chr2_3loci_distCT.npz"
PLATE_XLSX = SHARED / "schic" / "Table 2-Table 3-Cis_total annotation in plate-based single-cell Hi-C.xlsx"
PLATE_SHEET = "H1-platebased-HiC-batch1"
DROPLET_XLSX = SHARED / "schic" / "Table 3-Cis_total annotation in droplet-based single-cell Hi-C.xlsx"
DROPLET_SHEET = "H1-dropletHiC >0.3kb filter"

CELLS = ["h1", "DE", "HB", "mHEP"]
CELL_TITLE = {"h1": "ES", "DE": "DE", "HB": "HB", "mHEP": "HLC"}
AUTOS = [f"chr{i}" for i in range(1, 23)]
V_LABEL = {"mean": f"Mean(dist to CT{CT_UNIT})", "iqr": f"IQR(dist to CT{CT_UNIT})", "sd": f"SD(dist to CT{CT_UNIT})"}
LOC5 = ["spe", "nuc", "env", "act", "rep"]                          # exclusive classes, order of the P_* columns
P5 = [f"P_{k}" for k in LOC5]
LOC_COLORS = {"spe": "#8b1e3f", "nuc": "#f4a261", "env": "#7fa9c9", "act": "#e57373", "rep": "#9b59b6", "others": "#b0b0b0"}
LOC_DISP = {"spe": "spe", "nuc": "nuc", "env": "env", "act": "NL-act", "rep": "NL-rep", "others": "other"}


def locale_table(cell):
    """Full per-bin locale table of one stage model (15,453 bins incl. chrX/Y and centromeric rows)."""
    return pd.read_csv(LOC_DIR / f"{cell}_occ06_locale_N565.tsv", sep="\t")


def load_h1_200kb():
    """H1 model, autosomal non-cen 200 kb bins: chr start end mean sd iqr norm_rad."""
    ct = pd.read_csv(CT_TSV, sep="\t").rename(columns={"mean_signed_dist": "mean", "sd_signed_dist": "sd",
                                                       "iqr_signed_dist": "iqr"})
    rd = pd.read_csv(RAD_BED, sep="\t")
    lc = locale_table("h1")
    assert len(ct) == len(rd) == len(lc) == 15453, (len(ct), len(rd), len(lc))
    d = ct[["chr", "start", "end", "mean", "sd", "iqr"]].merge(rd[["chr", "start", "norm_rad"]], on=["chr", "start"],
                                                               validate="1:1") \
        .merge(lc[["chr", "start", "label"]], on=["chr", "start"], validate="1:1")
    d = d[d.chr.isin(AUTOS) & (d.label != "cen")].drop(columns="label").reset_index(drop=True)
    assert len(d) == 13315 and np.isfinite(d[["mean", "sd", "iqr", "norm_rad"]].to_numpy(float)).all(), len(d)
    return d


def to_1mb(d, cols):
    """Mean of the 200 kb bins (already cen-filtered) inside each 1 Mb bin (start // 1e6)."""
    g = d.assign(start=(d.start // 1_000_000) * 1_000_000)
    out = g.groupby(["chr", "start"])[cols].mean()
    out["n_200kb"] = g.groupby(["chr", "start"]).size()
    return out.reset_index()


def read_schic_table(xlsx, sheet, cols):
    t = pd.read_excel(xlsx, sheet_name=sheet).rename(columns={"chrom": "chr"})
    t["chr"] = t["chr"].astype(str)
    assert (t.start % 1_000_000 == 0).all()
    return t[["chr", "start", "end"] + cols]


def density_panel(ax, x, y, s=3):
    """Per-point gaussian_kde density, normalised to the panel maximum, jet colormap, dense points on top."""
    xy = np.vstack([x, y]); z = gaussian_kde(xy)(xy); o = np.argsort(z)
    ax.scatter(x[o], y[o], c=z[o] / z.max(), cmap="jet", vmin=0, vmax=1, s=s, edgecolor="none", rasterized=True)


def density_colorbar(fig, axs, fs=16):
    cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap="jet"), ax=axs, fraction=0.025, pad=0.02)
    cb.set_label("relative point density (KDE)", fontsize=fs); cb.ax.tick_params(labelsize=fs - 2)
    return cb


def spear_box(ax, x, y, fontsize=19, n=True, pos=(0.04, 0.95, "left", "top")):
    r, p = spearmanr(x, y); tx, ty, ha, va = pos
    ax.text(tx, ty, f"Spearman r = {r:.2f}" + (f"\nn = {len(x)}" if n else ""), transform=ax.transAxes,
            fontsize=fontsize, ha=ha, va=va, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    return r, p


def axis_range(v, lo=0.5, hi=99.5, pad=0.04):
    a, b = np.percentile(np.asarray(v, float), [lo, hi]); w = (b - a) * pad
    return a - w, b + w


def save(fig, name, dpi=300):
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote output/figures/{name}.{{pdf,png}}")


def save_table(df, name, float_format="%.6g"):
    df.to_csv(SRC_DIR / name, sep="\t", index=False, float_format=float_format)
    print(f"wrote output/source_data/{name}")
