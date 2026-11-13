"""Shared inputs, constants and helpers for the Extended Data Figure 6 panels (IGM population models of four
differentiation stages: ES (H1), DE, HB, HLC (mHEP); occupancy 0.6, 200 kb, 1000 diploid structures each).

Definitions
-----------
dist to CT   signed distance of a locus to its chromosome-territory core / core centroid-to-surface distance
             (dimensionless, > 0 inside); core = convex hull of the top 30 % beads of each chromosome copy by strict-cis
             contact ratio (2.5 bead radii); SD / IQR over 2 copies x 1000 structures per 200 kb bin
RAD          mean normalised radial position, sqrt((x/a)^2 + (y/b)^2 + (z/c)^2) of the model's ellipsoid (0 centre, 1 envelope)
locale       per bead and structure: speckle (spe), nucleolus (nuc), nuclear envelope (env) within 5 / 6 / 5 bead radii of
             exactly one of them; near none -> NL-act / NL-rep by A/B ratio > 0.5 or not; locale entropy = Shannon entropy
             (bits) of the five exclusive class frequencies; dominant locale = the most frequent of the five
locale strength  top1 / top2 of the five (renormalised) class frequencies
active / repressive  dominant locale spe or NL-act / env, NL-rep or nuc
Bins         autosomal non-centromeric 200 kb bins (13,315)
"""
import os
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
FIG_DIR, SRC_DIR = OUT / "figures", OUT / "source_data"
FIG_DIR.mkdir(parents=True, exist_ok=True); SRC_DIR.mkdir(parents=True, exist_ok=True)

CELLS = ["h1", "DE", "HB", "mHEP"]
CELL_TITLE = {"h1": "ES", "DE": "DE", "HB": "HB", "mHEP": "HLC"}
CELL_COLOR = {"h1": "tab:blue", "DE": "tab:orange", "HB": "tab:green", "mHEP": "tab:purple"}
AUTOS = [f"chr{i}" for i in range(1, 23)]
V_LABEL = {"mean": f"Mean(dist to CT{CT_UNIT})", "iqr": f"IQR(dist to CT{CT_UNIT})", "sd": f"SD(dist to CT{CT_UNIT})"}
LOC5 = ["spe", "nuc", "env", "act", "rep"]
P5 = [f"P_{k}" for k in LOC5]
LOC_COLORS = {"spe": "#8b1e3f", "nuc": "#f4a261", "env": "#7fa9c9", "act": "#e57373", "rep": "#9b59b6"}
LOC_DISP = {"spe": "spe", "nuc": "nuc", "env": "env", "act": "NL-act", "rep": "NL-rep"}
ACTIVE, REPRESSIVE = ("spe", "act"), ("env", "rep", "nuc")
CLASS_COLOR = {"all": "#555555", "active": "#c0392b", "repressive": "#4a5fa8"}


def load_stage(cell):
    """Autosomal non-centromeric 200 kb bins of one stage model: chr start mean sd iqr norm_rad entropy_bit P_* dom strength."""
    lc = pd.read_csv(LOC_DIR / f"{cell}_occ06_locale_N565.tsv", sep="\t")
    ct = pd.read_csv(CT_DIR / f"{cell}_occ06_distCT_top30_2p5r_stats.tsv", sep="\t").rename(
        columns={"mean_signed_dist": "mean", "sd_signed_dist": "sd", "iqr_signed_dist": "iqr"})
    rd = pd.read_csv(SHARED / "radial" / f"{cell}_occ06_radial_200kb.bed", sep="\t")
    assert len(lc) == len(ct) == len(rd) == 15453
    d = lc[["chr", "start", "label", "entropy_bit"] + P5].merge(ct[["chr", "start", "mean", "sd", "iqr"]], on=["chr", "start"],
                                                                 validate="1:1") \
        .merge(rd[["chr", "start", "norm_rad"]], on=["chr", "start"], validate="1:1")
    d = d[d.chr.isin(AUTOS) & (d.label != "cen")].drop(columns="label").reset_index(drop=True)
    assert len(d) == 13315 and d[P5 + ["entropy_bit", "mean", "sd", "iqr", "norm_rad"]].notna().all().all(), cell
    P = d[P5].to_numpy(float); Pn = P / P.sum(1, keepdims=True); srt = np.sort(Pn, axis=1)[:, ::-1]
    d["dom"] = np.array(LOC5)[np.argmax(Pn, axis=1)]
    d["strength"] = srt[:, 0] / np.clip(srt[:, 1], 1e-12, None)
    d["cls"] = np.where(d.dom.isin(ACTIVE), "active", "repressive")
    return d


def save(fig, name, dpi=300):
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote output/figures/{name}.{{pdf,png}}")


def save_table(df, name, float_format="%.6g"):
    df.to_csv(SRC_DIR / name, sep="\t", index=False, float_format=float_format)
    print(f"wrote output/source_data/{name}")


# ---- structure features (LaD, SpD, TransAB) and RNA-seq
SF_COLS = {"LaD": "LaD_nm", "SpD": "SpD_nm", "TransAB": "TransAB"}
STAGE_RNA = {"h1": "ESC", "DE": "DE", "HB": "HB", "mHEP": "mHEP"}       # DESeq2 file names / normcounts columns
NORM_COL = {"h1": "h1", "DE": "DE", "HB": "HB", "mHEP": "mHEP"}


def load_sf_z(cell):
    """Per autosomal non-centromeric 200 kb bin: LaD_nm, SpD_nm, TransAB and LaD_z, SpD_z (z-scored within the stage over
    these 13,315 bins, ddof = 1)."""
    sf = pd.read_csv(DATA / f"{cell}_occ06_sf_LaD_SpD_TransAB_200kb.tsv", sep="\t").drop(columns="label")
    lc = pd.read_csv(LOC_DIR / f"{cell}_occ06_locale_N565.tsv", sep="\t", usecols=["chr", "start", "label"])
    d = lc.merge(sf, on=["chr", "start"], validate="1:1")
    d = d[d.chr.isin(AUTOS) & (d.label != "cen")].drop(columns="label").reset_index(drop=True)
    assert len(d) == 13315 and d[["LaD_nm", "SpD_nm", "TransAB"]].notna().all().all(), cell
    for k in ("LaD", "SpD"):
        v = d[SF_COLS[k]]; d[f"{k}_z"] = (v - v.mean()) / v.std(ddof=1)
    return d
