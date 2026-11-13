"""Shared inputs, constants and helpers for the Figure 5 plots (IGM population models of four differentiation stages:
ES (H1), DE, HB, HLC (mHEP); occupancy 0.6, 200 kb, 1000 diploid structures each).
Bins: autosomal, non-centromeric 200 kb bins (13,315 per stage).
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
# colours of the original five-stage figure (ES, DE, HB, [iHLC red], HLC), kept fixed so HLC stays purple without iHLC
CELL_COLOR = {"h1": "tab:blue", "DE": "tab:orange", "HB": "tab:green", "mHEP": "tab:purple"}
AUTOS = [f"chr{i}" for i in range(1, 23)]
V_LABEL = {"mean": f"Mean(dist to CT{CT_UNIT})", "iqr": f"IQR(dist to CT{CT_UNIT})", "sd": f"SD(dist to CT{CT_UNIT})"}


def domain_bins(cell):
    """chr/start of the autosomal non-centromeric 200 kb bins (label from the locale table)."""
    lc = pd.read_csv(LOC_DIR / f"{cell}_occ06_locale_N565.tsv", sep="\t", usecols=["chr", "start", "label"])
    d = lc[lc.chr.isin(AUTOS) & (lc.label != "cen")][["chr", "start"]]
    assert len(d) == 13315, (cell, len(d))
    return d


def read_ct(path):
    return pd.read_csv(path, sep="\t").rename(columns={"mean_signed_dist": "mean", "sd_signed_dist": "sd",
                                                       "iqr_signed_dist": "iqr"})[["chr", "start", "mean", "sd", "iqr"]]


def read_radial(path):
    return pd.read_csv(path, sep="\t")[["chr", "start", "norm_rad"]]


def radial_deciles(d, n=10):
    """Equal-count shells of the per-bin mean radial position within one model (1 = interior ... n = periphery)."""
    return pd.qcut(d.norm_rad.rank(method="first"), q=n, labels=range(1, n + 1)).astype(int)


def shell_summary(d, col, shell_col="shell"):
    g = d.groupby(shell_col)[col]
    s = g.agg(["mean", "std", "count"]).reset_index()
    s["ci95"] = 1.96 * s["std"] / np.sqrt(s["count"])
    return s


def save(fig, name, dpi=300):
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote output/figures/{name}.{{pdf,png}}")


def save_table(df, name, float_format="%.6g"):
    df.to_csv(SRC_DIR / name, sep="\t", index=False, float_format=float_format)
    print(f"wrote output/source_data/{name}")


# ---- nuclear locale (N565 standard)
LOC5 = ["spe", "nuc", "env", "act", "rep"]                          # exclusive classes, order of the P_* columns
P5 = [f"P_{k}" for k in LOC5]
P_ALL = P5 + ["P_spe_nuc", "P_spe_env", "P_nuc_env", "P_triple"]
LOC_COLORS = {"spe": "#8b1e3f", "nuc": "#f4a261", "env": "#7fa9c9", "act": "#e57373", "rep": "#9b59b6", "others": "#b0b0b0"}
LOC_DISP = {"spe": "spe", "nuc": "nuc", "env": "env", "act": "NL-act", "rep": "NL-rep", "others": "other"}
# active / inactive grouping of the five exclusive locales (Fig. 5c, 5e)
GROUP_OF = {"spe": "active", "act": "active", "env": "inactive", "nuc": "inactive", "rep": "inactive"}
GROUP_COLORS = {"active": "#c0392b", "inactive": "#5b7fa6", "other": "#b0b0b0"}
GROUP_DISP = {"active": "Active (spe + NL-act)", "inactive": "Inactive (env + nuc + NL-rep)", "other": "other"}


def load_locale(cell):
    """Autosomal non-centromeric 200 kb bins of one model: chr start n_obs n_valid entropy_bit P_* dom5."""
    t = pd.read_csv(LOC_DIR / f"{cell}_occ06_locale_N565.tsv", sep="\t")
    t = t[t.chr.isin(AUTOS) & (t.label != "cen")].copy()
    assert len(t) == 13315, (cell, len(t))
    t["dom5"] = np.array(LOC5)[np.argmax(t[P5].to_numpy(float), axis=1)]   # dominant = argmax of the 5 exclusive classes
    return t.reset_index(drop=True)
