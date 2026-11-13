"""Shared inputs, constants and helpers for the Extended Data Figure 5 panels (IGM population models, occupancy 0.6, 200 kb,
1000 diploid structures each).

Models
------
a, b  four differentiation stages: ES (H1), DE, HB, HLC (mHEP)
c     H1 restraint ablation: the production H1 model ("full": intra + inter Hi-C, cis + trans) and 7 models built with
      subsets of the restraints (same geometry, seed and population size); "Unrestrained" = the RelaxInit snapshot
d     H1 Hi-C restraints in four nuclear envelopes: ES or HB ellipsoid shape x ES or HB nuclear volume
      (ES shape & ES volume = the production H1 model)

Definitions
-----------
dist to CT  signed distance of a locus to its chromosome-territory core / core centroid-to-surface distance
            (dimensionless, > 0 inside); core = convex hull of the top 30 % beads of each chromosome copy ranked by
            strict-cis contact ratio (2.5 bead radii); SD / IQR over 2 copies x 1000 structures per 200 kb bin
RAD         mean normalised radial position, sqrt((x/a)^2 + (y/b)^2 + (z/c)^2) of the model's ellipsoid (0 centre, 1 envelope)
Bins        autosomal non-centromeric 200 kb bins (13,315)
"""
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
DATA, OUT = ROOT / "data", ROOT / "output"
FIG_DIR, SRC_DIR = OUT / "figures", OUT / "source_data"
FIG_DIR.mkdir(parents=True, exist_ok=True); SRC_DIR.mkdir(parents=True, exist_ok=True)

CELLS = ["h1", "DE", "HB", "mHEP"]
CELL_TITLE = {"h1": "ES", "DE": "DE", "HB": "HB", "mHEP": "HLC"}
CELL_COLOR = {"h1": "tab:blue", "DE": "tab:orange", "HB": "tab:green", "mHEP": "tab:purple"}
AUTOS = [f"chr{i}" for i in range(1, 23)]
V_LABEL = {"mean": "Mean(dist to CT)", "iqr": "IQR(dist to CT)", "sd": "SD(dist to CT)"}
ABL = [("null_relaxinit", "Unrestrained"), ("damid_only", "DamID"), ("trans_only", "Trans"), ("cis_only", "Cis"),
       ("cis_damid", "Cis+DamID"), ("trans_damid", "Trans+DamID"), ("cis_trans_damid", "Cis+Trans+DamID")]
ENV = [("h1shape_h1V", "ES shape & ES volume"), ("h1shape_HBV", "ES shape & HB volume"),
       ("HBshape_h1V", "HB shape & ES volume"), ("HBshape_HBV", "HB shape & HB volume")]


def domain_bins(cell="h1"):
    """chr/start of the autosomal non-centromeric 200 kb bins (label from the stage's locale table)."""
    lc = pd.read_csv(DATA / f"{cell}_occ06_locale_N565.tsv", sep="\t", usecols=["chr", "start", "label"])
    d = lc[lc.chr.isin(AUTOS) & (lc.label != "cen")][["chr", "start"]]
    assert len(d) == 13315, (cell, len(d))
    return d


def read_ct(path):
    t = pd.read_csv(path, sep="\t"); assert len(t) == 15453, (path, len(t))
    return t.rename(columns={"mean_signed_dist": "mean", "sd_signed_dist": "sd",
                             "iqr_signed_dist": "iqr"})[["chr", "start", "mean", "sd", "iqr"]]


def read_radial(path):
    t = pd.read_csv(path, sep="\t"); assert len(t) == 15453, (path, len(t))
    return t[["chr", "start", "norm_rad"]]


def model_table(ct_path, rad_path):
    """Per-bin mean / sd / iqr / norm_rad of one H1-based model on the 13,315 bins."""
    d = domain_bins("h1").merge(read_ct(ct_path), on=["chr", "start"], validate="1:1") \
        .merge(read_radial(rad_path), on=["chr", "start"], validate="1:1")
    assert len(d) == 13315 and np.isfinite(d[["mean", "sd", "iqr", "norm_rad"]].to_numpy(float)).all()
    return d


def abl_table(cond):
    return model_table(DATA / "ablation" / f"h1_abl_{cond}_distCT_top30_2p5r_stats.tsv",
                       DATA / "ablation" / f"h1_abl_{cond}_radial_200kb.bed")


def env_table(run):
    if run == "h1shape_h1V":                                         # production H1 model = ablation reference "full"
        return abl_table("full")
    return model_table(DATA / "envswap" / f"es_{run}_distCT_top30_2p5r_stats.tsv",
                       DATA / "envswap" / f"es_{run}_radial_200kb.bed")


def radial_deciles(d, n=10):
    """Equal-count shells of the per-bin mean radial position within one model (1 = interior ... n = periphery)."""
    return pd.qcut(d.norm_rad.rank(method="first"), q=n, labels=range(1, n + 1)).astype(int)


def shell_summary(d, col, shell_col="shell"):
    s = d.groupby(shell_col)[col].agg(["mean", "std", "count"]).reset_index()
    s["ci95"] = 1.96 * s["std"] / np.sqrt(s["count"])
    return s


def density_panel(ax, x, y, s=3):
    """Per-point gaussian_kde density, normalised to the panel maximum, jet colormap, dense points on top."""
    xy = np.vstack([x, y]); z = gaussian_kde(xy)(xy); o = np.argsort(z)
    ax.scatter(x[o], y[o], c=z[o] / z.max(), cmap="jet", vmin=0, vmax=1, s=s, edgecolor="none", rasterized=True)


def density_colorbar(fig, axs, fs=16, **kw):
    cb = fig.colorbar(ScalarMappable(norm=Normalize(0, 1), cmap="jet"), ax=axs, **kw)
    cb.set_label("relative point density (KDE)", fontsize=fs); cb.ax.tick_params(labelsize=fs - 2)
    return cb


def spear_box(ax, x, y, fontsize=19, pos=(0.04, 0.95, "left", "top")):
    r, p = spearmanr(x, y); tx, ty, ha, va = pos
    ax.text(tx, ty, f"Spearman r = {r:.2f}", transform=ax.transAxes, fontsize=fontsize, ha=ha, va=va,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    return r, p


def save(fig, name, dpi=300):
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote output/figures/{name}.{{pdf,png}}")


def save_table(df, name, float_format="%.6g"):
    df.to_csv(SRC_DIR / name, sep="\t", index=False, float_format=float_format)
    print(f"wrote output/source_data/{name}")
