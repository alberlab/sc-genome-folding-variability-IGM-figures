"""Shared inputs, constants and helpers for the Figure 2 panels (H1 ESC IGM model, occupancy 0.6, 200 kb).

Definitions used throughout
---------------------------
dist to CT   Signed distance of a locus to its chromosome-territory (CT) core, normalised by the core centroid-to-
             surface distance (dimensionless; positive = inside the CT core). CT core = convex hull of the top 30 % beads
             of each chromosome copy ranked by strict-cis contact ratio (contacts within 2.5 bead radii, same chromosome
             and same copy). Per 200 kb bin, 2 copies x 1000 structures = 2000 values are pooled; the input table holds
             their mean, SD and IQR (75th - 25th percentile).
locale       Nuclear locale of each bead in each structure: speckle (spe), nucleolus (nuc), nuclear envelope (env) within
             5 / 6 / 5 bead radii, otherwise NL-act / NL-rep by A/B ratio. Entropy (bits) over these 5 classes per bin;
             dominant locale = the class with the highest probability.
Eig1 (1 Mb)  Bulk Hi-C compartment eigenvector, column E1 of Supplementary Table 2 (sheet H1-platebased-HiC-batch1).
Eig1 (200kb) Bulk Hi-C compartment eigenvector at 200 kb (cooltools eigs_cis, GC-phased so that A > 0).
Bins         Autosomal 200 kb bins that are not centromeric (label != "cen"): 13,315 bins.
1 Mb values  Mean of the non-centromeric 200 kb bin statistics inside each 1 Mb bin (start // 1 Mb).
SD(radial)   SD of the ellipsoid-normalised radial position over 2 copies x 1000 structures = 2000 values per 200 kb bin
             (same pooling as SD(dist to CT)), from the per-copy mean / SD over structures (load_sd_radial).
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
LOC_TSV = LOC_DIR / "h1_occ06_locale_N565.tsv"
E1_200K = DATA / "h1_bulk_E1_200kb.tsv"
RAD_BED = SHARED / "radial" / "h1_occ06_radial_200kb.bed"
RAD_COPY_BED = DATA / "h1_occ06_radial_200kb_copysep.bed"
LOCI_NPZ = CT_DIR / "h1_occ06_chr11_5loci_distCT.npz"
BULK_TSV = DATA / "h1_bulk_cis_total_1Mb.tsv"
PLATE_XLSX = SHARED / "schic" / "Table 2-Table 3-Cis_total annotation in plate-based single-cell Hi-C.xlsx"
PLATE_SHEET = "H1-platebased-HiC-batch1"
DROPLET_XLSX = SHARED / "schic" / "Table 3-Cis_total annotation in droplet-based single-cell Hi-C.xlsx"
DROPLET_SHEET = "H1-dropletHiC >0.3kb filter"

P5 = ["P_spe", "P_nuc", "P_env", "P_act", "P_rep"]
LOC_COLORS = {"spe": "#8b1e3f", "nuc": "#f4a261", "env": "#7fa9c9", "act": "#e57373", "rep": "#9b59b6"}
LOC_DISP = {"spe": "spe", "nuc": "nuc", "env": "env", "act": "NL-act", "rep": "NL-rep"}
BWR = mcolors.LinearSegmentedColormap.from_list("blue_white_red", ["blue", "white", "red"])
V_LABEL = {"mean": f"Mean(dist to CT{CT_UNIT})", "iqr": f"IQR(dist to CT{CT_UNIT})", "sd": f"SD(dist to CT{CT_UNIT})"}
AUTOS = [f"chr{i}" for i in range(1, 23)]


def load_200kb():
    """Per-200kb-bin table: chr start end mean sd iqr entropy_bit P_* dom5 E1_200kb (autosomal non-cen bins only)."""
    ct = pd.read_csv(CT_TSV, sep="\t").rename(columns={"mean_signed_dist": "mean", "sd_signed_dist": "sd",
                                                       "iqr_signed_dist": "iqr"})
    lc = pd.read_csv(LOC_TSV, sep="\t")
    e1 = pd.read_csv(E1_200K, sep="\t").rename(columns={"chrom": "chr", "E1": "E1_200kb"})
    assert len(ct) == len(lc) == len(e1) == 15453, (len(ct), len(lc), len(e1))
    d = ct[["chr", "start", "end", "mean", "sd", "iqr", "n_samples"]].merge(
        lc[["chr", "start", "label", "entropy_bit"] + P5], on=["chr", "start"], how="inner", validate="1:1")
    d = d.merge(e1[["chr", "start", "E1_200kb"]], on=["chr", "start"], how="inner", validate="1:1")
    assert len(d) == 15453, "dist-to-CT / locale / E1 200kb bin grids do not align"
    d = d[d.chr.isin(AUTOS) & (d.label != "cen")].copy()
    assert len(d) == 13315, len(d)
    d["dom5"] = np.array([p[2:] for p in P5])[np.argmax(d[P5].to_numpy(float), axis=1)]
    return d.drop(columns=["label"]).reset_index(drop=True)


def read_schic_table(xlsx, sheet, cols):
    t = pd.read_excel(xlsx, sheet_name=sheet).rename(columns={"chrom": "chr"})
    t["chr"] = t["chr"].astype(str)
    assert (t.start % 1_000_000 == 0).all()
    return t[["chr", "start", "end"] + cols]


def load_e1_1mb():
    return read_schic_table(PLATE_XLSX, PLATE_SHEET, ["E1"]).rename(columns={"E1": "E1_1mb"})


def to_1mb(d, cols):
    """Mean of the 200 kb bins (already cen-filtered) inside each 1 Mb bin (start // 1e6)."""
    g = d.assign(start=(d.start // 1_000_000) * 1_000_000)
    out = g.groupby(["chr", "start"])[cols].mean()
    out["n_200kb"] = g.groupby(["chr", "start"]).size()
    return out.reset_index()


def load_sd_radial():
    """Pooled SD (ddof = 1) of the normalised radial position over 2 copies x 1000 structures per haploid 200 kb bin.
    The copy-separated bed holds, per copy, the mean and population SD over structures; pooled population variance =
    mean of the per-copy variances + variance of the per-copy means (equal n per copy). Returns chr start sd_rad."""
    c = pd.read_csv(RAD_COPY_BED, sep="\t")
    c["chr"] = c["chr"].replace({"chr23": "chrX", "chr24": "chrY"})     # HSS chromstr stores chrX / chrY as chr23 / 24
    c["v"] = c.norm_rad_std ** 2
    g = c.groupby(["chr", "start"], sort=False)
    out = pd.DataFrame({"v_within": g.v.mean(), "v_between": g.norm_rad.var(ddof=0), "mean_rad": g.norm_rad.mean(),
                        "n_copy": g.size()}).reset_index()
    n = out.n_copy * 1000
    out["sd_rad"] = np.sqrt((out.v_within + out.v_between) * n / (n - 1))
    rd = pd.read_csv(RAD_BED, sep="\t").merge(out, on=["chr", "start"], validate="1:1")
    assert len(rd) == 15453 and np.allclose(rd.norm_rad, rd.mean_rad, atol=2e-6), "copysep and haploid radial differ"
    return rd[["chr", "start", "sd_rad"]]


def chr_tracks_1mb(chrom):
    """1 Mb model tracks for one chromosome on the scHi-C table grid. A 1 Mb bin is masked (NaN, grey span) when it
    contains any centromeric 200 kb bin or its 1 Mb E1 is NaN. Returns (table indexed by start, list of masked spans)."""
    ct = pd.read_csv(CT_TSV, sep="\t").rename(columns={"mean_signed_dist": "mean", "sd_signed_dist": "sd",
                                                       "iqr_signed_dist": "iqr"})
    lc = pd.read_csv(LOC_TSV, sep="\t"); rd = pd.read_csv(RAD_BED, sep="\t")
    d = ct[["chr", "start", "mean", "sd", "iqr"]].merge(lc[["chr", "start", "label", "entropy_bit"]],
                                                         on=["chr", "start"], validate="1:1") \
        .merge(rd[["chr", "start", "norm_rad"]], on=["chr", "start"], validate="1:1") \
        .merge(load_sd_radial(), on=["chr", "start"], validate="1:1")
    assert len(d) == 15453
    d = d[d.chr == chrom].copy(); d["start_1mb"] = (d.start // 1_000_000) * 1_000_000
    cols = ["mean", "sd", "iqr", "entropy_bit", "norm_rad", "sd_rad"]
    m = d[d.label != "cen"].groupby("start_1mb")[cols].mean()
    cen_1mb = set(d.loc[d.label == "cen", "start_1mb"])
    e1 = load_e1_1mb(); e1 = e1[e1.chr == chrom].set_index("start")
    t = e1[["end", "E1_1mb"]].join(m, how="left")
    t["masked"] = t.index.isin(cen_1mb) | ~np.isfinite(t.E1_1mb)
    t.loc[t.masked, cols] = np.nan
    x = t.index.to_numpy() / 1e6
    spans, run = [], None
    for xi, mk in zip(x, t.masked):
        if mk and run is None: run = [xi, xi + 1]
        elif mk: run[1] = xi + 1
        elif run is not None: spans.append(run); run = None
    if run is not None: spans.append(run)
    print(f"{chrom}: {len(t)} 1 Mb bins, masked {int(t.masked.sum())} "
          f"(contain cen {len(cen_1mb)}, E1 NaN {int((~np.isfinite(t.E1_1mb)).sum())})")
    return t, spans


def plot_tracks(t, spans, tracks, chrom, title, name, xlim=(-5, 140)):
    """Stacked genome tracks (1 Mb). tracks = [(column, title, y label, colour, ylim or None), ...];
    ylim = (bottom, top) — give bottom > top to flip the axis; "flip" = auto range, flipped. x axis padded to xlim, ticks only on real coordinates."""
    TITLE_FS, LABEL_FS, TICK_FS = 28, 24, 20
    x = t.index.to_numpy() / 1e6
    chrom_end = t["end"].max() / 1e6
    fig, axes = plt.subplots(len(tracks), 1, figsize=(20, 3.5 * len(tracks)), sharex=True, squeeze=False)
    for ax, (c, ttl, ylab, col, yl) in zip(axes[:, 0], tracks):
        for a, b in spans:
            ax.axvspan(a, b, color="lightgray", alpha=0.5, lw=0, zorder=0)
        ax.plot(x, t[c].to_numpy(float), lw=2.5, color=col)
        if yl == "flip":                                              # auto range, large values at the bottom
            ax.invert_yaxis()
        elif yl is not None:
            ax.set_ylim(*yl)
        ax.set_title(ttl, fontsize=TITLE_FS, fontweight="bold"); ax.set_ylabel(ylab, fontsize=LABEL_FS)
        ax.grid(True, linestyle=":", lw=0.9, color="gray", alpha=0.7)
        ax.tick_params(labelsize=TICK_FS)
    axes[-1, 0].set_xlabel(f"{chrom} genomic position (Mb)", fontsize=LABEL_FS + 4)
    axes[-1, 0].set_xlim(*xlim)
    axes[-1, 0].set_xticks(np.arange(0, chrom_end + 1e-9, 20))       # no tick outside the chromosome
    H = 3.5 * len(tracks)                                              # figure height (inch)
    fig.suptitle(title, fontsize=TITLE_FS + 4, fontweight="bold", y=1 - 0.1 / H, va="top")
    fig.tight_layout(rect=[0, 0, 1, 1 - 0.45 / H])                     # keep 0.45 inch for the title
    save(fig, name)


TEXT_POS = {"ul": (0.05, 0.95, "left", "top"), "ur": (0.95, 0.95, "right", "top"),
            "ll": (0.05, 0.05, "left", "bottom"), "ml": (0.05, 0.52, "left", "top"), "lr": (0.95, 0.05, "right", "bottom")}


def spear_text(ax, x, y, fontsize=10, n=True, loc="ul"):
    r, p = spearmanr(x, y); tx, ty, ha, va = TEXT_POS[loc]
    ax.text(tx, ty, f"Spearman r = {r:.2f}\n" + (f"p-value = {p:.2e}" if p > 0 else "p < 1e-300") + (f"\nn = {len(x)}" if n else ""),
            transform=ax.transAxes, fontsize=fontsize, ha=ha, va=va,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=2))
    return r, p


def prange(v, lo=2, hi=98):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return tuple(np.percentile(v, [lo, hi]))


def lim(v, lo=0, hi=100, pad=0.03):
    a, b = prange(v, lo, hi); w = (b - a) * pad
    return a - w, b + w


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote output/figures/{name}.{{pdf,png}}")


def save_table(df, name, **kw):
    df.to_csv(SRC_DIR / name, sep="\t", index=False, float_format=kw.get("float_format", "%.6g"))
    print(f"wrote output/source_data/{name}")
