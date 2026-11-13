#!/usr/bin/env python3
"""Extended Data Fig. 4d (IGM row): chr11, 5 representative loci (imaging-defined; 119-120, 107-109, 11-12, 41-42, 88.0-88.4 Mb),
change of the distance to the CT core across the four stage models (ES, DE, HB, HLC). One panel per locus, two y axes:
left (blue) Mean(dist to CT), right (red) SD(dist to CT) (IQR version dropped, 2026-11). Layout follows
the original Fig. 3G / S5C notebook (cis-ratio mean and CV per region): the value of a locus = mean of the per-bin
statistics of its 200 kb bins (each over 2 copies x 1000 structures). Axis ranges are shared by the five panels.
The plate / droplet scHi-C rows of the manuscript panel are not made here.
Outputs: output/figures/ExtFig4d_IGM_chr11_5loci_Mean_SDdistCT.*, output/source_data/ExtFig4d_chr11_5loci.tsv
"""
import os
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.family": "sans-serif",
                            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
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
OUT = ROOT / ("output" if DISTCT_NORM == DEFAULT_DISTCT else f"output_distCT-{DISTCT_NORM}")
FIG_DIR, SRC_DIR = OUT / "figures", OUT / "source_data"
FIG_DIR.mkdir(parents=True, exist_ok=True); SRC_DIR.mkdir(parents=True, exist_ok=True)

CELLS = ["h1", "DE", "HB", "mHEP"]
CELL_TITLE = {"h1": "ES", "DE": "DE", "HB": "HB", "mHEP": "HLC"}
LOCI = [("chr11", 119_000_000, 120_000_000, "chr11:119–120Mb"), ("chr11", 107_000_000, 109_000_000, "chr11:107–109Mb"),
        ("chr11", 11_000_000, 12_000_000, "chr11:11–12Mb"), ("chr11", 41_000_000, 42_000_000, "chr11:41–42Mb"),
        ("chr11", 88_000_000, 88_400_000, "chr11:88.0–88.4Mb")]
V_LABEL = {"mean": f"Mean(dist to CT{CT_UNIT})", "sd": f"SD(dist to CT{CT_UNIT})", "iqr": f"IQR(dist to CT{CT_UNIT})"}

rows = []
for c in CELLS:
    t = pd.read_csv(CT_DIR / f"{c}_occ06_distCT_top30_2p5r_stats.tsv", sep="\t")
    assert len(t) == 15453, (c, len(t))
    for ch, s, e, lab in LOCI:
        b = t[(t.chr == ch) & (t.start >= s) & (t.end <= e)]
        assert len(b) == (e - s) // 200_000 and np.isfinite(b[["mean_signed_dist", "sd_signed_dist", "iqr_signed_dist"]]).all().all()
        rows.append(dict(locus=lab, chr=ch, start=s, end=e, n_200kb_bins=len(b), cell=c, stage=CELL_TITLE[c],
                         mean=b.mean_signed_dist.mean(), sd=b.sd_signed_dist.mean(), iqr=b.iqr_signed_dist.mean()))
T = pd.DataFrame(rows)
T.to_csv(SRC_DIR / "ExtFig4d_chr11_5loci.tsv", sep="\t", index=False, float_format="%.6g")
print("wrote output/source_data/ExtFig4d_chr11_5loci.tsv")
print(T.pivot(index="locus", columns="stage", values=["mean", "sd", "iqr"]).round(3).to_string())


def shared_lim(v, pad=0.12):
    a, b = np.nanmin(v), np.nanmax(v); w = (b - a) * pad
    return a - w, b + w


FS = 15
x = np.arange(len(CELLS))
for v in ("sd",):                                                       # SD only (2026-11)
    LY, RY = shared_lim(T["mean"]), shared_lim(T[v])
    fig, axes = plt.subplots(1, len(LOCI), figsize=(4.6 * len(LOCI), 4.4))
    for i, (ax, (_, _, _, lab)) in enumerate(zip(axes, LOCI)):
        d = T[T.locus == lab].set_index("cell").loc[CELLS]
        ln1, = ax.plot(x, d["mean"], lw=2.5, color="blue", marker="o", ms=5, label=V_LABEL["mean"])
        ax.set_ylim(*LY); ax.tick_params(axis="y", colors="blue", labelsize=FS - 3)
        ax.set_ylabel(V_LABEL["mean"], color="blue", fontsize=FS)
        ax2 = ax.twinx()
        ln2, = ax2.plot(x, d[v], lw=2.5, color="red", marker="o", ms=5, label=V_LABEL[v])
        ax2.set_ylim(*RY); ax2.tick_params(axis="y", colors="red", labelsize=FS - 3)
        ax2.set_ylabel(V_LABEL[v], color="red", fontsize=FS)
        ax.set_xticks(x, [CELL_TITLE[c] for c in CELLS], fontsize=FS - 2); ax.set_xlim(-0.3, len(CELLS) - 0.7)
        ax.set_title(lab, fontsize=FS, loc="left")
    axes[len(LOCI) // 2].set_xlabel("Cell stage", fontsize=FS)
    fig.tight_layout(w_pad=3.0)
    fig.legend([ln1, ln2], [ln1.get_label(), ln2.get_label()], loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2,
               frameon=False, fontsize=FS - 1)
    name = f"ExtFig4d_IGM_chr11_5loci_Mean_{v.upper()}distCT"
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig); print(f"wrote output/figures/{name}.{{pdf,png}}")
