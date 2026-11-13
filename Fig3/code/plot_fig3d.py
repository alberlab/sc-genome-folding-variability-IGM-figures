#!/usr/bin/env python3
"""Figure 3d (IGM model columns): per-bin distributions in the four stage models, one 4 x 1 figure per metric —
SD(dist to CT), IQR(dist to CT) and nuclear locale entropy. 200 kb autosomal non-centromeric bins (13,315 per stage).
Bin edges shared by the four stages: 50 bins over a fixed range — SD(dist to CT) x 0-1.4 (y 0-2000 bins), locale entropy
x 0.6-2.4 (y 0-1000 bins); IQR(dist to CT) between the 0.5th and 99.5th percentiles of the pooled values. Values outside the
range are counted in the edge bins. The four stages of a metric share the y range. SD / IQR: light-grey bars with black edges (as the entropy panel), red dotted line = median. Locale entropy: two-component Gaussian
mixture for DE / HB / mHEP (low component red, high component blue, mixture black; dotted lines = component means)
and a single Gaussian for ES, as in the original figure (note: by BIC the ES distribution also prefers two components).
Mean(dist to CT): same layout, x -1 to 1 (values outside counted in the edge bins), y automatic; red dotted line = median.
GMM curves drawn thick (components lw 8, mixture lw 6.5; 2026-11-09). Locale = ENV v2.
Outputs: output/figures/Fig3d_IGM_SDdistCT_hist.*, Fig3d_IGM_entropy_hist_GMM.*, output/source_data/Fig3d_{summary,entropy_GMM}.tsv
(only the SD and entropy metrics are drawn here; the IQR / Mean branches of the code are unused)
"""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from fig3_common import CT_SCALE, CT_UNIT, CELLS, CELL_TITLE, load_200kb, save, save_table

D = {c: load_200kb(c) for c in CELLS}
METRICS = [("SDdistCT", "sd", f"SD(dist to CT{CT_UNIT})"),                       # manuscript panels only (IQR / Mean dropped, 2026-11)
           ("entropy", "entropy_bit", "Nuclear locale entropy (bits)")]
OUT_NAME = {"SDdistCT": "Fig3d_IGM_SDdistCT_hist", "entropy": "Fig3d_IGM_entropy_hist_GMM"}
RANGE = {"SDdistCT": ((0.0, 1.4), (0, 2000)), "entropy": ((0.6, 2.4), (0, 1000)), "MeandistCT": ((-1.0, 1.0), None)}    # (xlim, ylim); IQR: automatic
if CT_SCALE != 1.0:                     # DISTCT_NORM r / Req: rescaled x range, automatic y
    RANGE = {**RANGE, "SDdistCT": ((0.0, 1.4 * CT_SCALE), None), "MeandistCT": ((-CT_SCALE, CT_SCALE), None)}
FS_T, FS_L, FS_K, FS_M = 40, 36, 30, 30
N_COMP = {"h1": 1, "DE": 2, "HB": 2, "mHEP": 2}
rows, gmm_rows = [], []


def set_ylim_checked(ax, cnt, yl):
    if yl is None:
        ax.set_ylim(0, cnt.max() * 1.6)
    else:
        assert cnt.max() < yl[1], f"tallest bin {cnt.max():.0f} exceeds ylim {yl}"
        ax.set_ylim(*yl)

for tag, col, lab in METRICS:
    vals = {c: D[c][col].to_numpy(float) for c in CELLS}
    vals = {c: v[np.isfinite(v)] for c, v in vals.items()}
    if tag in RANGE:
        (lo, hi), YL = RANGE[tag]
    else:
        (lo, hi), YL = np.percentile(np.concatenate(list(vals.values())), [0.5, 99.5]), None
    edges = np.linspace(lo, hi, 51)
    bw = edges[1] - edges[0]
    fig, axs = plt.subplots(len(CELLS), 1, figsize=(13, 17 * len(CELLS) / 4), sharex=True)   # HB removed (2026-11-04)
    for ax, c in zip(axs, CELLS):
        v = vals[c]; med = float(np.median(v))
        if tag != "entropy":
            cnt, _, _ = ax.hist(np.clip(v, lo, hi), bins=edges, color="lightgrey", edgecolor="black", lw=0.5)
            set_ylim_checked(ax, cnt, YL)
            ax.axvline(med, color="#d62728", ls=":", lw=3)
            ax.text(0.02, 0.95, CELL_TITLE[c], transform=ax.transAxes, ha="left", va="top", fontsize=FS_M, fontweight="bold")
            ax.text(0.98, 0.95, f"median = {med:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=FS_M - 8,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
        else:
            cnt, _, _ = ax.hist(np.clip(v, lo, hi), bins=edges, color="lightgrey", edgecolor="black", lw=0.5)
            k = N_COMP[c]
            gm = GaussianMixture(n_components=k, random_state=0, n_init=5).fit(v.reshape(-1, 1))
            mu, sd, w = gm.means_.ravel(), np.sqrt(gm.covariances_.ravel()), gm.weights_
            o = np.argsort(mu); mu, sd, w = mu[o], sd[o], w[o]
            xs = np.linspace(lo, hi, 400); scale = len(v) * bw
            if k == 1:
                ax.plot(xs, w[0] * norm.pdf(xs, mu[0], sd[0]) * scale, "k--", lw=8)
                ax.axvline(mu[0], color="black", ls=":", lw=3)
                ttl = f"μ = {mu[0]:.2f}"
            else:
                for i, (colr, nm) in enumerate((("#d62728", "low"), ("steelblue", "high"))):
                    ax.plot(xs, w[i] * norm.pdf(xs, mu[i], sd[i]) * scale, color=colr, ls="--", lw=8)
                    ax.axvline(mu[i], color=colr, ls=":", lw=3)
                ax.plot(xs, sum(w[i] * norm.pdf(xs, mu[i], sd[i]) for i in range(2)) * scale, color="black", lw=6.5, alpha=0.6)
                ttl = f"low μ={mu[0]:.2f} w={w[0]:.2f} | high μ={mu[1]:.2f} w={w[1]:.2f}"
            set_ylim_checked(ax, cnt, YL)
            ax.text(0.02, 0.95, CELL_TITLE[c], transform=ax.transAxes, ha="left", va="top", fontsize=FS_M, fontweight="bold")
            ax.text(0.98, 0.95, ttl, transform=ax.transAxes, ha="right", va="top", fontsize=FS_M - 8,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
            bic = {kk: GaussianMixture(n_components=kk, random_state=0, n_init=5).fit(v.reshape(-1, 1)).bic(v.reshape(-1, 1))
                   for kk in (1, 2)}
            for i in range(k):
                gmm_rows.append(dict(cell=c, n_components=k, component=["low", "high"][i] if k == 2 else "single",
                                     mean=mu[i], sd=sd[i], weight=w[i], BIC_1comp=bic[1], BIC_2comp=bic[2]))
        ax.set_ylabel("bins", fontsize=FS_L - 4); ax.tick_params(labelsize=FS_K, length=9, width=1.8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        q1, q3 = np.percentile(v, [25, 75])
        rows.append(dict(metric=lab, cell=c, n=len(v), median=med, q25=q1, q75=q3, mean=float(v.mean()),
                         n_below_range=int((v < lo).sum()), n_above_range=int((v > hi).sum())))
    axs[-1].set_xlabel(lab, fontsize=FS_L); axs[0].set_xlim(lo, hi)
    fig.suptitle(f"IGM model: {lab}", fontsize=FS_T, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, OUT_NAME[tag], dpi=150)

S = pd.DataFrame(rows); save_table(S, "Fig3d_summary.tsv", float_format="%.4f")
G = pd.DataFrame(gmm_rows); save_table(G, "Fig3d_entropy_GMM.tsv", float_format="%.4f")
pd.set_option("display.width", 200)
print(S.pivot(index="metric", columns="cell", values="median")[CELLS].round(2).to_string())
print(G.round(3).to_string(index=False))
