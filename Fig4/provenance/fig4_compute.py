#!/usr/bin/env python
"""Figure 4 structure-level computations for one IGM occ0.6 stage model (run on the cluster; needs the HSS).

Usage: fig4_compute.py <cell> <hss> <a> <b> <c> <seedmask_dir> <outdir>
  a b c         ellipsoidal nucleus semi-axes (nm)
  seedmask_dir  directory with chunk_<s>_<e>_seedmask.npy (bool, n_beads x 100): CT-core seed beads per structure
                (top 30 % strict-cis beads at 2.5 r of each chromosome copy, the definition used for dist to CT)

Outputs (in outdir):
  cmm_pdb_chrom/Struct<k>/<cell>_Struct<k>_<chrom><A|B>.{cmm,pdb}   Fig 4a/4e: one pair per chromosome copy, tab20 colour
                                                                     per chromosome, bead radius from the HSS (nm in cmm,
                                                                     um in pdb)
  cmm_pdb_shannon/Struct<k>/...                                      Fig 4i: same structures, colour silver->red by
                                                                     clip(H**6, 0, 10), B-factor = clip(H**6, 0, 10)
  <cell>_nucleus_ellipsoid.cxc                                       ChimeraX commands drawing the nuclear envelope
  <cell>_shannon_per_structure.npy                                   float32 (n_struct, n_beads), H per bead per structure
  <cell>_Shannon_haploid_200kb.tsv                                   chr start end shannon_mean shannon_sd n_neighbors_mean
                                                                     frac_no_neighbor (2 copies x n_struct pooled)
  <cell>_chr11_ct_volfrac_compactness.tsv                            per structure x copy: CT hull volume / nuclear volume,
                                                                     Rg, bead radius, compactness = r / Rg
Shannon index (chromosome intermingling): neighbours of bead i = beads j != i with centre distance < 2.5 x bead radius;
counts by chromosome identity (the two homologous copies pooled); H = -sum p ln p; H = 0 when there is no neighbour.
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy.spatial import cKDTree, ConvexHull
from xml.sax.saxutils import escape
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from matplotlib.colors import ListedColormap
from alabtools import HssFile

cell, hss_path = sys.argv[1], sys.argv[2]
A, B, C = map(float, sys.argv[3:6])
seed_dir, outdir = sys.argv[6], sys.argv[7]
STRUCTS = [82, 399, 427, 602, 833]
CUT_MULT = 2.5
os.makedirs(outdir, exist_ok=True)

hss = HssFile(hss_path, "r")
crd = np.asarray(hss.coordinates, dtype=np.float32)            # (n_beads, n_struct, 3)
radii = np.asarray(hss.radii, dtype=np.float64)
idx = hss.index
chrom = np.asarray(idx.chrom).astype(int); copy = np.asarray(idx.copy).astype(int)
start = np.asarray(idx.start).astype(np.int64); end = np.asarray(idx.end).astype(np.int64)
chromstr = np.asarray(idx.chromstr).astype(str)
chroms = [str(c) for c in hss.genome.chroms]
n_beads, n_struct = crd.shape[0], crd.shape[1]
r0 = float(np.median(radii)); cut = CUT_MULT * r0
print(f"[{cell}] beads {n_beads}, structures {n_struct}, chroms {len(chroms)}; radius min/median/max "
      f"{radii.min():.4f}/{r0:.4f}/{radii.max():.4f} nm; Shannon cutoff {cut:.2f} nm; copies {np.unique(copy)}", flush=True)


# ---------------- cmm / pdb writers ----------------
def write_cmm(path, xyz, r, marker, rgb, link_radius=20):
    with open(path, "w") as f:
        f.write(f'<marker_set name="marker set {escape(str(marker))}">\n')
        for i in range(len(xyz)):
            cr, cg, cb = rgb[i]
            f.write(f'<marker id="{i + 1}" x="{xyz[i, 0]:.3f}" y="{xyz[i, 1]:.3f}" z="{xyz[i, 2]:.3f}" r="{cr:.3f}" '
                    f'g="{cg:.3f}" b="{cb:.3f}" radius="{r[i]:.3f}" note="" nr="{cr:.3f}" ng="{cg:.3f}" nb="{cb:.3f}"/>\n')
        for i in range(len(xyz) - 1):
            cr, cg, cb = rgb[i]
            f.write(f'<link id1="{i + 1}" id2="{i + 2}" r="{cr:.3f}" g="{cg:.3f}" b="{cb:.3f}" radius="{link_radius}" />\n')
        f.write("</marker_set>\n")


def write_pdb(path, xyz, bfactor):
    with open(path, "w") as f:
        for i in range(len(xyz)):
            x, y, z = xyz[i] / 1000.0
            f.write(f"ATOM  {i + 1:5d}  CA  ALA A{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{bfactor[i]:6.2f}           C\n")
        for i in range(1, len(xyz)):
            f.write(f"CONECT{i:5d}{i + 1:5d}\n")
        f.write("END\n")


base = plt.cm.get_cmap("tab20")
chrom_rgb = np.array(ListedColormap(base(np.linspace(0, 1, len(chroms)))).colors)[:, :3]
silver_red = mcol.LinearSegmentedColormap.from_list("silver_to_red", [mcol.to_rgba("silver"), mcol.to_rgba("red")])


def copy_tag(c):
    return "AB"[c]


# ---------------- Shannon index per structure ----------------
n_chrom = len(chroms)
H = np.zeros((n_struct, n_beads), dtype=np.float32)
NN = np.zeros((n_struct, n_beads), dtype=np.int32)
for s in range(n_struct):
    xyz = crd[:, s, :].astype(np.float64)
    pairs = cKDTree(xyz).query_pairs(cut, output_type="ndarray")
    cnt = np.zeros((n_beads, n_chrom), dtype=np.int32)
    np.add.at(cnt, (pairs[:, 0], chrom[pairs[:, 1]]), 1)
    np.add.at(cnt, (pairs[:, 1], chrom[pairs[:, 0]]), 1)
    tot = cnt.sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = cnt / np.maximum(tot, 1)[:, None]
        h = -np.where(p > 0, p * np.log(p), 0.0).sum(1)
    H[s] = h; NN[s] = tot
    if s % 100 == 0:
        print(f"  shannon structure {s}: pairs {len(pairs)}, mean neighbours {tot.mean():.2f}, mean H {h.mean():.3f}", flush=True)

# brute-force check on structure STRUCTS[0] for 3000 beads
s = STRUCTS[0]; xyz = crd[:, s, :].astype(np.float64); sub = np.arange(0, n_beads, 10)
err = 0.0
for i in sub:
    d = np.sqrt(((xyz - xyz[i]) ** 2).sum(1)); nb = np.flatnonzero(d < cut); nb = nb[nb != i]
    if len(nb):
        cc = np.bincount(chrom[nb], minlength=n_chrom) / len(nb); hb = -(cc[cc > 0] * np.log(cc[cc > 0])).sum()
    else:
        hb = 0.0
    err = max(err, abs(hb - float(H[s, i])))
print(f"[{cell}] brute-force check on structure {s}, {len(sub)} beads: max |dH| = {err:.2e}", flush=True)
assert err < 1e-5
np.save(f"{outdir}/{cell}_shannon_per_structure.npy", H)

# haploid 200 kb table (mean over both copies and all structures)
key = pd.DataFrame({"chr": chromstr, "start": start, "end": end})
g = key.assign(i=np.arange(n_beads)).groupby(["chr", "start", "end"], sort=False)["i"].apply(list)
rows = []
for (ch, st, en), ii in g.items():
    v = H[:, ii].ravel(); nn = NN[:, ii].ravel()
    rows.append((ch, st, en, float(v.mean()), float(v.std(ddof=1)), float(nn.mean()), float((nn == 0).mean())))
T = pd.DataFrame(rows, columns=["chr", "start", "end", "shannon_mean", "shannon_sd", "n_neighbors_mean", "frac_no_neighbor"])
T.to_csv(f"{outdir}/{cell}_Shannon_haploid_200kb.tsv", sep="\t", index=False, float_format="%.6g")
print(f"[{cell}] haploid bins {len(T)}; H mean {T.shannon_mean.mean():.3f}, beads with no neighbour {(NN == 0).mean():.4f}; "
      f"H quantiles 5/50/95 % {np.percentile(H, [5, 50, 95]).round(3)}; H**6 > 10 fraction {(H.astype(float) ** 6 > 10).mean():.4f}",
      flush=True)

# ---------------- cmm / pdb for Fig 4a/4e (chromosome colour) and Fig 4i (Shannon colour) ----------------
for s in STRUCTS:
    xyz_all = crd[:, s, :].astype(np.float64)
    bf_all = np.clip(H[s].astype(float) ** 6, 0, 10)
    col_all = silver_red(plt.Normalize(0, 10)(bf_all))[:, :3]
    for sub in ("cmm_pdb_chrom", "cmm_pdb_shannon"):
        os.makedirs(f"{outdir}/{sub}/Struct{s}", exist_ok=True)
    for ci, chname in enumerate(chroms):
        for c in np.unique(copy[chrom == ci]):
            m = (chrom == ci) & (copy == c)
            tag = f"{cell}_Struct{s}_{chname}{copy_tag(c)}"
            xyz, r = xyz_all[m], radii[m]
            rgb = np.repeat(chrom_rgb[ci][None, :], m.sum(), 0)
            write_cmm(f"{outdir}/cmm_pdb_chrom/Struct{s}/{tag}.cmm", xyz, r, f"{s}{cell}{chname}{copy_tag(c)}", rgb)
            write_pdb(f"{outdir}/cmm_pdb_chrom/Struct{s}/{tag}.pdb", xyz, np.zeros(m.sum()))
            write_cmm(f"{outdir}/cmm_pdb_shannon/Struct{s}/{tag}.cmm", xyz, r, f"{s}{cell}{chname}{copy_tag(c)}", col_all[m])
            write_pdb(f"{outdir}/cmm_pdb_shannon/Struct{s}/{tag}.pdb", xyz, bf_all[m])
    print(f"[{cell}] structure {s}: cmm/pdb written", flush=True)
with open(f"{outdir}/{cell}_nucleus_ellipsoid.cxc", "w") as f:
    f.write(f"# nuclear envelope of the {cell} model: ellipsoid centred at the origin, semi-axes {A:g} / {B:g} / {C:g} nm\n")
    f.write(f"# use with the .cmm files (nm):\nshape ellipsoid radius {A:g},{B:g},{C:g} center 0,0,0 color #b0b0b0 "
            f"mesh true name {cell}_nucleus\n")
    f.write(f"# or with the .pdb files (um):\n# shape ellipsoid radius {A / 1000:g},{B / 1000:g},{C / 1000:g} center 0,0,0 "
            f"color #b0b0b0 mesh true name {cell}_nucleus_um\n")

# ---------------- chr11 CT volume fraction and compactness ----------------
seed = np.zeros((n_beads, n_struct), dtype=bool)
files = sorted(glob.glob(f"{seed_dir}/chunk_*_seedmask.npy"), key=lambda p: int(os.path.basename(p).split("_")[1]))
for p in files:
    s0, s1 = map(int, os.path.basename(p).split("_")[1:3])
    sm = np.load(p); assert sm.shape == (n_beads, s1 - s0), (p, sm.shape)
    seed[:, s0:s1] = sm
assert len(files) * 100 == n_struct, len(files)
ci = chroms.index("chr11"); nuc_vol = 4.0 / 3.0 * np.pi * A * B * C
rows = []
for s in range(n_struct):
    for c in (0, 1):
        m = (chrom == ci) & (copy == c)
        xyz = crd[m, s, :].astype(np.float64); sd = seed[m, s]
        vol = ConvexHull(xyz[sd]).volume if sd.sum() >= 4 else np.nan
        rg = float(np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(1).mean())); r = float(np.median(radii[m]))
        rows.append((s, c, int(m.sum()), int(sd.sum()), vol, vol / nuc_vol, rg, r, r / rg))
V = pd.DataFrame(rows, columns=["structure", "copy", "n_beads", "n_seed", "ct_hull_volume_nm3", "ct_volume_fraction",
                                "rg_nm", "bead_radius_nm", "compactness_r_over_rg"])
V.to_csv(f"{outdir}/{cell}_chr11_ct_volfrac_compactness.tsv", sep="\t", index=False, float_format="%.6g")
print(f"[{cell}] chr11: n {len(V)}, seed fraction {V.n_seed.sum() / V.n_beads.sum():.3f}, "
      f"median volfrac {V.ct_volume_fraction.median():.4f}, median r/Rg {V.compactness_r_over_rg.median():.4f}", flush=True)
hss.close()
print(f"[{cell}] DONE", flush=True)
