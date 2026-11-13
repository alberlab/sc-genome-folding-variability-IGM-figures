#!/usr/bin/env python
"""Figure 5a/5b ChimeraX files for one IGM occ0.6 stage model (run on the cluster; needs the HSS).

Usage: fig5_structures.py <cell> <hss> <a> <b> <c> <ct_stats.tsv> <locale.tsv> <speckle_dir> <outdir>
For structures 82, 399, 427, 602, 833 writes one cmm/pdb pair per chromosome copy, coloured by a per-200 kb-bin value of
the model population (every bead gets the value of its bin; both copies the same):
  cmm_pdb_IQR/      IQR(dist to CT)   blue-white-red over [0.42, 1.32]   (shared by the four stages: pooled 5-95 % range)
  cmm_pdb_SD/       SD(dist to CT)    blue-white-red over [0.41, 1.16]
  cmm_pdb_entropy/  locale entropy    cividis over [0.8, 2.3]
B-factor = the raw value (%6.3f); centromeric / gap bins: B-factor = -1, grey.
Speckles (Fig 5a): <cell>_Struct<k>_speckle_centroids.{cmm,pdb}, one marker per speckle (MCL clusters of speckle-
associated beads with >= 5 beads on >= 2 chromosomes, the speckles used for the nuclear-locale definition), centroid of
the cluster, B-factor = cluster size. <cell>_fig5.cxc: nuclear envelope + the colour palettes for ChimeraX.
Coordinates: nm in cmm, um in pdb; bead radius from the HSS.
"""
import os, sys
import numpy as np, pandas as pd
from xml.sax.saxutils import escape
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
from alabtools import HssFile

cell, hss_path = sys.argv[1], sys.argv[2]
A, B, C = map(float, sys.argv[3:6])
ct_path, loc_path, spk_dir, outdir = sys.argv[6:10]
STRUCTS = [82, 399, 427, 602, 833]
BWR = mcol.LinearSegmentedColormap.from_list("bwr3", ["blue", "white", "red"])
SIGNALS = {"IQR": ("iqr_signed_dist", BWR, 0.42, 1.32), "SD": ("sd_signed_dist", BWR, 0.41, 1.16),
           "entropy": ("entropy_bit", plt.get_cmap("cividis"), 0.8, 2.3)}
os.makedirs(outdir, exist_ok=True)

hss = HssFile(hss_path, "r")
radii = np.asarray(hss.radii, dtype=float)
idx = hss.index
chrom = np.asarray(idx.chrom).astype(int); copy = np.asarray(idx.copy).astype(int)
start = np.asarray(idx.start).astype(np.int64); chromstr = np.asarray(idx.chromstr).astype(str)
chroms = [str(c) for c in hss.genome.chroms]
n_beads = len(radii)

# per-bin values mapped to every bead (by chr/start); cen / gap bins -> NaN
ct = pd.read_csv(ct_path, sep="\t"); lc = pd.read_csv(loc_path, sep="\t")
tab = ct.merge(lc[["chr", "start", "label", "entropy_bit"]], on=["chr", "start"], how="inner", validate="1:1")
assert len(tab) == 15453
tab.loc[tab.label == "cen", ["iqr_signed_dist", "sd_signed_dist", "entropy_bit"]] = np.nan
key = pd.DataFrame({"chr": chromstr, "start": start})
bead = key.merge(tab, on=["chr", "start"], how="left")
assert len(bead) == n_beads and bead.iqr_signed_dist.notna().sum() > 0
print(f"[{cell}] beads {n_beads}; beads without value (cen/gap): {int(bead.entropy_bit.isna().sum())}", flush=True)


def write_cmm(path, xyz, r, marker, rgb, link_radius=20, links=True):
    with open(path, "w") as f:
        f.write(f'<marker_set name="marker set {escape(str(marker))}">\n')
        for i in range(len(xyz)):
            cr, cg, cb = rgb[i]
            f.write(f'<marker id="{i + 1}" x="{xyz[i, 0]:.3f}" y="{xyz[i, 1]:.3f}" z="{xyz[i, 2]:.3f}" r="{cr:.3f}" '
                    f'g="{cg:.3f}" b="{cb:.3f}" radius="{r[i]:.3f}" note="" nr="{cr:.3f}" ng="{cg:.3f}" nb="{cb:.3f}"/>\n')
        if links:
            for i in range(len(xyz) - 1):
                cr, cg, cb = rgb[i]
                f.write(f'<link id1="{i + 1}" id2="{i + 2}" r="{cr:.3f}" g="{cg:.3f}" b="{cb:.3f}" radius="{link_radius}" />\n')
        f.write("</marker_set>\n")


def write_pdb(path, xyz, bfactor, resname="ALA", chain="A", conect=True):
    with open(path, "w") as f:
        for i in range(len(xyz)):
            x, y, z = xyz[i] / 1000.0
            f.write(f"ATOM  {i + 1:5d}  CA  {resname:<3s} {chain}{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{bfactor[i]:6.3f}"
                    f"           C\n")
        if conect:
            for i in range(1, len(xyz)):
                f.write(f"CONECT{i:5d}{i + 1:5d}\n")
        f.write("END\n")


for s in STRUCTS:
    crd = np.asarray(hss.get_struct_crd(s), dtype=float)
    for tag, (col, cmap, vmin, vmax) in SIGNALS.items():
        val = bead[col].to_numpy(float)
        rgb = cmap(plt.Normalize(vmin, vmax)(np.clip(val, vmin, vmax)))[:, :3]
        rgb[~np.isfinite(val)] = (0.6, 0.6, 0.6)
        bf = np.where(np.isfinite(val), val, -1.0)
        d = f"{outdir}/cmm_pdb_{tag}/Struct{s}"; os.makedirs(d, exist_ok=True)
        for ci, chname in enumerate(chroms):
            for c in np.unique(copy[chrom == ci]):
                m = (chrom == ci) & (copy == c); name = f"{cell}_Struct{s}_{chname}{'AB'[c]}"
                write_cmm(f"{d}/{name}.cmm", crd[m], radii[m], f"{s}{cell}{chname}{'AB'[c]}", rgb[m])
                write_pdb(f"{d}/{name}.pdb", crd[m], bf[m])
    # speckle centroids (N565 speckles)
    z = np.load(f"{spk_dir}/struct{s}.npz", allow_pickle=True)
    assert int(z["struct_id"]) == s
    cen = np.asarray(z["kept_cluster_centroids"], dtype=float); size = np.asarray(z["cluster_sizes"])[np.asarray(z["kept_mask"], bool)]
    assert len(cen) == len(size)
    d = f"{outdir}/speckles"; os.makedirs(d, exist_ok=True)
    write_cmm(f"{d}/{cell}_Struct{s}_speckle_centroids.cmm", cen, np.full(len(cen), 150.0), f"{s}{cell}speckles",
              np.tile([0.55, 0.27, 0.07], (len(cen), 1)), links=False)
    write_pdb(f"{d}/{cell}_Struct{s}_speckle_centroids.pdb", cen, size.astype(float), resname="SPK", chain="S", conect=False)
    print(f"[{cell}] structure {s}: cmm/pdb for IQR / SD / entropy + {len(cen)} speckle centroids", flush=True)

with open(f"{outdir}/{cell}_fig5.cxc", "w") as f:
    f.write(f"# {cell} model: nuclear envelope (ellipsoid at the origin, semi-axes {A:g} / {B:g} / {C:g} nm)\n")
    f.write("# for the .cmm files (nm):\n")
    f.write(f"shape ellipsoid radius {A:g},{B:g},{C:g} center 0,0,0 color #b0b0b0 mesh true name {cell}_nucleus\n")
    f.write("# for the .pdb files (um) use instead:\n")
    f.write(f"# shape ellipsoid radius {A / 1000:g},{B / 1000:g},{C / 1000:g} center 0,0,0 color #b0b0b0 mesh true name {cell}_nucleus_um\n")
    f.write("# colour the loaded pdb files by B-factor (choose the line matching the files you opened):\n")
    f.write("# Fig 5a IQR(dist to CT):  color byattribute bfactor palette 0.42,blue:0.87,white:1.32,red\n")
    f.write("# Fig 5a SD(dist to CT):   color byattribute bfactor palette 0.41,blue:0.785,white:1.16,red\n")
    f.write("# Fig 5b locale entropy:   color byattribute bfactor palette 0.8,#00204d:1.175,#414d6b:1.55,#7c7b78:1.925,#bcaf6f:2.3,#ffe838\n")
    f.write("# hide centromere / gap beads:  hide @@bfactor<0\n")
    f.write("# speckle centroids: open <cell>_Struct<k>_speckle_centroids.cmm (brown spheres, radius 150 nm)\n")
hss.close()
print(f"[{cell}] DONE", flush=True)
