#!/usr/bin/env python
"""Extended Fig 2c ChimeraX files: occ0.6 stage model coloured by nuclear locale (run on the cluster; needs the HSS).

Usage: extfig2c_locale_structures.py <cell> <hss> <a> <b> <c> <labels.npz> <locale.tsv> <speckle_dir> <outdir>
For structures 82, 399, 427, 602, 833 writes one cmm/pdb pair per chromosome copy, every bead coloured by its OWN locale
in that structure (N565 labels: spe / nuc / env = within 5 / 6 / 5 bead radii of a speckle / nucleolus / the envelope;
NL-act / NL-rep = near none, A/B ratio > 0.5 or not; overlap classes = near two or three landmarks):
  cmm  bead colour = locale colour, bead radius from the HSS, chain links; coordinates nm
  pdb  residue name encodes the locale (spe SER, nuc ASN, env GLU, spe_nuc HIS, spe_env LYS, nuc_env ARG, act ALA,
       rep PRO, triple TRP, missing GLY), B-factor = locale code (0-9, order of `classes`); coordinates um
Speckles: <cell>_Struct<k>_speckle_centroids.{cmm,pdb} (N565 speckles: MCL clusters >= 5 beads on >= 2 chromosomes,
centroid per cluster, B-factor = cluster size). <cell>_extfig2c.cxc: envelope + colour-by-residue-name commands.
Check: per-bin locale frequencies rebuilt from the per-bead labels must equal the per-bin table (bead order sanity).
"""
import os, sys
import numpy as np, pandas as pd
from xml.sax.saxutils import escape
from alabtools import HssFile

cell, hss_path = sys.argv[1], sys.argv[2]
A, B, C = map(float, sys.argv[3:6])
lab_path, loc_path, spk_dir, outdir = sys.argv[6:10]
STRUCTS = [82, 399, 427, 602, 833]
CLASSES = ["spe", "nuc", "env", "spe_nuc", "spe_env", "nuc_env", "act", "rep", "triple", "missing"]
RESNAME = {"spe": "SER", "nuc": "ASN", "env": "GLU", "spe_nuc": "HIS", "spe_env": "LYS", "nuc_env": "ARG",
           "act": "ALA", "rep": "PRO", "triple": "TRP", "missing": "GLY"}
HEX = {"spe": "#8b1e3f", "nuc": "#f4a261", "env": "#7fa9c9", "act": "#e57373", "rep": "#9b59b6",
       "spe_nuc": "#c06050", "spe_env": "#856484", "nuc_env": "#baa695", "triple": "#aa7878", "missing": "#cccccc"}
RGB = np.array([[int(HEX[k][i:i + 2], 16) / 255 for i in (1, 3, 5)] for k in CLASSES])
os.makedirs(outdir, exist_ok=True)

hss = HssFile(hss_path, "r")
radii = np.asarray(hss.radii, dtype=float)
idx = hss.index
chrom = np.asarray(idx.chrom).astype(int); copy = np.asarray(idx.copy).astype(int)
start = np.asarray(idx.start).astype(np.int64); chromstr = np.asarray(idx.chromstr).astype(str)
chroms = [str(c) for c in hss.genome.chroms]
n_beads = len(radii)

z = np.load(lab_path, allow_pickle=True)
labels = np.asarray(z["labels"])                                     # (n_beads, n_struct), int8
assert [str(c) for c in z["classes"]] == CLASSES, z["classes"]
assert labels.shape[0] == n_beads, (labels.shape, n_beads)
print(f"[{cell}] beads {n_beads}, labels {labels.shape}", flush=True)

# sanity: rebuild per-bin P_* (fraction of valid observations over 2 copies x all structures) and compare to the table
lc = pd.read_csv(loc_path, sep="\t")
key = pd.DataFrame({"chr": chromstr, "start": start})
valid = labels != CLASSES.index("missing")
rows = {"chr": chromstr, "start": start, "n_valid": valid.sum(1)}
for k in ["spe", "nuc", "env", "act", "rep"]:
    rows[k] = (labels == CLASSES.index(k)).sum(1)
g = pd.DataFrame(rows).groupby(["chr", "start"], sort=False).sum().reset_index()
m = g.merge(lc, on=["chr", "start"], how="inner", validate="1:1")
m = m[m.label != "cen"]
worst = 0.0
for k in ["spe", "nuc", "env", "act", "rep"]:
    worst = max(worst, float(np.nanmax(np.abs(m[k] / m.n_valid_x - m[f"P_{k}"]))))
print(f"[{cell}] per-bin check: {len(m)} bins, max |P_rebuilt - P_table| = {worst:.2e}", flush=True)
assert worst < 1e-3, "bead order of labels does not match the HSS / table"


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


def write_pdb(path, xyz, bfactor, resnames, chain="A", conect=True):
    with open(path, "w") as f:
        for i in range(len(xyz)):
            x, y, zz = xyz[i] / 1000.0
            f.write(f"ATOM  {i + 1:5d}  CA  {resnames[i]:<3s} {chain}{i + 1:4d}    {x:8.3f}{y:8.3f}{zz:8.3f}  1.00"
                    f"{bfactor[i]:6.2f}           C\n")
        if conect:
            for i in range(1, len(xyz)):
                f.write(f"CONECT{i:5d}{i + 1:5d}\n")
        f.write("END\n")


counts = []
for s in STRUCTS:
    crd = np.asarray(hss.get_struct_crd(s), dtype=float)
    code = labels[:, s].astype(int)
    rgb = RGB[code]; res = np.array([RESNAME[CLASSES[c]] for c in code])
    d = f"{outdir}/cmm_pdb_locale/Struct{s}"; os.makedirs(d, exist_ok=True)
    for ci, chname in enumerate(chroms):
        for c in np.unique(copy[chrom == ci]):
            mm = (chrom == ci) & (copy == c); name = f"{cell}_Struct{s}_{chname}{'AB'[c]}"
            write_cmm(f"{d}/{name}.cmm", crd[mm], radii[mm], f"{s}{cell}{chname}{'AB'[c]}", rgb[mm])
            write_pdb(f"{d}/{name}.pdb", crd[mm], code[mm].astype(float), res[mm])
    counts.append(dict(cell=cell, struct=s, **{k: int((code == i).sum()) for i, k in enumerate(CLASSES)}))
    zz = np.load(f"{spk_dir}/struct{s}.npz", allow_pickle=True)
    assert int(zz["struct_id"]) == s
    cen = np.asarray(zz["kept_cluster_centroids"], dtype=float)
    size = np.asarray(zz["cluster_sizes"])[np.asarray(zz["kept_mask"], bool)]
    assert len(cen) == len(size)
    d = f"{outdir}/speckles"; os.makedirs(d, exist_ok=True)
    write_cmm(f"{d}/{cell}_Struct{s}_speckle_centroids.cmm", cen, np.full(len(cen), 150.0), f"{s}{cell}speckles",
              np.tile([0.55, 0.27, 0.07], (len(cen), 1)), links=False)
    write_pdb(f"{d}/{cell}_Struct{s}_speckle_centroids.pdb", cen, size.astype(float), ["SPK"] * len(cen), chain="S",
              conect=False)
    print(f"[{cell}] structure {s}: locale cmm/pdb + {len(cen)} speckle centroids", flush=True)
pd.DataFrame(counts).to_csv(f"{outdir}/{cell}_locale_bead_counts_per_structure.tsv", sep="\t", index=False)

with open(f"{outdir}/{cell}_extfig2c.cxc", "w") as f:
    f.write(f"# {cell} model coloured by nuclear locale. Open the per-chromosome pdb files first, then this file.\n")
    f.write(f"# nuclear envelope (ellipsoid at the origin, semi-axes {A:g} / {B:g} / {C:g} nm); pdb coordinates are um\n")
    f.write(f"shape ellipsoid radius {A / 1000:g},{B / 1000:g},{C / 1000:g} center 0,0,0 color #b0b0b0 mesh true name {cell}_nucleus\n")
    f.write(f"# for the cmm files (nm) use instead:\n")
    f.write(f"# shape ellipsoid radius {A:g},{B:g},{C:g} center 0,0,0 color #b0b0b0 mesh true name {cell}_nucleus\n")
    f.write("# residue name = locale: spe SER, nuc ASN, env GLU, NL-act ALA, NL-rep PRO,\n")
    f.write("# overlaps spe_nuc HIS, spe_env LYS, nuc_env ARG, triple TRP, missing GLY\n")
    for k in CLASSES:
        f.write(f"color :{RESNAME[k]} {HEX[k]}\n")
    f.write("# speckle centroids: open <cell>_Struct<k>_speckle_centroids.pdb, then\n")
    f.write("# color :SPK #8c4512\n")
    f.write("# show only some locales, e.g. speckle beads:  hide ~:SPK atoms ; show :SER atoms\n")
hss.close()
print(f"[{cell}] DONE", flush=True)
