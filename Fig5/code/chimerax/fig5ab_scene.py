# Figure 5a / 5b ChimeraX scene for one stage model (run inside ChimeraX; driven by run_fig5{a,b}_chimerax.sh).
#   runscript fig5ab_scene.py <package_root> <cell> <SD|IQR|entropy> <structure>
# Opens the 46 per-chromosome pdb files of the structure (B-factor = SD / IQR of dist to CT, or locale entropy, of each
# bead's 200 kb bin), the nuclear envelope map (Fig4i/<cell>_nucleus_shell_pores.mrc, grey surface at a per-cell level)
# and, for 5a (SD / IQR), the speckle centroids. White background, worm cartoon coloured by B-factor (5a blue-white-red
# over 0.2-1.5; 5b blue-yellow over 0.8-2.3), view orient, near clipping (speckle centroids not clipped), speckle
# centroids as dark-red spheres, turn x -60, 2 um scale bar, zoom pixelsize 0.03. Saves "<cell> <metric>.cxs" and
# "<cell> <metric>.png" (2000 x 2000, supersample 3, transparent) to output/figures/Fig5a (SD, IQR) or Fig5b (entropy).
import os, sys, glob
from chimerax.core.commands import run

argv = [a for a in sys.argv[1:] if not a.endswith(".py")]
ROOT = argv[0]
CELL = argv[1] if len(argv) > 1 else "h1"
METRIC = argv[2] if len(argv) > 2 else "SD"
STRUCT = argv[3] if len(argv) > 3 else "82"
SETTINGS = {   # palette, range, draw speckles, output folder
    "SD":      ("blue:white:red", "0.2,1.5", True, "Fig5a"),
    "IQR":     ("blue:white:red", "0.2,1.5", True, "Fig5a"),
    "entropy": ("blue:yellow", "0.8,2.3", False, "Fig5b"),                   # as the earlier manual 5b recipe
}
PALETTE, RANGE, SPECKLES, PANEL = SETTINGS[METRIC]
SPK_COLOR, SPK_RADIUS = "darkred", 0.2
TURN_X, PIXELSIZE = -60, 0.03
SCALEBAR = "scalebar 2 color gray thickness 6 xpos .80 ypos .12"
IMG_PX, SUPERSAMPLE = 2000, 3
HIDE_CEN = False                                   # True: hide centromere / gap beads (B-factor = -1)
ENV_DIR = os.path.normpath(f"{ROOT}/../Fig4/output/figures/Fig4i")   # <cell>_nucleus_shell_pores.mrc
# nuclear-envelope surface level per cell (lower = thicker shell): h1 from the manual Fig 5a session (2026-10-17),
# DE / HB / mHEP from the manual Fig 4i sessions; colour = grey 178 (as in those sessions)
ENV_LEVEL = {"h1": 0.0777, "DE": 0.0479, "HB": 0.0965, "mHEP": 0.0245}
ENV_COLOR = "#b2b2b2"

src = f"{ROOT}/structures/{CELL}"
out = f"{ROOT}/output/figures/{PANEL}"; os.makedirs(out, exist_ok=True)
pdbs = sorted(glob.glob(f"{src}/cmm_pdb_{METRIC}/Struct{STRUCT}/*.pdb"))
spk = f"{src}/speckles/{CELL}_Struct{STRUCT}_speckle_centroids.pdb"
if len(pdbs) != 46 or (SPECKLES and not os.path.exists(spk)):
    raise FileNotFoundError(f"{CELL} {METRIC} Struct{STRUCT}: found {len(pdbs)} pdb files under {src} (run extract_structures.sh)")

run(session, "close session")
for p in pdbs:                                      # chromosome copies -> models #1-#46
    run(session, f"open {p}")
C = f"#1-{len(pdbs)}"
S = None
if SPECKLES:
    run(session, f"open {spk}")
    S = "#" + [m for m in session.models.list() if "speckle_centroids" in (m.name or "")][-1].id_string
env = f"{ENV_DIR}/{CELL}_nucleus_shell_pores.mrc"
if os.path.exists(env):                             # nuclear envelope (clipped together with the chromatin)
    run(session, f"open {env}")
    E = "#" + [m for m in session.models.list() if (m.name or "").endswith("_nucleus_shell_pores.mrc")][-1].id_string
    run(session, f"volume {E} style surface level {ENV_LEVEL.get(CELL, 0.08)} color {ENV_COLOR} step 1")
else:
    print(f"WARNING: no nuclear envelope map {env}")

run(session, "set bgcolor white")
run(session, "cartoon byattr worm_radius min:0.04 max:0.1")
run(session, f"color bfactor {C} palette {PALETTE} range {RANGE}")
if HIDE_CEN:
    run(session, f"hide {C}@@bfactor<0 cartoons")
run(session, "view orient")
run(session, "clip")
if S:
    run(session, f"clip model {S} false")
    run(session, f"show {S} atoms")
    run(session, f"style {S} sphere")
    run(session, f"color {S} {SPK_COLOR}")
    run(session, f"size {S} atomRadius {SPK_RADIUS}")
run(session, "~select")
run(session, f"turn x {TURN_X}")
run(session, SCALEBAR)
run(session, f"zoom pixelsize {PIXELSIZE}")

name = f"{CELL} {METRIC}" + ("" if STRUCT == "82" else f" Struct{STRUCT}")
run(session, f'save "{out}/{name}.cxs"')
run(session, f'save "{out}/{name}.png" width {IMG_PX} height {IMG_PX} supersample {SUPERSAMPLE} transparentBackground true')
print(f"FIG5_DONE {name}")
