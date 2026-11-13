# Extended Data Fig. 2c ChimeraX scenes for one stage model: one image per nuclear locale (run inside ChimeraX).
#   runscript extfig2c_locale_scene.py <package_root> <cell> [structure=82]
# Opens the 46 per-chromosome pdb files of the structure (residue name = locale of each bead in that structure) and the
# nuclear envelope map (../Fig4/output/figures/Fig4i/<cell>_nucleus_shell_pores.mrc, grey surface at a per-cell level).
# Beads drawn as spheres of the model's bead radius; white background, view orient, near clipping, turn x -60, 2 um scale
# bar, zoom pixelsize 0.03 (one view for all five locales). For each locale only its beads are shown in its colour;
# clipping: only the nuclear envelope is clipped; for every locale (env included, changed 2026-10-21) the whole bead
# cloud is shown above the half shell. Saves "<cell> <locale>.cxs" and
# "<cell> <locale>.png" (2000 x 2000, supersample 3, transparent) to output/figures/ExtFig2c/.
import os, sys, glob, time
from chimerax.core.commands import run as _run

def run(session, cmd):                               # every command timed into output/figures/ExtFig2c/steps_<cell>.log
    t = time.time(); _run(session, cmd)
    with open(STEPLOG, "a") as f:
        f.write(f"{time.time() - t:7.2f}s  {cmd[:120]}\n")

argv = [a for a in sys.argv[1:] if not a.endswith(".py")]
ROOT = argv[0]; CELL = argv[1] if len(argv) > 1 else "h1"; STRUCT = argv[2] if len(argv) > 2 else "82"
LOCALES = [("spe", "SER", "#B23A48", False), ("nuc", "ASN", "#f4a261", False), ("env", "GLU", "#7AAED1", False),
           ("NL-act", "ALA", "#e57373", False), ("NL-rep", "PRO", "#9b59b6", False)]   # name, residue, colour, clip beads
BEAD_R = {"h1": 0.0949, "DE": 0.1062, "HB": 0.0990, "mHEP": 0.1081}                 # model bead radius (um)
ENV_DIR = os.path.normpath(f"{ROOT}/../Fig4/output/figures/Fig4i")
ENV_LEVEL = {"h1": 0.0777, "DE": 0.0479, "HB": 0.0965, "mHEP": 0.0245}             # as Fig. 5a/5b (lower = thicker)
ENV_COLOR = "#b2b2b2"
TURN_X, PIXELSIZE = -60, 0.03
SCALEBAR = "scalebar 2 color gray thickness 6 xpos .80 ypos .12"
IMG_PX, SUPERSAMPLE = 2000, 3

src = f"{ROOT}/structures/{CELL}/cmm_pdb_locale/Struct{STRUCT}"
out = f"{ROOT}/output/figures/ExtFig2c"; os.makedirs(out, exist_ok=True)
STEPLOG = f"{out}/steps_{CELL}.log"; open(STEPLOG, "w").close()
pdbs = sorted(glob.glob(f"{src}/*.pdb"))
if len(pdbs) != 46:
    raise FileNotFoundError(f"{CELL} Struct{STRUCT}: found {len(pdbs)} pdb files under {src} (run extract_structures.sh)")

run(session, "close session")
for p in pdbs:                                      # chromosome copies -> models #1-#46
    run(session, f"open {p}")
C = f"#1-{len(pdbs)}"
env = f"{ENV_DIR}/{CELL}_nucleus_shell_pores.mrc"
if os.path.exists(env):
    run(session, f"open {env}")
    E = "#" + [m for m in session.models.list() if (m.name or "").endswith("_nucleus_shell_pores.mrc")][-1].id_string
    run(session, f"volume {E} style surface level {ENV_LEVEL.get(CELL, 0.08)} color {ENV_COLOR} step 1")
else:
    print(f"WARNING: no nuclear envelope map {env}")
run(session, "set bgcolor white")
run(session, "lighting soft")
run(session, f"cartoon hide {C}")
run(session, f"style {C} sphere")
run(session, f"size {C} atomRadius {BEAD_R[CELL]}")
run(session, "view orient")
run(session, "clip")
run(session, f"turn x {TURN_X}")
run(session, SCALEBAR)
run(session, f"zoom pixelsize {PIXELSIZE}")
for name, res, col, clip_beads in LOCALES:
    run(session, f"hide {C} atoms")
    run(session, f"show {C}&:{res} atoms")
    run(session, f"color {C}&:{res} {col}")
    run(session, f"clip model {C} {'true' if clip_beads else 'false'}")
    run(session, "~select")
    run(session, f'save "{out}/{CELL} {name}.cxs"')
    run(session, f'save "{out}/{CELL} {name}.png" width {IMG_PX} height {IMG_PX} supersample {SUPERSAMPLE} transparentBackground true')
print(f"EXTFIG2C_DONE {CELL}")
