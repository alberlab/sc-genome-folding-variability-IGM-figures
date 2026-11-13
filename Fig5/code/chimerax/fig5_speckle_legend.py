# Speckle-centroid legend image for Fig. 5a: the speckle centroids of one model only, dark-red spheres, transparent.
#   runscript fig5_speckle_legend.py <package_root> [cell=h1] [structure=82]
import os, sys
from chimerax.core.commands import run
argv = [a for a in sys.argv[1:] if not a.endswith(".py")]
ROOT = argv[0]; CELL = argv[1] if len(argv) > 1 else "h1"; STRUCT = argv[2] if len(argv) > 2 else "82"
out = f"{ROOT}/output/figures/Fig5a"; os.makedirs(out, exist_ok=True)
run(session, "close session")
run(session, f"open {ROOT}/structures/{CELL}/speckles/{CELL}_Struct{STRUCT}_speckle_centroids.pdb")
run(session, "set bgcolor white")
run(session, "show atoms"); run(session, "style sphere"); run(session, "color darkred"); run(session, "size atomRadius 0.2")
run(session, "view orient"); run(session, "turn x -60"); run(session, "view")
run(session, f'save "{out}/speckle legend {CELL}.png" width 1000 height 1000 supersample 3 transparentBackground true')
print("FIG5_SPECKLE_LEGEND_DONE")
