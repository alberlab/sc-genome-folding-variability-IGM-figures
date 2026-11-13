# Figure 4a / 4e / 4i: re-render the PNGs from the saved ChimeraX sessions (inside ChimeraX).
#   runscript fig4_render_sessions.py <figure_root>
# The panels were composed interactively in ChimeraX (views, clipping, lighting); the .cxs sessions in
# output/figures/Fig4{a,e,i}/ are the render recipe (the envelope maps <cell>_nucleus_shell_pores.mrc of Fig4i are read by
# the sessions and by Fig5 / ExtFig2c). Each session is opened and saved as 2000 x 2000 px, supersample 3, transparent
# background (the settings of the published PNGs) to output/figures/rerender/<panel>/<session name>.png, so the published
# PNGs are not overwritten. Run log: output/figures/rerender/fig4_rerender.log
import glob, os, sys, time, traceback
from chimerax.core.commands import run

argv = [a for a in sys.argv[1:] if not a.endswith(".py")]
ROOT = argv[0]
out_root = f"{ROOT}/output/figures/rerender"; os.makedirs(out_root, exist_ok=True)
with open(f"{out_root}/fig4_rerender.log", "w") as log:
    log.write(f"{time.ctime()}  start\n"); log.flush()
    for panel in ("Fig4a", "Fig4e", "Fig4i"):
        os.makedirs(f"{out_root}/{panel}", exist_ok=True)
        for cxs in sorted(glob.glob(f"{ROOT}/output/figures/{panel}/*.cxs")):
            name = os.path.splitext(os.path.basename(cxs))[0]
            try:
                run(session, "close session")
                run(session, f'open "{cxs}"')
                run(session, f'save "{out_root}/{panel}/{name}.png" width 2000 height 2000 supersample 3 transparentBackground true')
                log.write(f"OK    {panel}/{name}\n")
            except Exception:
                log.write(f"FAIL  {panel}/{name}\n{traceback.format_exc()}\n")
            log.flush()
    log.write(f"{time.ctime()}  done\n")
