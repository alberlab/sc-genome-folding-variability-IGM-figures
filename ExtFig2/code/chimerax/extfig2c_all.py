# Extended Data Fig. 2c: run extfig2c_locale_scene.py for the requested stage models (inside ChimeraX).
#   runscript extfig2c_all.py <package_root> [cell ...]      (default: h1 DE HB mHEP)
# Run log: output/figures/ExtFig2c/extfig2c_run.log (one line per cell, traceback if it fails).
import os, sys, time, traceback
from chimerax.core.commands import run
argv = [a for a in sys.argv[1:] if not a.endswith(".py")]
ROOT, CELLS = argv[0], (argv[1:] or ["h1", "DE", "HB", "mHEP"])
out = f"{ROOT}/output/figures/ExtFig2c"; os.makedirs(out, exist_ok=True)
with open(f"{out}/extfig2c_run.log", "w") as log:
    log.write(f"{time.ctime()}  cells {CELLS}\n"); log.flush()
    for cell in CELLS:
        try:
            run(session, f"runscript {ROOT}/code/chimerax/extfig2c_locale_scene.py {ROOT} {cell} 82")
            log.write(f"OK    {cell} (spe, nuc, env, NL-act, NL-rep)\n")
        except Exception:
            log.write(f"FAIL  {cell}\n{traceback.format_exc()}\n")
        log.flush()
    log.write(f"{time.ctime()}  done\n")
