# Figures 5a / 5b: run fig5ab_scene.py for the four stage models and the requested metrics (inside ChimeraX);
# with SD or IQR also renders the speckle legend. Run log: output/figures/Fig5{a,b}/fig5_run.log.
#   runscript fig5ab_all.py <package_root> <metric> [metric ...]        (SD IQR -> Fig5a; entropy -> Fig5b)
import os, sys, time, traceback
from chimerax.core.commands import run

argv = [a for a in sys.argv[1:] if not a.endswith(".py")]
ROOT, METRICS = argv[0], argv[1:]
CELLS, STRUCT = ["h1", "DE", "HB", "mHEP"], "82"
PANEL = "Fig5b" if METRICS == ["entropy"] else "Fig5a"
out = f"{ROOT}/output/figures/{PANEL}"; os.makedirs(out, exist_ok=True)
with open(f"{out}/fig5_run.log", "w") as log:
    log.write(f"{time.ctime()}  metrics {METRICS}\n"); log.flush()
    jobs = [(f"{ROOT}/code/chimerax/fig5ab_scene.py {ROOT} {c} {m} {STRUCT}", f"{c} {m}") for m in METRICS for c in CELLS]
    if any(m in ("SD", "IQR") for m in METRICS):
        jobs.append((f"{ROOT}/code/chimerax/fig5_speckle_legend.py {ROOT} h1 {STRUCT}", "speckle legend h1"))
    for cmd, lab in jobs:
        try:
            run(session, f"runscript {cmd}")
            log.write(f"OK    {lab}\n")
        except Exception:
            log.write(f"FAIL  {lab}\n{traceback.format_exc()}\n")
        log.flush()
    log.write(f"{time.ctime()}  done\n")
