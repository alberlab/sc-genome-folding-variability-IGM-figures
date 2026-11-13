#!/usr/bin/env bash
# Reproduce every plotted model panel of the manuscript figures: runs <Figure>/run.sh for each figure (~1 min).
# Usage:  conda activate sc-folding-figures && bash run_all.sh        (from any directory)
#         PYTHON=/path/to/python bash run_all.sh                     (any interpreter with the packages of environment.yml)
# Outputs: <Figure>/output/figures/ (PNG + PDF) and <Figure>/output/source_data/ (TSV tables behind each panel).
# Optional switches (alternative input tables are NOT shipped, see README "Switches"):
#   LOCALE_VERSION=<v> bash run_all.sh   -> data/locale/<v>/,      outputs <Figure>/output_<v>/
#   DISTCT_NORM=<r|Req> bash run_all.sh  -> data/distCT/<norm>/,   outputs <Figure>/output_distCT-<norm>/
# ChimeraX renders are NOT part of this script (need the Zenodo download, see fetch_data.sh, and UCSF ChimeraX):
#   Fig4/run_chimerax.sh, Fig5/run_chimerax.sh, ExtFig2/run_chimerax.sh.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"
[[ -n "$PY" ]] || { echo "python not found: conda activate sc-folding-figures (see environment.yml) or set PYTHON=" >&2; exit 1; }
export PYTHON="$PY"                       # every <Figure>/run.sh uses the same interpreter
echo "python: $PY  ($("$PY" --version 2>&1))"
"$PY" - <<'PYEOF'
import importlib
for m in ["numpy", "pandas", "scipy", "matplotlib", "seaborn", "sklearn", "openpyxl"]:
    print(f"  {m:<11}{importlib.import_module(m).__version__}")
PYEOF
for f in Fig2 Fig3 Fig4 Fig5 ExtFig2 ExtFig4 ExtFig5 ExtFig6; do
  echo "################ $f"
  bash "$ROOT/$f/run.sh"
done
echo "All figures done."
