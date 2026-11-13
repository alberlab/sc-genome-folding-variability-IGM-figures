#!/usr/bin/env bash
# Reproduce the Extended Data Figure 4d IGM row.
# Usage:  bash run.sh   (from any directory). Python: $PYTHON if set, else `python` of the active environment
# (create it with ../environment.yml: conda env create -f environment.yml; conda activate sc-folding-figures).
# ChimeraX renders are not made here.
set -euo pipefail
cd "$(dirname "$0")/code"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"   # $PYTHON, else python of the active env
[[ -n "$PY" ]] || { echo "python not found: conda activate sc-folding-figures (see ../environment.yml) or set PYTHON=" >&2; exit 1; }
echo "[ExtFig4] $("$PY" --version 2>&1) ($PY)"
echo "[ExtFig4] python: $PY"
for s in plot_extfig4.py; do
  echo "=== ExtFig4/code/$s"
  "$PY" "$s"
done
echo "[ExtFig4] done. Figures in ExtFig4/output/figures/, source data in ExtFig4/output/source_data/."
