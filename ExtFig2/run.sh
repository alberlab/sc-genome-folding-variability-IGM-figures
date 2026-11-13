#!/usr/bin/env bash
# Reproduce the Extended Data Figure 2 model panels (a, b, c bars, d, e, h). Panels d and e read the shared four-stage tables in ../data.
# Usage:  bash run.sh   (from any directory). Python: $PYTHON if set, else `python` of the active environment
# (create it with ../environment.yml: conda env create -f environment.yml; conda activate sc-folding-figures).
# Locale tables: data/locale/${LOCALE_VERSION:-v2_ab4r}/ (outputs in output_<version>/ for a non-default version).
# ChimeraX renders are not made here (see run_chimerax.sh).
set -euo pipefail
cd "$(dirname "$0")/code"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"   # $PYTHON, else python of the active env
[[ -n "$PY" ]] || { echo "python not found: conda activate sc-folding-figures (see ../environment.yml) or set PYTHON=" >&2; exit 1; }
echo "[ExtFig2] $("$PY" --version 2>&1) ($PY)"
OUTD=output; [[ "${LOCALE_VERSION:-v2_ab4r}" == "v2_ab4r" ]] || OUTD="output_$LOCALE_VERSION"
[[ "${DISTCT_NORM:-e2c}" == "e2c" ]] || OUTD="${OUTD}_distCT-$DISTCT_NORM"   # dist-to-CT normalisation (set_distct_norm.sh)
echo "[ExtFig2] locale version: ${LOCALE_VERSION:-v2_ab4r} (data/locale/<version>/), dist-to-CT norm: ${DISTCT_NORM:-e2c} -> $OUTD/"
for s in plot_extfig2a.py plot_extfig2b.py plot_extfig2b_right.py plot_extfig2c.py plot_extfig2d.py plot_extfig2e.py plot_extfig2h.py; do
  echo "=== ExtFig2/code/$s"
  "$PY" "$s"
done
echo "[ExtFig2] done. Figures in ExtFig2/$OUTD/figures/, source data in ExtFig2/$OUTD/source_data/."
