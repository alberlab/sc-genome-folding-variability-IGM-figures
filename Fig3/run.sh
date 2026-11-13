#!/usr/bin/env bash
# Reproduce the Figure 3 IGM panels (c, d, e: model rows / columns).
# Usage:  bash run.sh   (from any directory). Python: $PYTHON if set, else `python` of the active environment
# (create it with ../environment.yml: conda env create -f environment.yml; conda activate sc-folding-figures).
# Locale tables: data/locale/${LOCALE_VERSION:-v2_ab4r}/ (outputs in output_<version>/ for a non-default version).
# ChimeraX renders are not made here.
set -euo pipefail
cd "$(dirname "$0")/code"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"   # $PYTHON, else python of the active env
[[ -n "$PY" ]] || { echo "python not found: conda activate sc-folding-figures (see ../environment.yml) or set PYTHON=" >&2; exit 1; }
echo "[Fig3] $("$PY" --version 2>&1) ($PY)"
OUTD=output; [[ "${LOCALE_VERSION:-v2_ab4r}" == "v2_ab4r" ]] || OUTD="output_$LOCALE_VERSION"
[[ "${DISTCT_NORM:-e2c}" == "e2c" ]] || OUTD="${OUTD}_distCT-$DISTCT_NORM"   # dist-to-CT normalisation (set_distct_norm.sh)
echo "[Fig3] locale version: ${LOCALE_VERSION:-v2_ab4r} (data/locale/<version>/), dist-to-CT norm: ${DISTCT_NORM:-e2c} -> $OUTD/"
for s in plot_fig3c.py plot_fig3d.py plot_fig3e.py; do
  echo "=== Fig3/code/$s"
  "$PY" "$s"
done
echo "[Fig3] done. Figures in Fig3/$OUTD/figures/, source data in Fig3/$OUTD/source_data/."
