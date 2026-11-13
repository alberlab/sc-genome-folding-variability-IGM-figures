#!/usr/bin/env bash
# Reproduce the Figure 5 plotted panels (c, d, e) and the colour bars of the 5a / 5b renders.
# Usage:  bash run.sh   (from any directory). Python: $PYTHON if set, else `python` of the active environment
# (create it with ../environment.yml: conda env create -f environment.yml; conda activate sc-folding-figures).
# Locale tables: data/locale/${LOCALE_VERSION:-v2_ab4r}/ (outputs in output_<version>/ for a non-default version).
# ChimeraX renders are not made here (see run_chimerax.sh).
set -euo pipefail
cd "$(dirname "$0")/code"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"   # $PYTHON, else python of the active env
[[ -n "$PY" ]] || { echo "python not found: conda activate sc-folding-figures (see ../environment.yml) or set PYTHON=" >&2; exit 1; }
echo "[Fig5] $("$PY" --version 2>&1) ($PY)"
OUTD=output; [[ "${LOCALE_VERSION:-v2_ab4r}" == "v2_ab4r" ]] || OUTD="output_$LOCALE_VERSION"
[[ "${DISTCT_NORM:-e2c}" == "e2c" ]] || OUTD="${OUTD}_distCT-$DISTCT_NORM"   # dist-to-CT normalisation (set_distct_norm.sh)
echo "[Fig5] locale version: ${LOCALE_VERSION:-v2_ab4r} (data/locale/<version>/), dist-to-CT norm: ${DISTCT_NORM:-e2c} -> $OUTD/"
for s in plot_fig5ab_colorbars.py plot_fig5ce_per_locale.py plot_fig5d.py; do
  echo "=== Fig5/code/$s"
  "$PY" "$s"
done
echo "[Fig5] done. Figures in Fig5/$OUTD/figures/, source data in Fig5/$OUTD/source_data/."
