#!/usr/bin/env bash
# Render Extended Data Fig. 2c (ES model coloured by nuclear locale, ENV v2) with ChimeraX from the terminal (not part of run.sh).
# Usage:  bash run_chimerax.sh
# Output: output/figures/ExtFig2c/"h1 <locale>.png" (2000 x 2000, supersample 3, transparent) + "h1 <locale>.cxs"
# for locale = spe, nuc, env, NL-act, NL-rep (structure 82). Needs UCSF ChimeraX (auto-detected in /Applications, or set
# CHIMERAX=/path/to/ChimeraX) and the nuclear envelope map ../Fig4/output/figures/Fig4i/h1_nucleus_shell_pores.mrc.
# On macOS ChimeraX cannot render without a window: a ChimeraX window opens, draws the scenes and quits by itself (close
# any open ChimeraX first). The folder path must not contain spaces. Overwrites the PNG / CXS shipped in output/figures/ExtFig2c.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
case "$ROOT" in *" "*) echo "path contains a space: $ROOT" >&2; exit 1;; esac
CX="${CHIMERAX:-}"
if [[ -z "$CX" ]]; then
  for c in /Applications/ChimeraX*.app/Contents/MacOS/ChimeraX "$HOME"/Applications/ChimeraX*.app/Contents/MacOS/ChimeraX \
           "$(command -v chimerax 2>/dev/null || true)" "$(command -v ChimeraX 2>/dev/null || true)"; do
    [[ -n "$c" && -x "$c" ]] && CX="$c"
  done
fi
[[ -n "$CX" && -x "$CX" ]] || { echo "ChimeraX not found; install it or set CHIMERAX=/path/to/ChimeraX" >&2; exit 1; }
echo "ChimeraX: $CX   cell: h1"
[[ -d "$ROOT/structures/h1/cmm_pdb_locale" ]] || bash "$ROOT/extract_structures.sh"
mkdir -p "$ROOT/output/figures/ExtFig2c"; rm -f "$ROOT/output/figures/ExtFig2c/extfig2c_run.log"
"$CX" --exit --script "$ROOT/code/chimerax/extfig2c_all.py $ROOT h1" > /dev/null 2>&1 || true
echo "run log:"; cat "$ROOT/output/figures/ExtFig2c/extfig2c_run.log" 2>/dev/null || echo "(no run log - ChimeraX did not run the script; close any open ChimeraX and retry)"
echo "files in output/figures/ExtFig2c:"; ls -1 "$ROOT/output/figures/ExtFig2c"
