#!/usr/bin/env bash
# Fig. 4a / 4e / 4i (ChimeraX, not part of run.sh): re-render the PNGs from the saved sessions in output/figures/Fig4{a,e,i}/
# into output/figures/rerender/<panel>/ (the published PNGs next to the sessions are left untouched).
# Needs UCSF ChimeraX (auto-detected in /Applications, or set CHIMERAX=/path/to/ChimeraX). On macOS a ChimeraX window opens,
# renders and quits by itself (close any open ChimeraX first). The folder path must not contain spaces.
# The cmm/pdb inputs of the sessions are in data/structures/*.tar.gz (unpack with extract_structures.sh to inspect or
# rebuild a scene by hand: open the .cmm files of one structure, then <cell>_nucleus_ellipsoid.cxc).
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
echo "ChimeraX: $CX"
"$CX" --exit --script "$ROOT/code/chimerax/fig4_render_sessions.py $ROOT" > /dev/null 2>&1 || true
echo "run log:"; cat "$ROOT/output/figures/rerender/fig4_rerender.log" 2>/dev/null || echo "(no run log - ChimeraX did not run the script)"
