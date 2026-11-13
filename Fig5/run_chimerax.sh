#!/usr/bin/env bash
# Render Figure 5a / 5b structures with ChimeraX from the terminal (not part of run.sh).
# Usage:  bash run_chimerax.sh [SD] [entropy]   (default: SD, then entropy)     metrics: SD (-> output/figures/Fig5a) | entropy (-> Fig5b)
# For each stage model ES (h1), DE, HB, HLC (mHEP): "<cell> <metric>.png" (2000 x 2000, supersample 3, transparent) and
# "<cell> <metric>.cxs"; with SD also "speckle legend h1.png". Needs UCSF ChimeraX (auto-detected in /Applications,
# or set CHIMERAX=/path/to/ChimeraX). On macOS ChimeraX cannot render without a window: a ChimeraX window opens, draws
# the scenes and quits by itself (close any open ChimeraX first). The package path must not contain spaces.
set -euo pipefail
# drop a pasted trailing comment ("# ...": zsh passes it as arguments) and check the arguments
ARGS=(); for a in "$@"; do [[ "$a" == \#* ]] && break; ARGS+=("$a"); done
ALLOWED=" SD entropy "
for a in ${ARGS[@]+"${ARGS[@]}"}; do
  [[ "$ALLOWED" == *" $a "* ]] || { echo "unknown argument '$a' (allowed:$ALLOWED)" >&2; exit 1; }
done
set -- ${ARGS[@]+"${ARGS[@]}"}
ROOT="$(cd "$(dirname "$0")" && pwd)"
case "$ROOT" in *" "*) echo "package path contains a space: $ROOT — move the package" >&2; exit 1;; esac
[[ $# -ge 1 ]] || set -- SD entropy                                   # default: both panels
CX="${CHIMERAX:-}"
if [[ -z "$CX" ]]; then
  for c in /Applications/ChimeraX*.app/Contents/MacOS/ChimeraX "$HOME"/Applications/ChimeraX*.app/Contents/MacOS/ChimeraX \
           "$(command -v chimerax 2>/dev/null || true)" "$(command -v ChimeraX 2>/dev/null || true)"; do
    [[ -n "$c" && -x "$c" ]] && CX="$c"
  done
fi
[[ -n "$CX" && -x "$CX" ]] || { echo "ChimeraX not found; install it or set CHIMERAX=/path/to/ChimeraX" >&2; exit 1; }
echo "ChimeraX: $CX   metrics: $*"
[[ -d "$ROOT/structures/mHEP/cmm_pdb_entropy" ]] || bash "$ROOT/extract_structures.sh"
for M in "$@"; do                                                   # one ChimeraX run per metric (SD -> Fig5a, entropy -> Fig5b)
  PANEL=Fig5a; [[ "$M" == "entropy" ]] && PANEL=Fig5b
  mkdir -p "$ROOT/output/figures/$PANEL"; rm -f "$ROOT/output/figures/$PANEL/fig5_run.log"
  "$CX" --exit --script "$ROOT/code/chimerax/fig5ab_all.py $ROOT $M" > /dev/null 2>&1 || true
  echo "run log ($M):"; cat "$ROOT/output/figures/$PANEL/fig5_run.log" 2>/dev/null || echo "(no run log - ChimeraX did not run the script; close any open ChimeraX and retry)"
  echo "files in output/figures/$PANEL:"; ls -1 "$ROOT/output/figures/$PANEL"
done
