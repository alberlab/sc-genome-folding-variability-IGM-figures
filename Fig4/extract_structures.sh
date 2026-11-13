#!/usr/bin/env bash
# Unpack the ChimeraX cmm/pdb files of Fig. 4a/e/i into structures/<cell>/{by_chromosome,by_shannon}/ (~420 MB).
# Usage:  bash extract_structures.sh   (from any directory)
set -euo pipefail
cd "$(dirname "$0")"
for cell in h1 DE HB mHEP; do
  for kind in by_chromosome by_shannon; do
    out="structures/${cell}/${kind}"; mkdir -p "$out"
    tar xzf "data/structures/${cell}_occ06_cmm_pdb_${kind}.tar.gz" -C "$out"
  done
done
echo "Done: structures/<cell>/by_chromosome (Fig. 4a/e) and structures/<cell>/by_shannon (Fig. 4i)."
echo "In ChimeraX: open the .cmm files of one structure folder, then 'open <cell>_nucleus_ellipsoid.cxc' for the envelope."
