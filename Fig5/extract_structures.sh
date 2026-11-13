#!/usr/bin/env bash
# Unpack the ChimeraX cmm/pdb files of Fig. 5a/b into structures/<cell>/ (~360 MB).
# Usage:  bash extract_structures.sh   (from any directory)
# cmm_pdb_SD/ (Fig. 5a; B-factor = SD(dist to CT)), cmm_pdb_entropy/ (Fig. 5b; B-factor = locale entropy, ENV v2 locale),
# speckles/ (Fig. 5a speckle centroids), <cell>_fig5.cxc. Structures 82, 399, 427, 602, 833 (renders use 82).
set -euo pipefail
cd "$(dirname "$0")"
for cell in h1 DE HB mHEP; do
  mkdir -p "structures/${cell}"
  tar xzf "data/structures/${cell}_occ06_fig5_SD_entropy_ENVv2_cmm_pdb.tar.gz" -C "structures/${cell}"
done
echo "Done: structures/<cell>/cmm_pdb_{SD,entropy}/, speckles/, <cell>_fig5.cxc"
