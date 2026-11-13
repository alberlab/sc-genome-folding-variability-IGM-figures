#!/usr/bin/env bash
# Unpack the ChimeraX cmm/pdb files of Ext. Fig. 2c (ES model, beads coloured by their nuclear locale, ENV v2) into structures/h1/.
# Usage:  bash extract_structures.sh   (from any directory)
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p structures/h1
tar xzf data/structures/h1_occ06_extfig2c_locale_ENVv2_cmm_pdb.tar.gz -C structures/h1
echo "Done: structures/h1/cmm_pdb_locale/Struct<k>/, speckles/, h1_extfig2c.cxc"
