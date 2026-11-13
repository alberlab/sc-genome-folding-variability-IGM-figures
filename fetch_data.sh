#!/usr/bin/env bash
# PLACEHOLDER - downloads the ChimeraX inputs and published renders from Zenodo (DOI to be added at publication).
# Not needed for the statistical plots (bash run_all.sh); needed only for Fig4/run_chimerax.sh, Fig5/run_chimerax.sh and
# ExtFig2/run_chimerax.sh (UCSF ChimeraX 1.11.1).
#
# What is downloaded (68 files, ~350 MB; full list with byte sizes and md5 in data/zenodo_manifest.tsv):
#   <Figure>/data/structures/*.tar.gz      cmm/pdb structures: Fig4 (8: <cell>_occ06_cmm_pdb_{by_chromosome,by_shannon}),
#                                          Fig5 (4: <cell>_occ06_fig5_SD_entropy_ENVv2_cmm_pdb), ExtFig2 (1: h1 ExtFig2c locale)
#   Fig4/output/figures/Fig4i/<cell>_nucleus_shell_pores.mrc   nuclear-envelope maps (opened by Fig4i, Fig5a/b, ExtFig2c scenes)
#   Fig4/output/figures/Fig4{a,e,i}/*.cxs  ChimeraX sessions = scene recipes re-rendered by Fig4/run_chimerax.sh
#   <Figure>/output/figures/{Fig4a,Fig4e,Fig4i,Fig5a,Fig5b,ExtFig2c}/*.{png,cxs}   published renders / sessions (reference)
# Each file goes to the repo-relative path in column 1 of the manifest; the Zenodo record stores it under its basename.
#
# Usage (after publication):  ZENODO_RECORD=<record id> bash fetch_data.sh      (or edit ZENODO_RECORD below)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
ZENODO_RECORD="${ZENODO_RECORD:-}"          # TODO: set the Zenodo record id (DOI 10.5281/zenodo.<id>) at publication
MAN="$ROOT/data/zenodo_manifest.tsv"
if [[ -z "$ZENODO_RECORD" ]]; then
  echo "Zenodo record not yet available (DOI to be added at publication). Files that will be fetched:" >&2
  tail -n +2 "$MAN" | cut -f1,2 | awk -F'\t' '{printf "  %-80s %8.1f MB\n", $1, $2/1e6}' >&2
  exit 1
fi
md5of() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
tail -n +2 "$MAN" | while IFS=$'\t' read -r path bytes md5 _; do
  dest="$ROOT/$path"; name="$(basename "$path")"
  if [[ -f "$dest" && "$(md5of "$dest")" == "$md5" ]]; then echo "ok    $path"; continue; fi
  mkdir -p "$(dirname "$dest")"
  url="https://zenodo.org/records/$ZENODO_RECORD/files/$(printf '%s' "$name" | sed 's/ /%20/g')?download=1"
  curl -fL --retry 3 -o "$dest" "$url"
  [[ "$(md5of "$dest")" == "$md5" ]] || { echo "md5 mismatch: $path" >&2; exit 1; }
  echo "got   $path"
done
echo "Done. Render with: bash Fig4/run_chimerax.sh; bash Fig5/run_chimerax.sh; bash ExtFig2/run_chimerax.sh"
