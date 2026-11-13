#!/usr/bin/env bash
# Install an alternative dist-to-CT normalisation as a selectable table set, then rebuild data/checksums.md5.
# Usage:  bash set_distct_norm.sh <norm> <src_dir>        norm = r | Req   (e2c = shipped default in data/distCT/e2c/)
#   <src_dir>/stats/<cell>_occ06_distCT_top30_2p5r_stats.tsv   (cell = h1 DE HB mHEP; 15,453 rows + header, same columns)
#   <src_dir>/npz/h1_occ06_chr11_5loci_distCT.npz, <src_dir>/npz/h1_occ06_chr2_3loci_distCT.npz   (per-locus samples)
#   <src_dir>/SDrad_vs_SDdistCT_4stage_perbin.tsv                                                 (ExtFig2b right)
# are copied to data/distCT/<norm>/ (one shared copy read by every figure).
# Then:  DISTCT_NORM=<norm> bash run_all.sh  -> outputs in <Figure>/output_distCT-<norm>/
#        (combined with LOCALE_VERSION=<v>: output_<v>_distCT-<norm>/). ExtFig5 is not switched.
# Only the default e2c tables (signed distance / core centroid-to-surface distance, used in the manuscript) are shipped in
# this repository; the r / Req tables are available from the authors on request.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
[[ $# -eq 2 ]] || { echo "usage: bash set_distct_norm.sh <r|Req> <src_dir>" >&2; exit 1; }
N="$1"; SRC="$2"
[[ "$N" == "r" || "$N" == "Req" ]] || { echo "norm must be r or Req (e2c = shipped default)" >&2; exit 1; }
REF="$ROOT/data/distCT/e2c"
for c in h1 DE HB mHEP; do
  f="$SRC/stats/${c}_occ06_distCT_top30_2p5r_stats.tsv"; [[ -f "$f" ]] || { echo "missing $f" >&2; exit 1; }
  [[ "$(head -1 "$f")" == "$(head -1 "$REF/${c}_occ06_distCT_top30_2p5r_stats.tsv")" ]] || { echo "$f: header differs" >&2; exit 1; }
  n=$(wc -l < "$f" | tr -d ' '); [[ "$n" == "15454" ]] || { echo "$f: $n lines, expected 15454" >&2; exit 1; }
done
for f in npz/h1_occ06_chr11_5loci_distCT.npz npz/h1_occ06_chr2_3loci_distCT.npz SDrad_vs_SDdistCT_4stage_perbin.tsv; do
  [[ -f "$SRC/$f" ]] || { echo "missing $SRC/$f" >&2; exit 1; }
done
D="$ROOT/data/distCT/$N"; mkdir -p "$D"
for c in h1 DE HB mHEP; do cp -L "$SRC/stats/${c}_occ06_distCT_top30_2p5r_stats.tsv" "$D/"; done
cp -L "$SRC/npz/h1_occ06_chr11_5loci_distCT.npz" "$SRC/npz/h1_occ06_chr2_3loci_distCT.npz" "$SRC/SDrad_vs_SDdistCT_4stage_perbin.tsv" "$D/"
echo "data/distCT/$N/ <- 4 stats tables + 2 npz + SDrad_vs_SDdistCT_4stage_perbin.tsv"
bash "$ROOT/set_locale_version.sh"          # no arguments: only rebuilds data/checksums.md5
