#!/usr/bin/env bash
# Install a set of nuclear-locale tables as a selectable locale version, then rebuild data/checksums.md5.
# Usage:  bash set_locale_version.sh <version> <src_dir>
#   <src_dir>/locale_<cell>.tsv (cell = h1 DE HB mHEP; same columns and 15,453 rows as data/locale/v2_ab4r/) are copied to
#   data/locale/<version>/<cell>_occ06_locale_N565.tsv (one shared copy read by every figure).
#   Then:  LOCALE_VERSION=<version> bash run_all.sh   (or <Figure>/run.sh)  ->  outputs in <Figure>/output_<version>/
# Only the default version v2_ab4r (ENV v2, used in the manuscript) is shipped in this repository; alternative tables
# (e.g. the authors' v2_ab2p5r = 2.5 r A/B neighbourhood) are available from the authors on request.
# Without arguments: only rebuild data/checksums.md5 (md5 of every input table; paths relative to the repo root).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
md5line() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1"; else echo "$(md5 -q "$1")  $1"; fi; }
checksums() {
  ( cd "$ROOT" && { find data -type f; find Fig2 Fig3 Fig4 Fig5 ExtFig2 ExtFig4 ExtFig5 ExtFig6 -path '*/data/*' -type f; } 2>/dev/null |
    grep -v -e '/checksums\.md5$' -e '/zenodo_manifest\.tsv$' -e '/\._' -e '\.DS_Store$' -e '/data/structures/' | LC_ALL=C sort |
    while IFS= read -r f; do md5line "$f"; done > data/checksums.md5 )
}
if [[ $# -ge 1 ]]; then
  [[ $# -eq 2 ]] || { echo "usage: bash set_locale_version.sh <version> <src_dir>" >&2; exit 1; }
  V="$1"; SRC="$2"
  [[ "$V" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "bad version name '$V'" >&2; exit 1; }
  [[ "$V" != "v2_ab4r" ]] || { echo "v2_ab4r is the shipped default; choose another name" >&2; exit 1; }
  REF="$ROOT/data/locale/v2_ab4r"
  for c in h1 DE HB mHEP; do                                     # check all four before copying anything
    f="$SRC/locale_$c.tsv"; [[ -f "$f" ]] || { echo "missing $f" >&2; exit 1; }
    [[ "$(head -1 "$f")" == "$(head -1 "$REF/${c}_occ06_locale_N565.tsv")" ]] || { echo "$f: header differs from v2_ab4r" >&2; exit 1; }
    n=$(wc -l < "$f" | tr -d ' '); [[ "$n" == "15454" ]] || { echo "$f: $n lines, expected 15454 (15,453 bins + header)" >&2; exit 1; }
  done
  mkdir -p "$ROOT/data/locale/$V"
  for c in h1 DE HB mHEP; do cp -L "$SRC/locale_$c.tsv" "$ROOT/data/locale/$V/${c}_occ06_locale_N565.tsv"; done
  echo "data/locale/$V/ <- h1 DE HB mHEP"
fi
checksums
echo "data/checksums.md5 rebuilt ($(wc -l < "$ROOT/data/checksums.md5" | tr -d ' ') files; check: md5sum -c data/checksums.md5 from the repo root)"
