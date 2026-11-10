#!/usr/bin/env bash
# Download the COMET example BAMs (and optionally the brain reference panel) from Zenodo.
#
# usage: bash download_example_data.sh [--ref-panel-dir DIR] [--panels-only]
#   example BAMs           -> example/data/assign/  (brain_01, brain_02, brain_03)
#   --ref-panel-dir DIR    also download the brain reference panel (613 MB) into DIR;
#                          then: export COMET_REF_PANEL_DIR=DIR
#   --panels-only          with --ref-panel-dir: skip the example BAMs
# Files already present with the right md5 are not downloaded again.
set -euo pipefail
RECORD=24064938
URL=https://zenodo.org/records/$RECORD/files

PANEL_DIR= PANELS=brain BAMS=1
while [ $# -gt 0 ]; do
  case $1 in
    --ref-panel-dir) PANEL_DIR=$2; shift 2 ;;
    --panels-only) BAMS=0; shift ;;
    -h|--help) sed -n '2,10p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done
[ "$BAMS" == 1 ] || [ -n "$PANEL_DIR" ] || { echo "--panels-only needs --ref-panel-dir" >&2; exit 1; }
# panel folder relative to where the script was called from
[ -z "$PANEL_DIR" ] || { mkdir -p "$PANEL_DIR"; PANEL_DIR=$(cd "$PANEL_DIR" && pwd); }
cd "$(dirname "$0")"

if command -v curl >/dev/null 2>&1; then fetch() { curl -fL --retry 3 -o "$2" "$1"; }
elif command -v wget >/dev/null 2>&1; then fetch() { wget -q --tries=3 -O "$2" "$1"; }
else echo "ERROR: needs curl or wget" >&2; exit 1; fi
md5() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }

mkdir -p data/assign
fetch "$URL/md5sum.txt?download=1" data/assign/.md5sum.txt

get() {   # get FILE DEST_DIR
  local f=$1 d=$2 want
  want=$(awk -v f="$f" '$2 == f {print $1}' data/assign/.md5sum.txt)
  if [ -s "$d/$f" ] && [ "$(md5 "$d/$f")" == "$want" ]; then echo "ok (present): $d/$f"; return; fi
  echo "downloading $f"
  fetch "$URL/$f?download=1" "$d/$f.part"
  [ "$(md5 "$d/$f.part")" == "$want" ] || { echo "ERROR: md5 mismatch for $f" >&2; exit 1; }
  mv "$d/$f.part" "$d/$f"
}

if [ "$BAMS" == 1 ]; then
  for s in brain_01 brain_02 brain_03; do
    get $s.bam data/assign; get $s.bam.csi data/assign
  done
fi
if [ -n "$PANEL_DIR" ]; then
  mkdir -p "$PANEL_DIR"
  for t in $PANELS; do get $t.atlas.bed.gz "$PANEL_DIR"; get $t.atlas.bed.gz.tbi "$PANEL_DIR"; done
  echo "reference panel in $PANEL_DIR: export COMET_REF_PANEL_DIR=$PANEL_DIR"
fi
[ "$BAMS" == 1 ] && echo "done: example BAMs in $(pwd)/data/assign" || echo "done"
