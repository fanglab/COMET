#!/usr/bin/env bash
# comet test: check the software environment, then run comet dmr and comet assign on the brain
# test data and compare the results with the expected ones. Called through the `comet` dispatcher.
set -euo pipefail
export LC_ALL=C

VERSION=${COMET_VERSION:?run through the comet command}
BIN=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)
ZENODO_RECORD=24064938
ZENODO=https://zenodo.org/records/$ZENODO_RECORD/files

usage() {
cat <<EOF
comet test (COMET v$VERSION) -- check the installation and reproduce the test results

Usage:
  comet test [--out DIR] [--threads N] [--ref-panel FILE] [--env-only]

Steps:
  1. environment  required tools and Python packages (with versions), Arial font
  2. comet dmr    AD-specific microglia DMRs on the bundled per-site BEDs (chr15:73.5-75.5 Mb, seconds)
  3. download     test BAMs (1.3 GB) and, unless given with --ref-panel, the brain reference panel (613 MB)
                  from Zenodo (DOI 10.5281/zenodo.24064938), md5-checked; files already present are not downloaded again
  4. comet assign three postmortem brain samples (PacBio HiFi, 0.5X; ~5 min with 16 threads)
Each step reports PASS or FAIL; the run stops at the first failure.

Options:
  --out DIR             working folder for the downloads and test outputs [comet_test]
  --threads INT         threads for comet assign / dmr [4]
  --ref-panel FILE      brain.atlas.bed.gz already downloaded (tabix-indexed); otherwise it is downloaded
                        into OUT/ref_panel
  --env-only            only step 1 (no download, no test runs)
  -h, --help            show this help
EOF
}

OUT=comet_test THREADS=4 REF= ENV_ONLY=0
while [ $# -gt 0 ]; do
  case $1 in  # an option given without its value
    --out|--threads|--ref-panel)
      [ $# -ge 2 ] && [[ $2 != --* ]] || { echo "ERROR: $1 requires a value (see comet test --help)" >&2; exit 1; } ;;
  esac
  case $1 in
    --out) OUT=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --ref-panel) REF=$2; shift 2 ;;
    --env-only) ENV_ONLY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    -v|--version) echo "comet $VERSION"; exit 0 ;;
    *) echo "ERROR: unknown option: $1 (see comet test --help)" >&2; exit 1 ;;
  esac
done
[[ "$THREADS" =~ ^[0-9]+$ ]] && [ "$THREADS" -ge 1 ] || { echo "ERROR: --threads must be a positive integer" >&2; exit 1; }

NFAIL=0
pass() { printf '  PASS  %s\n' "$*"; }
fail() { printf '  FAIL  %s\n' "$*"; NFAIL=$((NFAIL + 1)); }
warn() { printf '  WARN  %s\n' "$*"; }
stop_if_failed() {
  [ "$NFAIL" -eq 0 ] && return
  echo; echo "comet test: $NFAIL check(s) failed in $1; see the messages above${2:+ and $2}." >&2; exit 1
}

# ---------------------------------------------------------------- 1. environment
echo "COMET v$VERSION test"
echo; echo "[1/4] Environment"
tool() {   # tool NAME VERSION_COMMAND [optional-note]
  local name=$1 cmd=$2 note=${3:-} v
  if command -v "$name" >/dev/null 2>&1; then
    v=$(bash -c "$cmd" 2>&1 | grep -m1 -E '[0-9]+\.[0-9]+' || true)
    pass "$name ${v:-(version not reported)}"
  elif [ -n "$note" ]; then
    warn "$name not found ($note)"
  else
    fail "$name not found"
  fi
}
tool python3   'python3 --version'
tool samtools  'samtools --version | head -1'
tool bedtools  'bedtools --version'
tool bgzip     'bgzip --version | head -1'
tool tabix     'tabix --version | head -1'
tool metilene  'echo "(metilene reports no version; 0.2.9 is used by COMET)"'
tool perl      'perl -e "print \"perl \$^V\n\""'
tool modkit    'modkit --version' 'needed for ONT input only'
tool awk       'awk --version 2>/dev/null | head -1 || echo awk'
if command -v metilene >/dev/null 2>&1; then
  m=$(dirname "$(readlink -f "$(command -v metilene)")")/metilene_input.pl
  { [ -f "$m" ] || command -v metilene_input.pl >/dev/null 2>&1; } && pass "metilene_input.pl" || fail "metilene_input.pl not found next to metilene or on PATH"
fi
if command -v curl >/dev/null 2>&1 || command -v wget >/dev/null 2>&1; then pass "curl/wget"
else fail "curl or wget (needed to download the test data)"; fi
python3 - <<'PY' || NFAIL=$((NFAIL + 1))
import importlib, os, sys
bad = 0
for m in ('numpy', 'pandas', 'scipy', 'pysam', 'matplotlib'):
    try:
        print('  PASS  python %-10s %s' % (m, importlib.import_module(m).__version__))
    except Exception as e:
        print('  FAIL  python %-10s %s' % (m, e)); bad = 1
try:
    from matplotlib import font_manager
    for d in (os.environ.get('COMET_FONT_DIR', ''), os.path.join(sys.prefix, 'fonts'),
              os.path.join(os.environ.get('CONDA_PREFIX', ''), 'fonts'), os.path.expanduser('~/.fonts')):
        if d and os.path.isdir(d):
            for f in font_manager.findSystemFonts(fontpaths=[d]):
                try:
                    font_manager.fontManager.addfont(f)
                except Exception:
                    pass
    if 'Arial' in {f.name for f in font_manager.fontManager.ttflist}:
        print('  PASS  Arial font (plots)')
    else:
        print('  WARN  Arial font not found: plots use another font (install mscorefonts or set COMET_FONT_DIR)')
except Exception:
    pass
sys.exit(bad)
PY
stop_if_failed "the environment" "install the missing software (see README: Requirements and Test)"
[ "$ENV_ONLY" == 1 ] && { echo; echo "comet test: environment OK (--env-only)"; exit 0; }

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
[ -z "$REF" ] || [ ! -e "$REF" ] || REF=$(cd "$(dirname "$REF")" && pwd)/$(basename "$REF")

# ---------------------------------------------------------------- 2. comet dmr
echo; echo "[2/4] comet dmr (bundled test data)"
DMR_DATA=
for d in "$BIN/../share/comet/test/dmr" "$BIN/../example/data/dmr"; do
  [ -s "$d/brain_deep_AD_1.Microglia.CG.bed.gz" ] && { DMR_DATA=$(cd "$d" && pwd); break; }
done
if [ -z "$DMR_DATA" ]; then
  fail "test BEDs not found (share/comet/test/dmr or example/data/dmr next to the comet install)"
else
  D=$DMR_DATA/brain_deep
  rm -rf "$OUT/dmr"
  if comet dmr --target "$D"_AD_1.Microglia.CG.bed.gz,"$D"_AD_2.Microglia.CG.bed.gz --target-name AD_Microglia \
       --background "$D"_MCI_1.Microglia.CG.bed.gz,"$D"_MCI_2.Microglia.CG.bed.gz,"$D"_NCI_1.Microglia.CG.bed.gz,"$D"_NCI_2.Microglia.CG.bed.gz \
       --background-name MCI_NCI_Microglia --dmr-diff 0.2 --dmr-level individual --threads "$THREADS" --out "$OUT/dmr" \
       > "$OUT/dmr.log" 2>&1; then
    SIG=$OUT/dmr/01.dmr/AD_Microglia_vs_MCI_NCI_Microglia.metilene.sig.tsv
    n=$(($(wc -l < "$SIG") - 1))
    [ "$n" -eq 3 ] && pass "3 significant DMRs" || fail "$n significant DMRs (expected 3)"
    # the two AD-specific microglia DMRs at chr15:74,174,956-74,175,254, with the AD_1 level of each
    for e in 74174956:74175055:0.147059 74175090:74175254:0.263158; do
      IFS=: read -r s en lv <<< "$e"
      got=$(awk -F'\t' -v s="$s" -v e="$en" 'NR==1 {for (i=1; i<=NF; i++) if ($i=="mean_brain_deep_AD_1.Microglia.CG") c=i}
                                             NR>1 && $1=="chr15" && $2==s && $3==e {print $c}' "$SIG")
      if [ -n "$got" ] && awk -v a="$got" -v b="$lv" 'BEGIN{exit !((a-b)^2 < 1e-6)}'; then
        pass "AD microglia DMR chr15:$s-$en (AD_1 $got)"
      else
        fail "AD microglia DMR chr15:$s-$en (expected AD_1 $lv, got ${got:-no DMR})"
      fi
    done
  else
    fail "comet dmr stopped (see $OUT/dmr.log)"
  fi
fi
stop_if_failed "comet dmr" "$OUT/dmr.log"

# ---------------------------------------------------------------- 3. download
echo; echo "[3/4] Download test data (Zenodo)"
if command -v curl >/dev/null 2>&1; then fetch() { curl -fsSL --retry 3 -o "$2" "$1"; }
else fetch() { wget -q --tries=3 -O "$2" "$1"; }; fi
md5() { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
DATA=$OUT/data
mkdir -p "$DATA"
fetch "$ZENODO/md5sum.txt?download=1" "$DATA/.md5sum.txt" && grep -q 'brain_01.bam$' "$DATA/.md5sum.txt" \
  || { fail "could not read md5sum.txt from $ZENODO"; stop_if_failed "the download"; }
get() {   # get FILE DEST_DIR
  local f=$1 d=$2 want
  want=$(awk -v f="$f" '$2 == f {print $1}' "$DATA/.md5sum.txt")
  if [ -s "$d/$f" ] && [ "$(md5 "$d/$f")" == "$want" ]; then pass "$f (present)"; return; fi
  echo "        downloading $f"
  if fetch "$ZENODO/$f?download=1" "$d/$f.part" && [ "$(md5 "$d/$f.part")" == "$want" ]; then
    mv "$d/$f.part" "$d/$f"; pass "$f (md5 OK)"
  else
    rm -f "$d/$f.part"; fail "$f: download failed or md5 mismatch"
  fi
}
for s in brain_01 brain_02 brain_03; do get $s.bam "$DATA"; get $s.bam.csi "$DATA"; done
if [ -n "$REF" ]; then
  [ -s "$REF" ] && { [ -s "$REF.tbi" ] || [ -s "$REF.csi" ]; } && pass "brain reference panel $REF" \
    || fail "--ref-panel $REF: file or its tabix index not found"
else
  mkdir -p "$OUT/ref_panel"
  get brain.atlas.bed.gz "$OUT/ref_panel"; get brain.atlas.bed.gz.tbi "$OUT/ref_panel"
  REF=$OUT/ref_panel/brain.atlas.bed.gz
fi
stop_if_failed "the download"

# ---------------------------------------------------------------- 4. comet assign
echo; echo "[4/4] comet assign (3 PacBio HiFi brain samples; this takes a few minutes)"
if comet assign --platform pacbio --ref-panel "$REF" --assign-thres 0.65 \
     --bam "$DATA/brain_01.bam,$DATA/brain_02.bam,$DATA/brain_03.bam" \
     --threads "$THREADS" --out "$OUT/assign" \
     > "$OUT/assign.log" 2>&1; then
  S=$OUT/assign/summary
  # sample group expected-value: composition (frac_single) and genome-wide methylation
  check() {   # check TABLE COLUMN LABEL sample:group:expected ...
    local tab=$1 col=$2 lab=$3; shift 3
    for e in "$@"; do
      IFS=: read -r s g x <<< "$e"
      got=$(awk -F'\t' -v s="$s" -v g="$g" -v c="$col" 'NR==1 {for (i=1; i<=NF; i++) if ($i==c) k=i}
                                                        NR>1 && $1==s && $2==g {print $k}' "$tab")
      if [ -n "$got" ] && awk -v a="$got" -v b="$x" 'BEGIN{exit !((a-b)^2 < 1e-6)}'; then
        pass "$lab $s $g $got"
      else
        fail "$lab $s $g: expected $x, got ${got:-nothing}"
      fi
    done
  }
  check "$S/all_samples.composition.tsv" frac_single "composition" \
    brain_01:Exc_neuron:0.126544 brain_01:Inh_neuron:0.137372 brain_01:Astrocyte:0.101231 \
    brain_01:Microglia:0.288579 brain_01:Oligodendrocyte:0.273446 brain_01:OPC:0.072828 \
    brain_02:Exc_neuron:0.211836 brain_02:Inh_neuron:0.156506 brain_02:Astrocyte:0.085726 \
    brain_02:Microglia:0.311399 brain_02:Oligodendrocyte:0.172898 brain_02:OPC:0.061634 \
    brain_03:Exc_neuron:0.285822 brain_03:Inh_neuron:0.145150 brain_03:Astrocyte:0.077836 \
    brain_03:Microglia:0.309384 brain_03:Oligodendrocyte:0.131872 brain_03:OPC:0.049936
  check "$S/all_samples.gw_methylation.tsv" methylation "methylation" \
    brain_01:Bulk:0.794730 brain_01:Exc_neuron:0.758388 brain_01:Inh_neuron:0.838912 \
    brain_01:Astrocyte:0.752617 brain_01:Microglia:0.737356 brain_01:Oligodendrocyte:0.716113 \
    brain_01:OPC:0.738632 brain_02:Bulk:0.815050 brain_02:Exc_neuron:0.763355 \
    brain_02:Inh_neuron:0.842282 brain_02:Astrocyte:0.775735 brain_02:Microglia:0.752121 \
    brain_02:Oligodendrocyte:0.727852 brain_02:OPC:0.756629 brain_03:Bulk:0.795649 \
    brain_03:Exc_neuron:0.736563 brain_03:Inh_neuron:0.831866 brain_03:Astrocyte:0.767003 \
    brain_03:Microglia:0.745924 brain_03:Oligodendrocyte:0.711784 brain_03:OPC:0.741966
  for p in composition cell_type_methylation; do
    [ -s "$S/all_samples.$p.pdf" ] && [ -s "$S/all_samples.$p.png" ] && pass "plot all_samples.$p.pdf/png" \
      || fail "plot all_samples.$p.pdf/png not written (see $OUT/assign.log)"
  done
else
  fail "comet assign stopped (see $OUT/assign.log)"
fi
stop_if_failed "comet assign" "$OUT/assign.log"

echo; echo "comet test: all checks passed. Outputs in $OUT"
