#!/usr/bin/env bash
# comet dmr: differentially methylated regions between two groups of per-site BED files
# (metilene). Called through the `comet` dispatcher.
set -euo pipefail
export LC_ALL=C

VERSION=${COMET_VERSION:?run through the comet command}
BIN=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)

usage() {
cat <<EOF
comet dmr (COMET v$VERSION) -- differentially methylated regions between two groups

Usage:
  comet dmr --target FILES --background FILES --dmr-diff 0.3 --out OUT_DIR [options]

Compares the target group against the background group (differences are target - background).
Input: per-site BED files as written by 'comet assign' (01.cell_type_assignment/BED/...):
chrom, start (0-based), end, n_mod, n_total; plain or gzipped. Typical use: one cell type
across two donor groups (e.g. Microglia of AD donors vs Microglia of control donors).

Required:
  --target FILES        target group: comma-separated BED files, or a text file with one path per line
  --background FILES    background group (same format); the target is compared against it
  --dmr-diff FLOAT      minimum |methylation difference| between the groups (fraction, e.g. 0.3);
                        passed to metilene (-d) and required of the final mean difference
  --out DIR             output directory

Optional:
  --dmr-level STR       group | individual [group]
                          group       the group means differ by >= --dmr-diff
                          individual  also every target file differs from every background file
                                      by >= --dmr-diff, in the same direction
  --target-name STR     label of the target group in outputs [target]
  --background-name STR label of the background group in outputs [background]
  --min-cpg INT         minimum CpGs (sites) in a DMR [5]
  --max-gap INT         maximum distance (bp) between neighbouring CpGs within one DMR [300]
                        (metilene -M)
  --min-depth FLOAT     minimum mean depth over the DMR, required in every sample used
                        (summed n_total >= min-depth x n_cpg) [3]
  --aggregate           pool the files of each group into one (summing n_mod and n_total per
                        site) and run 1 vs 1; for shallow cell-type data
  --alpha FLOAT         significance cutoff: metilene BH-adjusted q-value [0.05]
  --name STR            output file prefix [TARGET_vs_BACKGROUND]
  --threads INT         threads [4]
  --metilene PATH       metilene executable (metilene_input.pl expected next to it) [metilene]
  --keep-tmp            keep intermediate files
  -h, --help            show this help

metilene settings (as in the COMET manuscript): de novo mode, -X/-Y = ceil(n/2) samples per
group, -M max-gap, -c 2 (BH), -m min-cpg, -d dmr-diff.

Results go to OUT_DIR/01.dmr and OUT_DIR/02.run_info; 02.run_info/OUTPUT_README.txt
describes every file.
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*" >&2; }

# ---------------------------------------------------------------- arguments
TARGET= BACKGROUND= OUT= TOOL=metilene TNAME=target BNAME=background MIN_CPG=5 MIN_DEPTH=3 MAX_GAP=300
AGGREGATE=0 MIN_DIFF= ALPHA=0.05 LEVEL=group NAME= THREADS=4 METILENE=metilene
KEEP_TMP=0

[ $# -eq 0 ] && { usage; exit 1; }
CMDLINE="comet dmr$(printf ' %q' "$@")"
while [ $# -gt 0 ]; do
  case $1 in
    --target) TARGET=$2; shift 2 ;;
    --background) BACKGROUND=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --dmr-level) LEVEL=$2; shift 2 ;;
    --target-name) TNAME=$2; shift 2 ;;
    --background-name) BNAME=$2; shift 2 ;;
    --min-cpg) MIN_CPG=$2; shift 2 ;;
    --max-gap) MAX_GAP=$2; shift 2 ;;
    --min-depth) MIN_DEPTH=$2; shift 2 ;;
    --aggregate) AGGREGATE=1; shift ;;
    --dmr-diff) MIN_DIFF=$2; shift 2 ;;
    --alpha) ALPHA=$2; shift 2 ;;
    --name) NAME=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --metilene) METILENE=$2; shift 2 ;;
    --keep-tmp) KEEP_TMP=1; shift ;;
    -h|--help) usage; exit 0 ;;
    -v|--version) echo "comet $VERSION"; exit 0 ;;
    *) die "unknown option: $1 (see comet dmr --help)" ;;
  esac
done

LEVEL=$(echo "$LEVEL" | tr '[:upper:]' '[:lower:]')
case $LEVEL in group|individual) ;; *) die "--dmr-level must be group or individual" ;; esac
[ -n "$TARGET" ] || die "--target is required"
[ -n "$BACKGROUND" ] || die "--background is required"
[ -n "$OUT" ] || die "--out is required"
[[ "$TNAME" =~ ^[A-Za-z0-9_.-]+$ ]] && [[ "$BNAME" =~ ^[A-Za-z0-9_.-]+$ ]] \
  || die "--target-name / --background-name may only contain letters, digits, '_', '.', '-'"
[ "$TNAME" != "$BNAME" ] || die "--target-name and --background-name must differ"
[[ "$MIN_CPG" =~ ^[0-9]+$ ]] && [ "$MIN_CPG" -ge 1 ] || die "--min-cpg must be a positive integer"
[[ "$MAX_GAP" =~ ^[0-9]+$ ]] && [ "$MAX_GAP" -ge 1 ] || die "--max-gap must be a positive integer (bp)"
awk -v x="$MIN_DEPTH" 'BEGIN{exit !(x+0 >= 0)}' || die "--min-depth must be >= 0"
awk -v x="$ALPHA" 'BEGIN{exit !(x+0 > 0 && x+0 <= 1)}' || die "--alpha must be in (0, 1]"
[ -n "$MIN_DIFF" ] || die "--dmr-diff is required (minimum |methylation difference|, e.g. 0.3)"
awk -v x="$MIN_DIFF" 'BEGIN{exit !(x ~ /^[0-9.]+$/ && x+0 > 0 && x+0 < 1)}' || die "--dmr-diff must be a fraction in (0, 1)"
[ -n "$NAME" ] || NAME=${TNAME}_vs_${BNAME}

# read a plain or gzipped file
rd() {
  if [ "$(head -c2 "$1" | od -An -tx1 | tr -d ' \n')" == "1f8b" ]; then gzip -cd "$1"; else cat "$1"; fi
}
export -f rd
# first line of a file (a reader stopped early by head is not an error)
first_line() { { rd "$1" | head -n1; } 2>/dev/null || true; }
is_site_bed() { first_line "$1" | awk -F'\t' '{exit !(NF >= 5 && $2 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/)}'; }

# file lists: comma-separated paths, or one text file listing paths
read_list() {
  local spec=$1
  if [[ "$spec" != *,* ]] && [ -f "$spec" ] && ! is_site_bed "$spec"; then
    grep -v '^\s*\(#\|$\)' "$spec"
  else
    tr ',' '\n' <<< "$spec" | grep -v '^\s*$'
  fi
}
mapfile -t T_FILES < <(read_list "$TARGET")
mapfile -t B_FILES < <(read_list "$BACKGROUND")
[ ${#T_FILES[@]} -ge 1 ] || die "no files in --target"
[ ${#B_FILES[@]} -ge 1 ] || die "no files in --background"

# sample labels from file names (without .bed / .bed.gz); must be unique
label() { local b; b=$(basename "$1"); b=${b%.gz}; echo "${b%.bed}"; }
declare -A SEEN=()
for f in "${T_FILES[@]}" "${B_FILES[@]}"; do
  [ -f "$f" ] || die "file not found: $f"
  is_site_bed "$f" || die "$f is not a per-site BED (chrom, start, end, n_mod, n_total)"
  l=$(label "$f")
  [ -z "${SEEN[$l]:-}" ] || die "two input files share the sample name '$l'; rename one"
  SEEN[$l]=1
done

check_deps() {
  local miss=()
  for t in python3 bedtools awk sort; do command -v "$t" >/dev/null 2>&1 || miss+=("$t"); done
  command -v "$METILENE" >/dev/null 2>&1 || miss+=("metilene")
  command -v perl >/dev/null 2>&1 || miss+=("perl")
  [ ${#miss[@]} -eq 0 ] || die "missing: ${miss[*]}"
}
check_deps

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
RES=$OUT/01.dmr
RUNINFO=$OUT/02.run_info
LOGS=$RUNINFO/logs
TMP=$OUT/tmp
rm -rf "$TMP"
mkdir -p "$RES" "$LOGS" "$TMP/in"

# ---------------------------------------------------------------- inputs
log "Preparing ${#T_FILES[@]} target and ${#B_FILES[@]} background files ($TOOL)"
{
  for f in "${T_FILES[@]}"; do printf '%s\ttarget\t%s\n' "$(label "$f")" "$(readlink -f "$f")"; done
  for f in "${B_FILES[@]}"; do printf '%s\tbackground\t%s\n' "$(label "$f")" "$(readlink -f "$f")"; done
} > "$TMP/inputs.tsv"
# sorted per-site files with coverage, one per sample
export TMP
# shellcheck disable=SC2016
cut -f1,3 "$TMP/inputs.tsv" | xargs -P "$THREADS" -L1 bash -c '
  set -o pipefail; rd "$1" | awk -F"\t" -v OFS="\t" "\$5 > 0 {print \$1, \$2, \$3, \$4, \$5}" |
    sort -k1,1 -k2,2n -S 1G > "$TMP/in/$0.bed"' || die "could not read the input BED files"

if [ "$AGGREGATE" == "1" ]; then
  log "Aggregating each group into one file"
  for g in target background; do
    n=$([ $g == target ] && echo "$TNAME" || echo "$BNAME")
    awk -F'\t' -v g=$g -v d="$TMP/in" '$2==g {printf "%s/%s.bed\n", d, $1}' "$TMP/inputs.tsv" |
      xargs sort -m -k1,1 -k2,2n |
      awk -F'\t' -v OFS='\t' '{k=$1 OFS $2 OFS $3; if (k!=p) {if (p!="") print p, m, t; p=k; m=0; t=0} m+=$4; t+=$5}
           END {if (p!="") print p, m, t}' > "$TMP/agg.$n.bed"
  done
  rm -f "$TMP"/in/*.bed
  mv "$TMP/agg.$TNAME.bed" "$TMP/in/$TNAME.bed"; mv "$TMP/agg.$BNAME.bed" "$TMP/in/$BNAME.bed"
  printf '%s\ttarget\t%s\n%s\tbackground\t%s\n' "$TNAME" "$TMP/in/$TNAME.bed" "$BNAME" "$TMP/in/$BNAME.bed" > "$TMP/sheet.tsv"
else
  awk -F'\t' -v OFS='\t' -v d="$TMP/in" '{print $1, $2, d "/" $1 ".bed"}' "$TMP/inputs.tsv" > "$TMP/sheet.tsv"
fi
NT=$(awk -F'\t' '$2=="target"' "$TMP/sheet.tsv" | wc -l)
NB=$(awk -F'\t' '$2=="background"' "$TMP/sheet.tsv" | wc -l)

# ---------------------------------------------------------------- DMR calling
RAW=$TMP/raw.tsv
case $TOOL in
  metilene)
    MIN_T=$(( (NT + 1) / 2 )); MIN_B=$(( (NB + 1) / 2 ))
    log "metilene: $NT vs $NB samples (-X $MIN_T -Y $MIN_B -m $MIN_CPG -M $MAX_GAP -d $MIN_DIFF -c 2)"
    mkdir -p "$TMP/bg"
    while IFS=$'\t' read -r s g f; do
      awk -F'\t' -v OFS='\t' '{printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, $4/$5}' "$f" > "$TMP/bg/$s.bg"
    done < "$TMP/sheet.tsv"
    in1=$(awk -F'\t' -v d="$TMP/bg" '$2=="target" {printf "%s%s/%s.bg", (n++ ? "," : ""), d, $1}' "$TMP/sheet.tsv")
    in2=$(awk -F'\t' -v d="$TMP/bg" '$2=="background" {printf "%s%s/%s.bg", (n++ ? "," : ""), d, $1}' "$TMP/sheet.tsv")
    MET_DIR=$(dirname "$(readlink -f "$(command -v "$METILENE")")")
    MET_IN=$MET_DIR/metilene_input.pl
    [ -f "$MET_IN" ] || MET_IN=$(command -v metilene_input.pl || true)
    [ -n "$MET_IN" ] && [ -f "$MET_IN" ] || die "metilene_input.pl not found next to metilene or on PATH"
    perl "$MET_IN" --in1 "$in1" --in2 "$in2" --h1 "$TNAME" --h2 "$BNAME" --out "$TMP/metilene.input" \
      -b "$(command -v bedtools)" > "$LOGS/metilene_input.log" 2>&1 || die "metilene_input.pl failed, see $LOGS/metilene_input.log"
    rm -rf "$TMP/bg"
    "$METILENE" -X "$MIN_T" -Y "$MIN_B" -m "$MIN_CPG" -M "$MAX_GAP" -d "$MIN_DIFF" -c 2 -t "$THREADS" \
      -a "$TNAME" -b "$BNAME" "$TMP/metilene.input" 2> "$LOGS/metilene.log" | sort -k1,1 -k2,2n -k3,3n > "$TMP/raw.body" \
      || die "metilene failed, see $LOGS/metilene.log"
    { printf 'chr\tstart\tend\tq_value\tmean_diff\tnCpG\tMWU_pval\tKS_pval\tmean_%s\tmean_%s\n' "$TNAME" "$BNAME"
      cat "$TMP/raw.body"; } > "$RAW"
    ;;
esac

# ---------------------------------------------------------------- standardise, map counts, filter
log "Mapping per-sample counts onto DMRs; filters: |diff| >= $MIN_DIFF ($LEVEL level), >= $MIN_CPG CpGs, mean depth >= $MIN_DEPTH in all $((NT + NB)) samples"
python3 "$BIN/comet_dmr_post.py" standardize "$TOOL" "$RAW" "$TMP/std.bed" --alpha "$ALPHA"
# union of covered sites (n_cpg), then summed n_mod / n_total of each sample
cut -f3 "$TMP/sheet.tsv" | xargs sort -m -k1,1 -k2,2n | cut -f1-3 | uniq > "$TMP/union.bed"
cp "$TMP/std.bed" "$TMP/mapped.bed"
if [ -s "$TMP/std.bed" ]; then
  bedtools map -a "$TMP/mapped.bed" -b "$TMP/union.bed" -c 1 -o count > "$TMP/m.tmp" && mv "$TMP/m.tmp" "$TMP/mapped.bed"
  while IFS=$'\t' read -r s g f; do
    bedtools map -a "$TMP/mapped.bed" -b "$f" -c 4,5 -o sum,sum > "$TMP/m.tmp" && mv "$TMP/m.tmp" "$TMP/mapped.bed"
  done < "$TMP/sheet.tsv"
fi
PRE=$RES/$NAME.$TOOL
python3 "$BIN/comet_dmr_post.py" finalize "$TMP/std.bed" "$TMP/mapped.bed" "$TMP/sheet.tsv" "$PRE" \
  --min-cpg "$MIN_CPG" --min-depth "$MIN_DEPTH" --min-diff "$MIN_DIFF" --dmr-level "$LEVEL" --names "$TNAME,$BNAME" | tee "$LOGS/summary.txt"
cp "$RAW" "$PRE.raw.tsv"

# ---------------------------------------------------------------- run info
tool_version() {
  echo "metilene (reports no version; COMET is tested with 0.2.9)"
}
{
  echo "COMET v$VERSION  $(date '+%F %T')"
  echo "command: $CMDLINE"
  echo
  printf '%-18s %s\n' dmr_tool "$TOOL" target_name "$TNAME" background_name "$BNAME" \
    dmr_level "$LEVEL" n_target_files "${#T_FILES[@]}" n_background_files "${#B_FILES[@]}" aggregate "$([ "$AGGREGATE" == 1 ] && echo yes || echo no)" \
    dmr_diff "$MIN_DIFF" min_cpg "$MIN_CPG" max_gap "$MAX_GAP" min_depth "$MIN_DEPTH" alpha "$ALPHA" threads "$THREADS"
  [ "$TOOL" == "metilene" ] && printf '%-18s %s\n' metilene_args "-X $MIN_T -Y $MIN_B -m $MIN_CPG -M $MAX_GAP -d $MIN_DIFF -c 2"
  echo
  echo "tools:"
  echo "  $(tool_version)"
  echo "  $(bedtools --version)"
  echo "  $(python3 --version 2>&1)"
} > "$RUNINFO/run_parameters.txt"
{ printf 'sample\tgroup\tfile\n'; cat "$TMP/inputs.tsv"; } > "$RUNINFO/input_files.tsv"

TOOLSIG="metilene BH q-value < $ALPHA (metilene also requires |diff| >= $MIN_DIFF)"
{
  cat <<EOF
COMET v$VERSION comet dmr output guide
$TOOL | target '$TNAME' ($NT sample$([ "$NT" -gt 1 ] && echo s)) vs background '$BNAME' ($NB sample$([ "$NB" -gt 1 ] && echo s))$([ "$AGGREGATE" == 1 ] && echo " | files of each group aggregated")
Parameters and tool version: 02.run_info/run_parameters.txt; inputs: 02.run_info/input_files.tsv

01.dmr/
  $NAME.$TOOL.all.tsv
      every region reported by $TOOL, one row each:
        chrom, start, end            BED coordinates (0-based start, end exclusive)
        n_cpg                        CpGs in the region as counted by the tool; used for
                                     --min-cpg and the depth filter
        n_cpg_any_sample             CpGs in the region covered in at least one sample
        mean_$TNAME, mean_$BNAME       mean methylation (fraction) of each group
        diff_${TNAME}_minus_$BNAME     mean_$TNAME - mean_$BNAME
        p_value, q_value, statistic  metilene MWU p-value, BH q-value, 2D-KS p-value
        tool_significant             1 if $TOOLSIG
        min_mean_depth               lowest per-sample mean depth (summed n_total / n_cpg)
        min_individual_diff          smallest gap between any target and any background sample,
                                     in the direction of the group difference (sample level =
                                     summed n_mod / n_total; NA if a sample has no coverage)
        pass_dmr_diff                1 if |diff| >= $MIN_DIFF (--dmr-diff)
        pass_individual_diff         1 if min_individual_diff >= $MIN_DIFF
        pass_min_cpg                 1 if n_cpg >= $MIN_CPG
        pass_min_depth               1 if every sample has mean depth >= $MIN_DEPTH
        significant                  1 if tool_significant, pass_dmr_diff, pass_min_cpg and
                                     pass_min_depth are 1$([ "$LEVEL" == individual ] && echo ", and pass_individual_diff
                                     (--dmr-level individual)" || echo " (--dmr-level group)")
        mean_<sample>                methylation level of each input file over the region
                                     (summed n_mod / summed n_total; NA if not covered)
        meth_<sample>, total_<sample>  summed n_mod and n_total of each sample over the region
  $NAME.$TOOL.sig.tsv     rows of all.tsv with significant = 1
  $NAME.$TOOL.sig.bed     the same regions as BED; name = $TNAME:<difference>
  $NAME.$TOOL.raw.tsv     unmodified $TOOL output

02.run_info/
  OUTPUT_README.txt      this guide
  run_parameters.txt     command line, parameters, tool versions
  input_files.tsv        sample name, group and path of every input file
  logs/                  tool logs and summary.txt (region counts)
EOF
} > "$RUNINFO/OUTPUT_README.txt"

[ "$KEEP_TMP" == "1" ] || rm -rf "$TMP"
log "Done. Results in $RES (file guide: 02.run_info/OUTPUT_README.txt)"
