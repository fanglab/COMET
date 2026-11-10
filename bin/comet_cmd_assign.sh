#!/usr/bin/env bash
# comet assign: read-level cell-type assignment for long-read (PacBio HiFi / ONT)
# methylation data. Called through the `comet` dispatcher.
set -euo pipefail

VERSION=${COMET_VERSION:?run through the comet command}
BIN=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)

usage() {
cat <<EOF
comet assign (COMET v$VERSION) -- assign long reads to their cell type of origin from CpG methylation

Usage:
  comet assign --platform {pacbio,ont} --bam SAMPLE.bam \\
        --assign-thres 0.65 --out OUT_DIR [options]

Reference panel: brain (K=6): Exc_neuron, Inh_neuron, Astrocyte, Microglia, Oligodendrocyte, OPC

Required:
  --platform STR        pacbio | ont
  --bam FILE[,FILE..]   aligned, coordinate-sorted, indexed modBAM (MM/ML tags; hg38, chr1..chr22);
                        several samples: comma-separated BAMs or a text file with one BAM per line
  --assign-thres FLOAT  posterior cutoff theta for assigning a read
                        (recommended from our simulations: 0.65 for brain)
  --out DIR             output directory (several BAMs: one sub-folder per sample + summary/)
  --fasta FILE          reference genome FASTA (required for --platform ont, used by modkit)

Optional:
  --output-format STR   cell-type outputs: bam | bed-call | both [both]
                          bam      reads of each cell type as an indexed BAM
                          bed-call per-site calls of each read set (bgzipped, tabix-indexed):
                                   chrom, start, end, n_mod, n_total
                                   pacbio: SAMPLE.TYPE.CG.bed.gz from COMET's read-level calls
                                   ont:    SAMPLE.TYPE.<mod>.bed.gz from modkit pileup
  --ont-mods LIST       ont: modifications to call in bed-call output, comma-separated from
                        5mCG, 5hmCG, 5mCH, 5hmCH [5mCG]. The BAM is checked first: 5hmC needs a
                        5mCG_5hmCG / 5mC_5hmC basecalling model, CH needs an all-context model
                        (5mC / 5mC_5hmC). Read assignment always uses CpG 5mC (+ 5hmC if present).
  --ml-cutoff INT       methylation call confidence cutoff tau0 on the 0-255 ML scale [180]
  --threads INT         parallel workers / threads [4]
  --sample STR          sample name used to prefix outputs [BAM basename]
                        (several BAMs: comma-separated, one per BAM)
  --no-plot             do not draw the cell-type composition and cell-type methylation plots
  --meth-plot-samples LIST  samples shown in the cell-type methylation plot [all]
                        (e.g. to leave out cell lines or sorted samples)
  --cell-types STR      comma-separated cell-type names in atlas column order
                        [brain cell-type names]
  --ref-panel FILE      reference panel (atlas) BED, e.g. brain.atlas.bed.gz; bgzipped and
                        tabix-indexed files are used directly, others are indexed once into OUT
                        [<--ref-panel-dir>/brain.atlas.bed.gz]
  --ref-panel-dir DIR   folder holding brain.atlas.bed.gz(.tbi), used when --ref-panel is not given
                        [\$COMET_REF_PANEL_DIR, else <install>/share/comet/ref_panel]
  --no-lineage          brain: turn off lineage rescue of unassigned reads (on by default:
                        Neuron = Exc + Inh, Macroglia = Astrocyte + Oligodendrocyte + OPC)
  --lineage-thres FLOAT summed-posterior cutoff for lineage rescue [0.65]
  --lineage-margin FLOAT margin over the best outside cell type for lineage rescue [0.10]
  --min-cpg INT         minimum informative CpGs per read [5]
  --min-read-len INT    minimum read length [7500 pacbio, 500 ont]
  --chunk-size INT      reads per chunk (pacbio) or calls per chunk (ont) [10000 / 500000]
  --no-pn-calibration   pacbio: use tau0 for every pass number instead of calibrating
  --pn-table FILE       pacbio: use an existing pass-number table instead of calibrating
  --hmc-cutoff INT      ont: a 5hmC call with ML >= this replaces the 5mC call [225]
  --pileup-threshold F  ont: modkit pileup --filter-threshold for CpG calls [0.70]
  --ch-pileup-threshold F ont: modkit pileup --filter-threshold for CH calls [0.95]
  --modkit PATH         modkit executable [modkit]
  --keep-tmp            keep intermediate chunk files
  --step STR            all | prepare | chunks | finalize [all]
  --chunk-id ID         with --step chunks: run only this chunk (for cluster array jobs)
  -h, --help            show this help
  -v, --version         show version

Results go to OUT_DIR/01.cell_type_assignment, 02.cell_type_quantification and 03.run_info;
03.run_info/OUTPUT_README.txt describes every file a run produced. With several BAMs, each sample
is run in turn into OUT_DIR/<sample>/ and OUT_DIR/summary/ holds the combined tables and plots.
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*" >&2; }

# ---------------------------------------------------------------- arguments
PLATFORM= TISSUE=brain BAM= REF= REF_DIR=${COMET_REF_PANEL_DIR:-} ASSIGN_THRES= OUT= FASTA=
ML_CUTOFF=180 THREADS=4 SAMPLE= CELL_TYPES= MIN_CPG=5 MIN_READ_LEN= CHUNK_SIZE=
NO_LINEAGE=0 LINEAGE_THRES=0.65 LINEAGE_MARGIN=0.10 NO_PN_CAL=0 USER_PN_TABLE= HMC_CUTOFF=225
PILEUP_THRES=0.70 CH_PILEUP_THRES=0.95 MODS=5mCG MODKIT=modkit OUTPUT_FORMAT=both KEEP_TMP=0 STEP=all CHUNK_ID=
NO_PLOT=0 METH_PLOT_SAMPLES=
PASS=()   # options handed on to the per-sample runs of a multi-BAM call

[ $# -eq 0 ] && { usage; exit 1; }
CMDLINE="comet assign$(printf ' %q' "$@")"
while [ $# -gt 0 ]; do
  case $1 in
    --bam|--out|--sample|--meth-plot-samples|--no-plot|-h|--help|-v|--version) ;;
    -*) PASS+=("$1"); case ${2:-} in ""|--*) ;; *) PASS+=("$2") ;; esac ;;
  esac
  case $1 in
    --platform) PLATFORM=$2; shift 2 ;;
    --bam) BAM=$2; shift 2 ;;
    --ref-panel) REF=$2; shift 2 ;;
    --ref-panel-dir) REF_DIR=$2; shift 2 ;;
    --assign-thres) ASSIGN_THRES=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --fasta) FASTA=$2; shift 2 ;;
    --output-format) OUTPUT_FORMAT=$2; shift 2 ;;
    --ml-cutoff) ML_CUTOFF=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --sample) SAMPLE=$2; shift 2 ;;
    --cell-types) CELL_TYPES=$2; shift 2 ;;
    --no-lineage) NO_LINEAGE=1; shift ;;
    --lineage-thres) LINEAGE_THRES=$2; shift 2 ;;
    --lineage-margin) LINEAGE_MARGIN=$2; shift 2 ;;
    --min-cpg) MIN_CPG=$2; shift 2 ;;
    --min-read-len) MIN_READ_LEN=$2; shift 2 ;;
    --chunk-size) CHUNK_SIZE=$2; shift 2 ;;
    --no-pn-calibration) NO_PN_CAL=1; shift ;;
    --pn-table) USER_PN_TABLE=$2; shift 2 ;;
    --hmc-cutoff) HMC_CUTOFF=$2; shift 2 ;;
    --pileup-threshold) PILEUP_THRES=$2; shift 2 ;;
    --ch-pileup-threshold) CH_PILEUP_THRES=$2; shift 2 ;;
    --ont-mods) MODS=$2; shift 2 ;;
    --modkit) MODKIT=$2; shift 2 ;;
    --keep-tmp) KEEP_TMP=1; shift ;;
    --step) STEP=$2; shift 2 ;;
    --chunk-id) CHUNK_ID=$2; shift 2 ;;
    --no-plot) NO_PLOT=1; shift ;;
    --meth-plot-samples) METH_PLOT_SAMPLES=$2; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    -v|--version) echo "comet $VERSION"; exit 0 ;;
    *) die "unknown option: $1 (see comet assign --help)" ;;
  esac
done

PLATFORM=$(echo "$PLATFORM" | tr '[:upper:]' '[:lower:]')
[ "$PLATFORM" == "pacbio" ] || [ "$PLATFORM" == "ont" ] || die "--platform must be pacbio or ont"
[ -n "$BAM" ] || die "--bam is required"
TISSUE=$(echo "$TISSUE" | tr '[:upper:]' '[:lower:]')
LINEAGE_LIST=
case $TISSUE in
  brain) DEF_TYPES=Exc_neuron,Inh_neuron,Astrocyte,Microglia,Oligodendrocyte,OPC
         [ "$NO_LINEAGE" == "1" ] || LINEAGE_LIST="Neuron=1,2 Macroglia=3,5,6" ;;
  *) die "--tissue must be brain" ;;
esac
if [ "$TISSUE" != "other" ]; then
  if [ -n "$REF" ]; then
    [ -s "$REF" ] || die "--ref-panel: file not found: $REF"
  else
    [ -n "$REF_DIR" ] || REF_DIR=$(cd "$BIN/.." && pwd)/share/comet/ref_panel
    REF=$REF_DIR/$TISSUE.atlas.bed.gz
    [ -e "$REF" ] || die "reference panel not found at $REF
       download it (see README) and pass --ref-panel FILE, --ref-panel-dir DIR or set COMET_REF_PANEL_DIR"
    [ -e "$REF.tbi" ] || [ -e "$REF.csi" ] || die "reference panel index missing: $REF.tbi"
  fi
fi
[ -n "$CELL_TYPES" ] || CELL_TYPES=$DEF_TYPES
[ -n "$ASSIGN_THRES" ] || die "--assign-thres is required"
[ -n "$OUT" ] || die "--out is required"
awk -v t="$ASSIGN_THRES" 'BEGIN{exit !(t+0 > 0 && t+0 <= 1)}' || die "--assign-thres must be in (0, 1]"
awk -v t="$ML_CUTOFF" 'BEGIN{exit !(t ~ /^[0-9]+$/ && t >= 128 && t <= 255)}' || die "--ml-cutoff must be an integer in 128..255"
case $OUTPUT_FORMAT in bam|bed-call|both) ;; *) die "--output-format must be bam, bed-call or both" ;; esac
case $STEP in all|prepare|chunks|finalize) ;; *) die "--step must be all, prepare, chunks or finalize" ;; esac

# ---------------------------------------------------------------- several BAMs
BAMS=()
if [[ "$BAM" == *,* ]]; then
  IFS=, read -r -a BAMS <<< "$BAM"
elif [ -f "$BAM" ] && [[ "$BAM" != *.bam ]] && [[ "$BAM" != *.cram ]]; then
  while read -r b; do [[ -z "$b" || "$b" == \#* ]] || BAMS+=("$b"); done < "$BAM"
fi
if [ ${#BAMS[@]} -gt 1 ]; then
  [ "$STEP" == "all" ] || die "--step works on one BAM at a time; run each sample separately for cluster arrays"
  SAMPLES=()
  if [ -n "$SAMPLE" ]; then
    IFS=, read -r -a SAMPLES <<< "$SAMPLE"
    [ ${#SAMPLES[@]} -eq ${#BAMS[@]} ] || die "--sample: ${#SAMPLES[@]} names for ${#BAMS[@]} BAMs"
  else
    for b in "${BAMS[@]}"; do SAMPLES+=("$(basename "$b" .bam)"); done
  fi
  dup=$(printf '%s\n' "${SAMPLES[@]}" | sort | uniq -d | head -1)
  [ -z "$dup" ] || die "duplicate sample name: $dup (use --sample to name the samples)"
  for b in "${BAMS[@]}"; do [ -f "$b" ] || die "BAM not found: $b"; done
  for s in ${METH_PLOT_SAMPLES//,/ }; do
    printf '%s\n' "${SAMPLES[@]}" | grep -qxF "$s" || die "--meth-plot-samples: unknown sample $s"
  done
  mkdir -p "$OUT/summary"
  export COMET_BATCH_CMD=$CMDLINE
  for i in "${!BAMS[@]}"; do
    child=("${PASS[@]}" --bam "${BAMS[$i]}" --sample "${SAMPLES[$i]}" --out "$OUT/${SAMPLES[$i]}" --no-plot)
    done_f=$OUT/${SAMPLES[$i]}/03.run_info/.complete
    if [ -f "$done_f" ] && [ "$(cat "$done_f")" == "comet assign$(printf ' %q' "${child[@]}")" ]; then
      log "=== sample $((i + 1))/${#BAMS[@]}: ${SAMPLES[$i]} already complete (same command), skipped"
      continue
    fi
    log "=== sample $((i + 1))/${#BAMS[@]}: ${SAMPLES[$i]}"
    bash "$BIN/comet_cmd_assign.sh" "${child[@]}" \
      || die "sample ${SAMPLES[$i]} failed (rerun the same command to resume; finished samples are skipped)"
  done
  args=()
  for s in "${SAMPLES[@]}"; do args+=(--sample "$s=$OUT/$s/02.cell_type_quantification"); done
  [ "$NO_PLOT" == "1" ] && args+=(--no-plot)
  python3 "$BIN/comet_plot_assign.py" "${args[@]}" --cell-types "$CELL_TYPES" --tables \
    --meth-samples "$METH_PLOT_SAMPLES" --out-prefix "$OUT/summary/all_samples" || die "summary failed"
  {
    echo "COMET v$VERSION summary of ${#BAMS[@]} samples ($(date '+%F %T'))"
    echo "command: $CMDLINE"
    echo
    echo "Per-sample results: <sample>/ (file guide: <sample>/03.run_info/OUTPUT_README.txt). Samples:"
    printf '  %s\n' "${SAMPLES[@]}"
    echo
    echo "summary/"
    echo "  all_samples.composition.tsv      every sample's composition table, with a leading sample column"
    echo "  all_samples.gw_methylation.tsv   every sample's genome-wide methylation table, same layout"
    if [ "$NO_PLOT" != "1" ]; then
      echo "  all_samples.composition.pdf/png  cell-type composition, one bar per sample: proportions of the"
      echo "                                   reads assigned to a single cell type (frac_single)"
      echo "  all_samples.cell_type_methylation.pdf/png"
      echo "                                   genome-wide CpG methylation of each cell type's reads, one row"
      echo "                                   per sample${METH_PLOT_SAMPLES:+ (samples: $METH_PLOT_SAMPLES)}"
    fi
  } > "$OUT/summary/OUTPUT_README.txt"
  log "Done. ${#BAMS[@]} samples in $OUT; combined tables and plots in $OUT/summary"
  exit 0
fi
[ ${#BAMS[@]} -eq 1 ] && BAM=${BAMS[0]}
[ -f "$BAM" ] || die "BAM not found: $BAM"
[ -f "$REF" ] || die "reference panel not found: $REF"
if [ "$PLATFORM" == "ont" ]; then
  [ -n "$FASTA" ] || die "--fasta is required for --platform ont"
  [ -f "$FASTA" ] || die "FASTA not found: $FASTA"
fi
# normalised, de-duplicated mod list in a fixed order; 5mCG is always called
mods_in=",${MODS// /},"
MODS=5mCG
for m in 5hmCG 5mCH 5hmCH; do
  case ${mods_in,,} in *",${m,,},"*) MODS+=",$m" ;; esac
done
for m in ${mods_in//,/ }; do
  case ${m,,} in 5mcg|5hmcg|5mch|5hmch) ;; *) die "--ont-mods: unknown modification '$m' (use 5mCG, 5hmCG, 5mCH, 5hmCH)" ;; esac
done
[ "$PLATFORM" == "ont" ] || [ "$MODS" == "5mCG" ] || die "--ont-mods beyond 5mCG is only available for --platform ont"

[ -n "$SAMPLE" ] || SAMPLE=$(basename "$BAM" .bam)
if [ "$PLATFORM" == "pacbio" ]; then
  : "${MIN_READ_LEN:=7500}" "${CHUNK_SIZE:=10000}"
else
  : "${MIN_READ_LEN:=500}" "${CHUNK_SIZE:=500000}"
fi

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
# absolute path without resolving symlinks, so an index next to a symlinked BAM is found
BAM=$(cd "$(dirname "$BAM")" && pwd)/$(basename "$BAM")
TMP=$OUT/tmp
RUNINFO=$OUT/03.run_info
LOGS=$RUNINFO/logs
mkdir -p "$TMP/chunks" "$LOGS"

# ---------------------------------------------------------------- dependencies
check_deps() {
  local miss=()
  for t in python3 samtools bedtools tabix bgzip awk sort; do
    command -v "$t" >/dev/null 2>&1 || miss+=("$t")
  done
  if [ "$PLATFORM" == "ont" ]; then
    command -v "$MODKIT" >/dev/null 2>&1 || miss+=("modkit")
  fi
  [ ${#miss[@]} -eq 0 ] || die "missing tools on PATH: ${miss[*]}"
  python3 -c 'import pysam, numpy, pandas, scipy' 2>/dev/null \
    || die "python3 needs pysam, numpy, pandas and scipy"
  [ -e "$BAM.bai" ] || [ -e "$BAM.csi" ] || [ -e "${BAM%.bam}.bai" ] \
    || die "BAM index not found; run: samtools index $BAM"
}

# ---------------------------------------------------------------- atlas
prepare_ref() {
  if [[ "$REF" == *.gz ]] && { [ -e "$REF.tbi" ] || [ -e "$REF.csi" ]; }; then
    REF_GZ=$(readlink -f "$REF")
    return
  fi
  mkdir -p "$RUNINFO/ref_panel"
  REF_GZ=$RUNINFO/ref_panel/$(basename "${REF%.gz}").gz
  if [ ! -e "$REF_GZ.tbi" ]; then
    log "Compressing and indexing reference panel (one-time) -> $REF_GZ"
    if [[ "$REF" == *.gz ]]; then
      { zcat "$REF" || true; } | bgzip -@ "$THREADS" -c > "$REF_GZ.part"
    else
      bgzip -@ "$THREADS" -c "$REF" > "$REF_GZ.part"
    fi
    mv "$REF_GZ.part" "$REF_GZ"
    tabix -p bed "$REF_GZ" || die "tabix failed: the reference panel must be sorted by chrom and start (sort -k1,1 -k2,2n)"
    log "Reuse it in later runs with: --ref-panel $REF_GZ"
  fi
}

# ---------------------------------------------------------------- run info
write_run_parameters() {
  local pn_mode=
  if [ "$PLATFORM" == "pacbio" ]; then
    if [ -n "$USER_PN_TABLE" ]; then pn_mode="user table ($USER_PN_TABLE)"
    elif [ "$NO_PN_CAL" == "1" ]; then pn_mode="off (tau0 for every pass number)"
    else pn_mode="calibrated per pass number"; fi
  fi
  {
    echo "COMET v$VERSION  $(date '+%F %T')"
    echo "command: $CMDLINE"
    [ -n "${COMET_BATCH_CMD:-}" ] && echo "run as part of: $COMET_BATCH_CMD"
    echo
    printf '%-22s %s\n' sample "$SAMPLE" platform "$PLATFORM" tissue "$TISSUE" bam "$BAM" \
      ref_panel "$REF_GZ" cell_types "$CELL_TYPES" assign_thres "$ASSIGN_THRES" ml_cutoff "$ML_CUTOFF" \
      min_cpg "$MIN_CPG" min_read_len "$MIN_READ_LEN" chunk_size "$CHUNK_SIZE" threads "$THREADS" \
      output_format "$OUTPUT_FORMAT" lineage_rescue "${LINEAGE_LIST:-off}"
    if [ -n "$LINEAGE_LIST" ]; then
      printf '%-22s %s\n' lineage_thres "$LINEAGE_THRES" lineage_margin "$LINEAGE_MARGIN"
    fi
    if [ "$PLATFORM" == "pacbio" ]; then
      printf '%-22s %s\n' pass_number_cutoff "$pn_mode"
    else
      printf '%-22s %s\n' fasta "$FASTA" modbam_type "$(cat "$TMP/modbam_type.txt")" ont_mods "$MODS" \
        hmc_cutoff "$HMC_CUTOFF" pileup_threshold_CpG "$PILEUP_THRES" pileup_threshold_CH "$CH_PILEUP_THRES"
    fi
    echo
    echo "tools:"
    echo "  $(python3 --version 2>&1)"
    echo "  $(samtools --version | head -1)"
    echo "  $(bedtools --version)"
    echo "  tabix $(tabix --version 2>&1 | head -1 | awk '{print $NF}')"
    [ "$PLATFORM" == "ont" ] && echo "  $("$MODKIT" --version 2>&1 | head -1)"
    python3 -c 'import pysam, numpy, pandas, scipy; print("  pysam %s, numpy %s, pandas %s, scipy %s" % (pysam.__version__, numpy.__version__, pandas.__version__, scipy.__version__))'
  } > "$RUNINFO/run_parameters.txt"
}

# describe the outputs this run actually produced
write_readme() {
  local A=01.cell_type_assignment Q=02.cell_type_quantification lin=
  [ -n "$LINEAGE_NAMES" ] && lin=1
  {
    cat <<EOF
COMET v$VERSION output guide
sample $SAMPLE | platform $PLATFORM | tissue $TISSUE | cell types $CELL_TYPES
Files are prefixed by the sample name. Parameters and tool versions: 03.run_info/run_parameters.txt

$A/  read-level cell-type assignment
  cell_type_prob/$SAMPLE.cell_type_prob.tsv
      one row per scored read (>= $MIN_CPG informative CpGs in the reference panel); columns:
      read_id, n_cpg (informative CpGs), label (panel column 1..$N_CT, 0 = unassigned),
      cell_type, P_<cell type> (posterior of each cell type), rescued_lineage ('.' if none).
      A read is assigned when its top posterior >= $ASSIGN_THRES and >= 2x the runner-up.
  read_list/cell_type/$SAMPLE.<cell type>.reads.txt
      IDs of the reads assigned to each cell type.
EOF
    [ -n "$lin" ] && cat <<EOF
  read_list/brain_lineage_rescued/$SAMPLE.<lineage>.reads.txt
      unassigned reads rescued at lineage level (${LINEAGE_LIST// /, }; summed posterior
      >= $LINEAGE_THRES and every member > best outside type + $LINEAGE_MARGIN).
      Neuron = Exc_neuron + Inh_neuron; Macroglia = Astrocyte + Oligodendrocyte + OPC.
EOF
    [ "$OUTPUT_FORMAT" != "bed-call" ] && cat <<EOF
  BAM/cell_type/$SAMPLE.<cell type>.bam (.csi)$([ -n "$lin" ] && echo "
  BAM/brain_lineage_rescued/$SAMPLE.<lineage>.bam (.csi)")
      the reads of each set taken from the input BAM with all tags (MM/ML) kept, indexed;
      usable like a sorted-cell library by any methylation tool.
EOF
    if [ "$OUTPUT_FORMAT" != "bam" ]; then
      echo "  BED/bulk/        all reads (bulk reference)"
      echo "  BED/cell_type/   each cell type"
      [ -n "$lin" ] && echo "  BED/brain_lineage_rescued/   each rescued lineage"
      echo "      per-site calls, bgzipped + tabix-indexed, columns: chrom, start, end, n_mod, n_total"
      echo "      (0-based start; methylation level = n_mod / n_total)"
      if [ "$PLATFORM" == "pacbio" ]; then
        echo "      *.CG.bed.gz    5mC at CpGs from COMET's read-level calls; a call counts as modified at"
        echo "                     ML >= tau and unmodified at ML <= 255 - tau (tau = pass-number cutoff);"
        echo "                     both strands are merged onto the C of the CpG"
      else
        local mc
        for mc in ${MODS//,/ }; do
          case $mc in
            5mCG)  echo "      *.mCG.bed.gz   5mC at CpGs (modkit pileup, filter $PILEUP_THRES; strands merged onto the C)" ;;
            5hmCG) echo "      *.hmCG.bed.gz  5hmC at CpGs (modkit pileup, filter $PILEUP_THRES; strands merged onto the C)" ;;
            5mCH)  echo "      *.mCH.bed.gz   5mC at CH cytosines (modkit pileup, filter $CH_PILEUP_THRES; per strand)" ;;
            5hmCH) echo "      *.hmCH.bed.gz  5hmC at CH cytosines (modkit pileup, filter $CH_PILEUP_THRES; per strand)" ;;
          esac
        done
      fi
    fi
    cat <<EOF

$Q/  cell-type quantification
  $SAMPLE.composition.tsv
      group, n_reads and three proportions:
        frac_single          among reads assigned to a single cell type
        frac_single_lineage  among single-type + lineage-rescued reads$([ -z "$lin" ] && echo " (= frac_single, no lineage rescue)")
        frac_all             among all scored reads, including unassigned
  $SAMPLE.gw_methylation.tsv
      genome-wide CpG methylation of each read set from all CpG calls of its reads
      (in or outside the reference panel): group, n_reads, n_meth_calls, n_unmeth_calls, methylation.
      Bulk = every read with CpG calls; Unassigned = scored reads not assigned or rescued.
EOF
    [ "$NO_PLOT" != "1" ] && cat <<EOF
  $SAMPLE.composition.pdf/png
      cell-type composition (proportions of the reads assigned to a single cell type)
  $SAMPLE.cell_type_methylation.pdf/png
      genome-wide CpG methylation of each cell type's reads
EOF
    cat <<EOF

03.run_info/
  OUTPUT_README.txt     this guide
  run_parameters.txt    command line, parameters, reference panel, tool versions
EOF
    if [ "$PLATFORM" == "pacbio" ] && [ "$NO_PN_CAL" != "1" ]; then cat <<EOF
  $SAMPLE.pass_number_table.tsv
      ML cutoff per CCS pass number: np, tau0, n_meth(tau0), n_unmeth(tau0), tau_np,
      n_meth(tau_np), n_unmeth(tau_np); tau_np = 256 means calls from that pass number are dropped
EOF
    fi
    echo "  logs/                 per-step logs"
    echo
    echo "Files in this run:"
    (cd "$OUT" && find "$A" "$Q" 03.run_info -type f ! -path '*/logs/*' ! -name '*.csi' ! -name '*.tbi' | sort | sed 's/^/  /')
  } > "$RUNINFO/OUTPUT_README.txt"
}

# ---------------------------------------------------------------- step 1: prepare
prepare() {
  rm -f "$RUNINFO/.complete"
  prepare_ref

  N_CT=$({ zcat "$REF_GZ" || true; } | head -1 | awk '{print NF-3}')
  [ "$N_CT" -ge 2 ] || die "reference panel must have >= 2 cell-type columns after chrom/start/end"
  if [ -z "$CELL_TYPES" ]; then
    CELL_TYPES=$(seq 1 "$N_CT" | awk '{printf "%sCellType_%d", (NR>1?",":""), $1}')
  fi
  [ "$(echo "$CELL_TYPES" | awk -F, '{print NF}')" -eq "$N_CT" ] \
    || die "--cell-types has $(echo "$CELL_TYPES" | awk -F, '{print NF}') names but the reference panel has $N_CT cell types"
  log "Reference panel ($TISSUE): $N_CT cell types ($CELL_TYPES)"
  [ -z "$LINEAGE_LIST" ] || log "Lineage rescue: $LINEAGE_LIST"

  rm -rf "$TMP"/chunks/*
  PN_TABLE=
  if [ "$PLATFORM" == "pacbio" ]; then
    log "Scanning BAM (length >= $MIN_READ_LEN, primary, autosomal)"
    python3 "$BIN/comet_pacbio_scan.py" "$BAM" "$TMP/scan" --min-len "$MIN_READ_LEN" --threads "$THREADS"
    PN_TABLE=$TMP/pass_number_table.tsv
    if [ -n "$USER_PN_TABLE" ]; then
      cp "$USER_PN_TABLE" "$PN_TABLE"
    elif [ "$NO_PN_CAL" == "1" ]; then
      python3 "$BIN/comet_calibrate.py" "$TMP/scan.ml_hist.tsv" "$PN_TABLE" --cutoff "$ML_CUTOFF" --no-adjust
    else
      log "Calibrating ML cutoff per pass number"
      python3 "$BIN/comet_calibrate.py" "$TMP/scan.ml_hist.tsv" "$PN_TABLE" --cutoff "$ML_CUTOFF"
    fi
    split -d -a 5 -l "$CHUNK_SIZE" "$TMP/scan.reads.tsv" "$TMP/chunks/chunk_" --additional-suffix=.reads
    [ "$KEEP_TMP" == "1" ] || rm -f "$TMP/scan.reads.tsv"
  else
    log "Checking modifications in the BAM (requested output: $MODS)"
    python3 "$BIN/comet_ont_modcheck.py" "$BAM" --mods "$MODS" > "$TMP/modbam_type.txt" \
      || die "the BAM does not support the requested --ont-mods (see messages above)"
    local mk_extract=(extract) mk_help
    mk_help=$("$MODKIT" extract --help 2>&1 || true)
    grep -qE '^\s+full\s' <<< "$mk_help" && mk_extract=(extract full)
    log "modkit ${mk_extract[*]} (CpG, $THREADS threads)"
    "$MODKIT" "${mk_extract[@]}" "$BAM" stdout --cpg --reference "$FASTA" --threads "$THREADS" \
      2> "$LOGS/modkit_extract.log" | bgzip -@ "$THREADS" -c > "$TMP/modkit_extract.tsv.gz"
    log "Collecting CpG calls (5hmC >= $HMC_CUTOFF counted as modified), chunking"
    # columns located by header name, so modkit 0.3.x and newer layouts both work
    # shellcheck disable=SC2016
    local hdr='FNR==1 && $1=="read_id" {for (i=1; i<=NF; i++) c[$i]=i; next}
               {code = (("mod_code" in c) ? $c["mod_code"] : $c["call_code"])
                p    = (("mod_qual" in c) ? $c["mod_qual"] : $c["call_prob"])
                chr=$c["chrom"]; pos=$c["ref_position"]; rid=$c["read_id"]
                if ($c["ref_strand"]=="-") pos=pos-1
                ml=int(p*255+0.5)}'
    zcat "$TMP/modkit_extract.tsv.gz" | awk -v minlen="$MIN_READ_LEN" -v hmc="$HMC_CUTOFF" 'BEGIN{FS=OFS="\t"}'"$hdr"'
      code=="h" && $c["read_length"]>=minlen && ml>=hmc && pos>=0 {print chr, pos, rid, ml}' > "$TMP/hmc.tsv"
    zcat "$TMP/modkit_extract.tsv.gz" | awk -v minlen="$MIN_READ_LEN" -v chunk="$CHUNK_SIZE" -v dir="$TMP/chunks" \
      -v hfile="$TMP/hmc.tsv" '
      BEGIN{FS=OFS="\t"; while ((getline l < hfile) > 0) {split(l, a, "\t"); h[a[1] FS a[2] FS a[3]]=a[4]}}'"$hdr"'
      code=="m" && $c["read_length"]>=minlen && pos>=0 && chr ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {
        k=chr FS pos FS rid
        if (k in h) ml=h[k]
        if (rid != last && nl >= chunk) {close(f); nl=0}
        if (nl == 0) {f=sprintf("%s/chunk_%05d.calls", dir, nc++)}
        print chr, pos, pos+1, ml, (("flag" in c) ? $c["flag"] : 0), $c["ref_strand"], rid > f
        nl++; last=rid}'
    [ "$KEEP_TMP" == "1" ] || rm -f "$TMP/modkit_extract.tsv.gz" "$TMP/hmc.tsv"
  fi

  ls "$TMP/chunks" | sed -E 's/\.(reads|calls)$//' | sort -u > "$TMP/chunks.list"
  [ -s "$TMP/chunks.list" ] || die "no reads/calls passed the filters"

  cat > "$TMP/run.config" <<EOF
PLATFORM=$PLATFORM
TISSUE=$TISSUE
BAM='$BAM'
REF='$REF_GZ'
N_CT=$N_CT
CELL_TYPES='$CELL_TYPES'
ML_CUTOFF=$ML_CUTOFF
ASSIGN_THRES=$ASSIGN_THRES
MIN_CPG=$MIN_CPG
PN_TABLE='$PN_TABLE'
LINEAGE_LIST='$LINEAGE_LIST'
LINEAGE_THRES=$LINEAGE_THRES
LINEAGE_MARGIN=$LINEAGE_MARGIN
OUTPUT_FORMAT=$OUTPUT_FORMAT
MODS=$MODS
NO_PN_CAL=$NO_PN_CAL
KEEP_TMP=$KEEP_TMP
BIN='$BIN'
EOF
  write_run_parameters
  touch "$TMP/prepare.done"
  log "Prepared $(wc -l < "$TMP/chunks.list") chunks"
}

# ---------------------------------------------------------------- step 2: chunks
run_chunks() {
  [ -e "$TMP/prepare.done" ] || die "run --step prepare first"
  if [ -n "$CHUNK_ID" ]; then
    # cluster array jobs may pass a 1-based task index instead of a chunk name
    if [[ "$CHUNK_ID" =~ ^[0-9]+$ ]]; then CHUNK_ID=$(sed -n "${CHUNK_ID}p" "$TMP/chunks.list"); fi
    grep -qx "$CHUNK_ID" "$TMP/chunks.list" || die "unknown chunk: $CHUNK_ID"
    bash "$BIN/comet_worker.sh" "$TMP" "$CHUNK_ID" > "$LOGS/$CHUNK_ID.log" 2>&1 \
      || die "chunk $CHUNK_ID failed, see $LOGS/$CHUNK_ID.log"
    return
  fi
  local n
  n=$(wc -l < "$TMP/chunks.list")
  log "Assigning reads: $n chunks on $THREADS workers"
  export TMP LOGS BIN
  # shellcheck disable=SC2016
  if ! xargs -P "$THREADS" -I{} bash -c \
      'bash "$BIN/comet_worker.sh" "$TMP" {} > "$LOGS/{}.log" 2>&1 || { echo "chunk {} failed, see $LOGS/{}.log" >&2; exit 1; }' \
      < "$TMP/chunks.list"; then
    die "some chunks failed; fix the cause and re-run the same command (finished chunks are skipped)"
  fi
}

# ---------------------------------------------------------------- step 3: finalize
# folder of a read set for one output kind (read_list, BAM, BED):
# 01.cell_type_assignment/<kind>/{cell_type, brain_lineage_rescued, bulk}
group_dir() {
  local d=$OUT/01.cell_type_assignment/$1
  case " $LINEAGE_NAMES " in
    *" $2 "*) d=$d/brain_lineage_rescued ;;
    *) if [ "$2" == "Bulk" ]; then d=$d/bulk; else d=$d/cell_type; fi ;;
  esac
  mkdir -p "$d"
  echo "$d"
}

# merge modkit CG pileup strands onto the + strand C: chrom, start, end, N_mod, N_valid_cov
# (out_m for 5mC, out_h for 5hmC); same rule as the paper's merge_cpg_strands.sh
merge_cpg_strands() {
  local src=$1 out_m=$2 out_h=$3
  awk -v om="$out_m" -v oh="$out_h" 'BEGIN{OFS="\t"}
    function flush(md) {if (pc[md] != "") print pc[md], pp[md], pp[md]+1, pn[md], pv[md] > (md=="m" ? om : oh)}
    $4=="m" || $4=="h" {
      md=$4; canon = ($6=="+") ? $2 : $2-1; nmod=int($12); ncov=int($10)
      if (pc[md]==$1 && pp[md]==canon) {print $1, canon, canon+1, pn[md]+nmod, pv[md]+ncov > (md=="m" ? om : oh); pc[md]=""}
      else {flush(md); pc[md]=$1; pp[md]=canon; pn[md]=nmod; pv[md]=ncov}}
    END {flush("m"); flush("h")}' "$src"
}

finalize() {
  [ -e "$TMP/prepare.done" ] || die "run --step prepare first"
  # shellcheck disable=SC1091
  source "$TMP/run.config"
  local missing
  missing=$(while read -r id; do [ -e "$TMP/chunks/$id.done" ] || echo "$id"; done < "$TMP/chunks.list" | wc -l)
  [ "$missing" -eq 0 ] || die "$missing chunks not finished; run --step chunks"

  LINEAGE_NAMES=$(for l in $LINEAGE_LIST; do printf '%s ' "${l%%=*}"; done)
  log "Merging chunks"
  while read -r id; do cat "$TMP/chunks/$id.prob"; done < "$TMP/chunks.list" > "$TMP/all.prob"
  while read -r id; do cat "$TMP/chunks/$id.gw"; done < "$TMP/chunks.list" > "$TMP/all.gw"
  rm -rf "$OUT/01.cell_type_assignment" "$OUT/02.cell_type_quantification"
  # shellcheck disable=SC2086
  python3 "$BIN/comet_summarize.py" "$TMP/all.prob" "$TMP/all.gw" --out-dir "$OUT" --sample "$SAMPLE" \
    --names "$CELL_TYPES" --lineages "$(echo $LINEAGE_NAMES | tr ' ' ',')"
  rm -f "$RUNINFO/$SAMPLE.pass_number_table.tsv"
  if [ "$PLATFORM" == "pacbio" ] && [ "$NO_PN_CAL" != "1" ]; then
    cp "$PN_TABLE" "$RUNINFO/$SAMPLE.pass_number_table.tsv"
  fi

  # read sets with at least one read: group, reads file
  local g f
  for g in ${CELL_TYPES//,/ } $LINEAGE_NAMES; do
    f=$(group_dir read_list "$g")/$SAMPLE.$g.reads.txt
    if [ -s "$f" ]; then printf '%s\t%s\n' "$g" "$f"; fi
  done > "$TMP/groups.tsv"

  # BAMs: kept for bam/both; ont bed-call also needs them (in tmp) for modkit pileup
  local bam_dir=
  if [ "$OUTPUT_FORMAT" != "bed-call" ]; then
    bam_dir=keep
  elif [ "$PLATFORM" == "ont" ]; then
    bam_dir=$TMP/bam; mkdir -p "$bam_dir"
  fi
  if [ -n "$bam_dir" ]; then
    log "Writing per-cell-type BAMs"
    while IFS=$'\t' read -r g f; do
      if [ "$bam_dir" == "keep" ]; then printf '%s\t%s\n' "$f" "$(group_dir BAM "$g")/$SAMPLE.$g.bam"
      else printf '%s\t%s\n' "$f" "$bam_dir/$SAMPLE.$g.bam"; fi
    done < "$TMP/groups.tsv" > "$TMP/bam_jobs.tsv"
    export BAM LOGS
    # shellcheck disable=SC2016
    xargs -P "$THREADS" -L1 bash -c \
      'samtools view -b --no-PG -F 0x904 -N "$0" -o "$1" --write-index "$BAM" 2>> "$LOGS/bam_split.log"' \
      < "$TMP/bam_jobs.tsv" || die "BAM splitting failed, see $LOGS/bam_split.log"
  fi

  if [ "$OUTPUT_FORMAT" != "bam" ]; then
    if [ "$PLATFORM" == "pacbio" ]; then
      log "Site-level CpG calls (read-level calls -> per-CpG counts)"
      { echo Bulk; cut -f1 "$TMP/groups.tsv"; } | while read -r g; do
        printf '%s\t%s\n' "$g" "$(group_dir BED "$g")/$SAMPLE.$g.CG.bed.gz"
      done > "$TMP/site_jobs.tsv"
      export TMP
      # merge the pre-sorted chunk partials and sum counts of the same CpG
      # shellcheck disable=SC2016
      xargs -P "$THREADS" -L1 bash -c '
        set -euo pipefail; export LC_ALL=C
        find "$TMP/chunks" -path "*.sites/$0.tsv" | sort > "$1.list"
        xargs -a "$1.list" sort -m -k1,1V -k2,2n |
          awk "BEGIN{OFS=\"\t\"} {k=\$1 OFS \$2 OFS \$3; if (k!=p) {if (p!=\"\") print p, m, t; p=k; m=0; t=0} m+=\$4; t+=\$5}
               END {if (p!=\"\") print p, m, t}" | bgzip -c > "$1"
        tabix -p bed "$1"; rm "$1.list"' \
        < "$TMP/site_jobs.tsv" || die "site-level merging failed"
    else
      log "Site-level calls ($MODS; modkit pileup, CpG threshold $PILEUP_THRES, CH $CH_PILEUP_THRES)"
      local jobs=("Bulk	$BAM") b d x
      while IFS=$'\t' read -r g f; do
        if [ "$bam_dir" == "keep" ]; then jobs+=("$g	$(group_dir BAM "$g")/$SAMPLE.$g.bam")
        else jobs+=("$g	$bam_dir/$SAMPLE.$g.bam"); fi
      done < "$TMP/groups.tsv"
      for x in "${jobs[@]}"; do
        g=${x%%	*}; b=${x#*	}; d=$(group_dir BED "$g")
        # CpG: both strands merged onto the C of the CpG (5mCG always; 5hmCG on request)
        "$MODKIT" pileup --threads "$THREADS" -r "$FASTA" --filter-threshold "$PILEUP_THRES" \
          --motif CG 0 "$b" "$TMP/$g.CG.bed" 2>> "$LOGS/modkit_pileup.log" \
          || die "modkit pileup failed for $g, see $LOGS/modkit_pileup.log"
        merge_cpg_strands "$TMP/$g.CG.bed" "$TMP/$g.5mCG.bed" "$TMP/$g.5hmCG.bed"
        # CH: strand-specific, one row per cytosine
        if [[ ",$MODS," == *CH,* ]]; then
          "$MODKIT" pileup --threads "$THREADS" -r "$FASTA" --filter-threshold "$CH_PILEUP_THRES" \
            --motif CH 0 "$b" "$TMP/$g.CH.bed" 2>> "$LOGS/modkit_pileup.log" \
            || die "modkit pileup (CH) failed for $g, see $LOGS/modkit_pileup.log"
          awk -v om="$TMP/$g.5mCH.bed" -v oh="$TMP/$g.5hmCH.bed" 'BEGIN{OFS="\t"}
            $4=="m" {print $1, $2, $3, int($12), int($10) > om}
            $4=="h" {print $1, $2, $3, int($12), int($10) > oh}' "$TMP/$g.CH.bed"
        fi
        for mc in ${MODS//,/ }; do
          touch "$TMP/$g.$mc.bed"
          bgzip -c "$TMP/$g.$mc.bed" > "$d/$SAMPLE.$g.${mc#5}.bed.gz"
          tabix -p bed "$d/$SAMPLE.$g.${mc#5}.bed.gz"
        done
        rm -f "$TMP/$g".*.bed
      done
    fi
  fi

  if [ "$NO_PLOT" != "1" ]; then
    log "Plots"
    python3 "$BIN/comet_plot_assign.py" --sample "$SAMPLE=$OUT/02.cell_type_quantification" \
      --cell-types "$CELL_TYPES" --out-prefix "$OUT/02.cell_type_quantification/$SAMPLE" 2>> "$LOGS/plot.log" \
      || log "WARNING: plotting failed, see $LOGS/plot.log (results are unaffected)"
  fi
  write_readme
  printf '%s\n' "$CMDLINE" > "$RUNINFO/.complete"   # lets a multi-BAM rerun skip this sample
  [ "$KEEP_TMP" == "1" ] || rm -rf "$TMP"
  log "Done. Results in $OUT (file guide: 03.run_info/OUTPUT_README.txt)"
}

check_deps
case $STEP in
  all)
    if [ -e "$TMP/prepare.done" ]; then
      log "Resuming from $TMP (delete it to start over)"
    else
      prepare
    fi
    run_chunks; finalize ;;
  prepare) prepare ;;
  chunks) run_chunks ;;
  finalize) finalize ;;
esac
