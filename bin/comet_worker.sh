#!/usr/bin/env bash
# COMET per-chunk worker: [extract] -> match atlas -> assign -> per-read genome-wide counts
# -> (pacbio, bed-call output) per-site partial counts for each read set.
# Usage: comet_worker.sh <tmp_dir> <chunk_id>
# Safe to re-run: a finished chunk (chunk_id.done) is skipped.
set -euo pipefail
export LC_ALL=C

tmp=$1
id=$2
source "$tmp/run.config"
c=$tmp/chunks/$id
[ -e "$c.done" ] && exit 0

# 1. per-read CpG calls (ONT calls were already written during prepare)
if [ "$PLATFORM" == "pacbio" ]; then
  python3 "$BIN/comet_pacbio_extract.py" "$BAM" "$c.reads" "$PN_TABLE" "$c.calls.part"
  mv "$c.calls.part" "$c.calls"
fi

# 2. match calls to the atlas, fetching only the chunk's genomic span from the indexed atlas
awk 'BEGIN{OFS="\t"}
     !($1 in mn) {mn[$1]=$2; mx[$1]=$3; o[++n]=$1}
     {if ($2<mn[$1]) mn[$1]=$2; if ($3>mx[$1]) mx[$1]=$3}
     END {for (i=1; i<=n; i++) print o[i], mn[o[i]], mx[o[i]]}' "$c.calls" > "$c.regions.bed"
if [ -s "$c.regions.bed" ]; then
  tabix -R "$c.regions.bed" "$REF" > "$c.atlas.bed"
else
  : > "$c.atlas.bed"
fi
if [ -s "$c.atlas.bed" ]; then
  bedtools intersect -a "$c.calls" -b "$c.atlas.bed" -wo > "$c.matched"
else
  : > "$c.matched"
fi

# 3. read-level assignment (+ lineage rescue for brain)
opt=()
[ "$PLATFORM" == "pacbio" ] && opt+=(--pn-table "$PN_TABLE")
for l in $LINEAGE_LIST; do opt+=(--lineage "$l"); done
python3 "$BIN/comet_assign.py" "$c.matched" "$c.prob.part" --n-ct "$N_CT" --cutoff "$ML_CUTOFF" \
  --assign-thres "$ASSIGN_THRES" --min-cpg "$MIN_CPG" \
  --lineage-thres "$LINEAGE_THRES" --lineage-margin "$LINEAGE_MARGIN" "${opt[@]+"${opt[@]}"}"
mv "$c.prob.part" "$c.prob"

# 4. per-read genome-wide call counts (all CpG calls, in or out of the atlas)
awk -v fixed="$ML_CUTOFF" -v pn="${PN_TABLE:-}" 'BEGIN{FS=OFS="\t"
       if (pn != "") while ((getline l < pn) > 0) {split(l, f, "\t"); t[f[1]]=f[5]}}
     {tau = (pn == "") ? fixed : (($5 in t) ? t[$5] : 256)
      if (!($7 in npc)) {npc[$7]=$5; o[++n]=$7}
      if ($4 >= tau) m[$7]++; else if ($4 <= 255-tau) u[$7]++}
     END {for (i=1; i<=n; i++) {r=o[i]; print r, npc[r], m[r]+0, u[r]+0}}' "$c.calls" > "$c.gw.part"
mv "$c.gw.part" "$c.gw"

# 5. pacbio site-level partial counts per read set: chrom, start, end, n_meth, n_total
#    (Bulk = every read with calls; cell types and rescued lineages from step 3)
if [ "$PLATFORM" == "pacbio" ] && [ "$OUTPUT_FORMAT" != "bam" ]; then
  mkdir -p "$c.sites"
  awk -v pn="$PN_TABLE" -v prob="$c.prob" -v names="$CELL_TYPES" -v k="$N_CT" -v dir="$c.sites" '
    BEGIN {FS=OFS="\t"; split(names, nm, ",")
      while ((getline l < pn) > 0) {split(l, f, "\t"); t[f[1]]=f[5]}
      while ((getline l < prob) > 0) {
        split(l, f, "\t")
        if (f[3] > 0) g[f[1]] = nm[f[3]]
        else if (f[4+k] != ".") g[f[1]] = f[4+k]}}
    {tau = ($5 in t) ? t[$5] : 256
     if ($4 >= tau) s=1; else if ($4 <= 255-tau) s=0; else next
     key=$1 OFS $2 OFS $3
     M["Bulk" SUBSEP key]+=s; T["Bulk" SUBSEP key]++
     if ($7 in g) {M[g[$7] SUBSEP key]+=s; T[g[$7] SUBSEP key]++}}
    END {for (x in T) {split(x, y, SUBSEP); print y[2], M[x], T[x] > (dir "/" y[1] ".tsv")}}' "$c.calls"
  for f in "$c.sites"/*.tsv; do
    [ -e "$f" ] || continue
    sort -k1,1V -k2,2n "$f" > "$f.part" && mv "$f.part" "$f"
  done
fi

if [ "$KEEP_TMP" != "1" ]; then
  rm -f "$c.calls" "$c.reads" "$c.regions.bed" "$c.atlas.bed" "$c.matched"
fi
touch "$c.done"
