#!/usr/bin/env bash
# COMET example: postmortem brain PacBio HiFi data, 0.5X each (autosomes; base
# qualities and unused tags removed).
#   brain_01, brain_02, brain_03 : three postmortem brain samples
#   data/dmr/ : per-site CpG BEDs (chr15:73,500,000-75,500,000) of the six brain cell types of six
#               deep-sequenced postmortem brain samples, two per group (AD, MCI, NCI)
#
# What it shows:
#   - comet assign: cell-type composition and cell-type methylation of the three samples
#   - comet dmr: AD-specific microglia DMRs, including chr15:74,174,956-74,175,254 (16 CpGs,
#     unmethylated in AD microglia and methylated in MCI and NCI microglia);
#     output/dmr/AD_microglia.tracks.pdf shows the microglia methylation of the six samples there
#
# usage: bash download_example_data.sh --ref-panel-dir PANELS   (once: example BAMs + brain panel)
#        export COMET_REF_PANEL_DIR=PANELS
#        bash run_example.sh [THREADS]
set -euo pipefail
T=${1:-8}
cd "$(dirname "$0")"
for s in brain_01 brain_02 brain_03; do
  [ -s data/assign/$s.bam ] || { echo "data/assign/$s.bam missing: run 'bash download_example_data.sh' first" >&2; exit 1; }
done

# 1. cell-type assignment of the three samples (brain panel: Exc_neuron, Inh_neuron, Astrocyte,
#    Microglia, Oligodendrocyte, OPC; --assign-thres 0.65)
comet assign --platform pacbio --assign-thres 0.65 \
      --bam data/assign/brain_01.bam,data/assign/brain_02.bam,data/assign/brain_03.bam \
      --threads $T --out output/assign
#    -> output/assign/<sample>/ per sample; output/assign/summary/ combined tables and plots:
#       all_samples.composition.pdf/png, all_samples.cell_type_methylation.pdf/png

# 2. cell-type composition and cell-type methylation
cut -f1-5 output/assign/summary/all_samples.composition.tsv
awk -F'\t' 'NR==1 || $2!="Unassigned"' output/assign/summary/all_samples.gw_methylation.tsv

# 3. AD-specific microglia DMRs: AD microglia (2 samples) vs MCI and NCI microglia (4 samples),
#    required against each of the four individually (--dmr-level individual)
D=data/dmr/brain_deep
comet dmr --target ${D}_AD_1.Microglia.CG.bed.gz,${D}_AD_2.Microglia.CG.bed.gz --target-name AD_Microglia \
      --background ${D}_MCI_1.Microglia.CG.bed.gz,${D}_MCI_2.Microglia.CG.bed.gz,${D}_NCI_1.Microglia.CG.bed.gz,${D}_NCI_2.Microglia.CG.bed.gz \
      --background-name MCI_NCI_Microglia --dmr-diff 0.2 --dmr-level individual --threads $T --out output/dmr
cut -f1-4,20-25 output/dmr/01.dmr/AD_Microglia_vs_MCI_NCI_Microglia.metilene.sig.tsv

# 4. methylation tracks of the six samples' microglia around the region
#    (example-only plotting script; needs matplotlib and pysam)
python3 plot_dmr_region.py --region chr15:74,174,700-74,175,500 --min-depth 1 \
      $(for s in AD_1 AD_2 MCI_1 MCI_2 NCI_1 NCI_2; do echo "--track $s=${D}_$s.Microglia.CG.bed.gz"; done) \
      --out output/dmr/AD_microglia.tracks
