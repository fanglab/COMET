COMET v1.1.0 summary of 3 samples (2026-11-10 13:13:31)
command: comet assign --platform pacbio --ref-panel ref_panel/brain.atlas.bed.gz --assign-thres 0.65 --bam data/assign/brain_01.bam\,data/assign/brain_02.bam\,data/assign/brain_03.bam --threads 16 --out output/assign

Per-sample results: <sample>/ (file guide: <sample>/03.run_info/OUTPUT_README.txt). Samples:
  brain_01
  brain_02
  brain_03

summary/
  all_samples.composition.tsv      every sample's composition table, with a leading sample column
  all_samples.gw_methylation.tsv   every sample's genome-wide methylation table, same layout
  all_samples.composition.pdf/png  cell-type composition, one bar per sample: proportions of the
                                   reads assigned to a single cell type (frac_single)
  all_samples.cell_type_methylation.pdf/png
                                   genome-wide CpG methylation of each cell type's reads, one row
                                   per sample
