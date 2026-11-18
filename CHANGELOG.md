# Changelog

All notable changes to COMET (package `comet-meth`) are recorded here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [1.1.1] - 2026-11-18

- `comet assign`, `comet dmr` and `comet test`: an option given without its value (e.g.
  `comet assign --bam`) now stops with `ERROR: --bam requires a value`, instead of a shell error.
- README: the reference panel is downloaded with `curl` from Zenodo, so a Bioconda installation
  does not need the repository.
- Installation with Bioconda or from source only; the Docker image is no longer built.
- `comet_plot_dmr.py`: methylation tracks of a region (`--region`) or of the N-th significant DMR
  (`--dmr N`) from a `comet dmr` output folder, one track per input file or per group (`--pool`).
  Replaces `example/plot_dmr_region.py`.

## [1.1.0] - 2026-11-18

COMET 1.1.0 with the brain reference panel.

- `comet assign --ref-panel FILE` (required): the cell-type reference panel of the tissue, e.g.
  `brain.atlas.bed.gz` (6 cell types).
- Example data and `comet test`: three postmortem brain PacBio HiFi samples (`brain_01`,
  `brain_02`, `brain_03`, 0.5X) and the brain reference panel, downloaded from Zenodo
  (DOI 10.5281/zenodo.24064938).
- `comet dmr` example and `comet test` step: an AD-specific microglia DMR on per-site CpG BEDs of
  six deeply sequenced postmortem brain samples (AD, MCI, NCI; chr15:73.5-75.5 Mb), bundled in
  `example/data/dmr/`. `plot_dmr_region.py` pools several BEDs given in one track.
- Installation with Bioconda (`comet-meth`), the Docker image `ghcr.io/fanglab/comet-meth` or from
  source.
