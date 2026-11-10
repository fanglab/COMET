# Changelog

All notable changes to COMET (package `comet-meth`) are recorded here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [1.1.0] - 2026-11-10

COMET 1.1.0 with the brain reference panel.

- `comet assign` uses the brain reference panel (6 cell types).
- Example data and `comet test`: three postmortem brain PacBio HiFi samples (`brain_01`,
  `brain_02`, `brain_03`, 0.5X) and the brain reference panel, downloaded from Zenodo
  (DOI 10.5281/zenodo.24064938).
- `comet dmr` example and `comet test` step: an AD-specific microglia DMR on per-site CpG BEDs of
  six deeply sequenced postmortem brain samples (AD, MCI, NCI; chr15:73.5-75.5 Mb), bundled in
  `example/data/dmr/`. `plot_dmr_region.py` pools several BEDs given in one track.
- Installation with Bioconda (`comet-meth`), the Docker image `ghcr.io/fanglab/comet-meth` or from
  source.
