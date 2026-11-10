# COMET example data

Three small postmortem brain PacBio HiFi data sets that run `comet assign` in ~5 min with 16 threads,
and per-site CpG BEDs of six deeply sequenced postmortem brain samples for `comet dmr` (seconds).

```bash
bash example/download_example_data.sh --panel-dir ref_panel   # example BAMs + brain panel from Zenodo
bash example/run_example.sh ref_panel/brain.atlas.bed.gz 16     # reference panel, threads
```

| Path | Content | Walk-through |
|---|---|---|
| `data/assign/` | three postmortem brain BAMs (`brain_01`, `brain_02`, `brain_03`), 0.5X (downloaded from Zenodo, [10.5281/zenodo.24064938](https://doi.org/10.5281/zenodo.24064938)) | [comet assign example](../docs/assign.md#example-three-postmortem-brain-samples) |
| `data/dmr/` | per-site CpG BEDs of the six brain cell types of six deep postmortem brain samples (AD, MCI, NCI; two each), chr15:73.5-75.5 Mb (in the repository) | [comet dmr example](../docs/dmr.md#example-an-ad-specific-microglia-dmr) |
| `output/` | expected results (tables, plots, run records) | |
| `plot_dmr_region.py` | methylation track plot of a region from per-site BEDs (example only) | |

Paths in the shipped run records are written relative to the repository (`example/...`,
`ref_panel/...`).
