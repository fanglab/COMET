# COMET example data

Three small postmortem brain PacBio HiFi data sets that run `comet assign` in ~5 min with 16 threads,
and per-site CpG BEDs of six deeply sequenced postmortem brain samples for `comet dmr` (seconds).
The blood and colon examples will be added upon publication.

```bash
bash example/download_example_data.sh --ref-panel-dir ref_panel   # example BAMs + brain panel from Zenodo
export COMET_REF_PANEL_DIR=$PWD/ref_panel                          # skip if the panel is already set up
bash example/run_example.sh 16                                     # threads
```

| Path | Content | Walk-through |
|---|---|---|
| `data/assign/` | three postmortem brain BAMs (`brain_01`, `brain_02`, `brain_03`), 0.5X (downloaded from Zenodo, [10.5281/zenodo.24064938](https://doi.org/10.5281/zenodo.24064938)) | [comet assign example](../docs/assign.md#example-three-postmortem-brain-samples) |
| `data/dmr/` | per-site CpG BEDs of the six brain cell types of six deep postmortem brain samples (AD, MCI, NCI; two each), chr15:73.5-75.5 Mb (in the repository) | [comet dmr example](../docs/dmr.md#example-an-ad-specific-microglia-dmr) |
| `output/` | expected results (tables, plots, run records) | |
| `plot_dmr_region.py` | methylation track plot of a region from per-site BEDs (example only) | |

Paths in the shipped run records are written relative to the repository (`example/...`,
`$COMET_REF_PANEL_DIR/...`).
