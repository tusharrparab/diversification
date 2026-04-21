# Diversification

This repository contains a reproducible ancestral climatic core analysis for the
pruned bird phylogeny and climate table.

## Inputs

- `data/raw/pruned_phylogeny.nex`
- `data/raw/climate.csv`

## Run

From the repository root:

```sh
Rscript scripts/run_ancestral_climate_reconstruction.R
```

The script reconstructs the crown/root climate on the original climate
variables, projects that crown point into PCA climate space, and ranks extant
species by their Mahalanobis distance to the projected crown climate. Because
the exact `fastAnc(..., vars = TRUE, CI = TRUE)` reconstruction is slow on this
large tree, the script reuses
`results/ancestral_climate_reconstruction/crown_ancestral_original_climate.csv`
when it is already present. To force a fresh full reconstruction:

```sh
FORCE_RECONSTRUCT_CROWN=true Rscript scripts/run_ancestral_climate_reconstruction.R
```

## Key Outputs

- `results/ancestral_climate_reconstruction/crown_ancestral_original_climate.csv`
- `results/ancestral_climate_reconstruction/crown_ancestral_projected_PCA_point.csv`
- `results/ancestral_climate_reconstruction/ancestral_core_species_membership.csv`
- `results/ancestral_climate_reconstruction/ancestral_core_summary.csv`
- `results/ancestral_climate_reconstruction/ancestral_climatic_core.png`
- `results/ancestral_climate_reconstruction/ancestral_climate_reconstruction_outputs.xlsx`

Current run summary:

- Species used: 9,621
- Crown node: 9,622
- Thermal axis: PC1
- Precipitation axis: PC4
- Species inside 50% ancestral core: 5,365
- Species inside 80% ancestral core: 7,990
- Species inside 95% ancestral core: 8,843
