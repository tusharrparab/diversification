# Diversification

This repository contains a reproducible ancestral climatic reconstruction
workflow for birds.

Result Block 1 asks: where did the radiation begin in climatic space?

The pipeline is intentionally strict: ancestral climate is reconstructed in the
original environmental variables, and PCA is used only afterward as a
visualization space.

## Inputs

- `data/raw/pruned_phylogeny.nex`
- `data/raw/climate.csv`

## Run

From the repository root:

```sh
Rscript scripts/run_ancestral_climate_reconstruction.R
```

`geiger::fitContinuous()` can be slow on the full tree. To allow geiger to use
more than one worker where supported:

```sh
GEIGER_NCORES=4 Rscript scripts/run_ancestral_climate_reconstruction.R
```

For smoke tests or alternate inputs, the paths can be overridden without
editing the script:

```sh
TREE_FILE=/path/to/pruned_phylogeny.nex \
CLIMATE_FILE=/path/to/climate.csv \
OUTPUT_DIR=/tmp/strict_ancestral_climate_test \
Rscript scripts/run_ancestral_climate_reconstruction.R
```

For each climate variable, the script fits BM, OU, EB, lambda, delta, and
mean_trend with `geiger::fitContinuous()`. Primary ancestral reconstruction is
restricted to BM, lambda, and delta. OU, EB, and mean_trend are sensitivity-only
models: they can be the best AICc model, but they are not silently forced into
the primary root ellipse.

The output note in `SUMMARY.md` explicitly lists the remaining approximations,
including diagonal uncertainty propagation from original variables into PCA
space.

## Key Outputs

Outputs are written to `results/result_block_1_strict_ancestral_climate/`.

- `strict_model_fits_by_variable.csv`
- `strict_ancestral_reconstruction_by_variable.csv`
- `strict_root_vector_for_PCA_projection.csv`
- `strict_uncertainty_assumptions.csv`
- `visualization_PCA_species_scores.csv`
- `visualization_PCA_loadings.csv`
- `visualization_PCA_eigenvalues.csv`
- `projected_primary_root_PCA_point.csv`
- `projected_primary_root_PCA_ellipses.csv`
- `strict_ancestral_climatic_core.png`
- `strict_ancestral_climatic_core.pdf`
- `strict_ancestral_climate_reconstruction_outputs.xlsx`
