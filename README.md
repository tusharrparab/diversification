# Diversification

This repository contains the Result Block 1 ancestral climatic reconstruction
workflow for an avian macroevolution analysis.

## Biological Question

Where did the radiation begin in climatic space?

The workflow is intentionally conservative:

> Ancestral climate was reconstructed in the original environmental variables
> after per-variable evolutionary model fitting. PCA was used only to visualize
> extant climatic structure and the projected ancestral root. Variables whose
> best-supported model fell outside the primary reconstructable model set were
> retained as sensitivity cases and explicitly flagged.

## Inputs

- `data/raw/pruned_phylogeny.nex`
- `data/raw/climate.csv`

The script preserves the upstream workflow for reading the pruned phylogeny,
cleaning species names, matching tree tips to climate rows, pruning unmatched
tips, and building the final climate matrix.

## Method Logic

PCA is visualization only. It is not the inferential basis for ancestral
reconstruction.

For each original climate variable, the script fits the following evolutionary
models with `geiger::fitContinuous()`:

- `BM`
- `OU`
- `EB`
- `lambda`
- `delta`
- `mean_trend`

Models are ranked by AICc separately for each climate variable.

Primary ancestral reconstruction is intentionally restricted to models that are
implemented transparently in this workflow:

- `BM`: `ape::ace(..., method = "ML", model = "BM", CI = TRUE)`
- `lambda`: transform tree with `geiger::lambdaTree()`, then reconstruct with
  `ape::ace(..., model = "BM")`
- `delta`: transform tree with `geiger::deltaTree()`, then reconstruct with
  `ape::ace(..., model = "BM")`

`OU`, `EB`, and `mean_trend` are treated as sensitivity-only models. They may be
the best AICc model, but they are not silently forced into the primary root
estimate or the primary root ellipse. This is deliberate: routine ancestral
reconstruction support is not equally straightforward across these models, and
model-specific ancestral reconstructions should be treated cautiously.

## Important Approximation

The projected root ellipse uses a diagonal covariance approximation in the
original climate variables, then projects that uncertainty into PCA space. This
ellipse is an uncertainty projection, not a fully jointly inferred multivariate
ancestral climatic niche.

Interpret the outputs in this order:

1. Per-variable model fits.
2. Ancestral root estimates in original climate variables.
3. Sensitivity and coherence diagnostics.
4. PCA projection of the already reconstructed root.

## Run

From the repository root:

```sh
Rscript scripts/run_ancestral_climate_reconstruction.R
```

The full production run is expensive because it fits many `fitContinuous`
models. Use `GEIGER_NCORES` where supported:

```sh
GEIGER_NCORES=4 Rscript scripts/run_ancestral_climate_reconstruction.R
```

For alternate inputs or smoke tests:

```sh
TREE_FILE=/path/to/pruned_phylogeny.nex \
CLIMATE_FILE=/path/to/climate.csv \
OUTPUT_DIR=/tmp/result_block_1_test \
Rscript scripts/run_ancestral_climate_reconstruction.R
```

## Production Controls

Model fits are cached under:

```text
results/result_block_1_strict_ancestral_climate/cache/
```

Useful controls:

```sh
# Fit/cache models only, then stop
Rscript scripts/run_ancestral_climate_reconstruction.R --stage=model_fit

# Reuse cached fits and rerun root reconstruction only
Rscript scripts/run_ancestral_climate_reconstruction.R --stage=root_reconstruction

# Reuse cached fits and rerun PCA projection outputs without plots
Rscript scripts/run_ancestral_climate_reconstruction.R --stage=pca_projection --skip-plots

# Force fresh model fitting
Rscript scripts/run_ancestral_climate_reconstruction.R --force-refit-models

# Disable model cache resume
Rscript scripts/run_ancestral_climate_reconstruction.R --no-resume-model-cache

# Save intermediate RDS objects
Rscript scripts/run_ancestral_climate_reconstruction.R --save-intermediate-rds

# Skip plot generation
SKIP_PLOTS=true Rscript scripts/run_ancestral_climate_reconstruction.R
```

Equivalent environment variables:

- `PIPELINE_STAGE`
- `GEIGER_NCORES`
- `SKIP_PLOTS`
- `FORCE_REFIT_MODELS`
- `RESUME_MODEL_CACHE`
- `SAVE_INTERMEDIATE_RDS`
- `MIN_VALID_VALUES`

## Cloud Runs With GitHub Actions

The production pipeline can be run manually in GitHub Actions:

1. Open the repository on GitHub.
2. Go to **Actions**.
3. Select **Run ancestral climate pipeline**.
4. Click **Run workflow**.
5. Choose the stage and core count.

Run the cloud stages in this order:

1. `model_fit`
2. `root_reconstruction`
3. `pca_projection`
4. `plot`

The `model_fit` stage is the main computational bottleneck because it runs
`geiger::fitContinuous()` repeatedly across climate variables and candidate
models. The workflow restores and saves the pipeline cache directory:

```text
results/result_block_1_strict_ancestral_climate/cache/
```

Cache keys include the R script, the pruned phylogeny, and the climate table, so
cached model fits are invalidated when relevant inputs change.

Each cloud run uploads a stage-specific artifact, for example:

- `ancestral-climate-model-fit`
- `ancestral-climate-root-reconstruction`
- `ancestral-climate-pca-projection`
- `ancestral-climate-plot`

After cloud runs, inspect the biological outputs first:

- `per_variable_model_fits.csv`
- `ancestral_root_original_climate_strict.csv`
- `biological_interpretation_summary.txt`

Then inspect the PCA outputs as visualization and diagnostic products. They are
not the primary inferential basis.

## Key Outputs

Outputs are written to:

```text
results/result_block_1_strict_ancestral_climate/
```

Primary audit files:

- `per_variable_model_fits.csv`
- `per_variable_model_fit_attempts_detailed.csv`
- `ancestral_root_original_climate_strict.csv`
- `sensitivity_comparison_summary.csv`
- `cross_variable_climate_coherence_diagnostics.csv`
- `ancestral_root_projected_PCA_point_strict.csv`
- `ancestral_root_projected_PCA_ellipses_strict.csv`
- `ancestral_root_interpretability_metrics.csv`
- `biological_interpretation_summary.txt`
- `pipeline_warnings_and_safeguards.csv`

Visualization files:

- `visualization_only_PCA_species_scores.csv`
- `visualization_only_PCA_loadings.csv`
- `visualization_only_PCA_eigenvalues.csv`
- `species_membership_in_projected_root_ellipses_strict.csv`
- `ancestral_root_projected_PCA_strict_plot.png`
- `ancestral_root_projected_PCA_strict_plot.pdf`

Combined workbook:

- `result_block_1_strict_ancestral_climate_outputs.xlsx`

## How To Read The Main Outputs

`per_variable_model_fits.csv` tells you whether each variable is primarily
reconstructable and whether the primary model is much worse than the best
overall model.

`ancestral_root_original_climate_strict.csv` is the biological core of the
analysis. It reports the root estimate in original climate units, extant ranges,
standardized root positions, CI width diagnostics, and interpretability flags.

`sensitivity_comparison_summary.csv` reports cases where `OU`, `EB`, or
`mean_trend` won AICc and compares any available sensitivity root estimate to
the primary estimate.

`ancestral_root_projected_PCA_point_strict.csv` projects the reconstructed root
into PCA space and reports whether it falls inside the extant convex hull and
how far it lies from the extant climatic cloud.

`ancestral_root_projected_PCA_ellipses_strict.csv` contains projected
uncertainty ellipses and classifies their area relative to the extant occupied
PCA space.

`ancestral_root_interpretability_metrics.csv` gives the final conservative
classification: `ROBUST`, `TENTATIVE`, or `WEAK`.

`biological_interpretation_summary.txt` is the manuscript-facing diagnostic
summary. It names robust, tentative, and fragile variables and reports whether
the projected PCA root should be interpreted as stable, tentative, or weak.

## Safeguards

The script fails or warns loudly when:

- tree tips and climate rows do not align after matching;
- a variable has too few valid or unique values;
- model fits return malformed objects;
- too many variables have sensitivity-only best models;
- a primary model is much worse than the best overall model;
- root confidence intervals are extremely wide;
- the projected root is outside or far beyond the extant PCA cloud;
- projected uncertainty covariance is near-singular or degenerate.

The goal is not to force a neat-looking result. Fragile variables remain visible
in the outputs and should be discussed directly.
