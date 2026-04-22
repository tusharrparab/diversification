# Diversification

This repository contains the Result Block 1 ancestral climatic reconstruction
workflow for an avian macroevolution analysis.

Result Block 1 estimates the ancestral climatic core of birds as the baseline
for later frontier expansion and packing analyses. It is not a standalone PCA
story.

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

The PCA projection is a visualization and diagnostic layer for the
original-variable root estimate. It is not the primary biological inference.

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

The cloud workflow is:

```text
.github/workflows/run_ancestral_climate.yml
```

It is manual-only because its only trigger is `workflow_dispatch`. Pushing code
does not run the models. GitHub documents `workflow_dispatch` as the manual
workflow trigger, with manual starts supported from the Actions tab, GitHub CLI,
or REST API; GitHub also notes that the workflow file must exist on the default
branch for manual dispatch to be available.

Operationally, a committed workflow file and a completed workflow run are
different things. A commit can make the run button available, but only a real
manual dispatch creates a GitHub Actions run ID, logs, verified outputs, and an
uploaded artifact.

- [GitHub: manually running a workflow](https://docs.github.com/actions/how-tos/manage-workflow-runs/manually-run-a-workflow)
- [GitHub: `workflow_dispatch` event](https://docs.github.com/actions/reference/workflows-and-actions/events-that-trigger-workflows#workflow_dispatch)

To start the production pipeline from the GitHub UI:

1. Open the repository on GitHub.
2. Go to **Actions**.
3. Select **Run ancestral climate pipeline**.
4. Click **Run workflow**.
5. Choose the stage and core count.

You can also dispatch it with GitHub CLI:

```sh
gh workflow run run_ancestral_climate.yml \
  -f stage=model_fit \
  -f ncores=4 \
  -f use_fit_cache=true \
  -f upload_results=true
```

Or use the GitHub REST API workflow-dispatch endpoint for
`run_ancestral_climate.yml`.

Run the cloud stages manually in this order:

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

Before installing R packages or starting model runtime, the workflow checks that
these required paths exist:

- `scripts/run_ancestral_climate_reconstruction.R`
- `data/raw/pruned_phylogeny.nex`
- `data/raw/climate.csv`

The workflow log also prints the selected stage, `GEIGER_NCORES`, fit-cache
setting, exact input paths, result directory, cache directory, expected key
output for the selected stage, all files verified for that stage, and whether a
previous fit cache was restored or the run is starting cold.

Each run writes:

```text
results/result_block_1_strict_ancestral_climate/cloud_run_trace.txt
```

That trace records the workflow name, run ID, run attempt, stage, commit SHA,
UTC timestamp, whether fit cache was requested, whether a cache restore was
detected, the expected key output, and the artifact name.

The workflow verifies these real script outputs by stage:

| Stage | Files required by the workflow |
| --- | --- |
| `model_fit` | `per_variable_model_fit_attempts_detailed.csv`; `per_variable_model_selection_preliminary.csv` |
| `root_reconstruction` | `per_variable_model_fits.csv`; `ancestral_root_original_climate_strict.csv`; `sensitivity_comparison_summary.csv`; `cross_variable_climate_coherence_diagnostics.csv`; `ancestral_root_vector_used_for_PCA_projection_strict.csv` |
| `pca_projection` | `ancestral_root_projected_PCA_point_strict.csv`; `ancestral_root_projected_PCA_ellipses_strict.csv`; `ancestral_root_interpretability_metrics.csv`; `biological_interpretation_summary.txt`; `result_block_1_strict_ancestral_climate_outputs.xlsx` |
| `plot` | `ancestral_root_projected_PCA_strict_plot.png`; `ancestral_root_projected_PCA_strict_plot.pdf` |

To confirm that `fitContinuous()` actually ran on real data, inspect a real
workflow run log for the model-fit loop:

- `Fitting or loading per-variable geiger::fitContinuous model fits.`
- `fitting all candidate models for <variable>`
- `<variable> under <model>`

If the workflow restored previous fits, the log instead reports
`using cached fitContinuous fits for <variable>`. That confirms the stage reused
the fit cache rather than refitting that variable from scratch.

For a `model_fit`-only run, the stage completion proof is the uploaded artifact
for that GitHub Actions run containing non-empty:

- `per_variable_model_fit_attempts_detailed.csv`
- `per_variable_model_selection_preliminary.csv`

The verifier fails the run if either file is missing or empty. The first file is
the detailed per-model attempt table; the second is written only when the
`model_fit` stage reaches its intended stopping point.

`per_variable_model_fits.csv` is produced by the next stage,
`root_reconstruction`, after the cached/raw model fits are assembled with the
root-reconstruction diagnostics. Treat it as the first compact biological
model-fit summary, not as a file emitted by a `model_fit`-only run.

Each cloud run uploads a stage-specific artifact, for example:

- `ancestral-climate-model-fit`
- `ancestral-climate-root-reconstruction`
- `ancestral-climate-pca-projection`
- `ancestral-climate-plot`

After cloud runs, inspect the biological outputs first:

- `per_variable_model_fits.csv` if the run reached `root_reconstruction` or a
  later stage
- `ancestral_root_original_climate_strict.csv`
- `biological_interpretation_summary.txt`

Then inspect the PCA outputs as visualization and diagnostic products. They are
not the primary inferential basis and should not be framed as an independent PCA
result.

## Key Outputs

Outputs are written to:

```text
results/result_block_1_strict_ancestral_climate/
```

Primary audit files:

- `cloud_run_trace.txt`
- `per_variable_model_fits.csv`
- `per_variable_model_fit_attempts_detailed.csv`
- `per_variable_model_selection_preliminary.csv`
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
