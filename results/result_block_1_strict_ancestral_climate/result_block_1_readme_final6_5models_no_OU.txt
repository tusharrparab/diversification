Result Block 1 downstream readme (final6_5models_no_OU)
=======================================================

Climate set: final_6_variable_core
Model universe used:
BM, EB, lambda, delta, mean_trend.
OU status: excluded/deferred, not unsupported.

Why OU is absent here:
Full-tree OU fitting remained computationally intractable under the available hosted-runner/local compute setup.
The completed five-model/no-OU model-fit aggregation is therefore the explicit downstream model universe for this run.

How to interpret these outputs:
1. ancestral_root_original_climate_final6_5models_no_OU.csv is the per-variable ancestral reconstruction in original climate variables.
2. ancestral_root_projected_PCA_point_final6_5models_no_OU.csv and ancestral_root_projected_PCA_ellipses_final6_5models_no_OU.csv are visualization/diagnostic projections only.
3. ancestral_climatic_core_plot_final6_5models_no_OU.* is a figure derived from the projected root, not the inferential basis of reconstruction.
4. ancestral_core_diagnostics_final6_5models_no_OU.csv contains convex-hull, Mahalanobis, and density diagnostics.
5. biological_interpretation_summary_root_final6_5models_no_OU.txt is the manuscript-facing summary for this climate-set-specific Result Block 1 run.
6. In this completed run, all 6 variables selected lambda and downstream reconstruction used lambda_fastAnc.

What should not be claimed:
- Do not describe this as a six-model result.
- Do not describe PCA as ancestral inference.
- Do not treat OU as rejected; it was excluded/deferred for operational reasons.
- Do not overinterpret the projected ellipse as a true multivariate ancestral niche.

Upstream model-fit provenance:
stage=model_fit
mode=final_named_model_set_aggregation
model_set=5models_no_OU
OU_status=excluded_deferred
models_included=BM,EB,lambda,delta,mean_trend
variables_included=annual_mean_temperature,annual_precipitation,maximum_temperature_warmest_month,mean_temperature_coldest_quarter,mean_temperature_driest_quarter,minimum_temperature_coldest_month,precipitation_in_coldest_quarter,precipitation_of_coldest_quarter,precipitation_of_wettest_quarter,precipitation_wettest_month,temperature_annual_range,temperature_seasonality
attempt_rows=60
completed_at_utc=2026-04-23 10:07:47 UTC
scientific_note=This aggregation changes only model-set inclusion for operational readiness; it does not change geiger fitting logic or primary reconstruction design.

Upstream five-model biological summary:
Five-model no-OU model-fit interpretation summary
=================================================

Scope:
This summary compares BM, EB, lambda, delta, and mean_trend fits for Result Block 1.
OU was explicitly excluded/deferred from this aggregation because OU shards exceeded the hosted-runner execution limit.
This is not the final six-model universe and should be labeled as a five-model/no-OU comparison.

Biological model meanings:
BM = neutral drift baseline.
EB = early rapid climatic change.
lambda = strong phylogenetic structure / signal.
delta = temporal shift in rate concentration.
mean_trend = directional change.

Best-model counts:
BM: 0
EB: 0
lambda: 12
delta: 0
mean_trend: 0

Dominant best-fit mode in this five-model comparison: lambda (strong phylogenetic structure / signal).

Grouped summary by climate-variable type:
precipitation / lambda: 5 variable(s); annual_precipitation, precipitation_in_coldest_quarter, precipitation_of_coldest_quarter, precipitation_of_wettest_quarter, precipitation_wettest_month
thermal / lambda: 7 variable(s); annual_mean_temperature, maximum_temperature_warmest_month, mean_temperature_coldest_quarter, mean_temperature_driest_quarter, minimum_temperature_coldest_month, temperature_annual_range, temperature_seasonality

Parameter notes for lambda/delta/mean_trend winners:
annual_mean_temperature winner=lambda; parameters: lambda=0.933278; sigsq=158.34; z0=207.155
annual_precipitation winner=lambda; parameters: lambda=0.911117; sigsq=28593.2; z0=1340.15
maximum_temperature_warmest_month winner=lambda; parameters: lambda=0.919498; sigsq=83.2299; z0=287.22
mean_temperature_coldest_quarter winner=lambda; parameters: lambda=0.930474; sigsq=353.489; z0=168.018
mean_temperature_driest_quarter winner=lambda; parameters: lambda=0.914901; sigsq=269.866; z0=190.685
minimum_temperature_coldest_month winner=lambda; parameters: lambda=0.928897; sigsq=397.965; z0=124.442
precipitation_in_coldest_quarter winner=lambda; parameters: lambda=0.77366; sigsq=2120.61; z0=282.313
precipitation_of_coldest_quarter winner=lambda; parameters: lambda=0.773662; sigsq=2120.62; z0=282.313
precipitation_of_wettest_quarter winner=lambda; parameters: lambda=0.887992; sigsq=3867.47; z0=590.854
precipitation_wettest_month winner=lambda; parameters: lambda=0.883922; sigsq=466.155; z0=207.729
temperature_annual_range winner=lambda; parameters: lambda=0.916039; sigsq=316.046; z0=162.855
temperature_seasonality winner=lambda; parameters: lambda=0.921302; sigsq=298067; z0=2649.43

Caution flags:
none

Interpretation rule:
Treat weak AICc separations as suggestive rather than decisive, and do not over-interpret noisy transform or trend estimates.
