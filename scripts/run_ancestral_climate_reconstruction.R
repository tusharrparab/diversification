# ============================================================
# Result Block 1: strict ancestral climatic reconstruction
# ============================================================
#
# README-style method note
# ------------------------
# This script asks: "Where did the radiation begin in climatic space?"
#
# The inferential rule is deliberately strict:
#   1. Ancestral climate is reconstructed in the original environmental
#      variables, never on PCA axes.
#   2. PCA is computed only from extant, scaled climate variables and is used
#      only as a visualization space.
#   3. Each climate variable is first fit with candidate evolutionary models
#      using geiger::fitContinuous(): BM, OU, EB, lambda, delta, mean_trend.
#   4. Primary root reconstruction is limited to models with an explicit,
#      defensible reconstruction path here: BM, lambda, and delta.
#   5. OU, EB, and mean_trend are fit-and-flag sensitivity models. They can win
#      AICc, but they are not silently forced into the primary root ellipse.
#
# Approximate pieces that remain:
#   - Root uncertainty is propagated to PCA space with a diagonal covariance
#     matrix across original climate variables. Cross-variable covariance in
#     ancestral uncertainty is not estimated.
#   - ape::ace() confidence intervals are converted to an approximate variance
#     as (CI width / 2 / 1.96)^2 when no direct variance is returned.
#   - If a variable has no valid primary reconstruction, its extant mean is used
#     only to place the root point in PCA space; it contributes no uncertainty
#     to the root ellipse and is explicitly flagged.

req_pkgs <- c(
  "ape", "geiger", "phytools", "dplyr", "readr", "stringr", "ggplot2",
  "ggrepel", "tibble", "writexl", "scales", "Matrix"
)

to_install <- req_pkgs[!req_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(phytools)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(tibble)
  library(writexl)
  library(scales)
})

tree_file <- Sys.getenv("TREE_FILE", "data/raw/pruned_phylogeny.nex")
climate_file <- Sys.getenv("CLIMATE_FILE", "data/raw/climate.csv")
output_dir <- Sys.getenv("OUTPUT_DIR", "results/result_block_1_strict_ancestral_climate")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

clean_text <- function(x) {
  x |>
    as.character() |>
    str_trim() |>
    str_replace_all("[\t\r\n]", " ") |>
    str_replace_all("\\s+", " ") |>
    str_to_lower()
}

normalize_species <- function(x) {
  x |>
    clean_text() |>
    str_replace_all(" ", "_") |>
    str_replace_all("_+", "_")
}

standardize_colname <- function(x) {
  x |>
    clean_text() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")
}

pick_existing_column <- function(df, candidates, required = TRUE) {
  standardized_names <- standardize_colname(names(df))
  standardized_candidates <- standardize_colname(candidates)
  hit_idx <- match(standardized_candidates, standardized_names)
  hit_idx <- hit_idx[!is.na(hit_idx)]

  if (length(hit_idx) == 0) {
    if (required) {
      stop(
        paste0(
          "None of these column names were found in climate.csv: ",
          paste(candidates, collapse = ", ")
        )
      )
    }
    return(NA_character_)
  }

  names(df)[hit_idx[1]]
}

read_tree_robust <- function(file) {
  obj <- tryCatch(read.nexus(file), error = function(e) NULL)

  if (is.null(obj)) {
    message("read.nexus() failed. Trying read.tree() instead.")
    obj <- tryCatch(read.tree(file), error = function(e) NULL)
  }

  if (is.null(obj)) {
    stop("Could not read phylogeny from: ", file)
  }

  if (inherits(obj, "multiPhylo")) {
    message("Multiple trees detected. Using the first tree.")
    obj <- obj[[1]]
  }

  if (!inherits(obj, "phylo")) {
    stop("Tree file was read, but result is not a valid phylo object.")
  }

  if (is.null(obj$tip.label) || length(obj$tip.label) == 0) {
    stop("Tree has no tip labels.")
  }

  obj
}

make_ellipse <- function(center, cov_mat, level = 0.95, n = 300) {
  theta <- seq(0, 2 * pi, length.out = n)
  circle <- cbind(cos(theta), sin(theta))
  radius <- sqrt(stats::qchisq(level, df = 2))
  eig <- eigen(cov_mat, symmetric = TRUE)
  transform <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = 2)
  pts <- sweep(circle %*% t(transform) * radius, 2, center, FUN = "+")
  out <- as.data.frame(pts)
  colnames(out) <- c("x", "y")
  out$level <- paste0(round(level * 100), "%")
  out$level_probability <- level
  out
}

climate_vars <- c(
  "annual_mean_temperature",
  "annual_precipitation",
  "maximum_temperature_warmest_month",
  "mean_temperature_coldest_quarter",
  "mean_temperature_driest_quarter",
  "minimum_temperature_coldest_month",
  "precipitation_in_coldest_quarter",
  "precipitation_of_coldest_quarter",
  "precipitation_of_wettest_quarter",
  "precipitation_wettest_month",
  "temperature_annual_range",
  "temperature_seasonality"
)

thermal_vars <- c(
  "annual_mean_temperature",
  "mean_temperature_coldest_quarter",
  "mean_temperature_driest_quarter",
  "maximum_temperature_warmest_month",
  "minimum_temperature_coldest_month",
  "temperature_annual_range",
  "temperature_seasonality"
)

precip_vars <- c(
  "annual_precipitation",
  "precipitation_in_coldest_quarter",
  "precipitation_of_coldest_quarter",
  "precipitation_of_wettest_quarter",
  "precipitation_wettest_month"
)

tree <- read_tree_robust(tree_file)
write_csv(
  tibble(tree_tip_label_raw = tree$tip.label),
  file.path(output_dir, "tree_tip_labels_raw.csv")
)

clim <- read_csv(climate_file, show_col_types = FALSE)

species_col <- pick_existing_column(clim, c("species"))
col_map <- list(
  annual_mean_temperature = c("annual_mean_temperature", "Annual mean temperature_mean"),
  annual_precipitation = c("annual_precipitation", "Annual precipitation_mean"),
  maximum_temperature_warmest_month = c(
    "maximum_temperature_warmest_month",
    "Max Temperature_warmest month_mean"
  ),
  mean_temperature_coldest_quarter = c(
    "mean_temperature_coldest_quarter",
    "Mean Temp coldest Quarter_mean"
  ),
  mean_temperature_driest_quarter = c(
    "mean_temperature_driest_quarter",
    "Mean Temp Driest Quarter_mean"
  ),
  minimum_temperature_coldest_month = c(
    "minimum_temperature_coldest_month",
    "Min Temp_Coldest Month_mean"
  ),
  precipitation_in_coldest_quarter = c(
    "precipitation_in_coldest_quarter",
    "Precipitation in Coldest Quarter_mean"
  ),
  precipitation_of_coldest_quarter = c(
    "precipitation_of_coldest_quarter",
    "Precipitation of coldest Quarter_mean"
  ),
  precipitation_of_wettest_quarter = c(
    "precipitation_of_wettest_quarter",
    "Precipitation of Wettest Quarter_mean"
  ),
  precipitation_wettest_month = c(
    "precipitation_wetest_month",
    "precipitation_wettest_month",
    "Precipitation Wetest Month_mean",
    "Precipitation Wettest Month_mean"
  ),
  temperature_annual_range = c("temperature_annual_range", "Temp Annual Range_mean"),
  temperature_seasonality = c(
    "temperature_seasonality",
    "Temp Seasonality (Standard Deviation)_mean"
  )
)

selected_cols <- lapply(col_map, function(candidates) pick_existing_column(clim, candidates))

write_csv(
  tibble(species_raw = clim[[species_col]]),
  file.path(output_dir, "climate_species_raw.csv")
)

clim2 <- clim %>%
  transmute(
    species_original = .data[[species_col]],
    species_clean = clean_text(.data[[species_col]]),
    tree_label = normalize_species(.data[[species_col]]),
    annual_mean_temperature = as.numeric(.data[[selected_cols$annual_mean_temperature]]),
    annual_precipitation = as.numeric(.data[[selected_cols$annual_precipitation]]),
    maximum_temperature_warmest_month = as.numeric(.data[[selected_cols$maximum_temperature_warmest_month]]),
    mean_temperature_coldest_quarter = as.numeric(.data[[selected_cols$mean_temperature_coldest_quarter]]),
    mean_temperature_driest_quarter = as.numeric(.data[[selected_cols$mean_temperature_driest_quarter]]),
    minimum_temperature_coldest_month = as.numeric(.data[[selected_cols$minimum_temperature_coldest_month]]),
    precipitation_in_coldest_quarter = as.numeric(.data[[selected_cols$precipitation_in_coldest_quarter]]),
    precipitation_of_coldest_quarter = as.numeric(.data[[selected_cols$precipitation_of_coldest_quarter]]),
    precipitation_of_wettest_quarter = as.numeric(.data[[selected_cols$precipitation_of_wettest_quarter]]),
    precipitation_wettest_month = as.numeric(.data[[selected_cols$precipitation_wettest_month]]),
    temperature_annual_range = as.numeric(.data[[selected_cols$temperature_annual_range]]),
    temperature_seasonality = as.numeric(.data[[selected_cols$temperature_seasonality]])
  )

duplicate_species <- clim2 %>%
  count(tree_label, sort = TRUE) %>%
  filter(n > 1)

if (nrow(duplicate_species) > 0) {
  write_csv(duplicate_species, file.path(output_dir, "duplicate_climate_species.csv"))
  clim2 <- clim2 %>% distinct(tree_label, .keep_all = TRUE)
}

tree$tip.label <- normalize_species(tree$tip.label)

tree_tips <- tree$tip.label
clim_species <- clim2$tree_label
overlap <- intersect(tree_tips, clim_species)

match_summary <- tibble(
  tree_tips = length(tree_tips),
  climate_species = length(clim_species),
  matched_species_before_missing_filter = length(overlap),
  tree_is_rooted_before_pruning = is.rooted(tree)
)

write_csv(match_summary, file.path(output_dir, "match_summary.csv"))

if (length(overlap) == 0) {
  stop(
    "Zero overlap between tree tip labels and climate species names after normalization. ",
    "Inspect tree_tip_labels_raw.csv and climate_species_raw.csv."
  )
}

write_csv(
  clim2 %>% filter(!tree_label %in% tree_tips),
  file.path(output_dir, "species_removed_from_csv_not_in_tree.csv")
)

clim_complete <- clim2 %>%
  filter(tree_label %in% tree_tips) %>%
  filter(if_all(all_of(climate_vars), ~ !is.na(.)))

matched_species <- intersect(tree$tip.label, clim_complete$tree_label)

if (length(matched_species) < 3) {
  stop("Fewer than 3 species remain after matching tree and complete climate data.")
}

tips_to_drop <- setdiff(tree$tip.label, matched_species)
write_csv(
  tibble(tree_tip_missing_complete_climate = tips_to_drop),
  file.path(output_dir, "tree_tips_missing_climate.csv")
)

tree_pruned <- if (length(tips_to_drop) > 0) drop.tip(tree, tips_to_drop) else tree

if (is.null(tree_pruned) || !inherits(tree_pruned, "phylo")) {
  stop("Pruned tree is invalid.")
}

clim_final <- clim_complete %>%
  filter(tree_label %in% tree_pruned$tip.label) %>%
  slice(match(tree_pruned$tip.label, tree_label))

if (!all(clim_final$tree_label == tree_pruned$tip.label)) {
  stop("Ordering mismatch between climate data and tree tips.")
}

write_csv(clim_final, file.path(output_dir, "matched_climate_data_used.csv"))

X <- clim_final %>%
  select(all_of(climate_vars)) %>%
  as.data.frame()
rownames(X) <- clim_final$tree_label

candidate_models <- c("BM", "OU", "EB", "lambda", "delta", "mean_trend")
primary_models <- c("BM", "lambda", "delta")
sensitivity_models <- setdiff(candidate_models, primary_models)
crown_node <- as.character(Ntip(tree_pruned) + 1)
geiger_ncores <- suppressWarnings(as.integer(Sys.getenv("GEIGER_NCORES", "1")))
if (is.na(geiger_ncores) || geiger_ncores < 1) geiger_ncores <- 1

collapse_flags <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) "none" else paste(unique(x), collapse = "; ")
}

extract_fit_stat <- function(fit_obj, candidates) {
  for (container in list(fit_obj, fit_obj$opt)) {
    if (is.null(container)) next
    for (nm in candidates) {
      if (!is.null(container[[nm]])) {
        val <- suppressWarnings(as.numeric(container[[nm]]))
        if (length(val) == 1 && is.finite(val)) return(val)
      }
    }
  }
  NA_real_
}

extract_fit_parameter <- function(fit_obj, candidates) {
  containers <- list(fit_obj$opt, fit_obj$opt$par, fit_obj)
  for (container in containers) {
    if (is.null(container)) next
    for (nm in candidates) {
      if (!is.null(container[[nm]])) {
        val <- suppressWarnings(as.numeric(container[[nm]]))
        if (length(val) == 1 && is.finite(val)) return(val)
      }
    }
  }
  NA_real_
}

summarize_fit_parameters <- function(fit_obj) {
  if (is.null(fit_obj$opt)) return(NA_character_)
  opt <- fit_obj$opt
  drop_names <- c(
    "lnL", "logLik", "loglik", "aic", "AIC", "aicc", "AICc", "k",
    "convergence", "counts", "message", "method"
  )
  scalar_names <- names(opt)[vapply(opt, function(z) {
    is.numeric(z) && length(z) == 1 && is.finite(z)
  }, logical(1))]
  scalar_names <- setdiff(scalar_names, drop_names)
  if (length(scalar_names) == 0) return(NA_character_)
  paste(
    sprintf("%s=%.6g", scalar_names, vapply(opt[scalar_names], as.numeric, numeric(1))),
    collapse = "; "
  )
}

fit_one_model <- function(tree, trait, variable, model) {
  warnings_seen <- character()
  result <- tryCatch(
    {
      fit <- withCallingHandlers(
        geiger::fitContinuous(
          phy = tree,
          dat = trait,
          model = model,
          control = list(
            method = c("subplex", "L-BFGS-B"),
            niter = 100,
            hessian = FALSE,
            CI = 0.95
          ),
          ncores = geiger_ncores
        ),
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      aicc <- extract_fit_stat(fit, c("aicc", "AICc"))
      fit_status <- if (is.finite(aicc)) "ok" else "fit_returned_no_finite_AICc"
      list(
        row = tibble(
          variable = variable,
          model = model,
          model_role = if_else(model %in% primary_models, "primary", "sensitivity_only"),
          fit_status = fit_status,
          logLik = extract_fit_stat(fit, c("lnL", "logLik", "loglik")),
          AIC = extract_fit_stat(fit, c("aic", "AIC")),
          AICc = aicc,
          parameter_summary = summarize_fit_parameters(fit),
          fit_warning = collapse_flags(warnings_seen),
          fit_error = NA_character_
        ),
        fit = fit
      )
    },
    error = function(e) {
      list(
        row = tibble(
          variable = variable,
          model = model,
          model_role = if_else(model %in% primary_models, "primary", "sensitivity_only"),
          fit_status = "failed",
          logLik = NA_real_,
          AIC = NA_real_,
          AICc = NA_real_,
          parameter_summary = NA_character_,
          fit_warning = collapse_flags(warnings_seen),
          fit_error = conditionMessage(e)
        ),
        fit = NULL
      )
    }
  )
  result
}

add_model_deltas <- function(model_fit_table) {
  model_fit_table %>%
    group_by(variable) %>%
    mutate(
      best_AICc_for_variable = if (any(is.finite(AICc))) min(AICc, na.rm = TRUE) else NA_real_,
      delta_AICc = if_else(is.finite(AICc), AICc - best_AICc_for_variable, NA_real_)
    ) %>%
    ungroup()
}

extract_root_value <- function(ace_like_obj, tree) {
  root_node <- as.character(Ntip(tree) + 1)
  ace_values <- ace_like_obj$ace
  if (is.null(ace_values)) return(NA_real_)
  if (!is.null(names(ace_values)) && root_node %in% names(ace_values)) {
    return(as.numeric(ace_values[root_node]))
  }
  as.numeric(ace_values[1])
}

extract_root_ci <- function(ace_obj, tree) {
  root_node <- as.character(Ntip(tree) + 1)
  ci <- ace_obj$CI95
  if (is.null(ci) || !is.matrix(ci) || ncol(ci) < 2) return(c(NA_real_, NA_real_))
  if (!is.null(rownames(ci)) && root_node %in% rownames(ci)) {
    return(as.numeric(ci[root_node, 1:2]))
  }
  as.numeric(ci[1, 1:2])
}

variance_from_ci <- function(ci_low, ci_high) {
  if (!is.finite(ci_low) || !is.finite(ci_high) || ci_high <= ci_low) {
    return(NA_real_)
  }
  ((ci_high - ci_low) / (2 * stats::qnorm(0.975)))^2
}

transform_tree_for_primary_model <- function(tree, fit_obj, primary_model) {
  if (primary_model == "BM") {
    return(list(tree = tree, transform_parameter = NA_real_, transform_status = "not_applicable"))
  }

  if (primary_model == "lambda") {
    lambda <- extract_fit_parameter(fit_obj, c("lambda"))
    if (!is.finite(lambda)) stop("lambda model selected, but no finite lambda parameter was found")
    return(list(
      tree = geiger::lambdaTree(tree, lambda = lambda),
      transform_parameter = lambda,
      transform_status = "lambdaTree"
    ))
  }

  if (primary_model == "delta") {
    delta <- extract_fit_parameter(fit_obj, c("delta"))
    if (!is.finite(delta) || delta <= 0) {
      stop("delta model selected, but no finite positive delta parameter was found")
    }
    return(list(
      tree = geiger::deltaTree(tree, delta = delta, rescale = TRUE),
      transform_parameter = delta,
      transform_status = "deltaTree_rescaled"
    ))
  }

  stop("Unsupported primary reconstruction model: ", primary_model)
}

reconstruct_primary_root <- function(tree, trait, variable, primary_model, fit_obj) {
  tryCatch(
    {
      tree_info <- transform_tree_for_primary_model(tree, fit_obj, primary_model)
      recon_tree <- tree_info$tree
      trait_reordered <- trait[recon_tree$tip.label]

      if (any(!is.finite(trait_reordered))) {
        stop("Trait vector contains non-finite values after tree reordering")
      }

      ace_fit <- ape::ace(
        x = trait_reordered,
        phy = recon_tree,
        type = "continuous",
        method = "ML",
        model = "BM",
        CI = TRUE
      )

      ci <- extract_root_ci(ace_fit, recon_tree)
      root_variance <- variance_from_ci(ci[1], ci[2])

      tibble(
        variable = variable,
        primary_model_used = primary_model,
        transform_parameter = tree_info$transform_parameter,
        transform_status = tree_info$transform_status,
        root_estimate = extract_root_value(ace_fit, recon_tree),
        root_variance = root_variance,
        CI_low = ci[1],
        CI_high = ci[2],
        root_variance_source = if_else(
          is.finite(root_variance),
          "CI95_width_assuming_normality",
          "unavailable"
        ),
        reconstruction_status = if_else(
          is.finite(root_variance),
          "ok",
          "ok_root_variance_unavailable"
        ),
        reconstruction_error = NA_character_
      )
    },
    error = function(e) {
      tibble(
        variable = variable,
        primary_model_used = primary_model,
        transform_parameter = NA_real_,
        transform_status = "failed",
        root_estimate = NA_real_,
        root_variance = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        root_variance_source = "unavailable",
        reconstruction_status = "failed",
        reconstruction_error = conditionMessage(e)
      )
    }
  )
}

compute_sensitivity_root <- function(tree, trait, variable, sensitivity_model) {
  if (!sensitivity_model %in% sensitivity_models) {
    return(tibble(
      variable = variable,
      sensitivity_model = NA_character_,
      sensitivity_root_estimate = NA_real_,
      sensitivity_status = "not_requested",
      sensitivity_error = NA_character_
    ))
  }

  if (sensitivity_model == "OU") {
    return(tibble(
      variable = variable,
      sensitivity_model = "OU",
      sensitivity_root_estimate = NA_real_,
      sensitivity_status = "flag_only_OU_not_forced_into_primary_reconstruction",
      sensitivity_error = NA_character_
    ))
  }

  tryCatch(
    {
      if (sensitivity_model == "EB") {
        fit <- phytools::anc.ML(tree, trait[tree$tip.label], model = "EB")
        convergence <- if (!is.null(fit$convergence)) as.character(fit$convergence) else NA_character_
        status <- if (!is.na(convergence) && convergence == "0") {
          "ok"
        } else {
          paste0("estimate_returned_optimizer_code_", convergence)
        }
      } else if (sensitivity_model == "mean_trend") {
        fit <- phytools::anc.trend(tree, trait[tree$tip.label])
        convergence <- if (!is.null(fit$convergence)) as.character(fit$convergence) else NA_character_
        status <- if (!is.na(convergence) && convergence == "0") {
          "ok"
        } else {
          paste0("estimate_returned_optimizer_code_", convergence)
        }
      } else {
        stop("No sensitivity reconstruction implemented for model: ", sensitivity_model)
      }

      tibble(
        variable = variable,
        sensitivity_model = sensitivity_model,
        sensitivity_root_estimate = extract_root_value(fit, tree),
        sensitivity_status = status,
        sensitivity_error = NA_character_
      )
    },
    error = function(e) {
      tibble(
        variable = variable,
        sensitivity_model = sensitivity_model,
        sensitivity_root_estimate = NA_real_,
        sensitivity_status = "failed",
        sensitivity_error = conditionMessage(e)
      )
    }
  )
}

select_best_row <- function(rows) {
  rows %>%
    filter(is.finite(AICc)) %>%
    arrange(AICc, model) %>%
    slice(1)
}

make_positive_definite <- function(cov_mat) {
  if (any(!is.finite(cov_mat))) {
    stop("Covariance matrix contains non-finite values")
  }
  eig <- eigen(cov_mat, symmetric = TRUE)$values
  if (any(eig <= sqrt(.Machine$double.eps))) {
    cov_mat <- as.matrix(Matrix::nearPD(cov_mat, corr = FALSE)$mat)
  }
  cov_mat
}

# -----------------------------
# Primary inference: model fitting and original-variable root states
# -----------------------------
message("Fitting per-variable evolutionary models with geiger::fitContinuous().")

fit_objects <- setNames(vector("list", length(climate_vars)), climate_vars)
fit_rows <- list()

for (v in climate_vars) {
  trait <- X[[v]]
  names(trait) <- rownames(X)
  trait <- trait[tree_pruned$tip.label]
  fit_objects[[v]] <- list()

  for (m in candidate_models) {
    message("  fitting ", v, " under ", m)
    fit_result <- fit_one_model(tree_pruned, trait, variable = v, model = m)
    fit_rows[[length(fit_rows) + 1]] <- fit_result$row
    fit_objects[[v]][[m]] <- fit_result$fit
  }
}

model_fit_table <- bind_rows(fit_rows) %>%
  add_model_deltas()

reconstruction_rows <- list()

for (v in climate_vars) {
  trait <- X[[v]]
  names(trait) <- rownames(X)
  trait <- trait[tree_pruned$tip.label]

  rows_v <- model_fit_table %>% filter(variable == v)
  best_overall <- select_best_row(rows_v)
  best_primary <- select_best_row(rows_v %>% filter(model %in% primary_models))

  best_overall_model <- if (nrow(best_overall) == 1) best_overall$model else NA_character_
  best_overall_aicc <- if (nrow(best_overall) == 1) best_overall$AICc else NA_real_
  best_primary_model <- if (nrow(best_primary) == 1) best_primary$model else NA_character_
  best_primary_aicc <- if (nrow(best_primary) == 1) best_primary$AICc else NA_real_
  delta_primary_minus_overall <- best_primary_aicc - best_overall_aicc

  if (is.na(best_primary_model)) {
    primary_recon <- tibble(
      variable = v,
      primary_model_used = NA_character_,
      transform_parameter = NA_real_,
      transform_status = "not_available",
      root_estimate = NA_real_,
      root_variance = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      root_variance_source = "unavailable",
      reconstruction_status = "failed_no_primary_model_fit",
      reconstruction_error = "No finite AICc fit among BM, lambda, or delta"
    )
  } else {
    primary_recon <- reconstruct_primary_root(
      tree = tree_pruned,
      trait = trait,
      variable = v,
      primary_model = best_primary_model,
      fit_obj = fit_objects[[v]][[best_primary_model]]
    )
  }

  sensitivity_recon <- if (best_overall_model %in% sensitivity_models) {
    compute_sensitivity_root(tree_pruned, trait, variable = v, sensitivity_model = best_overall_model)
  } else {
    tibble(
      variable = v,
      sensitivity_model = NA_character_,
      sensitivity_root_estimate = NA_real_,
      sensitivity_status = "not_needed_best_overall_is_primary",
      sensitivity_error = NA_character_
    )
  }

  caution_flag <- collapse_flags(c(
    if (best_overall_model %in% sensitivity_models) {
      paste0("best_overall_model_is_sensitivity_only_", best_overall_model)
    },
    if (!is.finite(best_primary_aicc)) "no_supported_primary_model_with_finite_AICc",
    if (!is.finite(primary_recon$root_estimate)) "primary_root_estimate_unavailable",
    if (!is.finite(primary_recon$root_variance)) "primary_root_uncertainty_unavailable",
    if (is.finite(delta_primary_minus_overall) && delta_primary_minus_overall > 2) {
      "primary_model_substantially_worse_than_overall_AICc"
    }
  ))

  reconstruction_rows[[length(reconstruction_rows) + 1]] <- primary_recon %>%
    left_join(sensitivity_recon, by = "variable") %>%
    mutate(
      best_overall_model = best_overall_model,
      best_overall_AICc = best_overall_aicc,
      best_primary_model = best_primary_model,
      best_primary_AICc = best_primary_aicc,
      delta_AICc_primary_minus_overall = delta_primary_minus_overall,
      extant_mean = mean(trait, na.rm = TRUE),
      caution_flag = caution_flag,
      .after = variable
    )
}

ancestral_reconstruction_summary <- bind_rows(reconstruction_rows) %>%
  select(
    variable,
    best_overall_model,
    best_overall_AICc,
    best_primary_model,
    best_primary_AICc,
    delta_AICc_primary_minus_overall,
    primary_model_used,
    transform_parameter,
    transform_status,
    root_estimate,
    root_variance,
    root_variance_source,
    CI_low,
    CI_high,
    reconstruction_status,
    reconstruction_error,
    caution_flag,
    sensitivity_model,
    sensitivity_root_estimate,
    sensitivity_status,
    sensitivity_error,
    extant_mean
  )

# For unresolved variables, use extant means only to place a point in PCA
# space. They are not allowed to add unsupported uncertainty to the ellipse.
root_projection_table <- ancestral_reconstruction_summary %>%
  mutate(
    projection_value = if_else(is.finite(root_estimate), root_estimate, extant_mean),
    projection_value_source = if_else(
      is.finite(root_estimate),
      "primary_reconstruction",
      "extant_mean_fallback_for_projection_only"
    ),
    variance_used_for_uncertainty = if_else(
      is.finite(root_variance) & root_variance >= 0,
      root_variance,
      0
    ),
    uncertainty_source = if_else(
      is.finite(root_variance) & root_variance >= 0,
      root_variance_source,
      "omitted_from_uncertainty"
    )
  ) %>%
  select(
    variable,
    projection_value,
    projection_value_source,
    variance_used_for_uncertainty,
    uncertainty_source
  )

# -----------------------------
# Visualization-only PCA
# -----------------------------
X_scaled <- scale(X)
x_center <- attr(X_scaled, "scaled:center")
x_scale <- attr(X_scaled, "scaled:scale")

pca_vis <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca_vis$x)
colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
scores$tree_label <- rownames(X)

loadings <- as.data.frame(pca_vis$rotation)
colnames(loadings) <- paste0("PC", seq_len(ncol(loadings)))
loadings$variable <- rownames(loadings)

eigvals <- pca_vis$sdev^2
eigenvalues <- tibble(
  PC = paste0("PC", seq_along(eigvals)),
  eigenvalue = eigvals,
  variance_explained = eigvals / sum(eigvals),
  cumulative_variance = cumsum(eigvals / sum(eigvals))
)

pc_count <- length(eigvals)
candidate_pcs <- seq_len(min(4, pc_count))

thermal_strength <- sapply(candidate_pcs, function(i) {
  sum(abs(loadings[loadings$variable %in% thermal_vars, paste0("PC", i)]), na.rm = TRUE)
})

precip_strength <- sapply(candidate_pcs, function(i) {
  sum(abs(loadings[loadings$variable %in% precip_vars, paste0("PC", i)]), na.rm = TRUE)
})

pc1_idx <- candidate_pcs[which.max(thermal_strength)]
remaining <- setdiff(candidate_pcs, pc1_idx)
if (length(remaining) == 0) {
  stop("Could not choose a second PCA axis.")
}
pc2_idx <- remaining[which.max(precip_strength[match(remaining, candidate_pcs)])]

pc1_sign <- sign(sum(loadings[loadings$variable %in% thermal_vars, paste0("PC", pc1_idx)], na.rm = TRUE))
if (is.na(pc1_sign) || pc1_sign == 0) pc1_sign <- 1

pc2_sign <- sign(sum(loadings[loadings$variable %in% precip_vars, paste0("PC", pc2_idx)], na.rm = TRUE))
if (is.na(pc2_sign) || pc2_sign == 0) pc2_sign <- 1

scores$PC1_thermal <- scores[[paste0("PC", pc1_idx)]] * pc1_sign
scores$PC2_precipitation <- scores[[paste0("PC", pc2_idx)]] * pc2_sign

loadings$PC1_thermal <- loadings[[paste0("PC", pc1_idx)]] * pc1_sign
loadings$PC2_precipitation <- loadings[[paste0("PC", pc2_idx)]] * pc2_sign

root_projection_vector <- root_projection_table$projection_value
names(root_projection_vector) <- root_projection_table$variable

root_scaled <- (root_projection_vector[names(x_center)] - x_center) / x_scale
root_scores_all <- as.numeric(root_scaled %*% pca_vis$rotation)
names(root_scores_all) <- colnames(pca_vis$x)

root_center <- c(
  PC1_thermal = unname(root_scores_all[paste0("PC", pc1_idx)]) * pc1_sign,
  PC2_precipitation = unname(root_scores_all[paste0("PC", pc2_idx)]) * pc2_sign
)

projected_root_point <- tibble(
  crown_node = crown_node,
  tree_is_rooted_after_pruning = is.rooted(tree_pruned),
  n_species_used = length(tree_pruned$tip.label),
  PC1_thermal = as.numeric(root_center["PC1_thermal"]),
  PC2_precipitation = as.numeric(root_center["PC2_precipitation"]),
  projection_note = "root reconstructed in original variables; PCA used only for visualization"
)

# -----------------------------
# Approximate uncertainty propagation into PCA space
# -----------------------------
root_var_vec <- root_projection_table$variance_used_for_uncertainty
names(root_var_vec) <- root_projection_table$variable

Sigma_root_original <- diag(root_var_vec[names(x_center)])
rownames(Sigma_root_original) <- names(x_center)
colnames(Sigma_root_original) <- names(x_center)

S_inv <- diag(1 / x_scale)
Sigma_root_scaled <- S_inv %*% Sigma_root_original %*% S_inv
Sigma_root_pca <- t(pca_vis$rotation) %*% Sigma_root_scaled %*% pca_vis$rotation

sign_mat <- diag(c(pc1_sign, pc2_sign))
Sigma2 <- Sigma_root_pca[c(pc1_idx, pc2_idx), c(pc1_idx, pc2_idx), drop = FALSE]
Sigma2 <- sign_mat %*% Sigma2 %*% sign_mat
Sigma2 <- make_positive_definite(Sigma2)

ellipse_eigenvalues <- eigen(Sigma2, symmetric = TRUE)$values
ellipse_status <- if (all(ellipse_eigenvalues > 0)) "ok" else "nearPD_returned_nonpositive_eigenvalue"

ell <- bind_rows(
  make_ellipse(root_center, Sigma2, level = 0.50),
  make_ellipse(root_center, Sigma2, level = 0.80),
  make_ellipse(root_center, Sigma2, level = 0.95)
)

uncertainty_assumptions <- tibble(
  assumption = c(
    "original_variable_covariance",
    "root_variance_source",
    "unresolved_variables",
    "positive_definite_guard"
  ),
  status = c(
    "diagonal_only_cross_variable_covariance_not_estimated",
    "CI95_width_converted_to_variance_when_needed",
    collapse_flags(root_projection_table$variable[root_projection_table$uncertainty_source == "omitted_from_uncertainty"]),
    ellipse_status
  )
)

scores_out <- scores %>%
  left_join(
    clim_final %>% select(tree_label, species_original, species_clean),
    by = "tree_label"
  ) %>%
  relocate(species_original, species_clean, tree_label, PC1_thermal, PC2_precipitation)

d2 <- as.numeric(mahalanobis(
  x = scores_out[, c("PC1_thermal", "PC2_precipitation")],
  center = as.numeric(root_center),
  cov = Sigma2
))

core_membership <- scores_out %>%
  transmute(
    species_original,
    species_clean,
    tree_label,
    PC1_thermal,
    PC2_precipitation,
    mahalanobis_d2 = d2,
    inside_50_core = mahalanobis_d2 <= qchisq(0.50, df = 2),
    inside_80_core = mahalanobis_d2 <= qchisq(0.80, df = 2),
    inside_95_core = mahalanobis_d2 <= qchisq(0.95, df = 2),
    closest_core_level = case_when(
      inside_50_core ~ "50%",
      inside_80_core ~ "80%",
      inside_95_core ~ "95%",
      TRUE ~ "outside_95%"
    )
  ) %>%
  arrange(mahalanobis_d2)

core_summary <- tibble(
  result_block = "Result Block 1: Where did the radiation begin in climatic space?",
  n_species_used = nrow(core_membership),
  crown_node = crown_node,
  tree_is_rooted = is.rooted(tree_pruned),
  selected_thermal_axis = paste0("PC", pc1_idx),
  selected_precipitation_axis = paste0("PC", pc2_idx),
  thermal_axis_variance_explained = eigenvalues$variance_explained[pc1_idx],
  precipitation_axis_variance_explained = eigenvalues$variance_explained[pc2_idx],
  n_variables_with_sensitivity_best_model = sum(
    ancestral_reconstruction_summary$best_overall_model %in% sensitivity_models,
    na.rm = TRUE
  ),
  n_variables_with_primary_reconstruction = sum(is.finite(ancestral_reconstruction_summary$root_estimate)),
  n_variables_omitted_from_uncertainty = sum(root_projection_table$uncertainty_source == "omitted_from_uncertainty"),
  n_species_inside_50_core = sum(core_membership$inside_50_core),
  n_species_inside_80_core = sum(core_membership$inside_80_core),
  n_species_inside_95_core = sum(core_membership$inside_95_core),
  pct_species_inside_50_core = mean(core_membership$inside_50_core),
  pct_species_inside_80_core = mean(core_membership$inside_80_core),
  pct_species_inside_95_core = mean(core_membership$inside_95_core),
  uncertainty_assumption = "diagonal original-variable root uncertainty projected into PCA space"
)

write_csv(model_fit_table, file.path(output_dir, "strict_model_fits_by_variable.csv"))
write_csv(ancestral_reconstruction_summary, file.path(output_dir, "strict_ancestral_reconstruction_by_variable.csv"))
write_csv(root_projection_table, file.path(output_dir, "strict_root_vector_for_PCA_projection.csv"))
write_csv(uncertainty_assumptions, file.path(output_dir, "strict_uncertainty_assumptions.csv"))
write_csv(scores_out, file.path(output_dir, "visualization_PCA_species_scores.csv"))
write_csv(loadings, file.path(output_dir, "visualization_PCA_loadings.csv"))
write_csv(eigenvalues, file.path(output_dir, "visualization_PCA_eigenvalues.csv"))
write_csv(projected_root_point, file.path(output_dir, "projected_primary_root_PCA_point.csv"))
write_csv(ell, file.path(output_dir, "projected_primary_root_PCA_ellipses.csv"))
write_csv(core_membership, file.path(output_dir, "strict_ancestral_core_species_membership.csv"))
write_csv(core_summary, file.path(output_dir, "strict_ancestral_core_summary.csv"))

write_xlsx(
  list(
    summary = core_summary,
    model_fits_by_variable = model_fit_table,
    ancestral_reconstruction = ancestral_reconstruction_summary,
    root_vector_for_PCA_projection = root_projection_table,
    uncertainty_assumptions = uncertainty_assumptions,
    projected_root_PCA_point = projected_root_point,
    projected_root_PCA_ellipses = ell,
    species_core_membership = core_membership,
    visualization_PCA_scores = scores_out,
    visualization_PCA_loadings = loadings,
    visualization_PCA_eigenvalues = eigenvalues,
    data_used = clim_final
  ),
  path = file.path(output_dir, "strict_ancestral_climate_reconstruction_outputs.xlsx")
)

arrow_mult <- 2.5
thermal_arrows_df <- loadings %>%
  filter(variable %in% thermal_vars) %>%
  transmute(variable, x = PC1_thermal * arrow_mult, y = PC2_precipitation * arrow_mult)

precip_arrows_df <- loadings %>%
  filter(variable %in% precip_vars) %>%
  transmute(variable, x = PC1_thermal * arrow_mult, y = PC2_precipitation * arrow_mult)

label_df <- bind_rows(
  core_membership %>% filter(inside_50_core) %>% slice_head(n = 15),
  core_membership %>% slice_head(n = 10)
) %>%
  distinct(tree_label, .keep_all = TRUE)

var_pc1 <- percent(eigenvalues$variance_explained[pc1_idx], accuracy = 0.1)
var_pc2 <- percent(eigenvalues$variance_explained[pc2_idx], accuracy = 0.1)

xlab_txt <- paste0("Thermal climatic gradient (PC", pc1_idx, ", ", var_pc1, ")")
ylab_txt <- paste0("Precipitation climatic gradient (PC", pc2_idx, ", ", var_pc2, ")")

p_core <- ggplot(scores_out, aes(x = PC1_thermal, y = PC2_precipitation)) +
  geom_polygon(
    data = ell,
    aes(x = x, y = y, group = level, fill = level),
    alpha = 0.16,
    color = "black",
    linewidth = 0.3,
    inherit.aes = FALSE
  ) +
  geom_point(alpha = 0.28, size = 1.35, color = "grey35") +
  geom_point(
    data = core_membership %>% filter(inside_50_core),
    aes(x = PC1_thermal, y = PC2_precipitation),
    inherit.aes = FALSE,
    size = 1.9,
    color = "darkgreen",
    alpha = 0.8
  ) +
  geom_point(
    data = projected_root_point,
    aes(x = PC1_thermal, y = PC2_precipitation),
    inherit.aes = FALSE,
    size = 3.3,
    shape = 21,
    fill = "red",
    color = "black",
    stroke = 0.5
  ) +
  geom_segment(
    data = thermal_arrows_df,
    aes(x = 0, y = 0, xend = x, yend = y),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.18, "cm")),
    linewidth = 0.45,
    color = "firebrick"
  ) +
  geom_text(
    data = thermal_arrows_df,
    aes(x = x, y = y, label = variable),
    inherit.aes = FALSE,
    color = "firebrick",
    size = 2.8
  ) +
  geom_segment(
    data = precip_arrows_df,
    aes(x = 0, y = 0, xend = x, yend = y),
    inherit.aes = FALSE,
    arrow = grid::arrow(length = grid::unit(0.18, "cm")),
    linewidth = 0.45,
    color = "steelblue4"
  ) +
  geom_text(
    data = precip_arrows_df,
    aes(x = x, y = y, label = variable),
    inherit.aes = FALSE,
    color = "steelblue4",
    size = 2.8
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = species_clean),
    size = 2.5,
    max.overlaps = 100,
    box.padding = 0.18,
    point.padding = 0.1,
    min.segment.length = 0
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed", color = "grey55") +
  scale_fill_manual(values = c("50%" = "#1b9e77", "80%" = "#d95f02", "95%" = "#7570b3")) +
  labs(
    x = xlab_txt,
    y = ylab_txt,
    fill = "Projected root CI",
    title = "Result Block 1: ancestral climatic origin of the crown radiation",
    subtitle = "Root reconstructed in original climate variables; PCA used only for visualization"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(output_dir, "strict_ancestral_climatic_core.png"),
  plot = p_core,
  width = 10,
  height = 8,
  dpi = 400
)

ggsave(
  filename = file.path(output_dir, "strict_ancestral_climatic_core.pdf"),
  plot = p_core,
  width = 10,
  height = 8
)

write_csv(thermal_arrows_df, file.path(output_dir, "thermal_loading_arrows.csv"))
write_csv(precip_arrows_df, file.path(output_dir, "precipitation_loading_arrows.csv"))
write_csv(label_df, file.path(output_dir, "species_labels_used_in_plot.csv"))

summary_lines <- c(
  "# Strict Ancestral Climate Reconstruction Summary",
  "",
  "Result Block 1: Where did the radiation begin in climatic space?",
  "",
  "Core statement supported by this pipeline:",
  "\"Ancestral climate was reconstructed in the original environmental variables and then visualized in climatic PCA space.\"",
  "",
  paste0("Species used: ", core_summary$n_species_used),
  paste0("Crown node: ", crown_node),
  paste0("Tree rooted after pruning: ", core_summary$tree_is_rooted),
  paste0(
    "Selected visualization axes: thermal = ", core_summary$selected_thermal_axis,
    ", precipitation = ", core_summary$selected_precipitation_axis
  ),
  paste0("Variables where best overall model was sensitivity-only: ", core_summary$n_variables_with_sensitivity_best_model),
  paste0("Variables with primary reconstruction: ", core_summary$n_variables_with_primary_reconstruction),
  paste0("Variables omitted from root uncertainty: ", core_summary$n_variables_omitted_from_uncertainty),
  paste0("Species inside 50% projected ancestral core: ", core_summary$n_species_inside_50_core),
  paste0("Species inside 80% projected ancestral core: ", core_summary$n_species_inside_80_core),
  paste0("Species inside 95% projected ancestral core: ", core_summary$n_species_inside_95_core),
  "",
  "Primary inference models: BM, lambda, delta.",
  "Sensitivity-only models: OU, EB, mean_trend.",
  "OU is fit and flagged only; it is not forced into primary ancestral reconstruction.",
  "EB and mean_trend get optional sensitivity-only root estimates when they are the best overall AICc model.",
  "",
  "Approximation/assumption flags:",
  "- Root covariance is diagonal across original climate variables; cross-variable ancestral uncertainty is not estimated.",
  "- ace() confidence intervals are converted to approximate variances when no direct root variance is available.",
  "- Unresolved variables use extant means only for plotting the root point and are omitted from uncertainty propagation.",
  "",
  "Key outputs:",
  "- strict_model_fits_by_variable.csv",
  "- strict_ancestral_reconstruction_by_variable.csv",
  "- strict_root_vector_for_PCA_projection.csv",
  "- visualization_PCA_species_scores.csv",
  "- projected_primary_root_PCA_point.csv",
  "- projected_primary_root_PCA_ellipses.csv",
  "- strict_ancestral_climatic_core.png",
  "- strict_ancestral_climatic_core.pdf"
)

writeLines(summary_lines, con = file.path(output_dir, "SUMMARY.md"))

cat("\n==============================\n")
cat("Strict ancestral reconstruction summary\n")
cat("==============================\n")
print(ancestral_reconstruction_summary)

cat("\n==============================\n")
cat("Projected ancestral core summary\n")
cat("==============================\n")
print(core_summary)

cat("\nFiles written to: ", output_dir, "\n", sep = "")
