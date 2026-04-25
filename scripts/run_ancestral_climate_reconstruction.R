# ============================================================
# Result Block 1: strict ancestral climatic reconstruction
# ============================================================
#
# Scientific question
# -------------------
# Where did the radiation begin in climatic space?
#
# Core rule
# ---------
# Ancestral climate is reconstructed in the original climate variables.
# PCA is only a visualization layer used after reconstruction.
#
# Model policy
# ------------
# Each variable is fit with geiger::fitContinuous() under:
#   BM, OU, EB, lambda, delta, mean_trend
#
# Primary ancestral reconstruction is intentionally restricted to:
#   BM, lambda, delta
#
# OU, EB, and mean_trend are sensitivity-only models. They can be the
# best-supported AICc model, but they are never silently forced into the primary
# root estimate or root ellipse.
#
# Approximation policy
# --------------------
# Every approximation is written into code comments and output summaries:
#   - PCA is not inferential.
#   - The projected ellipse is not a fully jointly inferred ancestral niche.
#   - The ellipse uses a diagonal covariance approximation across original
#     climate variables, then projects that uncertainty into PCA space.
#   - ape::ace() CI95 widths are converted to approximate variances when needed.
#   - If a variable has no valid primary estimate, its extant mean can be used
#     only to place a root point in PCA space and does not add uncertainty.

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

# -----------------------------
# Runtime options
# -----------------------------
truthy <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "yes", "y")
}

get_cli_value <- function(args, key, default = NA_character_) {
  key_eq <- paste0("--", key, "=")
  hit_eq <- args[startsWith(args, key_eq)]
  if (length(hit_eq) > 0) {
    return(sub(key_eq, "", hit_eq[1], fixed = TRUE))
  }

  key_plain <- paste0("--", key)
  hit <- which(args == key_plain)
  if (length(hit) > 0 && hit[1] < length(args)) {
    return(args[hit[1] + 1])
  }

  default
}

has_cli_flag <- function(args, key) {
  paste0("--", key) %in% args
}

args <- commandArgs(trailingOnly = TRUE)

tree_file <- Sys.getenv("TREE_FILE", "data/raw/pruned_phylogeny.nex")
climate_file <- Sys.getenv("CLIMATE_FILE", "data/raw/climate.csv")
output_dir <- Sys.getenv("OUTPUT_DIR", "results/result_block_1_strict_ancestral_climate")
pipeline_stage <- get_cli_value(args, "stage", Sys.getenv("PIPELINE_STAGE", "all"))
pipeline_stage <- match.arg(
  pipeline_stage,
  choices = c("all", "model_fit", "root_reconstruction", "pca_projection", "plot")
)

geiger_ncores <- suppressWarnings(as.integer(Sys.getenv("GEIGER_NCORES", "1")))
if (is.na(geiger_ncores) || geiger_ncores < 1) geiger_ncores <- 1

skip_plots <- truthy(Sys.getenv("SKIP_PLOTS", "false")) || has_cli_flag(args, "skip-plots")
force_refit_models <- truthy(Sys.getenv("FORCE_REFIT_MODELS", "false")) ||
  has_cli_flag(args, "force-refit-models")
resume_model_cache <- truthy(Sys.getenv("RESUME_MODEL_CACHE", "true")) &&
  !has_cli_flag(args, "no-resume-model-cache")
save_intermediate_rds <- truthy(Sys.getenv("SAVE_INTERMEDIATE_RDS", "false")) ||
  has_cli_flag(args, "save-intermediate-rds")
model_fit_variable <- str_trim(get_cli_value(args, "model-fit-variable", Sys.getenv("MODEL_FIT_VARIABLE", "all")))
model_fit_model <- str_trim(get_cli_value(args, "model-fit-model", Sys.getenv("MODEL_FIT_MODEL", "all")))
downstream_model_set <- str_trim(get_cli_value(args, "model-set", Sys.getenv("MODEL_SET", "standard")))
downstream_model_set <- match.arg(
  downstream_model_set,
  choices = c("standard", "5models_no_OU")
)
using_named_downstream_model_set <- !identical(downstream_model_set, "standard")

min_valid_values <- suppressWarnings(as.integer(Sys.getenv("MIN_VALID_VALUES", "10")))
if (is.na(min_valid_values) || min_valid_values < 3) min_valid_values <- 10

primary_delta_warning <- suppressWarnings(as.numeric(Sys.getenv("PRIMARY_DELTA_WARNING", "4")))
primary_delta_severe <- suppressWarnings(as.numeric(Sys.getenv("PRIMARY_DELTA_SEVERE", "10")))
sensitivity_winner_prop_warning <- suppressWarnings(as.numeric(Sys.getenv("SENSITIVITY_WINNER_PROP_WARNING", "0.33")))
ci_nearly_whole_range_ratio <- suppressWarnings(as.numeric(Sys.getenv("CI_NEARLY_WHOLE_RANGE_RATIO", "0.8")))
ci_extremely_wide_ratio <- suppressWarnings(as.numeric(Sys.getenv("CI_EXTREMELY_WIDE_RATIO", "2")))
root_abs_z_warning <- suppressWarnings(as.numeric(Sys.getenv("ROOT_ABS_Z_WARNING", "2.5")))
cov_condition_warning <- suppressWarnings(as.numeric(Sys.getenv("COV_CONDITION_WARNING", "100000000")))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(output_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

pipeline_warnings <- character()
add_pipeline_warning <- function(message_text) {
  pipeline_warnings <<- c(pipeline_warnings, message_text)
  warning(message_text, immediate. = TRUE, call. = FALSE)
}

# -----------------------------
# General helpers
# -----------------------------
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

safe_file_component <- function(x) {
  x |>
    as.character() |>
    str_replace_all("[^A-Za-z0-9_]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")
}

model_set_file_path <- function(stem, ext = "csv", model_set_label = downstream_model_set) {
  file.path(output_dir, paste0(stem, "_", model_set_label, ".", ext))
}

collapse_flags <- function(x) {
  x <- x[!is.na(x) & nzchar(x) & x != "none"]
  if (length(x) == 0) "none" else paste(unique(x), collapse = "; ")
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
  if (!file.exists(file)) {
    stop("Tree file does not exist: ", file)
  }

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

shoelace_area <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  idx_next <- c(seq_along(x)[-1], 1)
  abs(sum(x * y[idx_next] - y * x[idx_next])) / 2
}

point_in_polygon <- function(px, py, poly_x, poly_y) {
  n <- length(poly_x)
  if (n < 3) return(FALSE)

  inside <- FALSE
  j <- n
  for (i in seq_len(n)) {
    intersects <- ((poly_y[i] > py) != (poly_y[j] > py)) &&
      (px < (poly_x[j] - poly_x[i]) * (py - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i])
    if (intersects) inside <- !inside
    j <- i
  }
  inside
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
  out$confidence_level <- level
  out$confidence_label <- paste0(round(level * 100), "%")
  out
}

ellipse_area <- function(cov_mat, level) {
  stats::qchisq(level, df = 2) * pi * sqrt(det(cov_mat))
}

read_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path)
  }
  read_csv(path, show_col_types = FALSE)
}

read_required_lines <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " does not exist: ", path)
  }
  readLines(path, warn = FALSE)
}

extract_parameter_from_summary <- function(parameter_summary, parameter_name) {
  if (!is.character(parameter_summary) || length(parameter_summary) == 0 || is.na(parameter_summary[1])) {
    return(NA_real_)
  }

  pieces <- str_split(parameter_summary[1], ";\\s*")[[1]]
  hit <- pieces[str_detect(pieces, paste0("^", parameter_name, "="))]
  if (length(hit) == 0) {
    return(NA_real_)
  }

  suppressWarnings(as.numeric(sub(paste0("^", parameter_name, "="), "", hit[1])))
}

load_named_downstream_model_set <- function(model_set_label) {
  list(
    fit_attempts = read_required_csv(
      model_set_file_path("per_variable_model_fit_attempts", model_set_label = model_set_label),
      paste0("Named model-fit attempts table for ", model_set_label)
    ),
    model_selection = read_required_csv(
      model_set_file_path("per_variable_model_selection", model_set_label = model_set_label),
      paste0("Named model-selection table for ", model_set_label)
    ),
    model_comparison = read_required_csv(
      model_set_file_path("per_variable_model_comparison", model_set_label = model_set_label),
      paste0("Named model-comparison table for ", model_set_label)
    ),
    grouped_summary = read_required_csv(
      model_set_file_path("grouped_model_summary", model_set_label = model_set_label),
      paste0("Named grouped model summary for ", model_set_label)
    ),
    metadata_lines = read_required_lines(
      model_set_file_path("model_fit_metadata", ext = "txt", model_set_label = model_set_label),
      paste0("Named model-fit metadata for ", model_set_label)
    ),
    biological_summary_lines = read_required_lines(
      model_set_file_path("biological_interpretation_summary", ext = "txt", model_set_label = model_set_label),
      paste0("Named biological interpretation summary for ", model_set_label)
    )
  )
}

classify_ratio_size <- function(ratio) {
  case_when(
    !is.finite(ratio) ~ "unknown",
    ratio < 0.05 ~ "small",
    ratio < 0.25 ~ "moderate",
    TRUE ~ "large"
  )
}

# -----------------------------
# Climate variable definitions
# -----------------------------
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

candidate_models <- c("BM", "OU", "EB", "lambda", "delta", "mean_trend")
primary_models <- c("BM", "lambda", "delta")
sensitivity_models <- setdiff(candidate_models, primary_models)
named_downstream_model_sets <- list(
  "5models_no_OU" = c("BM", "EB", "lambda", "delta", "mean_trend")
)

if (using_named_downstream_model_set && identical(pipeline_stage, "model_fit")) {
  stop(
    "--model-set / MODEL_SET is for downstream root_reconstruction, pca_projection, plot, or all. ",
    "Do not combine it with --stage=model_fit."
  )
}

if (!identical(model_fit_variable, "all")) {
  if (!identical(pipeline_stage, "model_fit")) {
    stop("--model-fit-variable / MODEL_FIT_VARIABLE can only be used with --stage=model_fit")
  }
  if (!model_fit_variable %in% climate_vars) {
    stop(
      "Unknown model-fit variable: ", model_fit_variable,
      ". Expected one of: all, ", paste(climate_vars, collapse = ", ")
    )
  }
}

if (!identical(model_fit_model, "all")) {
  if (!identical(pipeline_stage, "model_fit")) {
    stop("--model-fit-model / MODEL_FIT_MODEL can only be used with --stage=model_fit")
  }
  if (identical(model_fit_variable, "all")) {
    stop("--model-fit-model requires --model-fit-variable to name one climate variable, not 'all'")
  }
  if (!model_fit_model %in% candidate_models) {
    stop(
      "Unknown model-fit model: ", model_fit_model,
      ". Expected one of: all, ", paste(candidate_models, collapse = ", ")
    )
  }
}

model_fit_variables <- if (identical(model_fit_variable, "all")) climate_vars else model_fit_variable
model_fit_models <- if (identical(model_fit_model, "all")) candidate_models else model_fit_model
model_fit_is_shard <- identical(pipeline_stage, "model_fit") &&
  (!identical(model_fit_variable, "all") || !identical(model_fit_model, "all"))

# -----------------------------
# Model fitting helpers
# -----------------------------
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
  if (is.null(fit_obj$opt) || !is.list(fit_obj$opt)) return(NA_character_)
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

  tryCatch(
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

      if (is.null(fit$opt) || !is.list(fit$opt)) {
        stop("geiger::fitContinuous returned a malformed fit object without a list-valued opt slot")
      }

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
}

cache_path_for_variable <- function(variable) {
  file.path(cache_dir, paste0(safe_file_component(variable), "_fitContinuous_cache.rds"))
}

validate_cached_fit <- function(cache_obj, variable, tree, trait) {
  is.list(cache_obj) &&
    identical(cache_obj$variable, variable) &&
    identical(cache_obj$tip_labels, tree$tip.label) &&
    identical(cache_obj$trait_names, names(trait)) &&
    length(cache_obj$model_rows) > 0 &&
    is.list(cache_obj$fit_objects)
}

fit_models_for_variable <- function(tree, trait, variable, models_to_fit = candidate_models) {
  cache_path <- cache_path_for_variable(variable)
  fit_objects <- list()
  fit_rows <- tibble()
  cache_valid <- FALSE

  if (resume_model_cache && !force_refit_models && file.exists(cache_path)) {
    cache_obj <- readRDS(cache_path)
    if (validate_cached_fit(cache_obj, variable, tree, trait)) {
      message("  using cached fitContinuous cache for ", variable)
      fit_objects <- cache_obj$fit_objects
      fit_rows <- cache_obj$model_rows
      cache_valid <- TRUE
    } else {
      add_pipeline_warning(paste0("Ignoring stale or malformed model cache for ", variable))
    }
  }

  existing_models <- if ("model" %in% names(fit_rows)) unique(as.character(fit_rows$model)) else character()

  if (length(models_to_fit) == length(candidate_models) && setequal(models_to_fit, candidate_models)) {
    message("  fitting or loading all candidate models for ", variable)
  } else {
    message("  fitting or loading selected candidate models for ", variable, ": ", paste(models_to_fit, collapse = ", "))
  }

  for (m in models_to_fit) {
    if (cache_valid && !force_refit_models && m %in% existing_models) {
      message("    using cached fitContinuous fit for ", variable, " under ", m)
      next
    }

    message("    ", variable, " under ", m)
    fit_result <- fit_one_model(tree, trait, variable = variable, model = m)
    fit_objects[[m]] <- fit_result$fit
    if ("model" %in% names(fit_rows)) {
      fit_rows <- fit_rows %>% filter(.data$model != m)
    }
    fit_rows <- bind_rows(fit_rows, fit_result$row) %>%
      mutate(model = factor(.data$model, levels = candidate_models)) %>%
      arrange(.data$model) %>%
      mutate(model = as.character(.data$model))

    out <- list(
      variable = variable,
      tip_labels = tree$tip.label,
      trait_names = names(trait),
      model_rows = fit_rows,
      fit_objects = fit_objects
    )
    saveRDS(out, cache_path)
    existing_models <- unique(as.character(fit_rows$model))
    message("    saved fitContinuous cache for ", variable, " after ", m)
  }

  list(
    variable = variable,
    tip_labels = tree$tip.label,
    trait_names = names(trait),
    model_rows = fit_rows,
    fit_objects = fit_objects
  )
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

select_best_row <- function(rows) {
  rows %>%
    filter(is.finite(AICc)) %>%
    arrange(AICc, model) %>%
    slice(1)
}

summarize_model_selection <- function(model_fit_table) {
  model_fit_table %>%
    group_by(variable) %>%
    group_modify(function(.x, .y) {
      best_overall <- select_best_row(.x)
      best_primary <- select_best_row(.x %>% filter(model %in% primary_models))

      tibble(
        best_overall_model = if (nrow(best_overall) == 1) best_overall$model else NA_character_,
        best_overall_AICc = if (nrow(best_overall) == 1) best_overall$AICc else NA_real_,
        primary_model = if (nrow(best_primary) == 1) best_primary$model else NA_character_,
        primary_AICc = if (nrow(best_primary) == 1) best_primary$AICc else NA_real_,
        primary_delta_from_best = primary_AICc - best_overall_AICc
      )
    }) %>%
    ungroup()
}

# -----------------------------
# Primary and sensitivity root reconstruction
# -----------------------------
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
  if (is.null(ci) || !is.matrix(ci) || ncol(ci) < 2) {
    return(c(NA_real_, NA_real_))
  }
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

transform_tree_for_primary_model <- function(tree, fit_obj, primary_model, parameter_summary = NA_character_) {
  if (primary_model == "BM") {
    return(list(tree = tree, transform_parameter = NA_real_, transform_status = "not_applicable"))
  }

  if (primary_model == "lambda") {
    lambda <- extract_fit_parameter(fit_obj, c("lambda"))
    if (!is.finite(lambda)) {
      lambda <- extract_parameter_from_summary(parameter_summary, "lambda")
    }
    if (!is.finite(lambda)) stop("lambda model selected, but no finite lambda parameter was found")
    return(list(
      tree = geiger::lambdaTree(tree, lambda = lambda),
      transform_parameter = lambda,
      transform_status = "lambdaTree"
    ))
  }

  if (primary_model == "delta") {
    delta <- extract_fit_parameter(fit_obj, c("delta"))
    if (!is.finite(delta)) {
      delta <- extract_parameter_from_summary(parameter_summary, "delta")
    }
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

reconstruct_primary_root <- function(tree, trait, variable, primary_model, fit_obj = NULL, parameter_summary = NA_character_) {
  tryCatch(
    {
      tree_info <- transform_tree_for_primary_model(
        tree,
        fit_obj = fit_obj,
        primary_model = primary_model,
        parameter_summary = parameter_summary
      )
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
        reconstruction_model = primary_model,
        transform_parameter = tree_info$transform_parameter,
        transform_status = tree_info$transform_status,
        root_estimate = extract_root_value(ace_fit, recon_tree),
        root_variance = root_variance,
        root_variance_source = if_else(
          is.finite(root_variance),
          "CI95_width_assuming_normality",
          "unavailable"
        ),
        CI_low = ci[1],
        CI_high = ci[2],
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
        reconstruction_model = primary_model,
        transform_parameter = NA_real_,
        transform_status = "failed",
        root_estimate = NA_real_,
        root_variance = NA_real_,
        root_variance_source = "unavailable",
        CI_low = NA_real_,
        CI_high = NA_real_,
        reconstruction_status = "failed",
        reconstruction_error = conditionMessage(e)
      )
    }
  )
}

lookup_model_fit_row <- function(model_fit_table, target_variable, target_model) {
  hits <- model_fit_table %>%
    filter(.data$variable == .env$target_variable, .data$model == .env$target_model) %>%
    slice(1)

  if (nrow(hits) == 0) {
    return(NULL)
  }

  hits
}

compute_sensitivity_root <- function(tree, trait, variable, sensitivity_model) {
  if (!sensitivity_model %in% sensitivity_models) {
    return(tibble(
      variable = variable,
      sensitivity_model = NA_character_,
      sensitivity_root_estimate = NA_real_,
      sensitivity_status = "not_needed_best_overall_is_primary",
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

reconstruct_named_model_root <- function(tree, trait, variable, best_overall_model, primary_model, model_fit_table) {
  if (best_overall_model %in% primary_models) {
    fit_row <- lookup_model_fit_row(model_fit_table, variable, best_overall_model)
    if (is.null(fit_row)) {
      return(tibble(
        variable = variable,
        reconstruction_model = best_overall_model,
        transform_parameter = NA_real_,
        transform_status = "failed",
        root_estimate = NA_real_,
        root_variance = NA_real_,
        root_variance_source = "unavailable",
        CI_low = NA_real_,
        CI_high = NA_real_,
        reconstruction_status = "failed",
        reconstruction_error = paste0("Missing fit row for ", variable, "/", best_overall_model, " in named model-fit attempts table."),
        reconstruction_fallback_reason = NA_character_
      ))
    }

    return(
      reconstruct_primary_root(
        tree = tree,
        trait = trait,
        variable = variable,
        primary_model = best_overall_model,
        parameter_summary = fit_row$parameter_summary[[1]]
      ) %>%
        mutate(reconstruction_fallback_reason = NA_character_)
    )
  }

  if (best_overall_model %in% c("EB", "mean_trend")) {
    sensitivity_row <- compute_sensitivity_root(tree, trait, variable, best_overall_model)
    if (is.finite(sensitivity_row$sensitivity_root_estimate[[1]])) {
      fit_row <- lookup_model_fit_row(model_fit_table, variable, best_overall_model)
      fitted_parameter <- NA_real_
      if (!is.null(fit_row)) {
        fitted_parameter <- extract_parameter_from_summary(
          fit_row$parameter_summary[[1]],
          if (best_overall_model == "EB") "a" else "drift"
        )
      }

      return(tibble(
        variable = variable,
        reconstruction_model = best_overall_model,
        transform_parameter = fitted_parameter,
        transform_status = "direct_named_model_sensitivity_reconstruction",
        root_estimate = sensitivity_row$sensitivity_root_estimate[[1]],
        root_variance = NA_real_,
        root_variance_source = "unavailable_for_direct_named_model_sensitivity_reconstruction",
        CI_low = NA_real_,
        CI_high = NA_real_,
        reconstruction_status = if_else(
          identical(sensitivity_row$sensitivity_status[[1]], "ok"),
          "ok_root_variance_unavailable",
          paste0("ok_root_variance_unavailable_", sensitivity_row$sensitivity_status[[1]])
        ),
        reconstruction_error = sensitivity_row$sensitivity_error[[1]],
        reconstruction_fallback_reason = NA_character_
      ))
    }
  }

  fallback_reason <- paste0(
    "Best model ", best_overall_model,
    " was not directly usable for downstream root reconstruction; falling back to reconstructable primary model ",
    primary_model,
    "."
  )

  if (!is.na(primary_model) && primary_model %in% primary_models) {
    fit_row <- lookup_model_fit_row(model_fit_table, variable, primary_model)
    if (!is.null(fit_row)) {
      return(
        reconstruct_primary_root(
          tree = tree,
          trait = trait,
          variable = variable,
          primary_model = primary_model,
          parameter_summary = fit_row$parameter_summary[[1]]
        ) %>%
          mutate(reconstruction_fallback_reason = fallback_reason)
      )
    }
  }

  tibble(
    variable = variable,
    reconstruction_model = NA_character_,
    transform_parameter = NA_real_,
    transform_status = "failed",
    root_estimate = NA_real_,
    root_variance = NA_real_,
    root_variance_source = "unavailable",
    CI_low = NA_real_,
    CI_high = NA_real_,
    reconstruction_status = "failed",
    reconstruction_error = paste0(
      "No reconstructable downstream model was available for ",
      variable,
      " after named-set winner ",
      best_overall_model,
      "."
    ),
    reconstruction_fallback_reason = fallback_reason
  )
}

# -----------------------------
# Biological diagnostics
# -----------------------------
calculate_extant_stats <- function(X) {
  bind_rows(lapply(names(X), function(v) {
    x <- X[[v]]
    tibble(
      variable = v,
      n_valid_values = sum(is.finite(x)),
      n_unique_values = length(unique(x[is.finite(x)])),
      extant_mean = mean(x, na.rm = TRUE),
      extant_sd = sd(x, na.rm = TRUE),
      extant_min = min(x, na.rm = TRUE),
      extant_max = max(x, na.rm = TRUE),
      extant_range = extant_max - extant_min
    )
  }))
}

classify_variable_interpretability <- function(
  root_estimate, root_z, root_within_range, ci_width_ratio,
  primary_delta_from_best, best_overall_model, reconstruction_status, reconstruction_model
) {
  flags <- c()
  if (!identical(reconstruction_status, "ok")) flags <- c(flags, "primary_reconstruction_uncertain")
  if (!is.na(best_overall_model) && !is.na(reconstruction_model) && best_overall_model != reconstruction_model) {
    flags <- c(flags, "best_model_not_used_for_reconstruction")
  }
  if (is.finite(primary_delta_from_best) && primary_delta_from_best > primary_delta_warning) {
    flags <- c(flags, "primary_model_AICc_penalty")
  }
  if (!isTRUE(root_within_range)) flags <- c(flags, "root_outside_extant_range")
  if (is.finite(abs(root_z)) && abs(root_z) > root_abs_z_warning) flags <- c(flags, "root_extreme_relative_to_extants")
  if (is.finite(ci_width_ratio) && ci_width_ratio >= ci_nearly_whole_range_ratio) {
    flags <- c(flags, "CI_spans_most_extant_range")
  }
  collapse_flags(flags)
}

calculate_root_diagnostics <- function(root_table) {
  root_table %>%
    mutate(
      CI_width = CI_high - CI_low,
      CI_width_to_extant_range = if_else(extant_range > 0, CI_width / extant_range, NA_real_),
      CI_spans_nearly_whole_extant_range = is.finite(CI_width_to_extant_range) &
        CI_width_to_extant_range >= ci_nearly_whole_range_ratio,
      CI_extremely_wide = is.finite(CI_width_to_extant_range) &
        CI_width_to_extant_range >= ci_extremely_wide_ratio,
      root_lies_within_extant_observed_range = is.finite(root_estimate) &
        root_estimate >= extant_min &
        root_estimate <= extant_max,
      standardized_root_position_relative_to_extant_mean_SD = if_else(
        is.finite(extant_sd) & extant_sd > 0,
        (root_estimate - extant_mean) / extant_sd,
        NA_real_
      ),
      root_position_relative_to_extant_distribution = case_when(
        !is.finite(root_estimate) ~ "unknown_no_root_estimate",
        !root_lies_within_extant_observed_range ~ "outside_extant_distribution",
        abs(standardized_root_position_relative_to_extant_mean_SD) <= 1 ~ "central_within_extant_distribution",
        abs(standardized_root_position_relative_to_extant_mean_SD) <= root_abs_z_warning ~ "marginal_within_extant_distribution",
        TRUE ~ "extreme_within_extant_distribution"
      ),
      ecologically_plausible_rule = case_when(
        !is.finite(root_estimate) ~ "unknown_no_root_estimate",
        root_lies_within_extant_observed_range &
          abs(standardized_root_position_relative_to_extant_mean_SD) <= root_abs_z_warning ~ "plausible_within_extant_distribution",
        root_lies_within_extant_observed_range ~ "caution_extreme_within_extant_range",
        TRUE ~ "caution_outside_extant_observed_range"
      ),
      interpretability_flag = mapply(
        classify_variable_interpretability,
        root_estimate,
        standardized_root_position_relative_to_extant_mean_SD,
        root_lies_within_extant_observed_range,
        CI_width_to_extant_range,
        primary_delta_from_best,
        best_overall_model,
        reconstruction_status,
        reconstruction_model
      ),
      variable_interpretation_class = case_when(
        interpretability_flag == "none" ~ "robust",
        str_detect(interpretability_flag, "primary_model_AICc_penalty|root_outside|CI_spans|best_model_not_used_for_reconstruction") ~ "fragile",
        TRUE ~ "tentative"
      ),
      contributes_to_main_root_estimate = is.finite(root_estimate) & reconstruction_status %in% c("ok", "ok_root_variance_unavailable"),
      contributes_to_uncertainty_ellipse = is.finite(root_variance) & root_variance > 0
    )
}

coherence_level <- function(root_table, variable) {
  z <- root_table$standardized_root_position_relative_to_extant_mean_SD[root_table$variable == variable]
  if (length(z) == 0 || !is.finite(z)) return("unknown")
  if (z >= 1.5) return("high")
  if (z <= -1.5) return("low")
  "near_extant_mean"
}

cross_variable_coherence_diagnostics <- function(root_table) {
  level <- function(v) coherence_level(root_table, v)

  rows <- list(
    tibble(
      rule_id = "thermal_mean_vs_cold_quarter",
      variables_involved = "annual_mean_temperature; mean_temperature_coldest_quarter; temperature_annual_range",
      screen_result = if (
        level("annual_mean_temperature") == "high" &&
          level("mean_temperature_coldest_quarter") == "low" &&
          level("temperature_annual_range") != "high"
      ) {
        "potentially_contradictory"
      } else {
        "no_strong_contradiction_detected"
      },
      explanation = "Warm annual mean but cold quarter without high annual range is biologically difficult to interpret."
    ),
    tibble(
      rule_id = "thermal_max_vs_min",
      variables_involved = "maximum_temperature_warmest_month; minimum_temperature_coldest_month; temperature_annual_range",
      screen_result = if (
        level("maximum_temperature_warmest_month") == "high" &&
          level("minimum_temperature_coldest_month") == "low" &&
          level("temperature_annual_range") != "high"
      ) {
        "potentially_contradictory"
      } else {
        "no_strong_contradiction_detected"
      },
      explanation = "Hot warmest month and cold coldest month should usually coincide with high annual range."
    ),
    tibble(
      rule_id = "temperature_seasonality_vs_range",
      variables_involved = "temperature_seasonality; temperature_annual_range",
      screen_result = if (
        (level("temperature_seasonality") == "high" && level("temperature_annual_range") == "low") ||
          (level("temperature_seasonality") == "low" && level("temperature_annual_range") == "high")
      ) {
        "potentially_contradictory"
      } else {
        "no_strong_contradiction_detected"
      },
      explanation = "Temperature seasonality and annual range are expected to broadly agree."
    ),
    tibble(
      rule_id = "annual_precip_vs_wettest_quarter",
      variables_involved = "annual_precipitation; precipitation_of_wettest_quarter; precipitation_wettest_month",
      screen_result = if (
        level("annual_precipitation") == "high" &&
          level("precipitation_of_wettest_quarter") == "low" &&
          level("precipitation_wettest_month") == "low"
      ) {
        "potentially_contradictory"
      } else if (
        level("annual_precipitation") == "low" &&
          (level("precipitation_of_wettest_quarter") == "high" || level("precipitation_wettest_month") == "high")
      ) {
        "seasonal_precipitation_caution"
      } else {
        "no_strong_contradiction_detected"
      },
      explanation = "Annual precipitation should broadly agree with wettest-period precipitation unless precipitation is highly seasonal."
    ),
    tibble(
      rule_id = "cold_quarter_precip_duplicate_check",
      variables_involved = "precipitation_in_coldest_quarter; precipitation_of_coldest_quarter",
      screen_result = if (
        level("precipitation_in_coldest_quarter") != "unknown" &&
          level("precipitation_of_coldest_quarter") != "unknown" &&
          level("precipitation_in_coldest_quarter") != level("precipitation_of_coldest_quarter")
      ) {
        "potentially_contradictory"
      } else {
        "no_strong_contradiction_detected"
      },
      explanation = "These two cold-quarter precipitation variables should not strongly disagree if they represent the same biology."
    )
  )

  bind_rows(rows)
}

# -----------------------------
# Visualization-only PCA and projection
# -----------------------------
run_visualization_pca <- function(X, thermal_vars, precip_vars) {
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

  candidate_pcs <- seq_len(min(4, length(eigvals)))
  thermal_strength <- sapply(candidate_pcs, function(i) {
    sum(abs(loadings[loadings$variable %in% thermal_vars, paste0("PC", i)]), na.rm = TRUE)
  })
  precip_strength <- sapply(candidate_pcs, function(i) {
    sum(abs(loadings[loadings$variable %in% precip_vars, paste0("PC", i)]), na.rm = TRUE)
  })

  pc1_idx <- candidate_pcs[which.max(thermal_strength)]
  remaining <- setdiff(candidate_pcs, pc1_idx)
  if (length(remaining) == 0) {
    stop("Could not choose a second PCA visualization axis.")
  }
  pc2_idx <- remaining[which.max(precip_strength[match(remaining, candidate_pcs)])]

  pc1_sign <- sign(sum(loadings[loadings$variable %in% thermal_vars, paste0("PC", pc1_idx)], na.rm = TRUE))
  if (is.na(pc1_sign) || pc1_sign == 0) pc1_sign <- 1
  pc2_sign <- sign(sum(loadings[loadings$variable %in% precip_vars, paste0("PC", pc2_idx)], na.rm = TRUE))
  if (is.na(pc2_sign) || pc2_sign == 0) pc2_sign <- 1

  scores$display_PCx_oriented <- scores[[paste0("PC", pc1_idx)]] * pc1_sign
  scores$display_PCy_oriented <- scores[[paste0("PC", pc2_idx)]] * pc2_sign
  loadings$display_PCx_oriented <- loadings[[paste0("PC", pc1_idx)]] * pc1_sign
  loadings$display_PCy_oriented <- loadings[[paste0("PC", pc2_idx)]] * pc2_sign

  list(
    pca = pca_vis,
    x_center = x_center,
    x_scale = x_scale,
    scores = scores,
    loadings = loadings,
    eigenvalues = eigenvalues,
    display_pc_x = pc1_idx,
    display_pc_y = pc2_idx,
    display_sign_x = pc1_sign,
    display_sign_y = pc2_sign
  )
}

make_positive_definite <- function(cov_mat, label) {
  if (any(!is.finite(cov_mat))) {
    stop(label, " covariance matrix contains non-finite values.")
  }

  eig <- eigen(cov_mat, symmetric = TRUE)$values
  condition_number <- max(abs(eig)) / min(abs(eig))
  used_near_pd <- FALSE

  if (any(eig <= sqrt(.Machine$double.eps)) || !is.finite(condition_number)) {
    add_pipeline_warning(paste0(label, " covariance was near-singular; applying Matrix::nearPD()."))
    cov_mat <- as.matrix(Matrix::nearPD(cov_mat, corr = FALSE)$mat)
    used_near_pd <- TRUE
    eig <- eigen(cov_mat, symmetric = TRUE)$values
    condition_number <- max(abs(eig)) / min(abs(eig))
  }

  if (any(eig <= sqrt(.Machine$double.eps)) || det(cov_mat) <= sqrt(.Machine$double.eps)) {
    stop(label, " covariance remained numerically degenerate after nearPD; refusing to draw root ellipse.")
  }

  if (is.finite(condition_number) && condition_number > cov_condition_warning) {
    add_pipeline_warning(paste0(label, " covariance has high condition number: ", signif(condition_number, 4)))
  }

  list(
    cov = cov_mat,
    eigenvalues = eig,
    condition_number = condition_number,
    used_near_pd = used_near_pd
  )
}

project_root_to_pca <- function(root_projection_table, pca_obj, scores_out) {
  root_projection_vector <- root_projection_table$projection_value
  names(root_projection_vector) <- root_projection_table$variable

  root_scaled <- (root_projection_vector[names(pca_obj$x_center)] - pca_obj$x_center) / pca_obj$x_scale
  root_scores_all <- as.numeric(root_scaled %*% pca_obj$pca$rotation)
  names(root_scores_all) <- colnames(pca_obj$pca$x)

  root_display_x <- unname(root_scores_all[paste0("PC", pca_obj$display_pc_x)]) * pca_obj$display_sign_x
  root_display_y <- unname(root_scores_all[paste0("PC", pca_obj$display_pc_y)]) * pca_obj$display_sign_y

  hull_idx <- chull(scores_out$display_PCx_oriented, scores_out$display_PCy_oriented)
  hull_x <- scores_out$display_PCx_oriented[hull_idx]
  hull_y <- scores_out$display_PCy_oriented[hull_idx]
  extant_hull_area <- shoelace_area(hull_x, hull_y)
  if (!is.finite(extant_hull_area) || extant_hull_area <= sqrt(.Machine$double.eps)) {
    stop("Displayed extant PCA convex hull has zero or near-zero area; root cloud diagnostics are not interpretable.")
  }
  root_inside_hull <- point_in_polygon(root_display_x, root_display_y, hull_x, hull_y)

  cloud_cov <- cov(scores_out[, c("display_PCx_oriented", "display_PCy_oriented")])
  cloud_cov_info <- make_positive_definite(cloud_cov, "Extant PCA cloud")
  centroid <- colMeans(scores_out[, c("display_PCx_oriented", "display_PCy_oriented")])
  root_md2 <- as.numeric(mahalanobis(
    cbind(root_display_x, root_display_y),
    center = centroid,
    cov = cloud_cov_info$cov
  ))
  root_centroid_distance <- sqrt(sum((c(root_display_x, root_display_y) - centroid)^2))

  cloud_status <- case_when(
    !root_inside_hull ~ "outside_extant_climatic_cloud",
    root_md2 <= qchisq(0.50, df = 2) ~ "central",
    root_md2 <= qchisq(0.95, df = 2) ~ "marginal",
    TRUE ~ "outside_extant_climatic_cloud"
  )

  if (!root_inside_hull || root_md2 > qchisq(0.99, df = 2)) {
    add_pipeline_warning("Projected root is outside or far beyond the displayed extant PCA cloud.")
  }

  all_pc_cols <- as.list(root_scores_all)
  names(all_pc_cols) <- paste0("root_", names(root_scores_all))

  root_point <- as_tibble(all_pc_cols) %>%
    mutate(
      selected_display_PCx = paste0("PC", pca_obj$display_pc_x),
      selected_display_PCy = paste0("PC", pca_obj$display_pc_y),
      display_PCx_sign = pca_obj$display_sign_x,
      display_PCy_sign = pca_obj$display_sign_y,
      display_PCx_oriented = root_display_x,
      display_PCy_oriented = root_display_y,
      root_inside_extant_convex_hull_displayed_2D = root_inside_hull,
      extant_centroid_distance_displayed_2D = root_centroid_distance,
      mahalanobis_distance_squared_from_extant_cloud_displayed_2D = root_md2,
      mahalanobis_distance_from_extant_cloud_displayed_2D = sqrt(root_md2),
      extant_cloud_membership_status = cloud_status,
      extant_convex_hull_area_displayed_2D = extant_hull_area
    )

  list(
    root_point = root_point,
    root_center = c(display_PCx_oriented = root_display_x, display_PCy_oriented = root_display_y),
    extant_hull_area = extant_hull_area,
    root_inside_hull = root_inside_hull,
    root_md2 = root_md2,
    root_cloud_status = cloud_status,
    hull = tibble(x = hull_x, y = hull_y)
  )
}

# This function intentionally uses a diagonal covariance matrix in the original
# climate variables. That means the ellipse is an uncertainty projection of
# per-variable ancestral estimates, not a full multivariate ancestral niche.
project_root_uncertainty_to_pca <- function(root_projection_table, pca_obj, root_center, extant_hull_area) {
  root_var_vec <- root_projection_table$variance_used_for_uncertainty
  names(root_var_vec) <- root_projection_table$variable

  if (sum(is.finite(root_var_vec) & root_var_vec > 0) < 2) {
    stop("Fewer than two variables have usable root uncertainty; projected ellipse would be degenerate.")
  }

  Sigma_root_original <- diag(root_var_vec[names(pca_obj$x_center)])
  rownames(Sigma_root_original) <- names(pca_obj$x_center)
  colnames(Sigma_root_original) <- names(pca_obj$x_center)

  S_inv <- diag(1 / pca_obj$x_scale)
  Sigma_root_scaled <- S_inv %*% Sigma_root_original %*% S_inv
  Sigma_root_pca <- t(pca_obj$pca$rotation) %*% Sigma_root_scaled %*% pca_obj$pca$rotation

  sign_mat <- diag(c(pca_obj$display_sign_x, pca_obj$display_sign_y))
  Sigma2 <- Sigma_root_pca[
    c(pca_obj$display_pc_x, pca_obj$display_pc_y),
    c(pca_obj$display_pc_x, pca_obj$display_pc_y),
    drop = FALSE
  ]
  Sigma2 <- sign_mat %*% Sigma2 %*% sign_mat
  cov_info <- make_positive_definite(Sigma2, "Projected root uncertainty")
  Sigma2 <- cov_info$cov

  levels <- c(0.50, 0.80, 0.95)
  ell <- bind_rows(lapply(levels, function(level) {
    area <- ellipse_area(Sigma2, level)
    ratio <- area / extant_hull_area
    make_ellipse(root_center, Sigma2, level = level) %>%
      mutate(
        ellipse_area = area,
        extant_convex_hull_area = extant_hull_area,
        ellipse_area_to_extant_hull_area_ratio = ratio,
        ellipse_area_relative_to_extant_space = classify_ratio_size(ratio),
        covariance_condition_number = cov_info$condition_number,
        covariance_used_nearPD = cov_info$used_near_pd
      )
  }))

  if (any(ell$ellipse_area_to_extant_hull_area_ratio[ell$confidence_level == 0.95] > 1, na.rm = TRUE)) {
    add_pipeline_warning("The 95% projected root ellipse is larger than the extant convex hull in displayed PCA space.")
  }

  list(ellipse = ell, cov = Sigma2, cov_info = cov_info)
}

# -----------------------------
# Output and interpretation helpers
# -----------------------------
build_sensitivity_comparison <- function(root_table) {
  root_table %>%
    filter(best_overall_model %in% sensitivity_models, best_overall_model != reconstruction_model) %>%
    transmute(
      variable,
      best_overall_model,
      primary_model,
      primary_root_estimate = root_estimate,
      sensitivity_root_estimate,
      absolute_difference = abs(primary_root_estimate - sensitivity_root_estimate),
      standardized_difference = if_else(
        is.finite(extant_sd) & extant_sd > 0,
        absolute_difference / extant_sd,
        NA_real_
      ),
      sensitivity_flag = case_when(
        is.na(sensitivity_model) ~ "no_sensitivity_estimate_needed",
        sensitivity_model == "OU" ~ "OU_best_AICc_flag_only_no_forced_root",
        !is.finite(sensitivity_root_estimate) ~ "sensitivity_estimate_unavailable",
        standardized_difference > 1 ~ "large_sensitivity_difference",
        standardized_difference > 0.5 ~ "moderate_sensitivity_difference",
        TRUE ~ "low_sensitivity_difference"
      )
    )
}

classify_final_result <- function(metrics) {
  if (
    metrics$number_of_usable_primary_variables >= 0.8 * length(climate_vars) &&
      metrics$proportion_of_variables_with_caution_flags <= 0.25 &&
      metrics$median_primary_delta_from_best <= 2 &&
      metrics$number_of_variables_where_root_outside_extant_range <= 1 &&
      metrics$extant_cloud_membership %in% c("central", "marginal") &&
      metrics$ellipse_95_area_to_extant_hull_area_ratio <= 0.5
  ) {
    return("ROBUST")
  }

  if (
    metrics$number_of_usable_primary_variables >= 0.6 * length(climate_vars) &&
      metrics$proportion_of_variables_with_caution_flags <= 0.5 &&
      metrics$extant_cloud_membership != "outside_extant_climatic_cloud" &&
      metrics$ellipse_95_area_to_extant_hull_area_ratio <= 1
  ) {
    return("TENTATIVE")
  }

  "WEAK"
}

write_interpretation_summary <- function(
  path, metrics, root_table, coherence_table, sensitivity_table,
  pipeline_warnings
) {
  robust_vars <- root_table$variable[root_table$variable_interpretation_class == "robust"]
  fragile_vars <- root_table$variable[root_table$variable_interpretation_class == "fragile"]
  tentative_vars <- root_table$variable[root_table$variable_interpretation_class == "tentative"]
  contradictions <- coherence_table %>% filter(screen_result %in% c("potentially_contradictory", "seasonal_precipitation_caution"))

  lines <- c(
    "Result Block 1 biological interpretation summary",
    "=================================================",
    "",
    "Question: Where did the radiation begin in climatic space?",
    "",
    "Core method statement:",
    "Ancestral climate was reconstructed in the original environmental variables after per-variable evolutionary model fitting.",
    "PCA was used only to visualize extant climatic structure and the projected ancestral root.",
    "Variables whose best-supported model fell outside the primary reconstructable model set were retained as sensitivity cases and explicitly flagged.",
    "",
    paste0("Final classification: ", metrics$final_classification),
    paste0("Usable primary variables: ", metrics$number_of_usable_primary_variables, " of ", length(climate_vars)),
    paste0("Variables with OU/EB/mean_trend as best overall model: ", metrics$number_of_variables_with_non_primary_best_model),
    paste0("Variables with caution flags: ", metrics$number_of_variables_with_caution_flags),
    paste0("Median primary delta AICc from best model: ", signif(metrics$median_primary_delta_from_best, 4)),
    paste0("Variables where root is outside extant range: ", metrics$number_of_variables_where_root_outside_extant_range),
    paste0("Projected root cloud status: ", metrics$extant_cloud_membership),
    paste0("Projected root Mahalanobis distance squared: ", signif(metrics$root_mahalanobis_distance_squared_in_PCA_space, 4)),
    paste0("95% ellipse / extant convex hull area ratio: ", signif(metrics$ellipse_95_area_to_extant_hull_area_ratio, 4)),
    "",
    "Robust variables:",
    if (length(robust_vars) == 0) "none" else paste(robust_vars, collapse = ", "),
    "",
    "Tentative variables:",
    if (length(tentative_vars) == 0) "none" else paste(tentative_vars, collapse = ", "),
    "",
    "Fragile variables:",
    if (length(fragile_vars) == 0) "none" else paste(fragile_vars, collapse = ", "),
    "",
    "Cross-variable coherence screen:",
    if (nrow(contradictions) == 0) {
      "No strong rule-based contradictions detected."
    } else {
      paste(paste0(contradictions$rule_id, ": ", contradictions$screen_result), collapse = "\n")
    },
    "",
    "Sensitivity cases:",
    if (nrow(sensitivity_table) == 0) {
      "No sensitivity-only model won AICc."
    } else {
      paste(paste0(
        sensitivity_table$variable,
        " best=", sensitivity_table$best_overall_model,
        " flag=", sensitivity_table$sensitivity_flag
      ), collapse = "\n")
    },
    "",
    "Warnings raised during run:",
    if (length(pipeline_warnings) == 0) "none" else paste(unique(pipeline_warnings), collapse = "\n"),
    "",
    "Important limitation:",
    "The projected root ellipse is an uncertainty projection from per-variable estimates using a diagonal covariance approximation.",
    "It should not be described as a fully jointly inferred multivariate ancestral climatic niche."
  )

  writeLines(lines, path)
}

make_plot <- function(scores_out, root_point, ellipse_df, loadings, core_membership, eigenvalues, pca_obj) {
  arrow_mult <- 2.5
  thermal_arrows_df <- loadings %>%
    filter(variable %in% thermal_vars) %>%
    transmute(variable, x = display_PCx_oriented * arrow_mult, y = display_PCy_oriented * arrow_mult)

  precip_arrows_df <- loadings %>%
    filter(variable %in% precip_vars) %>%
    transmute(variable, x = display_PCx_oriented * arrow_mult, y = display_PCy_oriented * arrow_mult)

  label_df <- bind_rows(
    core_membership %>% filter(inside_50_projected_root_ellipse) %>% slice_head(n = 15),
    core_membership %>% slice_head(n = 10)
  ) %>%
    distinct(tree_label, .keep_all = TRUE)

  var_pc1 <- percent(eigenvalues$variance_explained[pca_obj$display_pc_x], accuracy = 0.1)
  var_pc2 <- percent(eigenvalues$variance_explained[pca_obj$display_pc_y], accuracy = 0.1)
  xlab_txt <- paste0("Thermal climatic gradient (PC", pca_obj$display_pc_x, ", ", var_pc1, ")")
  ylab_txt <- paste0("Precipitation climatic gradient (PC", pca_obj$display_pc_y, ", ", var_pc2, ")")

  p <- ggplot(scores_out, aes(x = display_PCx_oriented, y = display_PCy_oriented)) +
    geom_polygon(
      data = ellipse_df,
      aes(x = x, y = y, group = confidence_label, fill = confidence_label),
      alpha = 0.16,
      color = "black",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    geom_point(alpha = 0.28, size = 1.35, color = "grey35") +
    geom_point(
      data = core_membership %>% filter(inside_50_projected_root_ellipse),
      aes(x = display_PCx_oriented, y = display_PCy_oriented),
      inherit.aes = FALSE,
      size = 1.9,
      color = "darkgreen",
      alpha = 0.8
    ) +
    geom_point(
      data = root_point,
      aes(x = display_PCx_oriented, y = display_PCy_oriented),
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
    filename = file.path(output_dir, "ancestral_root_projected_PCA_strict_plot.png"),
    plot = p,
    width = 10,
    height = 8,
    dpi = 400
  )
  ggsave(
    filename = file.path(output_dir, "ancestral_root_projected_PCA_strict_plot.pdf"),
    plot = p,
    width = 10,
    height = 8
  )

  write_csv(thermal_arrows_df, file.path(output_dir, "visualization_thermal_loading_arrows.csv"))
  write_csv(precip_arrows_df, file.path(output_dir, "visualization_precipitation_loading_arrows.csv"))
  write_csv(label_df, file.path(output_dir, "visualization_species_labels_used_in_plot.csv"))
}

# -----------------------------
# Data reading, cleaning, matching
# -----------------------------
if (!file.exists(climate_file)) {
  stop("Climate file does not exist: ", climate_file)
}

tree <- read_tree_robust(tree_file)
write_csv(
  tibble(tree_tip_label_raw = tree$tip.label),
  file.path(output_dir, "input_tree_tip_labels_raw.csv")
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
  file.path(output_dir, "input_climate_species_raw.csv")
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
  write_csv(duplicate_species, file.path(output_dir, "input_duplicate_climate_species_after_normalization.csv"))
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
write_csv(match_summary, file.path(output_dir, "input_name_match_summary.csv"))

if (length(overlap) == 0) {
  stop(
    "Zero overlap between tree tip labels and climate species names after normalization. ",
    "Inspect input_tree_tip_labels_raw.csv and input_climate_species_raw.csv."
  )
}

write_csv(
  clim2 %>% filter(!tree_label %in% tree_tips),
  file.path(output_dir, "input_species_removed_from_csv_not_in_tree.csv")
)

clim_complete <- clim2 %>%
  filter(tree_label %in% tree_tips) %>%
  filter(if_all(all_of(climate_vars), ~ !is.na(.)))

matched_species <- intersect(tree$tip.label, clim_complete$tree_label)
if (length(matched_species) < min_valid_values) {
  stop("Too few species remain after matching tree and complete climate data: ", length(matched_species))
}

tips_to_drop <- setdiff(tree$tip.label, matched_species)
write_csv(
  tibble(tree_tip_missing_complete_climate = tips_to_drop),
  file.path(output_dir, "input_tree_tips_missing_complete_climate.csv")
)

tree_pruned <- if (length(tips_to_drop) > 0) drop.tip(tree, tips_to_drop) else tree
if (is.null(tree_pruned) || !inherits(tree_pruned, "phylo")) {
  stop("Pruned tree is invalid.")
}

clim_final <- clim_complete %>%
  filter(tree_label %in% tree_pruned$tip.label) %>%
  slice(match(tree_pruned$tip.label, tree_label))

if (nrow(clim_final) != length(tree_pruned$tip.label) ||
    !all(clim_final$tree_label == tree_pruned$tip.label)) {
  stop("Ordering mismatch between climate data and tree tips after pruning.")
}

write_csv(clim_final, file.path(output_dir, "input_matched_climate_data_used.csv"))

X <- clim_final %>%
  select(all_of(climate_vars)) %>%
  as.data.frame()
rownames(X) <- clim_final$tree_label

extant_stats <- calculate_extant_stats(X)
bad_variables <- extant_stats %>%
  filter(n_valid_values < min_valid_values | n_unique_values < 4)
if (nrow(bad_variables) > 0) {
  stop(
    "Variables have too few valid or unique values for model fitting: ",
    paste(bad_variables$variable, collapse = ", ")
  )
}

crown_node <- as.character(Ntip(tree_pruned) + 1)

# -----------------------------
# Stage: model fitting
# -----------------------------
fit_cache_objects <- list()
named_downstream_inputs <- NULL

if (using_named_downstream_model_set && !identical(pipeline_stage, "model_fit")) {
  message("Using precomputed downstream model universe: ", downstream_model_set)
  named_downstream_inputs <- load_named_downstream_model_set(downstream_model_set)
  expected_models <- named_downstream_model_sets[[downstream_model_set]]

  model_fit_attempts <- named_downstream_inputs$fit_attempts %>%
    filter(.data$variable %in% climate_vars, .data$model %in% expected_models) %>%
    add_model_deltas()

  model_selection <- named_downstream_inputs$model_selection %>%
    filter(.data$variable %in% climate_vars)

  if (nrow(model_selection) != length(climate_vars)) {
    stop(
      "Named downstream model-selection table for ",
      downstream_model_set,
      " must contain exactly ",
      length(climate_vars),
      " climate variables."
    )
  }

  if (!all(model_selection$variable %in% climate_vars)) {
    stop("Named downstream model-selection table contains unexpected climate variables.")
  }

  if (any(is.na(model_selection$best_overall_model))) {
    stop(
      "Named downstream model-selection table has missing best_overall_model entries for: ",
      paste(model_selection$variable[is.na(model_selection$best_overall_model)], collapse = ", ")
    )
  }

  if (any(!model_selection$best_overall_model %in% expected_models)) {
    stop(
      "Named downstream model-selection table contains models outside ",
      downstream_model_set,
      ": ",
      paste(unique(model_selection$best_overall_model[!model_selection$best_overall_model %in% expected_models]), collapse = ", ")
    )
  }
} else {
  message("Fitting or loading per-variable geiger::fitContinuous model fits.")
  message("Model-fit variable scope: ", paste(model_fit_variables, collapse = ", "))
  message("Model-fit model scope: ", paste(model_fit_models, collapse = ", "))

  for (v in model_fit_variables) {
    trait <- X[[v]]
    names(trait) <- rownames(X)
    trait <- trait[tree_pruned$tip.label]
    fit_cache_objects[[v]] <- fit_models_for_variable(
      tree_pruned,
      trait,
      v,
      models_to_fit = model_fit_models
    )
  }

  model_fit_attempts <- bind_rows(lapply(fit_cache_objects, `[[`, "model_rows")) %>%
    add_model_deltas()

  if (model_fit_is_shard) {
    write_csv(model_fit_attempts, file.path(output_dir, "model_fit_shard_attempts.csv"))
    if (setequal(model_fit_models, candidate_models)) {
      write_csv(
        summarize_model_selection(model_fit_attempts),
        file.path(output_dir, "model_fit_shard_selection.csv")
      )
    }

    completion_lines <- c(
      "stage=model_fit",
      "mode=shard",
      paste0("model_fit_variable=", model_fit_variable),
      paste0("model_fit_model=", model_fit_model),
      paste0("variables_completed=", paste(model_fit_variables, collapse = ",")),
      paste0("models_completed_or_loaded=", paste(model_fit_models, collapse = ",")),
      paste0("cache_files=", paste(vapply(model_fit_variables, cache_path_for_variable, character(1)), collapse = ",")),
      paste0("attempt_rows=", nrow(model_fit_attempts)),
      paste0("completed_at_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE))
    )
    writeLines(completion_lines, file.path(output_dir, "model_fit_shard_completion.txt"))
    message("Stage model_fit shard complete. Stopping before full model-fit aggregation.")
    quit(save = "no", status = 0)
  }

  write_csv(model_fit_attempts, file.path(output_dir, "per_variable_model_fit_attempts_detailed.csv"))

  if (any(!is.finite(model_fit_attempts$AICc))) {
    failed_fit_count <- sum(!is.finite(model_fit_attempts$AICc))
    add_pipeline_warning(paste0(failed_fit_count, " model fits did not return finite AICc; see per_variable_model_fit_attempts_detailed.csv."))
  }

  model_selection <- summarize_model_selection(model_fit_attempts)
  if (any(is.na(model_selection$best_overall_model))) {
    stop(
      "At least one variable has no finite AICc model fit: ",
      paste(model_selection$variable[is.na(model_selection$best_overall_model)], collapse = ", ")
    )
  }
  if (any(is.na(model_selection$primary_model))) {
    stop(
      "At least one variable has no finite primary model fit among BM/lambda/delta: ",
      paste(model_selection$variable[is.na(model_selection$primary_model)], collapse = ", ")
    )
  }
}

if (pipeline_stage == "model_fit") {
  write_csv(model_selection, file.path(output_dir, "per_variable_model_selection_preliminary.csv"))
  if (save_intermediate_rds) {
    saveRDS(fit_cache_objects, file.path(cache_dir, "all_fitContinuous_cache_objects.rds"))
  }
  message("Stage model_fit complete. Stopping before root reconstruction.")
  quit(save = "no", status = 0)
}

# -----------------------------
# Stage: primary root reconstruction and sensitivity
# -----------------------------
reuse_named_root_reconstruction <- using_named_downstream_model_set &&
  pipeline_stage %in% c("pca_projection", "plot") &&
  file.exists(file.path(output_dir, "ancestral_root_original_climate_strict.csv"))

if (reuse_named_root_reconstruction) {
  message("Reusing previously written named-set root reconstruction outputs.")
  ancestral_root_original <- read_required_csv(
    file.path(output_dir, "ancestral_root_original_climate_strict.csv"),
    "Existing root reconstruction table"
  )
  named_root_output <- if (file.exists(model_set_file_path("ancestral_root_original_climate"))) {
    read_required_csv(
      model_set_file_path("ancestral_root_original_climate"),
      paste0("Existing named root output table for ", downstream_model_set)
    )
  } else {
    NULL
  }
  per_variable_model_fits <- if (file.exists(file.path(output_dir, "per_variable_model_fits.csv"))) {
    read_required_csv(file.path(output_dir, "per_variable_model_fits.csv"), "Existing per-variable model fits table")
  } else {
    ancestral_root_original %>%
      transmute(
        variable,
        best_overall_model,
        best_overall_AICc,
        primary_model,
        primary_AICc,
        primary_delta_from_best,
        caution_flag = interpretability_flag,
        contributes_to_main_root_estimate,
        contributes_to_uncertainty_ellipse
      )
  }
  sensitivity_comparison <- if (file.exists(file.path(output_dir, "sensitivity_comparison_summary.csv"))) {
    read_required_csv(file.path(output_dir, "sensitivity_comparison_summary.csv"), "Existing sensitivity comparison table")
  } else {
    build_sensitivity_comparison(ancestral_root_original)
  }
  coherence_diagnostics <- if (file.exists(file.path(output_dir, "cross_variable_climate_coherence_diagnostics.csv"))) {
    read_required_csv(file.path(output_dir, "cross_variable_climate_coherence_diagnostics.csv"), "Existing coherence diagnostics table")
  } else {
    cross_variable_coherence_diagnostics(ancestral_root_original)
  }
  root_projection_table <- if (file.exists(file.path(output_dir, "ancestral_root_vector_used_for_PCA_projection_strict.csv"))) {
    read_required_csv(
      file.path(output_dir, "ancestral_root_vector_used_for_PCA_projection_strict.csv"),
      "Existing root projection table"
    )
  } else {
    ancestral_root_original %>%
      mutate(
        projection_value = if_else(is.finite(root_estimate), root_estimate, extant_mean),
        projection_value_source = if_else(
          is.finite(root_estimate),
          "primary_reconstruction",
          "extant_mean_fallback_for_projection_only"
        ),
        variance_used_for_uncertainty = if_else(
          is.finite(root_variance) & root_variance > 0,
          root_variance,
          0
        ),
        uncertainty_source = if_else(
          is.finite(root_variance) & root_variance > 0,
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
  }
} else {
  message("Reconstructing primary root states in original climate variables.")

  root_rows <- list()
  for (v in climate_vars) {
    trait <- X[[v]]
    names(trait) <- rownames(X)
    trait <- trait[tree_pruned$tip.label]

    selected <- model_selection %>% filter(variable == v)
    if (using_named_downstream_model_set) {
      root_rows[[length(root_rows) + 1]] <- reconstruct_named_model_root(
        tree = tree_pruned,
        trait = trait,
        variable = v,
        best_overall_model = selected$best_overall_model[[1]],
        primary_model = selected$primary_model[[1]],
        model_fit_table = model_fit_attempts
      )
    } else {
      root_rows[[length(root_rows) + 1]] <- reconstruct_primary_root(
        tree = tree_pruned,
        trait = trait,
        variable = v,
        primary_model = selected$primary_model,
        fit_obj = fit_cache_objects[[v]]$fit_objects[[selected$primary_model]]
      ) %>%
        mutate(reconstruction_fallback_reason = NA_character_)
    }
  }

  primary_root <- bind_rows(root_rows)

  sensitivity_rows <- list()
  for (v in climate_vars) {
    trait <- X[[v]]
    names(trait) <- rownames(X)
    trait <- trait[tree_pruned$tip.label]
    selected <- model_selection %>% filter(variable == v)
    sensitivity_rows[[length(sensitivity_rows) + 1]] <- compute_sensitivity_root(
      tree_pruned,
      trait,
      variable = v,
      sensitivity_model = selected$best_overall_model
    )
  }
  sensitivity_root <- bind_rows(sensitivity_rows)

  ancestral_root_original <- model_selection %>%
    left_join(primary_root, by = "variable") %>%
    left_join(sensitivity_root, by = "variable") %>%
    left_join(extant_stats, by = "variable") %>%
    calculate_root_diagnostics()

  named_root_output <- NULL
  if (using_named_downstream_model_set) {
    named_root_output <- ancestral_root_original %>%
      transmute(
        variable,
        best_model_5models_no_OU = best_overall_model,
        reconstruction_model_used = reconstruction_model,
        winning_model_used_for_reconstruction = best_overall_model == reconstruction_model,
        reconstruction_fallback_reason = if_else(
          is.na(reconstruction_fallback_reason) | reconstruction_fallback_reason == "",
          "none",
          reconstruction_fallback_reason
        ),
        root_estimate,
        CI_low,
        CI_high,
        extant_mean,
        extant_sd,
        extant_min,
        extant_max,
        root_position_relative_to_extant_distribution,
        caution_flag = interpretability_flag
      )
  }

  per_variable_model_fits <- ancestral_root_original %>%
    transmute(
      variable,
      best_overall_model,
      best_overall_AICc,
      primary_model,
      primary_AICc,
      primary_delta_from_best,
      caution_flag = interpretability_flag,
      contributes_to_main_root_estimate,
      contributes_to_uncertainty_ellipse
    )

  write_csv(per_variable_model_fits, file.path(output_dir, "per_variable_model_fits.csv"))
  write_csv(ancestral_root_original, file.path(output_dir, "ancestral_root_original_climate_strict.csv"))
  if (using_named_downstream_model_set) {
    write_csv(named_root_output, model_set_file_path("ancestral_root_original_climate"))
    write_csv(
      named_root_output %>% filter(!winning_model_used_for_reconstruction),
      model_set_file_path("reconstruction_model_differences")
    )
  }

  sensitivity_comparison <- build_sensitivity_comparison(ancestral_root_original)
  write_csv(sensitivity_comparison, file.path(output_dir, "sensitivity_comparison_summary.csv"))

  coherence_diagnostics <- cross_variable_coherence_diagnostics(ancestral_root_original)
  write_csv(coherence_diagnostics, file.path(output_dir, "cross_variable_climate_coherence_diagnostics.csv"))

  non_primary_best_prop <- mean(ancestral_root_original$best_overall_model %in% sensitivity_models, na.rm = TRUE)
  if (non_primary_best_prop > sensitivity_winner_prop_warning) {
    add_pipeline_warning(paste0(
      "High proportion of variables have sensitivity-only best models: ",
      signif(non_primary_best_prop, 3)
    ))
  }

  fallback_prop <- mean(
    ancestral_root_original$best_overall_model != ancestral_root_original$reconstruction_model,
    na.rm = TRUE
  )
  if (using_named_downstream_model_set && is.finite(fallback_prop) && fallback_prop > 0.25) {
    add_pipeline_warning(paste0(
      "More than 25% of variables required reconstruction fallback in the ",
      downstream_model_set,
      " downstream universe: ",
      signif(fallback_prop, 3)
    ))
  }

  bad_delta <- ancestral_root_original %>%
    filter(is.finite(primary_delta_from_best), primary_delta_from_best > primary_delta_warning)
  if (nrow(bad_delta) > 0) {
    add_pipeline_warning(paste0(
      "Primary model is worse than best overall model by > ",
      primary_delta_warning,
      " AICc for: ",
      paste(bad_delta$variable, collapse = ", ")
    ))
  }

  wide_ci <- ancestral_root_original %>%
    filter(CI_extremely_wide)
  if (nrow(wide_ci) > 0) {
    add_pipeline_warning(paste0(
      "Extremely wide root CIs relative to extant ranges for: ",
      paste(wide_ci$variable, collapse = ", ")
    ))
  }

  root_projection_table <- ancestral_root_original %>%
    mutate(
      projection_value = if_else(is.finite(root_estimate), root_estimate, extant_mean),
      projection_value_source = if_else(
        is.finite(root_estimate),
        "primary_reconstruction",
        "extant_mean_fallback_for_projection_only"
      ),
      variance_used_for_uncertainty = if_else(
        is.finite(root_variance) & root_variance > 0,
        root_variance,
        0
      ),
      uncertainty_source = if_else(
        is.finite(root_variance) & root_variance > 0,
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
  write_csv(root_projection_table, file.path(output_dir, "ancestral_root_vector_used_for_PCA_projection_strict.csv"))

  if (save_intermediate_rds) {
    saveRDS(
      list(
        model_fit_attempts = model_fit_attempts,
        model_selection = model_selection,
        ancestral_root_original = ancestral_root_original,
        root_projection_table = root_projection_table,
        sensitivity_comparison = sensitivity_comparison,
        coherence_diagnostics = coherence_diagnostics
      ),
      file.path(cache_dir, "root_reconstruction_intermediate.rds")
    )
  }
}

if (pipeline_stage == "root_reconstruction") {
  message("Stage root_reconstruction complete. Stopping before PCA projection.")
  quit(save = "no", status = 0)
}

# -----------------------------
# Stage: visualization-only PCA projection
# -----------------------------
message("Running visualization-only PCA and projecting reconstructed root.")
pca_obj <- run_visualization_pca(X, thermal_vars, precip_vars)

scores_out <- pca_obj$scores %>%
  left_join(
    clim_final %>% select(tree_label, species_original, species_clean),
    by = "tree_label"
  ) %>%
  relocate(species_original, species_clean, tree_label, display_PCx_oriented, display_PCy_oriented)

projection <- project_root_to_pca(root_projection_table, pca_obj, scores_out)
projected_root_point <- projection$root_point %>%
  mutate(
    crown_node = crown_node,
    n_species_used = nrow(scores_out),
    tree_is_rooted_after_pruning = is.rooted(tree_pruned),
    projection_note = "root reconstructed in original variables; PCA used only for visualization",
    .before = 1
  )

uncertainty_projection <- project_root_uncertainty_to_pca(
  root_projection_table,
  pca_obj,
  root_center = projection$root_center,
  extant_hull_area = projection$extant_hull_area
)
ellipse_df <- uncertainty_projection$ellipse

core_md2 <- as.numeric(mahalanobis(
  scores_out[, c("display_PCx_oriented", "display_PCy_oriented")],
  center = projection$root_center,
  cov = uncertainty_projection$cov
))

core_membership <- scores_out %>%
  transmute(
    species_original,
    species_clean,
    tree_label,
    display_PCx_oriented,
    display_PCy_oriented,
    mahalanobis_d2_to_projected_root_uncertainty = core_md2,
    inside_50_projected_root_ellipse = core_md2 <= qchisq(0.50, df = 2),
    inside_80_projected_root_ellipse = core_md2 <= qchisq(0.80, df = 2),
    inside_95_projected_root_ellipse = core_md2 <= qchisq(0.95, df = 2),
    closest_projected_root_ellipse_level = case_when(
      inside_50_projected_root_ellipse ~ "50%",
      inside_80_projected_root_ellipse ~ "80%",
      inside_95_projected_root_ellipse ~ "95%",
      TRUE ~ "outside_95%"
    )
  ) %>%
  arrange(mahalanobis_d2_to_projected_root_uncertainty)

ellipse_95_ratio <- ellipse_df$ellipse_area_to_extant_hull_area_ratio[ellipse_df$confidence_level == 0.95][1]

interpretability_metrics <- tibble(
  number_of_usable_primary_variables = sum(ancestral_root_original$contributes_to_main_root_estimate),
  proportion_of_variables_with_caution_flags = mean(ancestral_root_original$interpretability_flag != "none"),
  number_of_variables_with_caution_flags = sum(ancestral_root_original$interpretability_flag != "none"),
  median_primary_delta_from_best = median(ancestral_root_original$primary_delta_from_best, na.rm = TRUE),
  number_of_variables_with_non_primary_best_model = sum(ancestral_root_original$best_overall_model %in% sensitivity_models),
  number_of_variables_where_reconstruction_differs_from_best_model = sum(
    ancestral_root_original$best_overall_model != ancestral_root_original$reconstruction_model,
    na.rm = TRUE
  ),
  number_of_variables_where_root_outside_extant_range = sum(!ancestral_root_original$root_lies_within_extant_observed_range, na.rm = TRUE),
  root_mahalanobis_distance_squared_in_PCA_space = projection$root_md2,
  root_mahalanobis_distance_in_PCA_space = sqrt(projection$root_md2),
  extant_cloud_membership = projection$root_cloud_status,
  root_inside_extant_convex_hull_displayed_2D = projection$root_inside_hull,
  ellipse_95_area_to_extant_hull_area_ratio = ellipse_95_ratio,
  covariance_condition_number = uncertainty_projection$cov_info$condition_number,
  covariance_used_nearPD = uncertainty_projection$cov_info$used_near_pd,
  model_universe = if_else(using_named_downstream_model_set, downstream_model_set, "standard_full_universe"),
  OU_status = if_else(using_named_downstream_model_set, "excluded_deferred", "included_if_fit")
)
interpretability_metrics$final_classification <- classify_final_result(interpretability_metrics)

write_csv(pca_obj$scores, file.path(output_dir, "visualization_only_PCA_species_scores.csv"))
write_csv(pca_obj$loadings, file.path(output_dir, "visualization_only_PCA_loadings.csv"))
write_csv(pca_obj$eigenvalues, file.path(output_dir, "visualization_only_PCA_eigenvalues.csv"))
write_csv(projected_root_point, file.path(output_dir, "ancestral_root_projected_PCA_point_strict.csv"))
write_csv(ellipse_df, file.path(output_dir, "ancestral_root_projected_PCA_ellipses_strict.csv"))
write_csv(core_membership, file.path(output_dir, "species_membership_in_projected_root_ellipses_strict.csv"))
write_csv(interpretability_metrics, file.path(output_dir, "ancestral_root_interpretability_metrics.csv"))
write_csv(
  tibble(warning = if (length(pipeline_warnings) == 0) "none" else unique(pipeline_warnings)),
  file.path(output_dir, "pipeline_warnings_and_safeguards.csv")
)

write_interpretation_summary(
  file.path(output_dir, "biological_interpretation_summary.txt"),
  interpretability_metrics,
  ancestral_root_original,
  coherence_diagnostics,
  sensitivity_comparison,
  pipeline_warnings
)

if (using_named_downstream_model_set) {
  write_csv(projected_root_point, model_set_file_path("ancestral_root_projected_PCA_point"))
  write_csv(ellipse_df, model_set_file_path("ancestral_root_projected_PCA_ellipses"))
  write_csv(interpretability_metrics, model_set_file_path("ancestral_root_interpretability_metrics"))

  named_fallback_vars <- ancestral_root_original$variable[
    ancestral_root_original$best_overall_model != ancestral_root_original$reconstruction_model
  ]
  robust_root_vars <- ancestral_root_original$variable[ancestral_root_original$interpretability_flag == "none"]
  fragile_root_vars <- ancestral_root_original$variable[ancestral_root_original$interpretability_flag != "none"]
  thermal_positions <- ancestral_root_original %>%
    filter(.data$variable %in% thermal_vars) %>%
    transmute(line = paste0(.data$variable, "=", .data$root_position_relative_to_extant_distribution)) %>%
    pull(line)
  precipitation_positions <- ancestral_root_original %>%
    filter(.data$variable %in% precip_vars) %>%
    transmute(line = paste0(.data$variable, "=", .data$root_position_relative_to_extant_distribution)) %>%
    pull(line)

  writeLines(
    c(
      "Result Block 1 root interpretation summary (5models_no_OU)",
      "========================================================",
      "",
      "Model universe:",
      "BM, EB, lambda, delta, and mean_trend.",
      "OU was excluded/deferred because full-tree OU fitting remained computationally intractable under the available compute setup.",
      "",
      "Core method statement:",
      "Ancestral climate was reconstructed in the original climate variables.",
      "PCA was used only to project the reconstructed root for visualization and diagnostics.",
      "",
      paste0("Projected root status in displayed PCA space: ", projection$root_cloud_status),
      paste0("Projected root inside displayed extant convex hull: ", projection$root_inside_hull),
      paste0("Projected root Mahalanobis distance squared: ", signif(projection$root_md2, 4)),
      paste0("95% ellipse / extant convex hull area ratio: ", signif(ellipse_95_ratio, 4)),
      paste0("Variables where reconstruction differed from the winning model: ", interpretability_metrics$number_of_variables_where_reconstruction_differs_from_best_model),
      "",
      "Thermal variable root positions relative to extant distributions:",
      if (length(thermal_positions) == 0) "none" else paste(thermal_positions, collapse = "\n"),
      "",
      "Precipitation variable root positions relative to extant distributions:",
      if (length(precipitation_positions) == 0) "none" else paste(precipitation_positions, collapse = "\n"),
      "",
      "Variables with robust root estimates:",
      if (length(robust_root_vars) == 0) "none" else paste(robust_root_vars, collapse = ", "),
      "",
      "Variables with fragile or flagged root estimates:",
      if (length(fragile_root_vars) == 0) "none" else paste(fragile_root_vars, collapse = ", "),
      "",
      "Variables requiring fallback reconstruction:",
      if (length(named_fallback_vars) == 0) "none" else paste(named_fallback_vars, collapse = ", "),
      "",
      "Warnings raised during this downstream run:",
      if (length(pipeline_warnings) == 0) "none" else paste(unique(pipeline_warnings), collapse = "\n"),
      "",
      "Important limitation:",
      "This is a five-model/no-OU Result Block 1 analysis, not the final six-model universe.",
      "Do not describe PCA as the inferential basis of ancestral reconstruction."
    ),
    model_set_file_path("biological_interpretation_summary_root", ext = "txt")
  )

  readme_lines <- c(
    "Result Block 1 downstream readme (5models_no_OU)",
    "===============================================",
    "",
    "Model universe used:",
    "BM, EB, lambda, delta, mean_trend.",
    "OU status: excluded/deferred, not unsupported.",
    "",
    "Why OU is absent here:",
    "Full-tree OU fitting remained computationally intractable under the available hosted-runner/local compute setup.",
    "The completed five-model/no-OU model-fit aggregation is therefore the explicit downstream model universe for this run.",
    "",
    "How to interpret these outputs:",
    "1. ancestral_root_original_climate_5models_no_OU.csv is the per-variable ancestral reconstruction in original climate variables.",
    "2. ancestral_root_projected_PCA_point_5models_no_OU.csv and ancestral_root_projected_PCA_ellipses_5models_no_OU.csv are visualization/diagnostic projections only.",
    "3. ancestral_climatic_core_plot_5models_no_OU.* is a figure derived from the projected root, not the inferential basis of reconstruction.",
    "4. biological_interpretation_summary_root_5models_no_OU.txt is the manuscript-facing summary for this five-model/no-OU Result Block 1 run.",
    "",
    "What should not be claimed:",
    "- Do not describe this as a six-model result.",
    "- Do not describe PCA as ancestral inference.",
    "- Do not treat OU as rejected; it was excluded/deferred for operational reasons.",
    "",
    "Upstream model-fit provenance:",
    named_downstream_inputs$metadata_lines,
    "",
    "Upstream five-model biological summary:",
    named_downstream_inputs$biological_summary_lines
  )
  writeLines(readme_lines, model_set_file_path("result_block_1_readme", ext = "txt"))
}

writeLines(
  c(
    "# Result Block 1 Strict Ancestral Climate Output Summary",
    "",
    "Core statement supported by this pipeline:",
    "\"Ancestral climate was reconstructed in the original environmental variables after per-variable evolutionary model fitting. PCA was used only to visualize extant climatic structure and the projected ancestral root. Variables whose best-supported model fell outside the primary reconstructable model set were retained as sensitivity cases and explicitly flagged.\"",
    "",
    "Read outputs in this order:",
    "1. per_variable_model_fits.csv",
    "2. ancestral_root_original_climate_strict.csv",
    "3. sensitivity_comparison_summary.csv",
    "4. cross_variable_climate_coherence_diagnostics.csv",
    "5. ancestral_root_projected_PCA_point_strict.csv",
    "6. ancestral_root_projected_PCA_ellipses_strict.csv",
    "7. ancestral_root_interpretability_metrics.csv",
    "8. biological_interpretation_summary.txt",
    "",
    "Important limitation:",
    "The projected root ellipse uses diagonal per-variable root uncertainty and should not be described as a fully jointly inferred multivariate ancestral climatic niche.",
    "",
    paste0("Final classification: ", interpretability_metrics$final_classification)
  ),
  con = file.path(output_dir, "SUMMARY.md")
)

write_xlsx(
  list(
    interpretability_metrics = interpretability_metrics,
    per_variable_model_fits = per_variable_model_fits,
    model_fit_attempts = model_fit_attempts,
    ancestral_root_original = ancestral_root_original,
    sensitivity_comparison = sensitivity_comparison,
    coherence_diagnostics = coherence_diagnostics,
    root_projection_vector = root_projection_table,
    projected_root_PCA_point = projected_root_point,
    projected_root_PCA_ellipses = ellipse_df,
    species_root_ellipse_membership = core_membership,
    visualization_PCA_scores = pca_obj$scores,
    visualization_PCA_loadings = pca_obj$loadings,
    visualization_PCA_eigenvalues = pca_obj$eigenvalues,
    data_used = clim_final
  ),
  path = file.path(output_dir, "result_block_1_strict_ancestral_climate_outputs.xlsx")
)

if (save_intermediate_rds) {
  saveRDS(
    list(
      pca_obj = pca_obj,
      projection = projection,
      uncertainty_projection = uncertainty_projection,
      interpretability_metrics = interpretability_metrics
    ),
    file.path(cache_dir, "pca_projection_intermediate.rds")
  )
}

if (!skip_plots && pipeline_stage %in% c("all", "plot")) {
  make_plot(
    scores_out = scores_out,
    root_point = projected_root_point,
    ellipse_df = ellipse_df,
    loadings = pca_obj$loadings,
    core_membership = core_membership,
    eigenvalues = pca_obj$eigenvalues,
    pca_obj = pca_obj
  )

  if (using_named_downstream_model_set) {
    file.copy(
      file.path(output_dir, "ancestral_root_projected_PCA_strict_plot.png"),
      model_set_file_path("ancestral_climatic_core_plot", ext = "png"),
      overwrite = TRUE
    )
    file.copy(
      file.path(output_dir, "ancestral_root_projected_PCA_strict_plot.pdf"),
      model_set_file_path("ancestral_climatic_core_plot", ext = "pdf"),
      overwrite = TRUE
    )
  }
} else {
  message("Plot generation skipped.")
}

cat("\n==============================\n")
cat("Result Block 1 strict workflow complete\n")
cat("==============================\n")
print(interpretability_metrics)
cat("\nFiles written to: ", output_dir, "\n", sep = "")
