# ============================================================
# Ancestral crown climate reconstruction and climatic core test
# ============================================================

req_pkgs <- c(
  "ape", "phytools", "dplyr", "readr", "stringr", "ggplot2",
  "ggrepel", "tibble", "writexl", "scales", "Matrix"
)

to_install <- req_pkgs[!req_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(ape)
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

tree_file <- "data/raw/pruned_phylogeny.nex"
climate_file <- "data/raw/climate.csv"
output_dir <- "results/ancestral_climate_reconstruction"

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

crown_node <- as.character(Ntip(tree_pruned) + 1)

crown_recon_list <- lapply(climate_vars, function(v) {
  trait_vec <- clim_final[[v]]
  names(trait_vec) <- clim_final$tree_label

  anc <- fastAnc(tree_pruned, trait_vec, vars = TRUE, CI = TRUE)

  tibble(
    variable = v,
    crown_node = crown_node,
    crown_estimate = as.numeric(anc$ace[crown_node]),
    crown_variance = as.numeric(anc$var[crown_node]),
    CI_low = as.numeric(anc$CI95[crown_node, 1]),
    CI_high = as.numeric(anc$CI95[crown_node, 2])
  )
})

crown_climate_df <- bind_rows(crown_recon_list)

crown_climate <- crown_climate_df$crown_estimate
names(crown_climate) <- crown_climate_df$variable

crown_var_vec <- crown_climate_df$crown_variance
names(crown_var_vec) <- crown_climate_df$variable

if (any(!is.finite(crown_climate)) || any(!is.finite(crown_var_vec))) {
  stop("Non-finite crown estimates or variances were returned by fastAnc().")
}

crown_scaled <- (crown_climate[names(x_center)] - x_center) / x_scale
crown_scores_all <- as.numeric(crown_scaled %*% pca_vis$rotation)
names(crown_scores_all) <- colnames(pca_vis$x)

crown_center <- c(
  PC1_thermal = crown_scores_all[paste0("PC", pc1_idx)] * pc1_sign,
  PC2_precipitation = crown_scores_all[paste0("PC", pc2_idx)] * pc2_sign
)

crown_point <- tibble(
  crown_node = crown_node,
  tree_is_rooted_after_pruning = is.rooted(tree_pruned),
  n_species_used = length(tree_pruned$tip.label),
  PC1_thermal = as.numeric(crown_center["PC1_thermal"]),
  PC2_precipitation = as.numeric(crown_center["PC2_precipitation"])
)

Sigma_crown_original <- diag(crown_var_vec[names(x_center)])
rownames(Sigma_crown_original) <- names(x_center)
colnames(Sigma_crown_original) <- names(x_center)

S_inv <- diag(1 / x_scale)
Sigma_crown_scaled <- S_inv %*% Sigma_crown_original %*% S_inv
Sigma_crown_pca <- t(pca_vis$rotation) %*% Sigma_crown_scaled %*% pca_vis$rotation

sign_mat <- diag(c(pc1_sign, pc2_sign))
Sigma2 <- Sigma_crown_pca[c(pc1_idx, pc2_idx), c(pc1_idx, pc2_idx), drop = FALSE]
Sigma2 <- sign_mat %*% Sigma2 %*% sign_mat

if (any(!is.finite(Sigma2))) {
  stop("Non-finite values in projected crown covariance matrix.")
}

eig_check <- eigen(Sigma2, symmetric = TRUE)$values
if (any(eig_check <= 0)) {
  Sigma2 <- as.matrix(Matrix::nearPD(Sigma2)$mat)
}

ell <- bind_rows(
  make_ellipse(crown_center, Sigma2, level = 0.50),
  make_ellipse(crown_center, Sigma2, level = 0.80),
  make_ellipse(crown_center, Sigma2, level = 0.95)
)

scores_out <- scores %>%
  left_join(
    clim_final %>% select(tree_label, species_original, species_clean),
    by = "tree_label"
  ) %>%
  relocate(species_original, species_clean, tree_label, PC1_thermal, PC2_precipitation)

d2 <- mahalanobis(
  x = scores_out[, c("PC1_thermal", "PC2_precipitation")],
  center = as.numeric(crown_center),
  cov = Sigma2
)

core_membership <- scores_out %>%
  transmute(
    species_original,
    species_clean,
    tree_label,
    PC1_thermal,
    PC2_precipitation,
    mahalanobis_d2 = as.numeric(d2),
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
  n_species_used = nrow(core_membership),
  crown_node = crown_node,
  tree_is_rooted = is.rooted(tree_pruned),
  selected_thermal_axis = paste0("PC", pc1_idx),
  selected_precipitation_axis = paste0("PC", pc2_idx),
  thermal_axis_variance_explained = eigenvalues$variance_explained[pc1_idx],
  precipitation_axis_variance_explained = eigenvalues$variance_explained[pc2_idx],
  n_species_inside_50_core = sum(core_membership$inside_50_core),
  n_species_inside_80_core = sum(core_membership$inside_80_core),
  n_species_inside_95_core = sum(core_membership$inside_95_core),
  pct_species_inside_50_core = mean(core_membership$inside_50_core),
  pct_species_inside_80_core = mean(core_membership$inside_80_core),
  pct_species_inside_95_core = mean(core_membership$inside_95_core)
)

write_csv(scores_out, file.path(output_dir, "PCA_species_scores_visualization.csv"))
write_csv(loadings, file.path(output_dir, "PCA_loadings_visualization.csv"))
write_csv(eigenvalues, file.path(output_dir, "PCA_eigenvalues_visualization.csv"))
write_csv(crown_climate_df, file.path(output_dir, "crown_ancestral_original_climate.csv"))
write_csv(crown_point, file.path(output_dir, "crown_ancestral_projected_PCA_point.csv"))
write_csv(ell, file.path(output_dir, "crown_ancestral_projected_PCA_ellipses.csv"))
write_csv(core_membership, file.path(output_dir, "ancestral_core_species_membership.csv"))
write_csv(core_summary, file.path(output_dir, "ancestral_core_summary.csv"))

write_xlsx(
  list(
    summary = core_summary,
    crown_original_climate = crown_climate_df,
    crown_projected_PCA_point = crown_point,
    species_core_membership = core_membership,
    species_scores_visualization = scores_out,
    PCA_loadings_visualization = loadings,
    PCA_eigenvalues_visualization = eigenvalues,
    data_used = clim_final
  ),
  path = file.path(output_dir, "ancestral_climate_reconstruction_outputs.xlsx")
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
    data = crown_point,
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
    fill = "Crown root CI",
    title = "Ancestral climatic core of the crown group",
    subtitle = "Crown climate reconstructed on original variables, then projected into PCA space"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(output_dir, "ancestral_climatic_core.png"),
  plot = p_core,
  width = 10,
  height = 8,
  dpi = 400
)

ggsave(
  filename = file.path(output_dir, "ancestral_climatic_core.pdf"),
  plot = p_core,
  width = 10,
  height = 8
)

write_csv(thermal_arrows_df, file.path(output_dir, "thermal_loading_arrows.csv"))
write_csv(precip_arrows_df, file.path(output_dir, "precipitation_loading_arrows.csv"))
write_csv(label_df, file.path(output_dir, "species_labels_used_in_plot.csv"))

summary_lines <- c(
  "# Ancestral Climate Reconstruction Summary",
  "",
  paste0("Species used: ", core_summary$n_species_used),
  paste0("Crown node: ", crown_node),
  paste0("Tree rooted after pruning: ", core_summary$tree_is_rooted),
  paste0(
    "Selected PCA axes: thermal = ", core_summary$selected_thermal_axis,
    ", precipitation = ", core_summary$selected_precipitation_axis
  ),
  paste0("Species inside 50% ancestral core: ", core_summary$n_species_inside_50_core),
  paste0("Species inside 80% ancestral core: ", core_summary$n_species_inside_80_core),
  paste0("Species inside 95% ancestral core: ", core_summary$n_species_inside_95_core),
  "",
  "The crown/root climate estimates are in `crown_ancestral_original_climate.csv`.",
  "Species are ranked by climatic closeness to the projected crown root in `ancestral_core_species_membership.csv`."
)

writeLines(summary_lines, con = file.path(output_dir, "SUMMARY.md"))

cat("\n==============================\n")
cat("Crown ancestral climate summary\n")
cat("==============================\n")
print(crown_climate_df)

cat("\n==============================\n")
cat("Ancestral core membership summary\n")
cat("==============================\n")
print(core_summary)

cat("\nFiles written to: ", output_dir, "\n", sep = "")
