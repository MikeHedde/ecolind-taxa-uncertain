# ============================================================================
# 01_run_multitaxa_uncertainty_2024_v3_taxref_coherent_confusion.R
# Empirical multi-taxon taxonomic-uncertainty pipeline for RMQS 2024
#
# Scope
#   - alpha diversity: Hill q0, q1, q2 and gamma richness
#   - beta diversity: Bray-Curtis and Sørensen (Jaccard = optional sensitivity)
#   - ecological inference: GDM on turnover versus climate, soil and geography
#
# Design principles
#   1. Each assemblage = taxon group x protocol and is analysed independently.
#   2. Each scenario is compared with its appropriate taxonomic baseline.
#   3. Family-level coarsening is skipped when an assemblage is already one family
#      (e.g. Carabidae, Formicidae).
#   4. Blowes alpha-gamma-occupancy diagnostics are deliberately NOT included
#      here; they will be added afterwards from the harmonised outputs.
#
# Prerequisite
#   Run 00_prepare_multitaxa_inputs_2024.R first.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(vegan)
  library(gdm)
  library(sf)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# 0. USER SETTINGS
# -----------------------------------------------------------------------------

INPUT_DIR <- "outputs_multitaxa_2024"
ENV_FILE  <- "data/derived-data/all_env_variables.csv"
OUT_DIR   <- "outputs_multitaxa_2024/uncertainty_results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "per_assemblage"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

# Groups retained in the first full run. Secondary assemblages remain included
# by default but can be disabled if GDM diagnostics are too sparse.
MAIN_ASSEMBLAGES <- c(
  "collembola_soil_core",
  "araneae_pitfall",
  "formicidae_pitfall",
  "carabidae_pitfall",
  "isopoda_pitfall"
)
SECONDARY_ASSEMBLAGES <- c(
  "isopoda_hand_sorting",
  "diplopoda_pitfall",
  "diplopoda_hand_sorting"
)
RUN_SECONDARY_ASSEMBLAGES <- TRUE

# Stochastic ID-error scenarios.
# The main paper uses one transparent reference intensity (10%). The full
# 1–20% response curves are generated separately by
# 02_run_multitaxa_error_gradients_appendix.R.
N_SIM <- 30L
GDM_MAX_STOCHASTIC_ITERS <- 12L
MAIN_ERROR_RATE <- 0.10
RARE_WEIGHTED_ERROR <- 0.10
EXPERT_ERROR_RATE <- 0.10
RARE_MAX_TOTAL_ABUNDANCE <- 3L
RARE_WEIGHTED_MAX_ERROR <- 0.25

# Alpha-beta conventions.
RUN_JACCARD_SENSITIVITY <- TRUE # exported, but not used in main figures
MIN_BETA_SITES <- 8L
MIN_GDM_SITES <- 20L

# Environmental predictors retained a priori for all assemblages.
ENV_PREDICTORS <- c("t360_mean", "mos", "p_h")

# Regional TAXREF pools and optional expert maps are built by:
#   05_build_multitaxa_taxref_pools.R
#
# Regional pool required columns:
#   taxon_unit, genus
# Expert map required columns:
#   source_taxon_unit, target_taxon_unit, enabled
# Optional expert-map columns: weight, comment.
REGIONAL_POOL_DIR <- "regional_pools"
EXPERT_MAP_DIR <- "expert_confusion_maps"
RUN_REGIONAL_POOL_SCENARIOS <- TRUE
RUN_EXPERT_CONFUSION_SCENARIOS <- TRUE
# Store one auditable source→target map per assemblage × scenario by default.
# Full map sets can be produced by setting WRITE_ALL_COHERENT_CONFUSION_MAPS = TRUE.
WRITE_EXAMPLE_COHERENT_CONFUSION_MAPS <- TRUE
WRITE_ALL_COHERENT_CONFUSION_MAPS <- FALSE

set.seed(123)

# -----------------------------------------------------------------------------
# 1. VALIDATION AND HELPERS
# -----------------------------------------------------------------------------

required_files <- c(
  file.path(INPUT_DIR, "assemblage_manifest.csv"),
  file.path(INPUT_DIR, "assemblages")
)
if (!all(file.exists(required_files))) {
  stop(
    "Input files not found. Run 00_prepare_multitaxa_inputs_2024.R first.\n",
    "Expected: ", paste(required_files, collapse = ", ")
  )
}
if (!file.exists(ENV_FILE)) {
  stop("Environmental file not found: ", ENV_FILE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || dplyr::n_distinct(x[keep]) < 2L || dplyr::n_distinct(y[keep]) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

safe_rank_cor <- function(x, y) {
  safe_cor(x, y, method = "spearman")
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA else x[[1]]
}

make_matrix <- function(comm_long, station_frame) {
  stations <- as.character(station_frame$station)
  taxa <- sort(unique(comm_long$taxon_unit))
  
  if (length(taxa) == 0L) {
    return(matrix(0, nrow = length(stations), ncol = 0,
                  dimnames = list(stations, character(0))))
  }
  
  comm_long %>%
    filter(!is.na(station), !is.na(taxon_unit), abundance > 0) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    tidyr::complete(station = stations, taxon_unit = taxa, fill = list(abundance = 0)) %>%
    pivot_wider(names_from = taxon_unit, values_from = abundance, values_fill = 0) %>%
    arrange(match(station, stations)) %>%
    tibble::column_to_rownames("station") %>%
    as.matrix()
}

alpha_metrics <- function(mat) {
  n <- nrow(mat)
  if (n == 0L || ncol(mat) == 0L) {
    return(tibble(
      station = rownames(mat), total_abundance = 0,
      q0 = 0, q1 = 0, q2 = 0
    ))
  }
  
  total <- rowSums(mat)
  q0 <- rowSums(mat > 0)
  q1 <- exp(vegan::diversity(mat, index = "shannon"))
  q2 <- vegan::diversity(mat, index = "invsimpson")
  q1[!is.finite(q1)] <- 0
  q2[!is.finite(q2)] <- 0
  
  tibble(
    station = rownames(mat),
    total_abundance = total,
    q0 = q0,
    q1 = q1,
    q2 = q2
  )
}

community_gamma <- function(mat) {
  if (ncol(mat) == 0L) return(0L)
  sum(colSums(mat) > 0)
}

prepare_beta_matrix <- function(mat) {
  if (nrow(mat) == 0L || ncol(mat) == 0L) return(NULL)
  keep_sites <- rowSums(mat) > 0
  mat <- mat[keep_sites, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (nrow(mat) < MIN_BETA_SITES || ncol(mat) < 2L) return(NULL)
  mat
}

compute_distance <- function(mat, metric = c("bray", "sorensen", "jaccard")) {
  metric <- match.arg(metric)
  mat <- prepare_beta_matrix(mat)
  if (is.null(mat)) return(NULL)
  
  if (metric == "bray") return(vegan::vegdist(mat, method = "bray"))
  if (metric == "sorensen") return(vegan::vegdist(mat > 0, method = "bray", binary = TRUE))
  vegan::vegdist(mat > 0, method = "jaccard", binary = TRUE)
}

procrustes_r2 <- function(d1, d2) {
  # cmdscale(add = TRUE) returns a list in some R versions. Extract the
  # coordinates explicitly; treating that list as a matrix made all
  # Procrustes values NA in the previous pipeline.
  get_points <- function(d) {
    out <- try(stats::cmdscale(d, k = 2, eig = TRUE, add = TRUE), silent = TRUE)
    if (inherits(out, "try-error")) return(NULL)
    pts <- if (is.list(out) && !is.null(out$points)) out$points else out
    pts <- try(as.matrix(pts), silent = TRUE)
    if (inherits(pts, "try-error") || nrow(pts) < MIN_BETA_SITES || ncol(pts) < 2L) {
      return(NULL)
    }
    if (any(!is.finite(pts[, 1:2, drop = FALSE]))) return(NULL)
    pts[, 1:2, drop = FALSE]
  }
  
  if (is.null(d1) || is.null(d2) || length(d1) != length(d2)) return(NA_real_)
  sites_1 <- attr(d1, "Labels")
  sites_2 <- attr(d2, "Labels")
  if (is.null(sites_1) || is.null(sites_2) || !setequal(sites_1, sites_2)) return(NA_real_)
  
  p1 <- get_points(d1)
  p2 <- get_points(d2)
  if (is.null(p1) || is.null(p2)) return(NA_real_)
  
  common <- intersect(rownames(p1), rownames(p2))
  if (length(common) < MIN_BETA_SITES) return(NA_real_)
  p1 <- p1[common, , drop = FALSE]
  p2 <- p2[common, , drop = FALSE]
  
  fit <- try(vegan::procrustes(p1, p2, symmetric = TRUE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit$ss) || !is.finite(fit$ss)) return(NA_real_)
  max(0, 1 - fit$ss)
}

compare_beta <- function(mat_base, mat_scen) {
  sites <- intersect(
    rownames(mat_base)[rowSums(mat_base) > 0],
    rownames(mat_scen)[rowSums(mat_scen) > 0]
  )
  if (length(sites) < MIN_BETA_SITES) {
    return(tibble(
      n_common_sites = length(sites),
      bray_stability = NA_real_, sorensen_stability = NA_real_, jaccard_stability = NA_real_,
      bray_mean_change_pct = NA_real_, sorensen_mean_change_pct = NA_real_, jaccard_mean_change_pct = NA_real_,
      ordination_procrustes_r2 = NA_real_
    ))
  }
  
  base <- mat_base[sites, , drop = FALSE]
  scen <- mat_scen[sites, , drop = FALSE]
  
  db <- compute_distance(base, "bray")
  ds <- compute_distance(scen, "bray")
  sb <- compute_distance(base, "sorensen")
  ss <- compute_distance(scen, "sorensen")
  
  jb <- if (RUN_JACCARD_SENSITIVITY) compute_distance(base, "jaccard") else NULL
  js <- if (RUN_JACCARD_SENSITIVITY) compute_distance(scen, "jaccard") else NULL
  
  tibble(
    n_common_sites = length(sites),
    bray_stability = if (!is.null(db) && !is.null(ds)) safe_cor(as.vector(db), as.vector(ds)) else NA_real_,
    sorensen_stability = if (!is.null(sb) && !is.null(ss)) safe_cor(as.vector(sb), as.vector(ss)) else NA_real_,
    jaccard_stability = if (!is.null(jb) && !is.null(js)) safe_cor(as.vector(jb), as.vector(js)) else NA_real_,
    bray_mean_change_pct = if (!is.null(db) && !is.null(ds) && mean(db) > 0) 100 * (mean(ds) / mean(db) - 1) else NA_real_,
    sorensen_mean_change_pct = if (!is.null(sb) && !is.null(ss) && mean(sb) > 0) 100 * (mean(ss) / mean(sb) - 1) else NA_real_,
    jaccard_mean_change_pct = if (!is.null(jb) && !is.null(js) && mean(jb) > 0) 100 * (mean(js) / mean(jb) - 1) else NA_real_,
    ordination_procrustes_r2 = procrustes_r2(db, ds)
  )
}

# -----------------------------------------------------------------------------
# 2. ENVIRONMENTAL METADATA FOR GDM
# -----------------------------------------------------------------------------

env_raw <- readr::read_csv(ENV_FILE, show_col_types = FALSE, name_repair = "unique") %>%
  janitor::clean_names()

required_env <- c("station", "longitude", "latitude", ENV_PREDICTORS)
missing_env <- setdiff(required_env, names(env_raw))
if (length(missing_env) > 0L) {
  stop("Environmental file is missing: ", paste(missing_env, collapse = ", "))
}

env_station <- env_raw %>%
  mutate(
    station = as.character(station),
    across(all_of(c("longitude", "latitude", ENV_PREDICTORS)), as.numeric)
  ) %>%
  group_by(station) %>%
  summarise(
    longitude = first_non_missing(longitude),
    latitude = first_non_missing(latitude),
    across(all_of(ENV_PREDICTORS), first_non_missing),
    .groups = "drop"
  ) %>%
  filter(complete.cases(across(all_of(c("longitude", "latitude", ENV_PREDICTORS)))))

# GDM calculates geographic distances from X/Y coordinates. Lambert-93 ensures
# that the geographic predictor is expressed in a meaningful projected space.
env_sf <- sf::st_as_sf(env_station, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
  sf::st_transform(2154)
coords_l93 <- sf::st_coordinates(env_sf)
env_station <- env_station %>%
  mutate(x_l93 = coords_l93[, 1], y_l93 = coords_l93[, 2])

readr::write_csv(env_station, file.path(OUT_DIR, "environmental_metadata_for_gdm.csv"))

# -----------------------------------------------------------------------------
# 3. GDM HELPERS
# -----------------------------------------------------------------------------

make_gdm_sitepair <- function(mat, env_meta, metric = c("bray", "sorensen")) {
  metric <- match.arg(metric)
  
  common <- intersect(rownames(mat), env_meta$station)
  mat <- mat[common, , drop = FALSE]
  meta <- env_meta %>%
    filter(station %in% common) %>%
    arrange(match(station, common))
  mat <- mat[meta$station, , drop = FALSE]
  
  # GDM is fitted only on sites with biological records and complete predictors.
  keep_sites <- rowSums(mat) > 0
  mat <- mat[keep_sites, , drop = FALSE]
  meta <- meta[keep_sites, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  
  if (nrow(mat) < MIN_GDM_SITES || ncol(mat) < 2L) return(NULL)
  if (anyDuplicated(meta$station) > 0L || any(!complete.cases(meta))) return(NULL)
  
  # For bioFormat = 1, use plain base data.frames. The gdm input formatter can
  # fail with tibble subclasses under some package versions. Coordinates are
  # included explicitly in BOTH tables for an unambiguous site-pair specification.
  bio_mat <- if (metric == "bray") mat else (mat > 0) * 1L
  bio_tbl <- data.frame(
    station = as.character(rownames(bio_mat)),
    x_l93 = as.numeric(meta$x_l93),
    y_l93 = as.numeric(meta$y_l93),
    as.data.frame(bio_mat, check.names = FALSE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  pred_tbl <- as.data.frame(
    meta[, c("station", "x_l93", "y_l93", ENV_PREDICTORS), drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  pred_tbl$station <- as.character(pred_tbl$station)
  for (nm in c("x_l93", "y_l93", ENV_PREDICTORS)) {
    pred_tbl[[nm]] <- as.numeric(pred_tbl[[nm]])
  }
  
  if (!identical(bio_tbl$station, pred_tbl$station) ||
      anyDuplicated(bio_tbl$station) > 0L ||
      anyNA(pred_tbl)) return(NULL)
  
  sp <- tryCatch(
    gdm::formatsitepair(
      bioData = bio_tbl,
      bioFormat = 1,
      dist = "bray",
      abundance = identical(metric, "bray"),
      siteColumn = "station",
      XColumn = "x_l93",
      YColumn = "y_l93",
      predData = pred_tbl,
      verbose = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(sp, "error") || is.null(sp) || nrow(sp) < 1L) return(NULL)
  
  # Guard against accidental non-finite values in the site-pair table.
  numeric_sp <- as.data.frame(sp)[vapply(as.data.frame(sp), is.numeric, logical(1))]
  if (ncol(numeric_sp) == 0L || any(!is.finite(as.matrix(numeric_sp)))) return(NULL)
  
  list(sitepair = sp, meta = meta, n_sites = nrow(mat))
}

fit_gdm_model <- function(mat, env_meta, metric = c("bray", "sorensen")) {
  metric <- match.arg(metric)
  prep <- make_gdm_sitepair(mat, env_meta, metric)
  if (is.null(prep)) return(NULL)
  
  model <- tryCatch(
    gdm::gdm(prep$sitepair, geo = TRUE),
    error = function(e) e
  )
  
  if (inherits(model, "error") || is.null(model) ||
      is.null(model$explained) || !is.finite(model$explained)) {
    return(NULL)
  }
  
  list(model = model, sitepair = prep$sitepair, n_sites = prep$n_sites)
}

# Canonical schemas prevent empty GDM runs from dropping column names.
empty_gdm_weight_table <- function() {
  tibble(
    predictor = character(),
    ispline_weight = numeric()
  )
}

empty_gdm_comparison_terms <- function() {
  tibble(
    predictor = character(),
    baseline_ispline_weight = numeric(),
    scenario_ispline_weight = numeric(),
    delta_ispline_weight = numeric()
  )
}

empty_gdm_terms_summary <- function() {
  tibble(
    assemblage_id = character(),
    scenario = character(),
    scenario_family = character(),
    baseline_scenario = character(),
    unit_type = character(),
    distance = character(),
    predictor = character(),
    baseline_ispline_weight_median = numeric(),
    scenario_ispline_weight_median = numeric(),
    delta_ispline_weight_median = numeric(),
    n_iter = integer()
  )
}

extract_gdm_terms <- function(model) {
  # gdm stores coefficients and spline counts as vectors/lists depending on
  # package version. Convert explicitly and return a typed empty table whenever
  # a term-level decomposition cannot be reconstructed safely.
  coefficients <- suppressWarnings(as.numeric(unlist(model$coefficients, use.names = FALSE)))
  spline_counts <- suppressWarnings(as.integer(unlist(model$splines, use.names = FALSE)))
  predictors <- as.character(unlist(model$predictors, use.names = FALSE))
  
  if (length(coefficients) == 0L || length(spline_counts) == 0L ||
      any(!is.finite(coefficients)) || any(is.na(spline_counts)) ||
      any(spline_counts <= 0L) || sum(spline_counts) != length(coefficients)) {
    return(empty_gdm_weight_table())
  }
  
  # In gdm, geographic distance may be included in the spline vector but
  # omitted from model$predictors. Handle both object conventions explicitly.
  if (length(predictors) == length(spline_counts)) {
    term_names <- predictors
  } else if (isTRUE(model$geo) && length(predictors) + 1L == length(spline_counts)) {
    term_names <- c("Geographic", predictors)
  } else {
    term_names <- paste0("term_", seq_along(spline_counts))
  }
  
  if (length(term_names) != length(spline_counts)) {
    return(empty_gdm_weight_table())
  }
  
  tibble(
    predictor = rep(term_names, times = spline_counts),
    coefficient = coefficients
  ) %>%
    group_by(predictor) %>%
    summarise(ispline_weight = sum(abs(coefficient)), .groups = "drop")
}

compare_gdm <- function(mat_base, mat_scen, env_meta, metric = c("bray", "sorensen")) {
  metric <- match.arg(metric)
  
  common_sites <- intersect(
    rownames(mat_base)[rowSums(mat_base) > 0],
    rownames(mat_scen)[rowSums(mat_scen) > 0]
  )
  common_sites <- intersect(common_sites, env_meta$station)
  if (length(common_sites) < MIN_GDM_SITES) {
    return(list(summary = tibble(
      gdm_n_sites = length(common_sites),
      gdm_explained_baseline = NA_real_,
      gdm_explained_scenario = NA_real_,
      gdm_delta_explained = NA_real_,
      gdm_predicted_stability = NA_real_,
      gdm_term_rank_stability = NA_real_
    ), terms = tibble()))
  }
  
  base_fit <- fit_gdm_model(mat_base[common_sites, , drop = FALSE], env_meta, metric)
  scen_fit <- fit_gdm_model(mat_scen[common_sites, , drop = FALSE], env_meta, metric)
  if (is.null(base_fit) || is.null(scen_fit)) {
    return(list(summary = tibble(
      gdm_n_sites = length(common_sites),
      gdm_explained_baseline = NA_real_,
      gdm_explained_scenario = NA_real_,
      gdm_delta_explained = NA_real_,
      gdm_predicted_stability = NA_real_,
      gdm_term_rank_stability = NA_real_
    ), terms = tibble()))
  }
  
  tb <- extract_gdm_terms(base_fit$model) %>%
    rename(baseline_ispline_weight = ispline_weight)
  ts <- extract_gdm_terms(scen_fit$model) %>%
    rename(scenario_ispline_weight = ispline_weight)
  
  terms <- if (nrow(tb) == 0L && nrow(ts) == 0L) {
    empty_gdm_comparison_terms()
  } else {
    full_join(tb, ts, by = "predictor") %>%
      mutate(
        delta_ispline_weight = scenario_ispline_weight - baseline_ispline_weight
      )
  }
  
  rank_stability <- safe_rank_cor(terms$baseline_ispline_weight, terms$scenario_ispline_weight)
  
  list(
    summary = tibble(
      gdm_n_sites = min(base_fit$n_sites, scen_fit$n_sites),
      gdm_explained_baseline = base_fit$model$explained,
      gdm_explained_scenario = scen_fit$model$explained,
      gdm_delta_explained = scen_fit$model$explained - base_fit$model$explained,
      gdm_predicted_stability = safe_cor(base_fit$model$predicted, scen_fit$model$predicted),
      gdm_term_rank_stability = rank_stability
    ),
    terms = terms
  )
}

# -----------------------------------------------------------------------------
# 4. SCENARIO BUILDERS
# -----------------------------------------------------------------------------

read_assemblage <- function(assemblage_id) {
  base <- file.path(INPUT_DIR, "assemblages", assemblage_id)
  list(
    station_frame = readr::read_csv(paste0(base, "__station_frame.csv"), show_col_types = FALSE) %>%
      mutate(station = as.character(station)),
    rtu = readr::read_csv(paste0(base, "__rtu_long.csv"), show_col_types = FALSE) %>%
      mutate(station = as.character(station)),
    species = readr::read_csv(paste0(base, "__species_long.csv"), show_col_types = FALSE) %>%
      mutate(station = as.character(station)),
    lookup = readr::read_csv(paste0(base, "__taxon_lookup.csv"), show_col_types = FALSE),
    adult_only_path = paste0(base, "__adult_only_strict_rtu_long.csv")
  )
}

lookup_maps <- function(lookup) {
  lookup %>%
    transmute(
      taxon_unit = taxon_unit_rtu,
      species_unit = taxon_unit_species,
      genus_unit = taxon_unit_genus,
      family_unit = taxon_unit_family,
      taxo_level,
      genus,
      family
    ) %>%
    filter(!is.na(taxon_unit)) %>%
    distinct(taxon_unit, .keep_all = TRUE)
}

build_rare_to_genus <- function(rtu_long, maps) {
  x <- rtu_long %>%
    left_join(maps, by = "taxon_unit")
  
  rare_species <- x %>%
    filter(taxo_level == "species", !is.na(species_unit)) %>%
    group_by(species_unit) %>%
    summarise(total_abundance = sum(abundance), .groups = "drop") %>%
    filter(total_abundance <= RARE_MAX_TOTAL_ABUNDANCE) %>%
    pull(species_unit)
  
  x %>%
    mutate(taxon_unit = if_else(
      taxo_level == "species" & species_unit %in% rare_species & !is.na(genus_unit),
      genus_unit,
      taxon_unit
    )) %>%
    select(station, taxon_unit, abundance) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
}

build_genus_resolution <- function(rtu_long, maps) {
  x <- rtu_long %>% left_join(maps, by = "taxon_unit") %>% filter(!is.na(genus_unit))
  list(
    baseline = x %>% transmute(station, taxon_unit, abundance) %>% group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop"),
    scenario = x %>% transmute(station, taxon_unit = genus_unit, abundance) %>% group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop")
  )
}

build_family_resolution <- function(rtu_long, maps) {
  x <- rtu_long %>% left_join(maps, by = "taxon_unit") %>% filter(!is.na(family_unit))
  list(
    baseline = x %>% transmute(station, taxon_unit, abundance) %>% group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop"),
    scenario = x %>% transmute(station, taxon_unit = family_unit, abundance) %>% group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop")
  )
}

make_rare_error_probabilities <- function(species_long, rate) {
  sp <- species_long %>%
    group_by(taxon_unit) %>%
    summarise(total_abundance = sum(abundance), .groups = "drop") %>%
    mutate(raw_weight = 1 / sqrt(total_abundance))
  
  scale_factor <- rate / weighted.mean(sp$raw_weight, w = sp$total_abundance)
  sp %>%
    transmute(
      taxon_unit,
      p_error = pmin(raw_weight * scale_factor, RARE_WEIGHTED_MAX_ERROR)
    )
}

# A coherent confusion map is drawn ONCE per source species and iteration, then
# applied at every station. The number of misidentified individuals can still
# vary by site through binomial sampling, but source -> target identity does not.
make_coherent_confusion_map <- function(source_meta, candidate_pool, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # CD_REF is optional for observed-pool maps. It is used for the regional
  # TAXREF pool to prevent a source from being reassigned to the same accepted
  # species under a different synonym/identifier.
  if (!"cd_ref" %in% names(source_meta)) source_meta$cd_ref <- NA_character_
  if (!"cd_ref" %in% names(candidate_pool)) candidate_pool$cd_ref <- NA_character_

  source_meta <- source_meta %>%
    transmute(
      source_taxon_unit = taxon_unit,
      source_genus = genus,
      source_cd_ref = as.character(cd_ref)
    ) %>%
    filter(!is.na(source_taxon_unit), !is.na(source_genus)) %>%
    distinct()

  candidate_pool <- candidate_pool %>%
    transmute(
      target_taxon_unit = taxon_unit,
      target_genus = genus,
      target_cd_ref = as.character(cd_ref)
    ) %>%
    filter(!is.na(target_taxon_unit), !is.na(target_genus)) %>%
    distinct()

  purrr::map_dfr(seq_len(nrow(source_meta)), function(i) {
    src <- source_meta[i, , drop = FALSE]

    candidates <- candidate_pool %>%
      filter(
        target_genus == src$source_genus,
        target_taxon_unit != src$source_taxon_unit,
        is.na(src$source_cd_ref) | is.na(target_cd_ref) | target_cd_ref != src$source_cd_ref
      )

    if (nrow(candidates) == 0L) {
      return(tibble(
        source_taxon_unit = src$source_taxon_unit,
        source_genus = src$source_genus,
        source_cd_ref = src$source_cd_ref,
        target_taxon_unit = NA_character_,
        target_genus = NA_character_,
        target_cd_ref = NA_character_,
        mapping_type = "no_candidate"
      ))
    }

    target <- candidates %>% slice_sample(n = 1L)

    tibble(
      source_taxon_unit = src$source_taxon_unit,
      source_genus = src$source_genus,
      source_cd_ref = src$source_cd_ref,
      target_taxon_unit = target$target_taxon_unit,
      target_genus = target$target_genus,
      target_cd_ref = target$target_cd_ref,
      mapping_type = "coherent_random"
    )
  })
}

read_regional_pool <- function(assemblage_id) {
  path <- file.path(REGIONAL_POOL_DIR, paste0(assemblage_id, "__regional_pool.csv"))
  if (!RUN_REGIONAL_POOL_SCENARIOS || !file.exists(path)) return(NULL)

  pool <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  required <- c("taxon_unit", "genus")
  if (!all(required %in% names(pool))) {
    warning("Regional pool skipped (missing columns): ", path)
    return(NULL)
  }

  pool %>%
    transmute(
      taxon_unit = as.character(taxon_unit),
      genus = stringr::str_squish(as.character(genus)),
      cd_ref = if ("cd_ref" %in% names(pool)) as.character(cd_ref) else NA_character_
    ) %>%
    filter(!is.na(taxon_unit), !is.na(genus), genus != "") %>%
    distinct()
}

read_regional_source_audit <- function(assemblage_id) {
  path <- file.path(
    REGIONAL_POOL_DIR,
    paste0(assemblage_id, "__observed_taxref_match_audit.csv")
  )
  if (!file.exists(path)) {
    return(tibble(taxon_unit = character(), cd_ref = character()))
  }

  audit <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  if (!all(c("observed_taxon_unit", "cd_ref") %in% names(audit))) {
    return(tibble(taxon_unit = character(), cd_ref = character()))
  }

  audit %>%
    transmute(
      taxon_unit = as.character(observed_taxon_unit),
      cd_ref = as.character(cd_ref)
    ) %>%
    filter(!is.na(taxon_unit)) %>%
    distinct(taxon_unit, .keep_all = TRUE)
}

read_expert_confusion_map <- function(assemblage_id) {
  path <- file.path(EXPERT_MAP_DIR, paste0(assemblage_id, "__expert_confusions.csv"))
  if (!RUN_EXPERT_CONFUSION_SCENARIOS || !file.exists(path)) return(NULL)
  
  map <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  required <- c("source_taxon_unit", "target_taxon_unit")
  if (!all(required %in% names(map))) {
    warning("Expert map skipped (missing columns): ", path)
    return(NULL)
  }
  
  if (!"enabled" %in% names(map)) map$enabled <- TRUE
  if (!"weight" %in% names(map)) map$weight <- 1
  
  map %>%
    mutate(
      enabled = tolower(as.character(enabled)) %in% c("true", "t", "1", "yes", "y"),
      weight = suppressWarnings(as.numeric(weight)),
      weight = dplyr::coalesce(weight, 1)
    ) %>%
    filter(enabled, !is.na(source_taxon_unit), !is.na(target_taxon_unit),
           source_taxon_unit != target_taxon_unit) %>%
    select(source_taxon_unit, target_taxon_unit, weight, any_of(c("comment"))) %>%
    distinct()
}

make_expert_confusion_map <- function(source_meta, expert_map, seed = NULL) {
  if (is.null(expert_map) || nrow(expert_map) == 0L) return(NULL)
  if (!is.null(seed)) set.seed(seed)
  
  source_allowed <- source_meta %>% select(source_taxon_unit = taxon_unit) %>% distinct()
  
  expert_map %>%
    inner_join(source_allowed, by = "source_taxon_unit") %>%
    group_by(source_taxon_unit) %>%
    group_modify(~ {
      choices <- .x
      target <- choices %>%
        slice_sample(n = 1L, weight_by = weight)
      tibble(
        target_taxon_unit = target$target_taxon_unit,
        mapping_type = "expert",
        expert_comment = target$comment %||% NA_character_
      )
    }) %>%
    ungroup()
}

simulate_from_coherent_map <- function(
    species_long,
    confusion_map,
    p_error_tbl,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  if (is.null(confusion_map) || nrow(confusion_map) == 0L) {
    return(species_long %>%
             group_by(station, taxon_unit) %>%
             summarise(abundance = sum(abundance), .groups = "drop"))
  }
  
  dat <- species_long %>%
    left_join(
      confusion_map %>%
        select(source_taxon_unit, target_taxon_unit),
      by = c("taxon_unit" = "source_taxon_unit")
    ) %>%
    left_join(p_error_tbl, by = "taxon_unit") %>%
    mutate(
      p_error = coalesce(p_error, 0),
      abundance = as.integer(round(abundance))
    )
  
  out <- purrr::pmap_dfr(
    dat,
    function(station, taxon_unit, abundance, target_taxon_unit, p_error, ...) {
      n <- as.integer(abundance)
      can_swap <- n > 0L &&
        !is.na(target_taxon_unit) &&
        is.finite(p_error) &&
        p_error > 0
      
      if (!can_swap) {
        return(tibble(station = station, taxon_unit = taxon_unit, abundance = n))
      }
      
      n_swap <- stats::rbinom(1L, size = n, prob = pmin(p_error, 1))
      bind_rows(
        if (n - n_swap > 0L) {
          tibble(station = station, taxon_unit = taxon_unit, abundance = n - n_swap)
        },
        if (n_swap > 0L) {
          tibble(station = station, taxon_unit = target_taxon_unit, abundance = n_swap)
        }
      )
    }
  )
  
  out %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
}

annotate_confusion_map <- function(confusion_map, species_long, p_error_tbl) {
  abundance_by_source <- species_long %>%
    group_by(taxon_unit) %>%
    summarise(source_total_abundance = sum(abundance), .groups = "drop") %>%
    rename(source_taxon_unit = taxon_unit)

  confusion_map %>%
    left_join(abundance_by_source, by = "source_taxon_unit") %>%
    left_join(
      p_error_tbl %>% rename(source_taxon_unit = taxon_unit),
      by = "source_taxon_unit"
    ) %>%
    mutate(
      source_total_abundance = coalesce(source_total_abundance, 0),
      p_error = coalesce(p_error, 0),
      has_target = !is.na(target_taxon_unit),
      expected_reassigned_individuals = source_total_abundance * p_error * has_target
    )
}

write_confusion_map <- function(
    assemblage_id,
    scenario,
    iter,
    map,
    species_long,
    p_error_tbl
) {
  if (is.null(map) || nrow(map) == 0L) return(invisible(NULL))
  if (!WRITE_ALL_COHERENT_CONFUSION_MAPS &&
      !(WRITE_EXAMPLE_COHERENT_CONFUSION_MAPS && iter == 1L)) {
    return(invisible(NULL))
  }

  path <- file.path(
    OUT_DIR, "per_assemblage",
    paste0(assemblage_id, "__", scenario, "__iter_", sprintf("%02d", iter), "__confusion_map.csv")
  )

  readr::write_csv(
    annotate_confusion_map(map, species_long, p_error_tbl),
    path
  )
  invisible(path)
}

build_scenarios <- function(assemblage_id, data, family_applicable, stage_candidate) {
  maps <- lookup_maps(data$lookup)
  
  scenario_list <- list(
    list(
      scenario = "rtu_mixed_best_available",
      scenario_family = "workflow",
      baseline_scenario = "rtu_mixed_best_available",
      unit_type = "mixed_RTU",
      iter = 1L,
      comm_long = data$rtu
    ),
    list(
      scenario = "rtu_species_only_drop_unresolved",
      scenario_family = "workflow",
      baseline_scenario = "rtu_mixed_best_available",
      unit_type = "species_only",
      iter = 1L,
      comm_long = data$species
    ),
    list(
      scenario = "rtu_rare_species_to_genus",
      scenario_family = "workflow",
      baseline_scenario = "rtu_mixed_best_available",
      unit_type = "mixed_RTU",
      iter = 1L,
      comm_long = build_rare_to_genus(data$rtu, maps)
    )
  )
  
  if (stage_candidate && file.exists(data$adult_only_path)) {
    adult <- readr::read_csv(data$adult_only_path, show_col_types = FALSE) %>%
      mutate(station = as.character(station))
    scenario_list <- append(scenario_list, list(list(
      scenario = "adult_only_strict",
      scenario_family = "ontogenetic_filtering",
      baseline_scenario = "rtu_mixed_best_available",
      unit_type = "adult_RTU",
      iter = 1L,
      comm_long = adult
    )))
  }
  
  genus_pair <- build_genus_resolution(data$rtu, maps)
  scenario_list <- append(scenario_list, list(
    list(
      scenario = "rtu_mixed_genus_resolvable",
      scenario_family = "resolution_genus",
      baseline_scenario = "rtu_mixed_genus_resolvable",
      unit_type = "genus_resolvable_RTU",
      iter = 1L,
      comm_long = genus_pair$baseline
    ),
    list(
      scenario = "rtu_genus_level",
      scenario_family = "resolution_genus",
      baseline_scenario = "rtu_mixed_genus_resolvable",
      unit_type = "genus_units",
      iter = 1L,
      comm_long = genus_pair$scenario
    )
  ))
  
  if (family_applicable) {
    family_pair <- build_family_resolution(data$rtu, maps)
    scenario_list <- append(scenario_list, list(
      list(
        scenario = "rtu_mixed_family_resolvable",
        scenario_family = "resolution_family",
        baseline_scenario = "rtu_mixed_family_resolvable",
        unit_type = "family_resolvable_RTU",
        iter = 1L,
        comm_long = family_pair$baseline
      ),
      list(
        scenario = "rtu_family_level",
        scenario_family = "resolution_family",
        baseline_scenario = "rtu_mixed_family_resolvable",
        unit_type = "family_units",
        iter = 1L,
        comm_long = family_pair$scenario
      )
    ))
  }
  
  species_meta <- maps %>%
    filter(!is.na(species_unit), !is.na(genus)) %>%
    transmute(taxon_unit = species_unit, genus) %>%
    distinct()
  
  if (nrow(data$species) == 0L || nrow(species_meta) == 0L) return(scenario_list)
  
  scenario_list <- append(scenario_list, list(list(
    scenario = "species_baseline",
    scenario_family = "species_error",
    baseline_scenario = "species_baseline",
    unit_type = "species_units",
    iter = 1L,
    comm_long = data$species
  )))
  
  observed_pool <- species_meta
  regional_pool <- read_regional_pool(assemblage_id)
  regional_source_meta <- species_meta %>%
    left_join(read_regional_source_audit(assemblage_id), by = "taxon_unit")
  expert_map <- read_expert_confusion_map(assemblage_id)
  
  fixed_p <- species_meta %>% transmute(taxon_unit, p_error = MAIN_ERROR_RATE)
  rare_p <- make_rare_error_probabilities(data$species, RARE_WEIGHTED_ERROR)
  expert_p <- species_meta %>%
    transmute(taxon_unit, p_error = EXPERT_ERROR_RATE)
  
  for (iter_i in seq_len(N_SIM)) {
    # Observed-pool source -> target map, constant across all sites in this iteration.
    observed_map <- make_coherent_confusion_map(
      source_meta = species_meta,
      candidate_pool = observed_pool,
      seed = 100000L + iter_i
    )
    write_confusion_map(
      assemblage_id, "species_observed_congeneric_10pct", iter_i,
      observed_map, data$species, fixed_p
    )
    
    scenario_list <- append(scenario_list, list(list(
      scenario = "species_observed_congeneric_10pct",
      scenario_family = "species_error_observed_pool",
      baseline_scenario = "species_baseline",
      unit_type = "species_units",
      iter = iter_i,
      comm_long = simulate_from_coherent_map(
        data$species, observed_map, fixed_p,
        seed = 110000L + iter_i
      )
    )))
    
    # Same coherent mapping principle, but probabilities are abundance-weighted
    # towards rare source species.
    rare_map <- make_coherent_confusion_map(
      source_meta = species_meta,
      candidate_pool = observed_pool,
      seed = 120000L + iter_i
    )
    write_confusion_map(
      assemblage_id, "species_observed_rare_weighted_10pct", iter_i,
      rare_map, data$species, rare_p
    )
    
    scenario_list <- append(scenario_list, list(list(
      scenario = "species_observed_rare_weighted_10pct",
      scenario_family = "species_error_rare_weighted",
      baseline_scenario = "species_baseline",
      unit_type = "species_units",
      iter = iter_i,
      comm_long = simulate_from_coherent_map(
        data$species, rare_map, rare_p,
        seed = 130000L + iter_i
      )
    )))
    
    # The TAXREF scenario is a single 10% reference scenario in the main text.
    # Its full 1–20% response is handled separately in the Appendix script.
    if (!is.null(regional_pool) && nrow(regional_pool) > 0L) {
      regional_map <- make_coherent_confusion_map(
        source_meta = regional_source_meta,
        candidate_pool = regional_pool,
        seed = 140000L + iter_i
      )
      write_confusion_map(
        assemblage_id, "species_regional_congeneric_10pct", iter_i,
        regional_map, data$species, fixed_p
      )
      
      scenario_list <- append(scenario_list, list(list(
        scenario = "species_regional_congeneric_10pct",
        scenario_family = "species_error_regional_pool",
        baseline_scenario = "species_baseline",
        unit_type = "species_units",
        iter = iter_i,
        comm_long = simulate_from_coherent_map(
          data$species, regional_map, fixed_p,
          seed = 150000L + iter_i
        )
      )))
    }
    
    # Expert scenarios are only built when a curated source -> target map exists.
    if (!is.null(expert_map) && nrow(expert_map) > 0L) {
      coherent_expert_map <- make_expert_confusion_map(
        source_meta = species_meta,
        expert_map = expert_map,
        seed = 160000L + iter_i
      )
      write_confusion_map(
        assemblage_id, "species_expert_confusion_10pct", iter_i,
        coherent_expert_map, data$species, expert_p
      )
      
      scenario_list <- append(scenario_list, list(list(
        scenario = "species_expert_confusion_10pct",
        scenario_family = "species_error_expert",
        baseline_scenario = "species_baseline",
        unit_type = "species_units",
        iter = iter_i,
        comm_long = simulate_from_coherent_map(
          data$species, coherent_expert_map, expert_p,
          seed = 170000L + iter_i
        )
      )))
    }
  }
  
  scenario_list
}

# -----------------------------------------------------------------------------
# 5. PER-ASSEMBLAGE ANALYSIS
# -----------------------------------------------------------------------------

analyse_assemblage <- function(def, env_meta) {
  id <- def$assemblage_id
  message("\nAnalysing ", id)
  data <- read_assemblage(id)
  scenarios <- build_scenarios(
    assemblage_id = id,
    data = data,
    family_applicable = def$family_coarsening_applicable,
    stage_candidate = def$stage_scenario_candidate
  )
  
  # Construct every matrix once.
  scenario_index <- purrr::map_dfr(scenarios, function(s) {
    tibble(
      assemblage_id = id,
      group_label = def$group_label,
      method = def$method,
      analysis_tier = def$analysis_tier,
      scenario = s$scenario,
      scenario_family = s$scenario_family,
      baseline_scenario = s$baseline_scenario,
      unit_type = s$unit_type,
      iter = s$iter
    )
  })
  
  matrices <- purrr::map(scenarios, ~ make_matrix(.x$comm_long, data$station_frame))
  names(matrices) <- purrr::map_chr(scenarios, ~ paste(.x$scenario, .x$iter, sep = "__"))
  
  # Alpha metrics and state summary for each scenario x iteration.
  alpha_by_site <- purrr::map2_dfr(scenarios, matrices, function(s, mat) {
    alpha_metrics(mat) %>%
      mutate(
        assemblage_id = id,
        scenario = s$scenario,
        scenario_family = s$scenario_family,
        baseline_scenario = s$baseline_scenario,
        unit_type = s$unit_type,
        iter = s$iter
      )
  })
  
  state_by_iter <- alpha_by_site %>%
    group_by(assemblage_id, scenario, scenario_family, baseline_scenario, unit_type, iter) %>%
    summarise(
      n_sampled_stations = n(),
      n_positive_stations = sum(total_abundance > 0),
      mean_total_abundance = mean(total_abundance),
      mean_q0 = mean(q0),
      mean_q1 = mean(q1),
      mean_q2 = mean(q2),
      .groups = "drop"
    )
  
  gamma_by_iter <- purrr::map2_dfr(scenarios, matrices, function(s, mat) {
    tibble(
      assemblage_id = id,
      scenario = s$scenario,
      scenario_family = s$scenario_family,
      baseline_scenario = s$baseline_scenario,
      unit_type = s$unit_type,
      iter = s$iter,
      gamma = community_gamma(mat)
    )
  })
  state_by_iter <- state_by_iter %>%
    left_join(gamma_by_iter, by = c("assemblage_id", "scenario", "scenario_family", "baseline_scenario", "unit_type", "iter"))
  
  # Compare every scenario with its own baseline.
  comparison_objects <- purrr::map2(scenarios, matrices, function(s, mat_scen) {
    base_key <- paste(s$baseline_scenario, 1L, sep = "__")
    mat_base <- matrices[[base_key]]
    a_base <- alpha_metrics(mat_base)
    a_scen <- alpha_metrics(mat_scen)
    alpha_compare <- tibble(
      abundance_stability = safe_cor(a_base$total_abundance, a_scen$total_abundance),
      q0_stability = safe_cor(a_base$q0, a_scen$q0),
      q1_stability = safe_cor(a_base$q1, a_scen$q1),
      q2_stability = safe_cor(a_base$q2, a_scen$q2)
    )
    base_state <- state_by_iter %>% filter(scenario == s$baseline_scenario, iter == 1L) %>% slice(1)
    scen_state <- state_by_iter %>% filter(scenario == s$scenario, iter == s$iter) %>% slice(1)
    beta_compare <- compare_beta(mat_base, mat_scen)
    
    run_gdm <- s$scenario == s$baseline_scenario || s$iter <= GDM_MAX_STOCHASTIC_ITERS
    gdm_bray <- if (run_gdm) compare_gdm(mat_base, mat_scen, env_meta, "bray") else list(summary = tibble(), terms = tibble())
    gdm_sor <- if (run_gdm) compare_gdm(mat_base, mat_scen, env_meta, "sorensen") else list(summary = tibble(), terms = tibble())
    
    info <- tibble(
      assemblage_id = id,
      group_label = def$group_label,
      method = def$method,
      analysis_tier = def$analysis_tier,
      scenario = s$scenario,
      scenario_family = s$scenario_family,
      baseline_scenario = s$baseline_scenario,
      unit_type = s$unit_type,
      iter = s$iter,
      gamma_change_pct = if (base_state$gamma > 0) 100 * (scen_state$gamma / base_state$gamma - 1) else NA_real_,
      mean_q0_change_pct = if (base_state$mean_q0 > 0) 100 * (scen_state$mean_q0 / base_state$mean_q0 - 1) else NA_real_,
      mean_q1_change_pct = if (base_state$mean_q1 > 0) 100 * (scen_state$mean_q1 / base_state$mean_q1 - 1) else NA_real_,
      mean_q2_change_pct = if (base_state$mean_q2 > 0) 100 * (scen_state$mean_q2 / base_state$mean_q2 - 1) else NA_real_
    )
    
    gb <- if (nrow(gdm_bray$summary) > 0L) {
      gdm_bray$summary %>% rename(
        gdm_bray_n_sites = gdm_n_sites,
        gdm_bray_explained_baseline = gdm_explained_baseline,
        gdm_bray_explained_scenario = gdm_explained_scenario,
        gdm_bray_delta_explained = gdm_delta_explained,
        gdm_bray_predicted_stability = gdm_predicted_stability,
        gdm_bray_term_rank_stability = gdm_term_rank_stability
      )
    } else tibble(
      gdm_bray_n_sites = NA_real_, gdm_bray_explained_baseline = NA_real_,
      gdm_bray_explained_scenario = NA_real_, gdm_bray_delta_explained = NA_real_,
      gdm_bray_predicted_stability = NA_real_, gdm_bray_term_rank_stability = NA_real_
    )
    
    gs <- if (nrow(gdm_sor$summary) > 0L) {
      gdm_sor$summary %>% rename(
        gdm_sorensen_n_sites = gdm_n_sites,
        gdm_sorensen_explained_baseline = gdm_explained_baseline,
        gdm_sorensen_explained_scenario = gdm_explained_scenario,
        gdm_sorensen_delta_explained = gdm_delta_explained,
        gdm_sorensen_predicted_stability = gdm_predicted_stability,
        gdm_sorensen_term_rank_stability = gdm_term_rank_stability
      )
    } else tibble(
      gdm_sorensen_n_sites = NA_real_, gdm_sorensen_explained_baseline = NA_real_,
      gdm_sorensen_explained_scenario = NA_real_, gdm_sorensen_delta_explained = NA_real_,
      gdm_sorensen_predicted_stability = NA_real_, gdm_sorensen_term_rank_stability = NA_real_
    )
    
    terms <- bind_rows(
      if (nrow(gdm_bray$terms) > 0L) {
        gdm_bray$terms %>% mutate(distance = "Bray-Curtis")
      } else {
        empty_gdm_comparison_terms() %>% mutate(distance = character())
      },
      if (nrow(gdm_sor$terms) > 0L) {
        gdm_sor$terms %>% mutate(distance = "Sørensen")
      } else {
        empty_gdm_comparison_terms() %>% mutate(distance = character())
      }
    ) %>% mutate(
      assemblage_id = id,
      scenario = s$scenario,
      scenario_family = s$scenario_family,
      baseline_scenario = s$baseline_scenario,
      unit_type = s$unit_type,
      iter = s$iter
    )
    
    list(summary = bind_cols(info, alpha_compare, beta_compare, gb, gs), terms = terms)
  })
  
  comparison_by_iter <- purrr::map_dfr(comparison_objects, "summary")
  gdm_terms_by_iter <- purrr::map_dfr(comparison_objects, "terms")
  
  # Write per-assemblage outputs.
  prefix <- file.path(OUT_DIR, "per_assemblage", id)
  readr::write_csv(scenario_index, paste0(prefix, "__scenario_definitions.csv"))
  readr::write_csv(alpha_by_site, paste0(prefix, "__alpha_by_site.csv"))
  readr::write_csv(state_by_iter, paste0(prefix, "__state_by_iter.csv"))
  readr::write_csv(comparison_by_iter, paste0(prefix, "__results_by_iter.csv"))
  readr::write_csv(gdm_terms_by_iter, paste0(prefix, "__gdm_terms_by_iter.csv"))
  
  list(
    scenario_index = scenario_index,
    alpha_by_site = alpha_by_site,
    state_by_iter = state_by_iter,
    results_by_iter = comparison_by_iter,
    gdm_terms_by_iter = gdm_terms_by_iter
  )
}

# -----------------------------------------------------------------------------
# 6. RUN ALL ASSEMBLAGES
# -----------------------------------------------------------------------------

manifest <- readr::read_csv(file.path(INPUT_DIR, "assemblage_manifest.csv"), show_col_types = FALSE) %>%
  mutate(
    family_coarsening_applicable = as.logical(family_coarsening_applicable),
    stage_scenario_candidate = as.logical(stage_scenario_candidate)
  )

selected_ids <- MAIN_ASSEMBLAGES
if (RUN_SECONDARY_ASSEMBLAGES) selected_ids <- c(selected_ids, SECONDARY_ASSEMBLAGES)

selected_manifest <- manifest %>%
  filter(assemblage_id %in% selected_ids) %>%
  mutate(analysis_tier = if_else(assemblage_id %in% MAIN_ASSEMBLAGES, "main", "secondary")) %>%
  arrange(match(assemblage_id, selected_ids))

if (nrow(selected_manifest) == 0L) stop("No selected assemblages found in manifest.")
readr::write_csv(selected_manifest, file.path(OUT_DIR, "selected_assemblages.csv"))

regional_pool_availability <- selected_manifest %>%
  transmute(
    assemblage_id,
    regional_pool_file = file.path(REGIONAL_POOL_DIR, paste0(assemblage_id, "__regional_pool.csv")),
    regional_pool_available = file.exists(regional_pool_file),
    expert_map_file = file.path(EXPERT_MAP_DIR, paste0(assemblage_id, "__expert_confusions.csv")),
    expert_map_available = file.exists(expert_map_file)
  )
readr::write_csv(
  regional_pool_availability,
  file.path(OUT_DIR, "regional_pool_availability.csv")
)

if (RUN_REGIONAL_POOL_SCENARIOS && !any(regional_pool_availability$regional_pool_available)) {
  warning(
    "No regional TAXREF pools were found. The main run will omit regional-pool ",
    "error scenarios. Run 05_build_multitaxa_taxref_pools_v2.R first."
  )
}

all_out <- purrr::map(seq_len(nrow(selected_manifest)), function(i) {
  analyse_assemblage(selected_manifest[i, ], env_station)
})

scenario_definitions_all <- purrr::map_dfr(all_out, "scenario_index")
alpha_by_site_all <- purrr::map_dfr(all_out, "alpha_by_site")
state_by_iter_all <- purrr::map_dfr(all_out, "state_by_iter")
results_by_iter_all <- purrr::map_dfr(all_out, "results_by_iter")
gdm_terms_by_iter_all <- purrr::map_dfr(all_out, "gdm_terms_by_iter")

# Summaries: medians are robust for stochastic error scenarios; deterministic
# scenarios simply have n_iter = 1.
result_numeric <- c(
  "gamma_change_pct", "mean_q0_change_pct", "mean_q1_change_pct", "mean_q2_change_pct",
  "abundance_stability", "q0_stability", "q1_stability", "q2_stability",
  "bray_stability", "sorensen_stability", "ordination_procrustes_r2",
  "bray_mean_change_pct", "sorensen_mean_change_pct",
  "gdm_bray_delta_explained", "gdm_bray_predicted_stability", "gdm_bray_term_rank_stability",
  "gdm_sorensen_delta_explained", "gdm_sorensen_predicted_stability", "gdm_sorensen_term_rank_stability"
)
if (RUN_JACCARD_SENSITIVITY) {
  result_numeric <- c(result_numeric, "jaccard_stability", "jaccard_mean_change_pct")
}

results_summary <- results_by_iter_all %>%
  group_by(assemblage_id, group_label, method, analysis_tier,
           scenario, scenario_family, baseline_scenario, unit_type) %>%
  summarise(
    across(
      all_of(result_numeric),
      list(
        median = ~ median(.x, na.rm = TRUE),
        p10 = ~ quantile(.x, 0.10, na.rm = TRUE, names = FALSE),
        p90 = ~ quantile(.x, 0.90, na.rm = TRUE, names = FALSE)
      ),
      .names = "{.col}_{.fn}"
    ),
    n_iter = n(),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_median") | ends_with("_p10") | ends_with("_p90"), ~ ifelse(is.nan(.x), NA_real_, .x)))

state_summary <- state_by_iter_all %>%
  group_by(assemblage_id, scenario, scenario_family, baseline_scenario, unit_type) %>%
  summarise(
    across(c(mean_total_abundance, mean_q0, mean_q1, mean_q2, gamma),
           list(median = ~ median(.x, na.rm = TRUE), p10 = ~ quantile(.x, .10, na.rm = TRUE, names = FALSE), p90 = ~ quantile(.x, .90, na.rm = TRUE, names = FALSE)),
           .names = "{.col}_{.fn}"),
    n_iter = n(),
    .groups = "drop"
  )

# GDM model fits can legitimately yield no non-zero I-spline terms, particularly
# for sparse assemblages or degenerate scenario matrices. Preserve the pipeline
# and write a typed empty table rather than failing during group_by().
.required_gdm_term_cols <- c(
  "assemblage_id", "scenario", "scenario_family", "baseline_scenario",
  "unit_type", "distance", "predictor",
  "baseline_ispline_weight", "scenario_ispline_weight", "delta_ispline_weight"
)

gdm_terms_summary <- if (
  nrow(gdm_terms_by_iter_all) == 0L ||
  !all(.required_gdm_term_cols %in% names(gdm_terms_by_iter_all))
) {
  empty_gdm_terms_summary()
} else {
  gdm_terms_by_iter_all %>%
    group_by(assemblage_id, scenario, scenario_family, baseline_scenario, unit_type, distance, predictor) %>%
    summarise(
      baseline_ispline_weight_median = median(baseline_ispline_weight, na.rm = TRUE),
      scenario_ispline_weight_median = median(scenario_ispline_weight, na.rm = TRUE),
      delta_ispline_weight_median = median(delta_ispline_weight, na.rm = TRUE),
      n_iter = n(),
      .groups = "drop"
    ) %>%
    mutate(across(ends_with("_median"), ~ ifelse(is.nan(.x), NA_real_, .x)))
}

if (nrow(gdm_terms_summary) == 0L) {
  warning(
    "No valid GDM term weights were extracted. The alpha/beta pipeline completed, ",
    "but inspect results_by_iter_all.csv for GDM explained-deviance fields before ",
    "interpreting GDM terms."
  )
}

readr::write_csv(scenario_definitions_all, file.path(OUT_DIR, "scenario_definitions_all.csv"))
readr::write_csv(alpha_by_site_all, file.path(OUT_DIR, "alpha_by_site_all.csv"))
readr::write_csv(state_by_iter_all, file.path(OUT_DIR, "state_by_iter_all.csv"))
readr::write_csv(results_by_iter_all, file.path(OUT_DIR, "results_by_iter_all.csv"))
readr::write_csv(results_summary, file.path(OUT_DIR, "results_summary.csv"))
readr::write_csv(state_summary, file.path(OUT_DIR, "state_summary.csv"))
readr::write_csv(gdm_terms_by_iter_all, file.path(OUT_DIR, "gdm_terms_by_iter_all.csv"))
readr::write_csv(gdm_terms_summary, file.path(OUT_DIR, "gdm_terms_summary.csv"))

# -----------------------------------------------------------------------------
# 7. DIAGNOSTIC FIGURES
# -----------------------------------------------------------------------------

scenario_labels <- c(
  "rtu_rare_species_to_genus" = "Rare taxa\nreported at genus",
  "rtu_species_only_drop_unresolved" = "Drop unresolved\nrecords",
  "adult_only_strict" = "Adults only",
  "rtu_genus_level" = "Genus-level\nresolution",
  "rtu_family_level" = "Family-level\nresolution",
  "species_observed_congeneric_10pct" = "Observed-pool\nerror 10%",
  "species_observed_rare_weighted_10pct" = "Rare-weighted\nerror 10%",
  "species_regional_congeneric_10pct" = "Regional-pool\nerror 10%",
  "species_expert_confusion_10pct" = "Expert\nconfusion 10%"
)

# Main balanced heatmap: alpha, beta, and ecological inference appear together.
heatmap_main <- results_summary %>%
  filter(scenario != baseline_scenario) %>%
  select(assemblage_id, group_label, method, analysis_tier, scenario,
         q0_stability_median, q1_stability_median, q2_stability_median,
         bray_stability_median, sorensen_stability_median, ordination_procrustes_r2_median,
         gdm_bray_predicted_stability_median, gdm_sorensen_predicted_stability_median) %>%
  pivot_longer(
    cols = -c(assemblage_id, group_label, method, analysis_tier, scenario),
    names_to = "metric", values_to = "stability"
  ) %>%
  mutate(
    metric = recode(
      metric,
      q0_stability_median = "Alpha: q0",
      q1_stability_median = "Alpha: q1",
      q2_stability_median = "Alpha: q2",
      bray_stability_median = "Beta: Bray-Curtis",
      sorensen_stability_median = "Beta: Sørensen",
      ordination_procrustes_r2_median = "Beta: ordination",
      gdm_bray_predicted_stability_median = "Inference: GDM Bray",
      gdm_sorensen_predicted_stability_median = "Inference: GDM Sørensen"
    ),
    scenario_label = recode(scenario, !!!scenario_labels, .default = scenario),
    loss = 1 - stability,
    metric = factor(metric, levels = c(
      "Alpha: q0", "Alpha: q1", "Alpha: q2",
      "Beta: Bray-Curtis", "Beta: Sørensen", "Beta: ordination",
      "Inference: GDM Bray", "Inference: GDM Sørensen"
    ))
  )

p_heatmap <- ggplot(heatmap_main, aes(x = scenario_label, y = metric, fill = loss)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_wrap(~ assemblage_id, ncol = 2, scales = "free_x") +
  scale_fill_gradient(limits = c(0, 1), na.value = "grey85", name = "Loss of\nstability") +
  labs(
    title = "Taxonomic uncertainty across RMQS 2024 soil-invertebrate assemblages",
    subtitle = "Balanced assessment of alpha diversity, beta diversity and ecological turnover inference",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(OUT_DIR, "figures", "FigM1_multitaxa_robustness_heatmap.png"), p_heatmap, width = 15, height = 11, dpi = 320)
ggsave(file.path(OUT_DIR, "figures", "FigM1_multitaxa_robustness_heatmap.pdf"), p_heatmap, width = 15, height = 11)

# Directional alpha-beta effects: gamma and mean beta-distance changes.
effects_long <- results_summary %>%
  filter(scenario != baseline_scenario) %>%
  select(assemblage_id, group_label, method, analysis_tier, scenario,
         gamma_change_pct_median, bray_mean_change_pct_median, sorensen_mean_change_pct_median) %>%
  pivot_longer(
    cols = ends_with("_median"), names_to = "metric", values_to = "change_pct"
  ) %>%
  mutate(
    metric = recode(
      metric,
      gamma_change_pct_median = "Gamma richness",
      bray_mean_change_pct_median = "Mean Bray-Curtis dissimilarity",
      sorensen_mean_change_pct_median = "Mean Sørensen dissimilarity"
    ),
    scenario_label = recode(scenario, !!!scenario_labels, .default = scenario)
  )

p_effects <- ggplot(effects_long, aes(x = scenario_label, y = change_pct, colour = assemblage_id, group = assemblage_id)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(
    title = "Directional effects of taxonomic uncertainty on inventory and beta diversity",
    x = NULL, y = "Change relative to the appropriate baseline (%)", colour = "Assemblage"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave(file.path(OUT_DIR, "figures", "FigM2_multitaxa_alpha_beta_effects.png"), p_effects, width = 15, height = 11, dpi = 320)
ggsave(file.path(OUT_DIR, "figures", "FigM2_multitaxa_alpha_beta_effects.pdf"), p_effects, width = 15, height = 11)

# Ecological-inference diagnostic: change in explained turnover by GDM.
gdm_effects <- results_summary %>%
  filter(scenario != baseline_scenario) %>%
  select(assemblage_id, group_label, method, analysis_tier, scenario,
         gdm_bray_delta_explained_median, gdm_sorensen_delta_explained_median) %>%
  pivot_longer(
    cols = ends_with("_median"), names_to = "distance", values_to = "delta_explained"
  ) %>%
  mutate(
    distance = recode(
      distance,
      gdm_bray_delta_explained_median = "Bray-Curtis GDM",
      gdm_sorensen_delta_explained_median = "Sørensen GDM"
    ),
    scenario_label = recode(scenario, !!!scenario_labels, .default = scenario)
  )

p_gdm <- ggplot(gdm_effects, aes(x = scenario_label, y = delta_explained, colour = assemblage_id, group = assemblage_id)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3)) +
  facet_wrap(~ distance, ncol = 1) +
  labs(
    title = "Taxonomic uncertainty and GDM-based ecological inference",
    subtitle = "Change in explained environmental turnover relative to the appropriate baseline",
    x = NULL, y = "Delta explained deviance", colour = "Assemblage"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave(file.path(OUT_DIR, "figures", "FigM3_multitaxa_gdm_effects.png"), p_gdm, width = 15, height = 9, dpi = 320)
ggsave(file.path(OUT_DIR, "figures", "FigM3_multitaxa_gdm_effects.pdf"), p_gdm, width = 15, height = 9)

message("\nCompleted multi-taxon uncertainty pipeline.")
message("Key outputs:")
message("  - ", file.path(OUT_DIR, "results_summary.csv"))
message("  - ", file.path(OUT_DIR, "gdm_terms_summary.csv"))
message("  - ", file.path(OUT_DIR, "figures", "FigM1_multitaxa_robustness_heatmap.png"))
message("\nCoherent-map audit files (iteration 1 by default) are stored in:")
message("  ", file.path(OUT_DIR, "per_assemblage"))
message("Regional pool availability audit:")
message("  ", file.path(OUT_DIR, "regional_pool_availability.csv"))
message("\nRun the separate 1–20% appendix-gradient script after reviewing this main 10% analysis.")
