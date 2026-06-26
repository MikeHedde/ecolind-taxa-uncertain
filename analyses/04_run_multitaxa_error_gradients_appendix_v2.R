# =============================================================================
# 02_run_multitaxa_error_gradients_appendix_v3_fast.R
# -----------------------------------------------------------------------------
# Fast, checkpointed appendix analysis: 1–20% taxonomic-ID error gradients.
#
# Compared with v2:
#   - source -> target maps remain coherent within an iteration;
#   - community matrices are created once, then simulated directly as matrices;
#   - baseline Bray and Sørensen distances are calculated once per assemblage;
#   - progress is printed at every iteration;
#   - each assemblage × iteration is checkpointed and can be resumed safely.
#
# It does NOT fit GDMs. This appendix focuses on alpha and beta response curves.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(janitor)
  library(ggplot2)
  library(vegan)
})

# ---- 0. Settings -------------------------------------------------------------
INPUT_DIR <- "outputs_multitaxa_2024"
OUT_DIR <- file.path(INPUT_DIR, "uncertainty_results", "appendix_error_gradients_v3")
CHECKPOINT_DIR <- file.path(OUT_DIR, "checkpoints")
REGIONAL_POOL_DIR <- "regional_pools"
EXPERT_MAP_DIR <- "expert_confusion_maps"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

ASSEMBLAGES <- c(
  "collembola_soil_core",
  "araneae_pitfall",
  "carabidae_pitfall",
  "formicidae_pitfall",
  "isopoda_pitfall",
  "isopoda_hand_sorting",
  "diplopoda_pitfall",
  "diplopoda_hand_sorting"
)

ERROR_RATES <- seq(0.01, 0.20, by = 0.01)

# 20 iterations are enough for smooth median/10th–90th percentile curves and
# make this appendix tractable. Use 30L only for a final sensitivity rerun.
N_SIM <- 20L

RARE_WEIGHTED_MAX_ERROR <- 0.25
MIN_BETA_SITES <- 8L

# TRUE: reuse completed checkpoints from an interrupted run.
RESUME_FROM_CHECKPOINTS <- TRUE

set.seed(20260626)

# ---- 1. Lightweight helpers --------------------------------------------------
safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || dplyr::n_distinct(x[keep]) < 2L ||
      dplyr::n_distinct(y[keep]) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

read_assemblage <- function(assemblage_id) {
  base <- file.path(INPUT_DIR, "assemblages", assemblage_id)
  
  list(
    station_frame = readr::read_csv(
      paste0(base, "__station_frame.csv"), show_col_types = FALSE
    ) %>%
      mutate(station = as.character(station)),
    species = readr::read_csv(
      paste0(base, "__species_long.csv"), show_col_types = FALSE
    ) %>%
      transmute(
        station = as.character(station),
        taxon_unit = as.character(taxon_unit),
        abundance = as.integer(round(abundance))
      ) %>%
      filter(!is.na(station), !is.na(taxon_unit), abundance > 0),
    lookup = readr::read_csv(
      paste0(base, "__taxon_lookup.csv"), show_col_types = FALSE
    )
  )
}

species_meta_from_lookup <- function(lookup) {
  lookup %>%
    transmute(
      taxon_unit = as.character(taxon_unit_species),
      genus = as.character(genus)
    ) %>%
    filter(!is.na(taxon_unit), !is.na(genus)) %>%
    distinct()
}

read_regional_pool <- function(assemblage_id) {
  path <- file.path(
    REGIONAL_POOL_DIR,
    paste0(assemblage_id, "__regional_pool.csv")
  )
  if (!file.exists(path)) return(NULL)
  
  pool <- readr::read_csv(path, show_col_types = FALSE) %>%
    janitor::clean_names()
  
  if (!all(c("taxon_unit", "genus") %in% names(pool))) return(NULL)
  
  pool %>%
    transmute(
      taxon_unit = as.character(taxon_unit),
      genus = as.character(genus),
      cd_ref = if ("cd_ref" %in% names(pool)) as.character(cd_ref) else NA_character_
    ) %>%
    filter(!is.na(taxon_unit), !is.na(genus)) %>%
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
  
  audit <- readr::read_csv(path, show_col_types = FALSE) %>%
    janitor::clean_names()
  
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

read_expert_map <- function(assemblage_id) {
  path <- file.path(
    EXPERT_MAP_DIR,
    paste0(assemblage_id, "__expert_confusions.csv")
  )
  if (!file.exists(path)) return(NULL)
  
  dat <- readr::read_csv(path, show_col_types = FALSE) %>%
    janitor::clean_names()
  
  if (!all(c("source_taxon_unit", "target_taxon_unit") %in% names(dat))) {
    return(NULL)
  }
  if (!"enabled" %in% names(dat)) dat$enabled <- TRUE
  if (!"weight" %in% names(dat)) dat$weight <- 1
  
  dat %>%
    mutate(
      enabled = tolower(as.character(enabled)) %in% c("true", "t", "1", "yes", "y"),
      weight = suppressWarnings(as.numeric(weight)),
      weight = coalesce(weight, 1)
    ) %>%
    filter(
      enabled,
      !is.na(source_taxon_unit),
      !is.na(target_taxon_unit),
      source_taxon_unit != target_taxon_unit
    ) %>%
    select(source_taxon_unit, target_taxon_unit, weight) %>%
    distinct()
}

# Construct a matrix with a predeclared universe of taxa. This is much faster
# than repeatedly completing / pivoting long tables within the simulation loop.
make_count_matrix <- function(species_long, stations, taxa) {
  n_sites <- length(stations)
  n_taxa <- length(taxa)
  
  if (n_taxa == 0L) {
    return(matrix(0L, nrow = n_sites, ncol = 0L,
                  dimnames = list(stations, character())))
  }
  
  i <- match(species_long$station, stations)
  j <- match(species_long$taxon_unit, taxa)
  keep <- !is.na(i) & !is.na(j) & species_long$abundance > 0
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "The recommended package 'Matrix' is required for the fast matrix builder. ",
      "Install it with install.packages('Matrix')."
    )
  }
  
  # sparseMatrix adds duplicate station × taxon records by construction, then
  # converts to a regular matrix for vegan.
  mat <- as.matrix(
    Matrix::sparseMatrix(
      i = i[keep],
      j = j[keep],
      x = species_long$abundance[keep],
      dims = c(n_sites, n_taxa),
      dimnames = list(stations, taxa)
    )
  )
  
  storage.mode(mat) <- "integer"
  mat
}

alpha_vectors <- function(mat) {
  q0 <- rowSums(mat > 0L)
  q1 <- exp(vegan::diversity(mat, index = "shannon"))
  q2 <- vegan::diversity(mat, index = "invsimpson")
  q1[!is.finite(q1)] <- 0
  q2[!is.finite(q2)] <- 0
  list(q0 = q0, q1 = q1, q2 = q2)
}

prepare_beta_matrix <- function(mat) {
  mat <- mat[rowSums(mat) > 0L, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0L, drop = FALSE]
  
  if (nrow(mat) < MIN_BETA_SITES || ncol(mat) < 2L) {
    return(NULL)
  }
  mat
}

beta_distances <- function(mat) {
  mat <- prepare_beta_matrix(mat)
  if (is.null(mat)) return(list(bray = NULL, sorensen = NULL))
  
  list(
    bray = vegan::vegdist(mat, method = "bray"),
    sorensen = vegan::vegdist((mat > 0L) * 1L, method = "bray", binary = TRUE)
  )
}

make_rare_probability_vector <- function(base_mat, rate) {
  total_abundance <- colSums(base_mat)
  active <- total_abundance > 0L
  
  p <- numeric(ncol(base_mat))
  if (!any(active)) return(p)
  
  raw_weight <- 1 / sqrt(total_abundance[active])
  scale_factor <- rate / stats::weighted.mean(
    raw_weight,
    w = total_abundance[active]
  )
  
  p[active] <- pmin(raw_weight * scale_factor, RARE_WEIGHTED_MAX_ERROR)
  p
}

make_coherent_map <- function(source_meta, candidate_pool, seed) {
  set.seed(seed)
  
  if (!"cd_ref" %in% names(source_meta)) source_meta$cd_ref <- NA_character_
  if (!"cd_ref" %in% names(candidate_pool)) candidate_pool$cd_ref <- NA_character_
  
  src <- source_meta %>%
    transmute(
      source_taxon_unit = taxon_unit,
      source_genus = genus,
      source_cd_ref = as.character(cd_ref)
    ) %>%
    distinct()
  
  candidates <- candidate_pool %>%
    transmute(
      target_taxon_unit = taxon_unit,
      target_genus = genus,
      target_cd_ref = as.character(cd_ref)
    ) %>%
    distinct()
  
  purrr::map_dfr(seq_len(nrow(src)), function(i) {
    x <- src[i, , drop = FALSE]
    
    possible <- candidates %>%
      filter(
        target_genus == x$source_genus,
        target_taxon_unit != x$source_taxon_unit,
        is.na(x$source_cd_ref) | is.na(target_cd_ref) |
          target_cd_ref != x$source_cd_ref
      )
    
    if (!nrow(possible)) {
      return(tibble(
        source_taxon_unit = x$source_taxon_unit,
        target_taxon_unit = NA_character_
      ))
    }
    
    tibble(
      source_taxon_unit = x$source_taxon_unit,
      target_taxon_unit = sample(possible$target_taxon_unit, 1L)
    )
  })
}

make_expert_map <- function(source_meta, expert_map, seed) {
  if (is.null(expert_map) || !nrow(expert_map)) return(NULL)
  
  set.seed(seed)
  allowed_sources <- source_meta %>% pull(taxon_unit)
  
  expert_map %>%
    filter(source_taxon_unit %in% allowed_sources) %>%
    group_by(source_taxon_unit) %>%
    group_modify(~ {
      choice <- slice_sample(.x, n = 1L, weight_by = weight)
      tibble(target_taxon_unit = choice$target_taxon_unit)
    }) %>%
    ungroup()
}

map_to_indices <- function(map, taxa) {
  out <- rep(NA_integer_, length(taxa))
  names(out) <- taxa
  
  if (is.null(map) || !nrow(map)) return(out)
  
  src_i <- match(map$source_taxon_unit, taxa)
  tgt_i <- match(map$target_taxon_unit, taxa)
  
  keep <- !is.na(src_i) & !is.na(tgt_i) & src_i != tgt_i
  out[src_i[keep]] <- tgt_i[keep]
  out
}

simulate_matrix <- function(base_mat, target_index, p_error, seed) {
  set.seed(seed)
  
  n_sites <- nrow(base_mat)
  n_taxa <- ncol(base_mat)
  
  valid_source <- !is.na(target_index) & p_error > 0
  p_eff <- ifelse(valid_source, pmin(p_error, 1), 0)
  
  # R matrices are column-major: repeat each taxon's probability for all sites.
  swapped <- matrix(
    stats::rbinom(
      n = length(base_mat),
      size = as.integer(base_mat),
      prob = rep(p_eff, each = n_sites)
    ),
    nrow = n_sites,
    ncol = n_taxa,
    dimnames = dimnames(base_mat)
  )
  
  scenario <- base_mat - swapped
  
  # A source may share a target with other sources. Add reassigned individuals
  # once for each distinct target.
  active <- which(valid_source)
  by_target <- split(active, target_index[active])
  
  for (target in names(by_target)) {
    source_cols <- by_target[[target]]
    scenario[, as.integer(target)] <- scenario[, as.integer(target)] +
      rowSums(swapped[, source_cols, drop = FALSE])
  }
  
  storage.mode(scenario) <- "integer"
  scenario
}

map_diagnostics_matrix <- function(base_mat, target_index, p_error) {
  source_total <- colSums(base_mat)
  eligible <- !is.na(target_index) & p_error > 0 & source_total > 0
  total <- sum(source_total)
  
  tibble(
    eligible_individual_share = if (total > 0) {
      sum(source_total[eligible]) / total
    } else NA_real_,
    expected_reassigned_pct = if (total > 0) {
      100 * sum(source_total[eligible] * p_error[eligible]) / total
    } else NA_real_,
    n_eligible_source_species = sum(eligible),
    n_source_species = sum(source_total > 0)
  )
}

summarise_gradient <- function(data) {
  data %>%
    group_by(assemblage_id, mechanism, error_rate) %>%
    summarise(
      across(
        c(
          eligible_individual_share, expected_reassigned_pct,
          n_eligible_source_species, n_source_species,
          gamma_change_pct, mean_q0_change_pct, mean_q1_change_pct,
          mean_q2_change_pct, bray_mean_change_pct, sorensen_mean_change_pct,
          q0_stability, bray_stability, sorensen_stability
        ),
        list(
          median = ~ median(.x, na.rm = TRUE),
          p10 = ~ quantile(.x, .10, na.rm = TRUE, names = FALSE),
          p90 = ~ quantile(.x, .90, na.rm = TRUE, names = FALSE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n_iter = n(),
      .groups = "drop"
    ) %>%
    mutate(
      across(
        ends_with("_median") | ends_with("_p10") | ends_with("_p90"),
        ~ ifelse(is.nan(.x), NA_real_, .x)
      )
    )
}

# ---- 2. Run gradient simulations --------------------------------------------
if (!dir.exists(file.path(INPUT_DIR, "assemblages"))) {
  stop("Prepared assemblages are missing. Run 00_prepare_multitaxa_inputs_2024.R first.")
}

manifest <- readr::read_csv(
  file.path(INPUT_DIR, "assemblage_manifest.csv"),
  show_col_types = FALSE
) %>%
  filter(assemblage_id %in% ASSEMBLAGES)

all_results <- list()

for (assemblage_id in manifest$assemblage_id) {
  message("\n========================================================================")
  message("Gradient analysis: ", assemblage_id)
  
  dat <- read_assemblage(assemblage_id)
  meta <- species_meta_from_lookup(dat$lookup)
  
  if (!nrow(dat$species) || !nrow(meta)) {
    message("  Skipped: no species-level community.")
    next
  }
  
  regional_pool <- read_regional_pool(assemblage_id)
  regional_meta <- meta %>%
    left_join(read_regional_source_audit(assemblage_id), by = "taxon_unit")
  expert_map <- read_expert_map(assemblage_id)
  
  # The universe must contain every possible target taxon across mechanisms.
  taxa_universe <- unique(c(
    meta$taxon_unit,
    if (!is.null(regional_pool)) regional_pool$taxon_unit else character(),
    if (!is.null(expert_map)) expert_map$target_taxon_unit else character()
  ))
  taxa_universe <- sort(taxa_universe)
  
  stations <- as.character(dat$station_frame$station)
  base_mat <- make_count_matrix(dat$species, stations, taxa_universe)
  
  base_alpha <- alpha_vectors(base_mat)
  base_gamma <- sum(colSums(base_mat) > 0L)
  base_dist <- beta_distances(base_mat)
  
  if (is.null(base_dist$bray) || is.null(base_dist$sorensen)) {
    message("  Skipped: insufficient non-empty stations or taxa for beta metrics.")
    next
  }
  
  base_bray_mean <- mean(as.numeric(base_dist$bray), na.rm = TRUE)
  base_sorensen_mean <- mean(as.numeric(base_dist$sorensen), na.rm = TRUE)
  
  # Fixed rare-error probability vectors only change with rate, not iteration.
  rare_probabilities <- setNames(
    lapply(ERROR_RATES, function(rate) make_rare_probability_vector(base_mat, rate)),
    sprintf("%.2f", ERROR_RATES)
  )
  
  assemblage_rows <- list()
  row_index <- 1L
  
  for (iter_i in seq_len(N_SIM)) {
    checkpoint_path <- file.path(
      CHECKPOINT_DIR,
      sprintf("%s__iter_%02d.rds", assemblage_id, iter_i)
    )
    
    if (RESUME_FROM_CHECKPOINTS && file.exists(checkpoint_path)) {
      message(sprintf("  Iteration %02d/%02d: reusing checkpoint", iter_i, N_SIM))
      iter_rows <- readRDS(checkpoint_path)
      assemblage_rows[[row_index]] <- iter_rows
      row_index <- row_index + 1L
      next
    }
    
    observed_map <- make_coherent_map(
      meta, meta,
      seed = 100000L + iter_i
    )
    regional_map <- if (!is.null(regional_pool) && nrow(regional_pool)) {
      make_coherent_map(
        regional_meta, regional_pool,
        seed = 200000L + iter_i
      )
    } else NULL
    coherent_expert_map <- make_expert_map(
      meta, expert_map,
      seed = 300000L + iter_i
    )
    
    mechanisms <- list(
      observed_pool = observed_map,
      rare_weighted_observed_pool = observed_map
    )
    if (!is.null(regional_map)) mechanisms$regional_pool <- regional_map
    if (!is.null(coherent_expert_map) && nrow(coherent_expert_map)) {
      mechanisms$expert_map <- coherent_expert_map
    }
    
    mechanism_indices <- lapply(
      mechanisms,
      map_to_indices,
      taxa = taxa_universe
    )
    
    message(sprintf(
      "  Iteration %02d/%02d: %d mechanism(s) × %d error rates",
      iter_i, N_SIM, length(mechanisms), length(ERROR_RATES)
    ))
    
    iter_rows <- vector(
      "list",
      length(mechanisms) * length(ERROR_RATES)
    )
    k <- 1L
    
    for (mechanism in names(mechanisms)) {
      target_index <- mechanism_indices[[mechanism]]
      
      for (rate in ERROR_RATES) {
        rate_key <- sprintf("%.2f", rate)
        
        p_error <- if (mechanism == "rare_weighted_observed_pool") {
          rare_probabilities[[rate_key]]
        } else {
          # The species baseline is represented only in the observed species
          # columns. Target-only regional taxa receive zero source probability.
          p <- numeric(length(taxa_universe))
          p[match(meta$taxon_unit, taxa_universe)] <- rate
          p
        }
        
        diag <- map_diagnostics_matrix(
          base_mat, target_index, p_error
        )
        
        scenario_mat <- simulate_matrix(
          base_mat = base_mat,
          target_index = target_index,
          p_error = p_error,
          seed = 400000L +
            iter_i * 10000L +
            round(rate * 100) * 20L +
            match(mechanism, names(mechanisms))
        )
        
        scenario_alpha <- alpha_vectors(scenario_mat)
        scenario_gamma <- sum(colSums(scenario_mat) > 0L)
        scenario_dist <- beta_distances(scenario_mat)
        
        iter_rows[[k]] <- tibble(
          assemblage_id = assemblage_id,
          mechanism = mechanism,
          error_rate = rate,
          iter = iter_i,
          eligible_individual_share = diag$eligible_individual_share,
          expected_reassigned_pct = diag$expected_reassigned_pct,
          n_eligible_source_species = diag$n_eligible_source_species,
          n_source_species = diag$n_source_species,
          gamma_change_pct = 100 * (scenario_gamma / base_gamma - 1),
          mean_q0_change_pct = 100 * (
            mean(scenario_alpha$q0) / mean(base_alpha$q0) - 1
          ),
          mean_q1_change_pct = 100 * (
            mean(scenario_alpha$q1) / mean(base_alpha$q1) - 1
          ),
          mean_q2_change_pct = 100 * (
            mean(scenario_alpha$q2) / mean(base_alpha$q2) - 1
          ),
          q0_stability = safe_cor(base_alpha$q0, scenario_alpha$q0),
          bray_stability = if (!is.null(scenario_dist$bray)) {
            safe_cor(
              as.numeric(base_dist$bray),
              as.numeric(scenario_dist$bray)
            )
          } else NA_real_,
          sorensen_stability = if (!is.null(scenario_dist$sorensen)) {
            safe_cor(
              as.numeric(base_dist$sorensen),
              as.numeric(scenario_dist$sorensen)
            )
          } else NA_real_,
          bray_mean_change_pct = if (!is.null(scenario_dist$bray)) {
            100 * (
              mean(as.numeric(scenario_dist$bray), na.rm = TRUE) /
                base_bray_mean - 1
            )
          } else NA_real_,
          sorensen_mean_change_pct = if (!is.null(scenario_dist$sorensen)) {
            100 * (
              mean(as.numeric(scenario_dist$sorensen), na.rm = TRUE) /
                base_sorensen_mean - 1
            )
          } else NA_real_
        )
        
        k <- k + 1L
      }
    }
    
    iter_rows <- bind_rows(iter_rows)
    saveRDS(iter_rows, checkpoint_path)
    
    assemblage_rows[[row_index]] <- iter_rows
    row_index <- row_index + 1L
  }
  
  assemblage_results <- bind_rows(assemblage_rows)
  readr::write_csv(
    assemblage_results,
    file.path(
      OUT_DIR,
      paste0("error_gradient_by_iter__", assemblage_id, ".csv")
    )
  )
  
  all_results[[assemblage_id]] <- assemblage_results
  message("  Completed: ", assemblage_id)
}

gradient_by_iter <- bind_rows(all_results)
gradient_summary <- summarise_gradient(gradient_by_iter)

readr::write_csv(
  gradient_by_iter,
  file.path(OUT_DIR, "error_gradient_by_iter.csv")
)
readr::write_csv(
  gradient_summary,
  file.path(OUT_DIR, "error_gradient_summary.csv")
)

eligibility_summary <- gradient_summary %>%
  select(
    assemblage_id, mechanism, error_rate,
    eligible_individual_share_median,
    eligible_individual_share_p10,
    eligible_individual_share_p90,
    expected_reassigned_pct_median,
    expected_reassigned_pct_p10,
    expected_reassigned_pct_p90,
    n_eligible_source_species_median,
    n_source_species_median
  )

readr::write_csv(
  eligibility_summary,
  file.path(OUT_DIR, "error_gradient_eligibility_audit.csv")
)

# ---- 3. Appendix figures -----------------------------------------------------
assemblage_labels <- c(
  collembola_soil_core = "Collembola\nsoil cores",
  araneae_pitfall = "Araneae\npitfall traps",
  carabidae_pitfall = "Carabidae\npitfall traps",
  formicidae_pitfall = "Formicidae\npitfall traps",
  isopoda_pitfall = "Isopoda\npitfall traps",
  isopoda_hand_sorting = "Isopoda\nhand sorting",
  diplopoda_pitfall = "Diplopoda\npitfall traps",
  diplopoda_hand_sorting = "Diplopoda\nhand sorting"
)

mechanism_labels <- c(
  observed_pool = "Observed-pool congeneric",
  rare_weighted_observed_pool = "Rare-weighted observed-pool",
  regional_pool = "Regional-pool congeneric",
  expert_map = "Expert confusion map"
)

mechanism_colours <- c(
  "Observed-pool congeneric" = "#0072B2",
  "Rare-weighted observed-pool" = "#CC79A7",
  "Regional-pool congeneric" = "#D55E00",
  "Expert confusion map" = "#009E73"
)

theme_paper <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      legend.position = "bottom",
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.y = element_blank()
    )
}

gradient_long_state <- gradient_summary %>%
  select(
    assemblage_id, mechanism, error_rate,
    gamma_change_pct_median, gamma_change_pct_p10, gamma_change_pct_p90,
    bray_mean_change_pct_median, bray_mean_change_pct_p10, bray_mean_change_pct_p90,
    sorensen_mean_change_pct_median, sorensen_mean_change_pct_p10, sorensen_mean_change_pct_p90
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, mechanism, error_rate),
    names_to = c("metric", ".value"),
    names_pattern = "(.*)_(median|p10|p90)"
  ) %>%
  mutate(
    metric = recode(
      metric,
      gamma_change_pct = "Regional richness (γ)",
      bray_mean_change_pct = "Mean Bray–Curtis",
      sorensen_mean_change_pct = "Mean Sørensen"
    ),
    assemblage = recode(assemblage_id, !!!assemblage_labels),
    mechanism = recode(mechanism, !!!mechanism_labels)
  )

p_state <- ggplot(
  gradient_long_state,
  aes(x = 100 * error_rate, y = median, colour = mechanism, fill = mechanism)
) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.3) +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  facet_grid(metric ~ assemblage, scales = "free_y") +
  scale_colour_manual(values = mechanism_colours, name = "Error mechanism") +
  scale_fill_manual(values = mechanism_colours, name = "Error mechanism") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  labs(
    title = "Response of inventory and beta diversity to error rates from 1% to 20%",
    subtitle = "Bands show 10th–90th percentiles; each source→target map is fixed within an iteration across the full gradient.",
    x = "Identification-error rate (%)",
    y = "Change relative to the species baseline (%)"
  ) +
  theme_paper(base_size = 7.5)

ggsave(
  file.path(OUT_DIR, "figures", "FigS_error_gradient_directional_effects.pdf"),
  p_state, width = 250, height = 165, units = "mm", device = grDevices::pdf
)
ggsave(
  file.path(OUT_DIR, "figures", "FigS_error_gradient_directional_effects.png"),
  p_state, width = 250, height = 165, units = "mm", dpi = 450
)

gradient_long_stability <- gradient_summary %>%
  select(
    assemblage_id, mechanism, error_rate,
    q0_stability_median, q0_stability_p10, q0_stability_p90,
    bray_stability_median, bray_stability_p10, bray_stability_p90,
    sorensen_stability_median, sorensen_stability_p10, sorensen_stability_p90
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, mechanism, error_rate),
    names_to = c("metric", ".value"),
    names_pattern = "(.*)_(median|p10|p90)"
  ) %>%
  mutate(
    metric = recode(
      metric,
      q0_stability = "q0 richness stability",
      bray_stability = "Bray–Curtis stability",
      sorensen_stability = "Sørensen stability"
    ),
    assemblage = recode(assemblage_id, !!!assemblage_labels),
    mechanism = recode(mechanism, !!!mechanism_labels)
  )

p_stability <- ggplot(
  gradient_long_stability,
  aes(x = 100 * error_rate, y = median, colour = mechanism, fill = mechanism)
) +
  geom_hline(yintercept = 1, colour = "grey35", linewidth = 0.3) +
  geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  facet_grid(metric ~ assemblage, scales = "free_y") +
  scale_colour_manual(values = mechanism_colours, name = "Error mechanism") +
  scale_fill_manual(values = mechanism_colours, name = "Error mechanism") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  labs(
    title = "Stability of alpha and beta metrics along a 1%–20% identification-error gradient",
    subtitle = "Stability is the Spearman correlation with the species-level baseline.",
    x = "Identification-error rate (%)",
    y = "Stability"
  ) +
  theme_paper(base_size = 7.5)

ggsave(
  file.path(OUT_DIR, "figures", "FigS_error_gradient_stability.pdf"),
  p_stability, width = 250, height = 165, units = "mm", device = grDevices::pdf
)
ggsave(
  file.path(OUT_DIR, "figures", "FigS_error_gradient_stability.png"),
  p_stability, width = 250, height = 165, units = "mm", dpi = 450
)

message("\n========================================================================")
message("Completed appendix error-gradient analysis.")
message("Outputs: ", normalizePath(OUT_DIR))
message("Checkpoints: ", normalizePath(CHECKPOINT_DIR))
