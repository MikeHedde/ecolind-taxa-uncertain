# Error-gradient appendix (v6).
# Extends the v4/v5 gradient with rare-weighted observed-pool errors.
# The source-to-target map is drawn once per iteration and retained across the
# full 1–20% gradient, so curves isolate error intensity.

source("R/utils.R")
source_pipeline_files()
settings <- read_settings()
load_pipeline_packages(require_gdm = FALSE)

out_dir <- file.path(settings$outputs_dir, "error_gradient")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- readr::read_csv(
  file.path(settings$outputs_dir, "assemblage_manifest.csv"),
  show_col_types = FALSE
)

gradient_rows <- list()
k <- 1L

for (i in seq_len(nrow(manifest))) {
  a <- manifest[i, ]
  message("Gradient: ", a$assemblage_id)

  data <- read_assemblage(a$assemblage_id, settings)
  maps <- lookup_maps(data$lookup)

  species_meta <- maps %>%
    filter(!is.na(species_unit), !is.na(genus)) %>%
    transmute(taxon_unit = species_unit, genus) %>%
    distinct()

  if (!nrow(species_meta) || !nrow(data$species)) next

  pool_path <- file.path(settings$regional_pool_dir, paste0(a$assemblage_id, "__regional_pool.csv"))
  audit_path <- file.path(settings$regional_pool_dir, paste0(a$assemblage_id, "__observed_taxref_match_audit.csv"))

  pool <- if (file.exists(pool_path)) {
    readr::read_csv(pool_path, show_col_types = FALSE)
  } else {
    tibble(taxon_unit = character(), cd_ref = character(), genus = character())
  }

  source_audit <- if (file.exists(audit_path)) {
    readr::read_csv(audit_path, show_col_types = FALSE)
  } else {
    tibble(observed_taxon_unit = character(), cd_ref = character())
  }

  base <- make_matrix(data$species, data$station_frame)
  alpha_base <- alpha_metrics(base)
  gamma_base <- community_gamma(base)

  for (iter in seq_len(settings$n_sim)) {

    # Maps are drawn once here and reused across every rate below.
    observed_map <- observed_pool_map(species_meta, seed = 100000 + iter)
    regional_map <- regional_pool_map(species_meta, pool, source_audit, seed = 200000 + iter)

    mechanisms <- list(
      observed_pool = list(map = observed_map, rare_weighted = FALSE),
      rare_weighted_observed_pool = list(map = observed_map, rare_weighted = TRUE),
      regional_pool = list(map = regional_map, rare_weighted = FALSE)
    )

    for (mechanism_id in names(mechanisms)) {
      mechanism <- mechanisms[[mechanism_id]]

      for (rate in settings$error_gradient) {
        probability_map <- NULL

        if (isTRUE(mechanism$rare_weighted)) {
          probability_map <- rare_probabilities(
            data$species,
            rate = rate,
            cap = settings$rare_error_cap
          )
        }

        scenario_long <- simulate_reassignment(
          data$species,
          mechanism$map,
          rate = rate,
          seed = 300000 + iter * 1000 + round(rate * 100) +
            ifelse(mechanism_id == "rare_weighted_observed_pool", 100, 0),
          probabilities = probability_map
        )

        scenario_matrix <- make_matrix(scenario_long, data$station_frame)
        alpha_scenario <- alpha_metrics(scenario_matrix)
        beta <- compare_beta(base, scenario_matrix, settings)

        gradient_rows[[k]] <- bind_cols(
          tibble(
            assemblage_id = a$assemblage_id,
            taxon_key = a$taxon_key,
            mechanism = mechanism_id,
            error_rate = rate,
            iter = iter,
            gamma_change_pct = 100 * (community_gamma(scenario_matrix) / gamma_base - 1),
            mean_q0_change_pct = 100 * (mean(alpha_scenario$q0) / mean(alpha_base$q0) - 1),
            q0_stability = safe_cor(alpha_base$q0, alpha_scenario$q0)
          ),
          beta
        )
        k <- k + 1L
      }
    }
  }
}

by_iter <- bind_rows(gradient_rows)
if (!nrow(by_iter)) stop("No gradient result was generated. Check species-level input matrices.")

metric_columns <- setdiff(
  names(by_iter)[vapply(by_iter, is.numeric, logical(1))],
  c("iter", "error_rate")
)

summary <- by_iter %>%
  group_by(assemblage_id, taxon_key, mechanism, error_rate) %>%
  summarise(
    across(
      all_of(metric_columns),
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
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

readr::write_csv(by_iter, file.path(out_dir, "error_gradient_by_iter.csv"))
readr::write_csv(summary, file.path(out_dir, "error_gradient_summary.csv"))

message("Finished error gradients, including rare-weighted observed-pool errors: ", out_dir)
