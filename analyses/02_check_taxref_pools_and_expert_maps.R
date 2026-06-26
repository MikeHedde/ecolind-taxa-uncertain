# =============================================================================
# 05b_check_taxref_pools_and_expert_maps.R
# Pre-flight validation for multi-taxon TAXREF regional pools and curated
# expert-confusion maps. Run after 05_build_multitaxa_taxref_pools_v2.R and
# before the main v3 scenario analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

INPUT_DIR <- "outputs_multitaxa_2024"
REGIONAL_POOL_DIR <- "regional_pools"
EXPERT_MAP_DIR <- "expert_confusion_maps"
OUT_FILE <- file.path(REGIONAL_POOL_DIR, "regional_pool_preflight_check.csv")

manifest_path <- file.path(INPUT_DIR, "assemblage_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("Missing assemblage manifest: ", manifest_path)
}

manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)

read_lookup_species <- function(id) {
  path <- file.path(INPUT_DIR, "assemblages", paste0(id, "__taxon_lookup.csv"))
  readr::read_csv(path, show_col_types = FALSE) %>%
    transmute(taxon_unit = taxon_unit_species, genus = genus) %>%
    filter(!is.na(taxon_unit), !is.na(genus)) %>%
    distinct()
}

checks <- purrr::map_dfr(manifest$assemblage_id, function(id) {
  observed <- read_lookup_species(id)

  pool_path <- file.path(REGIONAL_POOL_DIR, paste0(id, "__regional_pool.csv"))
  audit_path <- file.path(REGIONAL_POOL_DIR, paste0(id, "__observed_taxref_match_audit.csv"))
  expert_path <- file.path(EXPERT_MAP_DIR, paste0(id, "__expert_confusions.csv"))

  if (!file.exists(pool_path) || !file.exists(audit_path)) {
    return(tibble(
      assemblage_id = id,
      pool_available = FALSE,
      n_observed_species = nrow(observed),
      observed_match_rate = NA_real_,
      share_sources_with_regional_alternative = NA_real_,
      n_regional_candidates = NA_integer_,
      expert_map_available = file.exists(expert_path),
      n_enabled_expert_pairs = NA_integer_,
      invalid_expert_sources = NA_integer_,
      invalid_expert_targets = NA_integer_,
      status = "MISSING_POOL_OR_AUDIT"
    ))
  }

  pool <- readr::read_csv(pool_path, show_col_types = FALSE) %>% janitor::clean_names()
  audit <- readr::read_csv(audit_path, show_col_types = FALSE) %>% janitor::clean_names()

  match_rate <- if ("in_regional_pool" %in% names(audit)) {
    mean(as.logical(audit$in_regional_pool), na.rm = TRUE)
  } else NA_real_

  alt_share <- if ("has_alternative_candidate" %in% names(audit)) {
    mean(as.logical(audit$has_alternative_candidate), na.rm = TRUE)
  } else NA_real_

  expert_summary <- tibble(
    expert_map_available = file.exists(expert_path),
    n_enabled_expert_pairs = 0L,
    invalid_expert_sources = 0L,
    invalid_expert_targets = 0L
  )

  if (file.exists(expert_path)) {
    exp_map <- readr::read_csv(expert_path, show_col_types = FALSE) %>% janitor::clean_names()
    if (all(c("source_taxon_unit", "target_taxon_unit") %in% names(exp_map))) {
      if (!"enabled" %in% names(exp_map)) exp_map$enabled <- TRUE
      exp_map <- exp_map %>%
        mutate(enabled = tolower(as.character(enabled)) %in% c("true", "t", "1", "yes", "y")) %>%
        filter(enabled)

      expert_summary <- tibble(
        expert_map_available = TRUE,
        n_enabled_expert_pairs = nrow(exp_map),
        invalid_expert_sources = sum(!exp_map$source_taxon_unit %in% observed$taxon_unit),
        invalid_expert_targets = sum(!exp_map$target_taxon_unit %in% pool$taxon_unit)
      )
    }
  }

  status <- case_when(
    !is.finite(match_rate) ~ "AUDIT_SCHEMA_UNRECOGNISED",
    match_rate < 0.8 ~ "REVIEW_MATCH_RATE",
    expert_summary$invalid_expert_sources > 0L ~ "INVALID_EXPERT_SOURCE",
    expert_summary$invalid_expert_targets > 0L ~ "INVALID_EXPERT_TARGET",
    TRUE ~ "OK"
  )

  tibble(
    assemblage_id = id,
    pool_available = TRUE,
    n_observed_species = nrow(observed),
    observed_match_rate = match_rate,
    share_sources_with_regional_alternative = alt_share,
    n_regional_candidates = nrow(pool),
    expert_map_available = expert_summary$expert_map_available,
    n_enabled_expert_pairs = expert_summary$n_enabled_expert_pairs,
    invalid_expert_sources = expert_summary$invalid_expert_sources,
    invalid_expert_targets = expert_summary$invalid_expert_targets,
    status = status
  )
})

readr::write_csv(checks, OUT_FILE)
print(checks)
message("\nSaved: ", OUT_FILE)
message("Resolve every status other than OK before interpreting regional or expert scenarios.")
