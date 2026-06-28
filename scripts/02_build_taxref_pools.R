source("R/utils.R")
source("R/taxref.R")
settings <- read_settings()
load_pipeline_packages(require_gdm = FALSE)
dir.create(settings$regional_pool_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(settings$outputs_dir, "assemblage_manifest.csv")
if (!file.exists(manifest_path)) stop("Run scripts/01_prepare_assemblages.R first.")
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
taxref <- make_taxref_objects(read_taxref_flexible(settings$taxref_file), settings)

summary <- purrr::map_dfr(seq_len(nrow(manifest)), function(i) {
  a <- manifest[i,]
  message("TAXREF pool: ", a$assemblage_id)
  out <- build_regional_pool_one(a, taxref, settings)
  readr::write_csv(out$pool, file.path(settings$regional_pool_dir, paste0(a$assemblage_id, "__regional_pool.csv")))
  readr::write_csv(out$audit, file.path(settings$regional_pool_dir, paste0(a$assemblage_id, "__observed_taxref_match_audit.csv")))
  tibble(
    assemblage_id = a$assemblage_id,
    taxon_key = a$taxon_key,
    n_regional_candidates = nrow(out$pool),
    n_observed_species = nrow(out$audit),
    match_rate = mean(out$audit$in_regional_pool, na.rm = TRUE),
    n_sources_with_alternative_candidate = sum(out$audit$n_alternative_candidates > 0, na.rm = TRUE)
  )
})
readr::write_csv(summary,file.path(settings$regional_pool_dir,"regional_pool_build_summary.csv"))
message("Finished TAXREF pools.")
