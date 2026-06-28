source("R/utils.R")
source("R/assemblages.R")
settings <- read_settings()
load_pipeline_packages(require_gdm = FALSE)
dir.create(settings$outputs_dir, recursive = TRUE, showWarnings = FALSE)

registry <- read_registry()
records <- read_and_standardise_raw(settings)
if (!nrow(records)) stop("No data left after project/year filtering.")
readr::write_csv(records, file.path(settings$outputs_dir, "standardised_records.csv"))

assemblages <- discover_assemblages(records, registry)
manifest <- export_assemblages(records, assemblages, settings)

message("\nPrepared ", nrow(manifest), " automatically discovered taxon × protocol assemblages.")
message("Review: ", file.path(settings$outputs_dir, "assemblage_manifest.csv"))
