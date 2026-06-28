# Creates folders and blank configuration templates without overwriting files.
source("R/utils.R")
source("config/analysis_settings.R")

dirs <- c(
  "data",
  settings$outputs_dir,
  settings$regional_pool_dir,
  settings$expert_input_dir,
  settings$expert_internal_dir,
  "config"
)

purrr::walk(dirs, ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE))
message("Project folders checked/created.")
message("Edit config/analysis_settings.R and config/taxon_registry.csv before running.")
