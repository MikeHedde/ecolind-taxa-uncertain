# Configuration shared by all scripts.
# Keep paths relative to the project root.

settings <- list(
  raw_input_file = "data/all_clean_sp.csv",
  env_input_file = "data/all_env_variables.csv",
  taxref_file = "data/TAXREFv17.txt",

  target_project = "RMQS",
  target_year = 2024L,

  outputs_dir = "outputs",
  regional_pool_dir = "regional_pools",
  expert_input_dir = "expert_input",
  expert_internal_dir = "expert_internal",

  # Analysis thresholds
  min_beta_sites = 8L,
  min_gdm_sites = 10L,
  min_main_stations = 30L,
  min_main_abundance = 100L,
  min_main_taxa = 10L,

  # Replication
  n_sim = 50L,
  gdm_max_stochastic_iters = 12L,

  # Taxonomic scenarios
  default_error_rate = 0.10,
  rare_error_cap = 0.25,
  rare_max_total_abundance = 3L,
  error_gradient = seq(0.01, 0.20, by = 0.01),

  # Environmental turnover model
  run_gdm = TRUE,
  env_predictors = c("t360_mean", "mos", "p_h"),

  # TAXREF
  taxref_fr_status_keep = c("P", "E", "C"),
  taxref_species_ranks = c("ES", "S"),

  # Outputs
  write_example_confusion_maps = TRUE,
  write_all_confusion_maps = FALSE
)
