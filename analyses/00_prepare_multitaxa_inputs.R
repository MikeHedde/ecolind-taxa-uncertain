# ============================================================================
# 00_prepare_multitaxa_inputs_2024.R
# Empirical multi-taxon foundation for the taxonomic-uncertainty paper
#
# Purpose
#   1. Restrict the database to RMQS 2024.
#   2. Define independent assemblages as taxon group × sampling protocol.
#   3. Audit sampling coverage, taxonomic resolution and life-stage information.
#   4. Export clean station-level community inputs for the later scenario engine.
#
# Important design choice
#   Each assemblage is analysed separately. Never merge protocols or taxa into
#   one community matrix: their capture processes and ecological compartments
#   differ fundamentally.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# -----------------------------------------------------------------------------
# 0. USER SETTINGS
# -----------------------------------------------------------------------------

INPUT_FILE <- "data/raw-data/1.faune/all_clean_sp.csv"  # change if needed
OUT_DIR <- "outputs_multitaxa_2024"

TARGET_PROJECT <- "RMQS"
TARGET_YEAR <- 2024L

# A taxon/protocol pair is considered suitable for a main analysis only after
# reviewing the audit. These thresholds are descriptive flags, not filters.
MIN_STATIONS_FLAG <- 20L
MIN_TOTAL_ABUNDANCE_FLAG <- 20L

# Taxonomic ranks retained in ecological matrices. Records at order/class level
# remain in audits, but are not meaningful community units for this paper.
SPECIES_RANKS <- c("S", "E", "SE", "sE")
GENUS_RANKS   <- c("G", "sG")
FAMILY_RANKS  <- c("F", "sF")
ORDER_RANKS   <- c("O", "sO")

# Independent assemblages used in the paper.
# `stage_scenario_candidate = TRUE` means that an adult-only sensitivity
# dataset will be exported. It remains an optional scenario, not a default.
ASSEMBLAGES <- tribble(
  ~assemblage_id,                 ~group_label,  ~selector_col, ~selector_value, ~method,          ~role,        ~stage_scenario_candidate,
  "collembola_soil_core",         "Collembola", "order",       "Collembola",    "soil-core",     "primary",     FALSE,
  "araneae_pitfall",              "Araneae",    "order",       "Araneae",       "pitfall-trap",  "primary",     TRUE,
  "formicidae_pitfall",           "Formicidae", "family",      "Formicidae",    "pitfall-trap",  "primary",     FALSE,
  "carabidae_pitfall",            "Carabidae",  "family",      "Carabidae",     "pitfall-trap",  "primary",     FALSE,
  "isopoda_pitfall",              "Isopoda",    "order",       "Isopoda",       "pitfall-trap",  "secondary",   FALSE,
  "isopoda_hand_sorting",         "Isopoda",    "order",       "Isopoda",       "hand-sorting",  "secondary",   FALSE,
  "diplopoda_pitfall",            "Diplopoda",  "class",       "Diplopoda",     "pitfall-trap",  "secondary",   TRUE,
  "diplopoda_hand_sorting",       "Diplopoda",  "class",       "Diplopoda",     "hand-sorting",  "secondary",   TRUE
)

# -----------------------------------------------------------------------------
# 1. HELPERS
# -----------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "assemblages"), showWarnings = FALSE, recursive = TRUE)

clean_text <- function(x) {
  x %>% as.character() %>% stringr::str_squish() %>% na_if("")
}

parse_count <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x), na = c("", "NA", "NaN")))
}

extract_binomial <- function(x) {
  x <- clean_text(x)
  # Genus + lower-case epithet; authorship and qualifiers are deliberately ignored.
  out <- stringr::str_extract(x, "^[A-Z][A-Za-zÀ-ÖØ-öø-ÿ-]+\\s+[a-z][A-Za-zÀ-ÖØ-öø-ÿ-]+")
  clean_text(out)
}

extract_genus <- function(x) {
  x <- clean_text(x)
  out <- stringr::str_extract(x, "^[A-Z][A-Za-zÀ-ÖØ-öø-ÿ-]+")
  clean_text(out)
}

rank_to_level <- function(x) {
  dplyr::case_when(
    x %in% SPECIES_RANKS ~ "species",
    x %in% GENUS_RANKS ~ "genus",
    x %in% FAMILY_RANKS ~ "family",
    x %in% ORDER_RANKS ~ "order",
    TRUE ~ "other"
  )
}

scenario_eligibility <- function(n_stations, total_abundance, n_taxa) {
  dplyr::case_when(
    n_stations >= 50 & total_abundance >= 500 & n_taxa >= 15 ~ "main_analysis_candidate",
    n_stations >= MIN_STATIONS_FLAG & total_abundance >= MIN_TOTAL_ABUNDANCE_FLAG & n_taxa >= 8 ~ "secondary_or_sensitivity_candidate",
    TRUE ~ "audit_only_or_insufficient"
  )
}

# Generate a station-level matrix only when later needed. This helper preserves
# sampled stations with zero records, which is essential for a monitoring design.
make_station_matrix <- function(comm_long, station_frame) {
  taxa <- sort(unique(comm_long$taxon_unit))
  if (length(taxa) == 0L) {
    return(matrix(numeric(0), nrow = nrow(station_frame), ncol = 0,
                  dimnames = list(station_frame$station, character(0))))
  }
  
  comm_long %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    tidyr::complete(station = station_frame$station, taxon_unit = taxa, fill = list(abundance = 0)) %>%
    tidyr::pivot_wider(names_from = taxon_unit, values_from = abundance, values_fill = 0) %>%
    arrange(match(station, station_frame$station)) %>%
    tibble::column_to_rownames("station") %>%
    as.matrix()
}

# -----------------------------------------------------------------------------
# 2. IMPORT AND STANDARDISATION
# -----------------------------------------------------------------------------

message("Reading input file: ", INPUT_FILE)

raw <- readr::read_delim(
  INPUT_FILE,
  delim = ";",
  locale = readr::locale(encoding = "Latin1"),
  show_col_types = FALSE,
  name_repair = "unique"
) %>%
  janitor::clean_names()

required_cols <- c(
  "code_ech", "project", "yr", "station", "method", "phylum", "class",
  "order", "family", "cd_nom", "name", "valid_name", "rank", "tot_abund"
)
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

rmqs_2024 <- raw %>%
  mutate(
    across(c(code_ech, project, station, method, phylum, subphylum, class, order,
             family, cd_nom, name, valid_name, rank, remarque_individu), clean_text),
    yr = suppressWarnings(as.integer(yr)),
    across(any_of(c("male", "female", "juv", "indet", "tot_abund")), parse_count),
    abundance = coalesce(tot_abund, 0),
    rank_raw = rank,
    taxo_level = rank_to_level(rank_raw),
    binomial = coalesce(extract_binomial(valid_name), extract_binomial(name)),
    genus = case_when(
      taxo_level == "species" ~ extract_genus(binomial),
      taxo_level == "genus" ~ coalesce(extract_genus(name), extract_genus(valid_name)),
      TRUE ~ NA_character_
    ),
    species_key = case_when(
      taxo_level == "species" & !is.na(cd_nom) ~ paste0("cdnom:", cd_nom),
      taxo_level == "species" & !is.na(binomial) ~ paste0("name:", binomial),
      TRUE ~ NA_character_
    ),
    taxon_unit_rtu = case_when(
      taxo_level == "species" & !is.na(species_key) ~ paste0("species:", species_key),
      taxo_level == "genus" & !is.na(genus) ~ paste0("genus:", genus),
      taxo_level == "family" & !is.na(family) ~ paste0("family:", family),
      TRUE ~ NA_character_
    ),
    taxon_unit_species = if_else(
      taxo_level == "species" & !is.na(species_key),
      paste0("species:", species_key),
      NA_character_
    ),
    taxon_unit_genus = if_else(
      taxo_level %in% c("species", "genus") & !is.na(genus),
      paste0("genus:", genus),
      NA_character_
    ),
    taxon_unit_family = if_else(
      taxo_level %in% c("species", "genus", "family") & !is.na(family),
      paste0("family:", family),
      NA_character_
    ),
    n_stage_fields = rowSums(!is.na(cbind(male, female, juv, indet))),
    stage_sum = rowSums(cbind(coalesce(male, 0), coalesce(female, 0), coalesce(juv, 0), coalesce(indet, 0))),
    adults_definite = coalesce(male, 0) + coalesce(female, 0),
    adults_plus_indet = coalesce(male, 0) + coalesce(female, 0) + coalesce(indet, 0),
    juveniles_explicit = coalesce(juv, 0),
    abundance_not_stage_resolved = case_when(
      n_stage_fields == 0 ~ abundance,
      TRUE ~ pmax(abundance - stage_sum, 0)
    ),
    stage_accounting_difference = abundance - stage_sum
  ) %>%
  filter(project == TARGET_PROJECT, yr == TARGET_YEAR)

if (nrow(rmqs_2024) == 0L) {
  stop("No records found after filtering project = '", TARGET_PROJECT,
       "' and year = ", TARGET_YEAR, ". Check project/year coding.")
}

readr::write_csv(rmqs_2024, file.path(OUT_DIR, "rmqs_2024_standardised_records.csv"))

# Frame of all sampled stations per method. It includes empty samples when they
# are represented in the raw database and is the reference for zero communities.
sampling_frame <- rmqs_2024 %>%
  distinct(code_ech, station, method) %>%
  filter(!is.na(code_ech), !is.na(station), !is.na(method))

readr::write_csv(sampling_frame, file.path(OUT_DIR, "rmqs_2024_sampling_frame.csv"))

# -----------------------------------------------------------------------------
# 3. PER-ASSEMBLAGE AUDIT AND EXPORT
# -----------------------------------------------------------------------------

assemblage_manifest <- vector("list", nrow(ASSEMBLAGES))
taxonomic_audits <- vector("list", nrow(ASSEMBLAGES))
stage_audits <- vector("list", nrow(ASSEMBLAGES))
station_audits <- vector("list", nrow(ASSEMBLAGES))

for (i in seq_len(nrow(ASSEMBLAGES))) {
  def <- ASSEMBLAGES[i, ]
  message("Processing: ", def$assemblage_id)
  
  station_frame <- sampling_frame %>%
    filter(method == def$method) %>%
    distinct(station) %>%
    arrange(station)
  
  x <- rmqs_2024 %>%
    filter(method == def$method, .data[[def$selector_col]] == def$selector_value) %>%
    filter(abundance > 0)
  
  # Ecological matrices only retain species/genus/family units.
  x_usable <- x %>%
    filter(taxo_level %in% c("species", "genus", "family"), !is.na(taxon_unit_rtu))
  
  # Audit of taxonomic structure, including unusable/coarse ranks.
  tax_audit <- x %>%
    mutate(
      resolution = case_when(
        taxo_level == "species" ~ "species",
        taxo_level == "genus" ~ "genus",
        taxo_level == "family" ~ "family",
        taxo_level == "order" ~ "order",
        TRUE ~ "other_or_unusable"
      )
    ) %>%
    group_by(assemblage_id = def$assemblage_id, group_label = def$group_label,
             method = def$method, resolution) %>%
    summarise(
      n_records = n(),
      total_abundance = sum(abundance),
      n_reported_names = n_distinct(coalesce(valid_name, name), na.rm = TRUE),
      n_stations_positive = n_distinct(station),
      .groups = "drop"
    ) %>%
    mutate(
      abundance_share = total_abundance / sum(total_abundance),
      record_share = n_records / sum(n_records)
    )
  
  # Station-level community vectors. Sampling stations with no occurrence are
  # retained separately in station_frame and must be reinstated as zeros later.
  comm_rtu <- x_usable %>%
    transmute(station, taxon_unit = taxon_unit_rtu, abundance, taxo_level, genus, family,
              species_key, binomial, valid_name, name) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
  
  comm_species <- x_usable %>%
    filter(taxo_level == "species", !is.na(taxon_unit_species)) %>%
    transmute(station, taxon_unit = taxon_unit_species, abundance) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
  
  comm_genus <- x_usable %>%
    filter(!is.na(taxon_unit_genus)) %>%
    transmute(station, taxon_unit = taxon_unit_genus, abundance) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
  
  comm_family <- x_usable %>%
    filter(!is.na(taxon_unit_family)) %>%
    transmute(station, taxon_unit = taxon_unit_family, abundance) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop")
  
  taxon_lookup <- x_usable %>%
    transmute(
      taxon_unit_rtu,
      taxon_unit_species,
      taxon_unit_genus,
      taxon_unit_family,
      taxo_level,
      species_key,
      binomial,
      genus,
      family,
      order,
      class,
      valid_name,
      name
    ) %>%
    distinct()
  
  # Life-stage audit. `indet` is intentionally kept separate: in Formicidae it
  # often represents workers, while in other groups it may mean unknown stage.
  stage_audit <- x %>%
    summarise(
      assemblage_id = def$assemblage_id,
      group_label = def$group_label,
      method = def$method,
      total_abundance = sum(abundance),
      male = sum(coalesce(male, 0)),
      female = sum(coalesce(female, 0)),
      juvenile_explicit = sum(coalesce(juv, 0)),
      indeterminate_stage_or_sex = sum(coalesce(indet, 0)),
      abundance_not_stage_resolved = sum(abundance_not_stage_resolved),
      stage_accounting_difference = sum(stage_accounting_difference),
      juvenile_share = juvenile_explicit / total_abundance,
      stage_candidate = def$stage_scenario_candidate
    )
  
  # Adult-only strict sensitivity input where biologically relevant. It only
  # retains explicitly sexed adults; unknown-stage individuals are excluded.
  if (isTRUE(def$stage_scenario_candidate)) {
    comm_adults_strict <- x_usable %>%
      mutate(adult_abundance = adults_definite) %>%
      filter(adult_abundance > 0) %>%
      transmute(station, taxon_unit = taxon_unit_rtu, abundance = adult_abundance) %>%
      group_by(station, taxon_unit) %>%
      summarise(abundance = sum(abundance), .groups = "drop")
    
    readr::write_csv(
      comm_adults_strict,
      file.path(OUT_DIR, "assemblages", paste0(def$assemblage_id, "__adult_only_strict_rtu_long.csv"))
    )
  }
  
  station_audit <- station_frame %>%
    left_join(
      comm_rtu %>%
        group_by(station) %>%
        summarise(total_abundance = sum(abundance), rtu_richness = n_distinct(taxon_unit), .groups = "drop"),
      by = "station"
    ) %>%
    mutate(
      assemblage_id = def$assemblage_id,
      group_label = def$group_label,
      method = def$method,
      total_abundance = replace_na(total_abundance, 0),
      rtu_richness = replace_na(rtu_richness, 0L),
      detected = total_abundance > 0
    )
  
  summary_one <- tibble(
    assemblage_id = def$assemblage_id,
    group_label = def$group_label,
    selector = paste0(def$selector_col, " = ", def$selector_value),
    method = def$method,
    role = def$role,
    n_sampled_stations = nrow(station_frame),
    n_positive_stations = sum(station_audit$detected),
    prevalence = mean(station_audit$detected),
    total_abundance = sum(comm_rtu$abundance),
    n_rtu = n_distinct(comm_rtu$taxon_unit),
    n_species = n_distinct(comm_species$taxon_unit),
    n_genera = n_distinct(comm_genus$taxon_unit),
    n_families = n_distinct(comm_family$taxon_unit),
    share_species_abundance = tax_audit %>% filter(resolution == "species") %>% pull(abundance_share) %>% dplyr::first(default = 0),
    share_genus_abundance = tax_audit %>% filter(resolution == "genus") %>% pull(abundance_share) %>% dplyr::first(default = 0),
    share_family_or_coarser_abundance = tax_audit %>% filter(resolution %in% c("family", "order", "other_or_unusable")) %>% summarise(x = sum(abundance_share)) %>% pull(x),
    juvenile_share = stage_audit$juvenile_share,
    stage_scenario_candidate = def$stage_scenario_candidate,
    family_coarsening_applicable = n_distinct(x_usable$family, na.rm = TRUE) >= 2,
    analysis_flag = scenario_eligibility(nrow(station_frame), sum(comm_rtu$abundance), n_distinct(comm_rtu$taxon_unit))
  )
  
  assemblage_dir <- file.path(OUT_DIR, "assemblages")
  prefix <- file.path(assemblage_dir, def$assemblage_id)
  
  readr::write_csv(station_frame %>% mutate(assemblage_id = def$assemblage_id), paste0(prefix, "__station_frame.csv"))
  readr::write_csv(comm_rtu, paste0(prefix, "__rtu_long.csv"))
  readr::write_csv(comm_species, paste0(prefix, "__species_long.csv"))
  readr::write_csv(comm_genus, paste0(prefix, "__genus_long.csv"))
  readr::write_csv(comm_family, paste0(prefix, "__family_long.csv"))
  readr::write_csv(taxon_lookup, paste0(prefix, "__taxon_lookup.csv"))
  readr::write_csv(tax_audit, paste0(prefix, "__taxonomic_resolution_audit.csv"))
  readr::write_csv(stage_audit, paste0(prefix, "__stage_audit.csv"))
  readr::write_csv(station_audit, paste0(prefix, "__station_audit.csv"))
  
  # Optional wide matrices are convenient for manual inspection and validation.
  saveRDS(make_station_matrix(comm_rtu, station_frame), paste0(prefix, "__rtu_matrix.rds"))
  
  assemblage_manifest[[i]] <- summary_one
  taxonomic_audits[[i]] <- tax_audit
  stage_audits[[i]] <- stage_audit
  station_audits[[i]] <- station_audit
}

# -----------------------------------------------------------------------------
# 4. CROSS-ASSEMBLAGE OUTPUTS
# -----------------------------------------------------------------------------

manifest <- bind_rows(assemblage_manifest) %>%
  arrange(match(role, c("primary", "secondary")), desc(n_sampled_stations), group_label, method)

taxonomic_audit_all <- bind_rows(taxonomic_audits)
stage_audit_all <- bind_rows(stage_audits)
station_audit_all <- bind_rows(station_audits)

readr::write_csv(ASSEMBLAGES, file.path(OUT_DIR, "assemblage_definitions.csv"))
readr::write_csv(manifest, file.path(OUT_DIR, "assemblage_manifest.csv"))
readr::write_csv(taxonomic_audit_all, file.path(OUT_DIR, "taxonomic_resolution_audit_all.csv"))
readr::write_csv(stage_audit_all, file.path(OUT_DIR, "stage_audit_all.csv"))
readr::write_csv(station_audit_all, file.path(OUT_DIR, "station_audit_all.csv"))

# Diagnostic figure: only an audit, not intended as a manuscript figure.
p_audit <- taxonomic_audit_all %>%
  mutate(
    assemblage_id = forcats::fct_reorder(assemblage_id, total_abundance, .fun = sum),
    resolution = factor(resolution, levels = c("species", "genus", "family", "order", "other_or_unusable"))
  ) %>%
  ggplot(aes(x = assemblage_id, y = abundance_share, fill = resolution)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "Taxonomic-resolution audit across RMQS 2024 assemblage × protocol datasets",
    x = NULL,
    y = "Share of total abundance",
    fill = "Reported resolution"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(OUT_DIR, "audit_taxonomic_resolution_2024.png"), p_audit, width = 10, height = 6, dpi = 300)

ggsave(file.path(OUT_DIR, "audit_taxonomic_resolution_2024.pdf"), p_audit, width = 10, height = 6)

message("\nDone. Key files:")
message("  - ", file.path(OUT_DIR, "assemblage_manifest.csv"))
message("  - ", file.path(OUT_DIR, "taxonomic_resolution_audit_all.csv"))
message("  - ", file.path(OUT_DIR, "stage_audit_all.csv"))
message("  - ", file.path(OUT_DIR, "assemblages/"))
message("\nReview the manifest before launching uncertainty scenarios.")
