# =============================================================================
# 05_build_multitaxa_taxref_pools_v2.R
# France-mainland TAXREF candidate pools and curated-expert templates
# for the RMQS 2024 multi-taxon taxonomic-uncertainty workflow.
#
# Purpose
# -------
# 1. Build one accepted-species candidate pool per assemblage × protocol.
# 2. Preserve the original RMQS species identifiers for taxa already observed.
# 3. Use accepted TAXREF CD_REF identifiers only for non-observed candidates.
# 4. Produce match audits and templates for biologically justified expert maps.
#
# Important
# ---------
# The national TAXREF pool is a deliberately broad taxonomic candidate pool.
# It is not interpreted as the local ecological species pool. The corresponding
# scenario is a bounded sensitivity analysis of potential taxonomic confusion.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# ---- 0. User settings --------------------------------------------------------
INPUT_DIR <- "outputs_multitaxa_2024"

# Point to the complete TAXREF export from INPN/MNHN.
# Both semicolon- and tab-delimited exports are supported.
TAXREF_FILE <- "data/raw-data/TAXREFv17.txt"

REGIONAL_POOL_DIR <- "regional_pools"
EXPERT_MAP_DIR <- "expert_confusion_maps"

# Include accepted / introduced / cryptogenic records in mainland France.
FR_STATUS_KEEP <- c("P", "E", "C")

# TAXREF uses ES for species in current exports. S is retained for compatibility
# with older exports.
SPECIES_RANKS <- c("ES", "S")

# Assemblages and TAXREF filters. Each row specifies the taxonomic scope of the
# candidate pool. Assemblage = group × protocol, even when the same taxon occurs
# under two protocols.
ASSEMBLAGE_FILTERS <- tribble(
  ~assemblage_id,              ~filter_rank, ~filter_value,
  "collembola_soil_core",      "ordre",      "Collembola",
  "araneae_pitfall",           "ordre",      "Araneae",
  "carabidae_pitfall",         "famille",    "Carabidae",
  "formicidae_pitfall",        "famille",    "Formicidae",
  "isopoda_pitfall",           "ordre",      "Isopoda",
  "isopoda_hand_sorting",      "ordre",      "Isopoda",
  "diplopoda_pitfall",         "classe",     "Diplopoda",
  "diplopoda_hand_sorting",    "classe",     "Diplopoda"
)

# Known Collembola groups for which expert confusion is plausible from the
# consultation already obtained. These rows create CANDIDATE templates only:
# no expert scenario is run until the user explicitly sets enabled = TRUE for
# a source-target pair in the final map.
#
# Cross-genus confusion is deliberately not generated automatically. Add an
# explicit source_genus / target_genus row here only when it is justified by an
# expert and an actual source-target pair can be reviewed.
EXPERT_GENUS_RULES <- tribble(
  ~assemblage_id,         ~source_genus,    ~target_genus,      ~review_note,
  "collembola_soil_core", "Folsomia",        "Folsomia",         "Difficult genus; expert review required",
  "collembola_soil_core", "Mesaphorura",     "Mesaphorura",      "Difficult genus; expert review required",
  "collembola_soil_core", "Isotoma",         "Isotoma",          "Difficult genus; expert review required",
  "collembola_soil_core", "Lepidocyrtus",    "Lepidocyrtus",     "Species complex; expert review required",
  "collembola_soil_core", "Entomobrya",      "Entomobrya",       "Species complex; expert review required",
  "collembola_soil_core", "Ceratophysella",  "Ceratophysella",   "Species complex; expert review required",
  "collembola_soil_core", "Sminthurides",    "Sminthurides",     "Species complex; expert review required"
)

# Optional manual overrides file. It is useful when an observed RMQS taxon is
# represented in TAXREF under a synonym or spelling not resolved automatically.
# Required columns if used:
#   assemblage_id, observed_taxon_unit, cd_ref, note
MANUAL_OVERRIDE_FILE <- file.path(REGIONAL_POOL_DIR, "observed_taxref_manual_overrides.csv")

# ---- 1. Generic helpers ------------------------------------------------------
clean_chr <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_squish() %>%
    na_if("")
}

normalise_taxon <- function(x) {
  x %>%
    clean_chr() %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    stringr::str_to_lower()
}

extract_binomial <- function(x) {
  x <- clean_chr(x)
  clean_chr(stringr::str_extract(x, "^[A-Z][A-Za-z-]+\\s+[a-z][A-Za-z-]+"))
}

extract_genus <- function(x) {
  x <- clean_chr(x)
  clean_chr(stringr::str_extract(x, "^[A-Z][A-Za-z-]+"))
}

first_existing_column <- function(data, candidates, required = TRUE) {
  hit <- intersect(candidates, names(data))
  if (length(hit) == 0L) {
    if (required) {
      stop(
        "Missing TAXREF column. Tried: ",
        paste(candidates, collapse = ", "),
        "\nAvailable columns: ", paste(names(data), collapse = ", ")
      )
    }
    return(rep(NA_character_, nrow(data)))
  }
  data[[hit[[1]]]]
}

read_taxref_flexible <- function(path) {
  if (!file.exists(path)) {
    stop(
      "TAXREF file not found: ", path, "\n",
      "Set TAXREF_FILE at the top of this script."
    )
  }

  first_line <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  delim <- if (stringr::str_count(first_line, ";") >= stringr::str_count(first_line, "\t")) ";" else "\t"

  readr::read_delim(
    path,
    delim = delim,
    locale = readr::locale(encoding = "Latin1"),
    show_col_types = FALSE,
    name_repair = "unique"
  ) %>%
    janitor::clean_names()
}

read_manual_overrides <- function(path) {
  if (!file.exists(path)) {
    return(tibble(
      assemblage_id = character(),
      observed_taxon_unit = character(),
      cd_ref_override = character(),
      note = character()
    ))
  }

  x <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  needed <- c("assemblage_id", "observed_taxon_unit", "cd_ref")
  if (!all(needed %in% names(x))) {
    stop(
      "Manual override file exists but is missing: ",
      paste(setdiff(needed, names(x)), collapse = ", ")
    )
  }

  note_vec <- if ("note" %in% names(x)) as.character(x$note) else rep(NA_character_, nrow(x))

  x %>%
    transmute(
      assemblage_id = as.character(assemblage_id),
      observed_taxon_unit = as.character(observed_taxon_unit),
      cd_ref_override = as.character(cd_ref),
      note = note_vec
    ) %>%
    distinct()
}

# ---- 2. TAXREF standardisation ----------------------------------------------
# We retain two linked objects:
# - taxref_alias: all CD_NOM -> accepted CD_REF aliases, including synonyms;
# - taxref_accepted: one accepted candidate record per CD_REF.
make_taxref_objects <- function(taxref_raw) {
  raw <- tibble(
    cd_nom = clean_chr(first_existing_column(taxref_raw, c("cd_nom"))),
    cd_ref = clean_chr(first_existing_column(taxref_raw, c("cd_ref"))),
    rang = toupper(clean_chr(first_existing_column(taxref_raw, c("rang", "rank")))),
    lb_nom = clean_chr(first_existing_column(taxref_raw, c("lb_nom", "nom_complet"), required = FALSE)),
    nom_valide = clean_chr(first_existing_column(taxref_raw, c("nom_valide"), required = FALSE)),
    classe = clean_chr(first_existing_column(taxref_raw, c("classe", "class"), required = FALSE)),
    ordre = clean_chr(first_existing_column(taxref_raw, c("ordre", "order"), required = FALSE)),
    famille = clean_chr(first_existing_column(taxref_raw, c("famille", "family"), required = FALSE)),
    fr = toupper(clean_chr(first_existing_column(taxref_raw, c("fr"), required = FALSE)))
  ) %>%
    mutate(
      cd_ref = coalesce(cd_ref, cd_nom),
      label_raw = coalesce(nom_valide, lb_nom),
      accepted_name = coalesce(extract_binomial(nom_valide), extract_binomial(lb_nom)),
      genus = extract_genus(accepted_name),
      is_accepted = !is.na(cd_nom) & !is.na(cd_ref) & cd_nom == cd_ref
    )

  # Candidate accepted species: continental France + species rank.
  accepted <- raw %>%
    filter(
      rang %in% SPECIES_RANKS,
      fr %in% FR_STATUS_KEEP,
      !is.na(cd_ref),
      !is.na(accepted_name),
      !is.na(genus)
    ) %>%
    arrange(desc(is_accepted)) %>%
    distinct(cd_ref, .keep_all = TRUE) %>%
    select(
      cd_ref, accepted_name, genus, famille, classe, ordre, fr
    )

  # Alias mapping includes every TAXREF row that resolves to an accepted
  # candidate CD_REF, not only the accepted row. This is critical for RMQS
  # records stored under synonyms.
  alias <- raw %>%
    filter(!is.na(cd_nom), !is.na(cd_ref)) %>%
    inner_join(accepted %>% select(cd_ref), by = "cd_ref") %>%
    transmute(
      cd_nom,
      cd_ref,
      alias_name = coalesce(extract_binomial(lb_nom), extract_binomial(nom_valide))
    ) %>%
    distinct()

  list(accepted = accepted, alias = alias)
}

filter_taxref_scope <- function(accepted, assemblage_id) {
  rule <- ASSEMBLAGE_FILTERS %>% filter(.data$assemblage_id == !!assemblage_id)
  if (nrow(rule) != 1L) stop("No unique TAXREF scope rule for: ", assemblage_id)

  rank_col <- rule$filter_rank[[1]]
  rank_value <- normalise_taxon(rule$filter_value[[1]])

  if (!rank_col %in% names(accepted)) {
    stop("TAXREF accepted table does not contain: ", rank_col)
  }

  accepted %>%
    filter(normalise_taxon(.data[[rank_col]]) == rank_value)
}

# ---- 3. Observed RMQS species ------------------------------------------------
read_observed_species <- function(assemblage_id) {
  path <- file.path(INPUT_DIR, "assemblages", paste0(assemblage_id, "__taxon_lookup.csv"))
  if (!file.exists(path)) stop("Missing taxon lookup: ", path)

  readr::read_csv(path, show_col_types = FALSE) %>%
    transmute(
      observed_taxon_unit = taxon_unit_species,
      observed_label = coalesce(binomial, valid_name, name),
      observed_genus = genus,
      observed_family = family
    ) %>%
    filter(!is.na(observed_taxon_unit)) %>%
    mutate(
      observed_cd_nom = stringr::str_match(observed_taxon_unit, "^species:cdnom:(.+)$")[, 2],
      observed_label = clean_chr(observed_label),
      observed_name_key = normalise_taxon(extract_binomial(observed_label)),
      observed_genus = clean_chr(observed_genus),
      observed_family = clean_chr(observed_family)
    ) %>%
    distinct()
}

# ---- 4. Pool construction ----------------------------------------------------
build_pool_one <- function(assemblage_id, taxref_objects, overrides) {
  scope <- filter_taxref_scope(taxref_objects$accepted, assemblage_id)
  observed <- read_observed_species(assemblage_id)

  # Limit aliases to the assemblage-specific accepted candidate pool.
  scope_alias <- taxref_objects$alias %>%
    inner_join(scope %>% select(cd_ref), by = "cd_ref")

  # First: direct CD_NOM alias mapping.
  direct <- observed %>%
    left_join(
      scope_alias %>% select(observed_cd_nom = cd_nom, cd_ref_direct = cd_ref),
      by = "observed_cd_nom"
    )

  # Second: exact accepted binomial fallback.
  name_match <- scope %>%
    transmute(
      observed_name_key = normalise_taxon(accepted_name),
      cd_ref_name = cd_ref
    )

  matches <- direct %>%
    left_join(name_match, by = "observed_name_key") %>%
    left_join(
      overrides %>%
        filter(.data$assemblage_id == !!assemblage_id) %>%
        select(observed_taxon_unit, cd_ref_override, override_note = note),
      by = "observed_taxon_unit"
    ) %>%
    mutate(
      cd_ref = coalesce(cd_ref_override, cd_ref_direct, cd_ref_name),
      match_method = case_when(
        !is.na(cd_ref_override) ~ "manual_override",
        !is.na(cd_ref_direct) ~ "cd_nom_alias",
        !is.na(cd_ref_name) ~ "accepted_binomial",
        TRUE ~ "unmatched"
      )
    ) %>%
    left_join(
      scope %>%
        select(
          cd_ref, accepted_name, genus_taxref = genus,
          family_taxref = famille
        ),
      by = "cd_ref"
    ) %>%
    mutate(
      genus = coalesce(genus_taxref, observed_genus),
      family = coalesce(family_taxref, observed_family),
      accepted_name = coalesce(accepted_name, observed_label),
      in_regional_pool = !is.na(cd_ref)
    ) %>%
    select(
      observed_taxon_unit, observed_cd_nom, observed_label,
      observed_genus, observed_family,
      cd_ref, accepted_name, genus, family,
      match_method, in_regional_pool, override_note
    )

  # Any CD_REF represented by an observed RMQS species keeps the original RMQS
  # taxon unit. TAXREF-only candidates use cdref ids.
  observed_by_ref <- matches %>%
    filter(!is.na(cd_ref)) %>%
    arrange(observed_taxon_unit) %>%
    distinct(cd_ref, .keep_all = TRUE) %>%
    select(
      cd_ref,
      observed_taxon_unit,
      observed_label
    )

  pool <- scope %>%
    left_join(observed_by_ref, by = "cd_ref") %>%
    transmute(
      taxon_unit = coalesce(observed_taxon_unit, paste0("species:cdref:", cd_ref)),
      cd_ref,
      accepted_name,
      genus,
      family = famille,
      observed_in_assemblage = !is.na(observed_taxon_unit),
      candidate_origin = if_else(observed_in_assemblage, "observed_RMQS", "TAXREF_mainland")
    ) %>%
    distinct(taxon_unit, .keep_all = TRUE) %>%
    arrange(genus, accepted_name)

  audit <- matches %>%
    left_join(
      pool %>%
        count(genus, name = "n_regional_candidates"),
      by = "genus"
    ) %>%
    mutate(
      n_regional_candidates = coalesce(n_regional_candidates, 0L),
      n_alternative_candidates = pmax(n_regional_candidates - 1L, 0L),
      has_alternative_candidate = n_alternative_candidates > 0L,
      assemblage_id = assemblage_id
    )

  list(pool = pool, audit = audit)
}

make_expert_template <- function(assemblage_id, pool, audit) {
  rules <- EXPERT_GENUS_RULES %>%
    filter(.data$assemblage_id == !!assemblage_id)

  if (nrow(rules) == 0L) {
    return(tibble(
      source_taxon_unit = character(),
      source_label = character(),
      source_genus = character(),
      target_taxon_unit = character(),
      target_label = character(),
      target_genus = character(),
      target_origin = character(),
      weight = numeric(),
      enabled = logical(),
      review_note = character(),
      comment = character()
    ))
  }

  sources <- audit %>%
    filter(!is.na(observed_taxon_unit), !is.na(genus)) %>%
    transmute(
      source_taxon_unit = observed_taxon_unit,
      source_label = observed_label,
      source_genus = genus
    ) %>%
    distinct()

  # Expert templates can span genus boundaries if an explicit rule is later
  # added. No cross-genus rule is enabled by default.
  rules %>%
    inner_join(sources, by = c("source_genus" = "source_genus")) %>%
    inner_join(
      pool %>%
        transmute(
          target_taxon_unit = taxon_unit,
          target_label = accepted_name,
          target_genus = genus,
          target_origin = candidate_origin
        ),
      by = c("target_genus" = "target_genus")
    ) %>%
    filter(source_taxon_unit != target_taxon_unit) %>%
    transmute(
      source_taxon_unit,
      source_label,
      source_genus,
      target_taxon_unit,
      target_label,
      target_genus,
      target_origin,
      weight = 1,
      enabled = FALSE,
      review_note,
      comment = NA_character_
    ) %>%
    arrange(source_genus, source_label, target_genus, target_label)
}

# ---- 5. Run ------------------------------------------------------------------
dir.create(REGIONAL_POOL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(EXPERT_MAP_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(file.path(INPUT_DIR, "assemblages"))) {
  stop("Prepared assemblages are missing. Run 00_prepare_multitaxa_inputs_2024.R first.")
}

# Write a blank override template once, if the user has not already created it.
if (!file.exists(MANUAL_OVERRIDE_FILE)) {
  readr::write_csv(
    tibble(
      assemblage_id = character(),
      observed_taxon_unit = character(),
      cd_ref = character(),
      note = character()
    ),
    MANUAL_OVERRIDE_FILE
  )
}

taxref_raw <- read_taxref_flexible(TAXREF_FILE)
taxref_objects <- make_taxref_objects(taxref_raw)
overrides <- read_manual_overrides(MANUAL_OVERRIDE_FILE)

message("Accepted mainland-France TAXREF species retained: ", nrow(taxref_objects$accepted))

summary_rows <- list()

for (assemblage_id in ASSEMBLAGE_FILTERS$assemblage_id) {
  message("Building TAXREF pool: ", assemblage_id)

  result <- build_pool_one(assemblage_id, taxref_objects, overrides)
  template <- make_expert_template(assemblage_id, result$pool, result$audit)

  readr::write_csv(
    result$pool,
    file.path(REGIONAL_POOL_DIR, paste0(assemblage_id, "__regional_pool.csv"))
  )
  readr::write_csv(
    result$audit,
    file.path(REGIONAL_POOL_DIR, paste0(assemblage_id, "__observed_taxref_match_audit.csv"))
  )
  readr::write_csv(
    template,
    file.path(EXPERT_MAP_DIR, paste0(assemblage_id, "__expert_confusions_TEMPLATE.csv"))
  )

  summary_rows[[assemblage_id]] <- tibble(
    assemblage_id = assemblage_id,
    n_regional_species = nrow(result$pool),
    n_observed_species = nrow(result$audit),
    n_observed_matched = sum(result$audit$in_regional_pool, na.rm = TRUE),
    observed_match_rate = mean(result$audit$in_regional_pool, na.rm = TRUE),
    n_observed_with_alternative_candidate = sum(result$audit$has_alternative_candidate, na.rm = TRUE),
    share_observed_with_alternative_candidate = mean(result$audit$has_alternative_candidate, na.rm = TRUE),
    median_n_alternative_candidates = median(result$audit$n_alternative_candidates, na.rm = TRUE),
    n_expert_template_rows = nrow(template)
  )
}

summary_tbl <- bind_rows(summary_rows)
readr::write_csv(summary_tbl, file.path(REGIONAL_POOL_DIR, "regional_pool_build_summary.csv"))

message("\nCompleted.")
message("Inspect before running scenarios:")
message("  - ", file.path(REGIONAL_POOL_DIR, "regional_pool_build_summary.csv"))
message("  - each *__observed_taxref_match_audit.csv")
message("  - ", MANUAL_OVERRIDE_FILE, " for unmatched observed names")
message("\nTo activate an expert scenario, copy a reviewed template to:")
message("  expert_confusion_maps/<assemblage_id>__expert_confusions.csv")
message("and set enabled = TRUE only for defensible source-target pairs.")
