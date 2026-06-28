empty_expert_rules <- function() {
  tibble(
    rule_id = character(), rule_type = character(),
    source_cd_ref = character(), source_unit = character(),
    target_cd_ref = character(), target_unit = character(),
    weight = numeric(), enabled = logical(), degree_confusion = numeric(), comment = character()
  )
}

read_expert_rules_for_taxon <- function(taxon_row, settings) {
  file_name <- taxon_row$expert_workbook[[1]]
  if (is.na(file_name) || !nzchar(file_name)) return(empty_expert_rules())
  path <- file.path(settings$expert_input_dir, file_name)
  if (!file.exists(path)) return(empty_expert_rules())

  ext <- tolower(tools::file_ext(path))
  x <- if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path, sheet = "rules") %>% janitor::clean_names()
  } else {
    readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  }

  required <- c("rule_id", "rule_type", "source_cd_ref", "source_unit",
                "target_cd_ref", "target_unit", "weight", "enabled")
  for (nm in setdiff(required, names(x))) x[[nm]] <- NA
  if (!"degree_confusion" %in% names(x)) x$degree_confusion <- NA_real_
  if (!"comment" %in% names(x)) x$comment <- NA_character_

  x %>%
    transmute(
      rule_id = clean_chr(rule_id),
      rule_type = clean_chr(rule_type),
      source_cd_ref = clean_chr(source_cd_ref),
      source_unit = clean_chr(source_unit),
      target_cd_ref = clean_chr(target_cd_ref),
      target_unit = clean_chr(target_unit),
      weight = coalesce(suppressWarnings(as.numeric(weight)), 1),
      enabled = parse_bool(enabled),
      degree_confusion = suppressWarnings(as.numeric(degree_confusion)),
      comment = clean_chr(comment)
    ) %>%
    filter(enabled, !is.na(rule_type)) %>%
    distinct()
}

read_assemblage_taxref_audit <- function(assemblage_id_value, settings) {
  p <- file.path(settings$regional_pool_dir, paste0(assemblage_id_value, "__observed_taxref_match_audit.csv"))
  if (!file.exists(p)) return(tibble(observed_taxon_unit = character(), cd_ref = character()))
  readr::read_csv(p, show_col_types = FALSE) %>%
    transmute(observed_taxon_unit, cd_ref) %>%
    mutate(across(everything(), as.character)) %>%
    filter(!is.na(observed_taxon_unit), !is.na(cd_ref)) %>%
    distinct(observed_taxon_unit, .keep_all = TRUE)
}

read_assemblage_pool <- function(assemblage_id_value, settings) {
  p <- file.path(settings$regional_pool_dir, paste0(assemblage_id_value, "__regional_pool.csv"))
  if (!file.exists(p)) return(tibble(taxon_unit = character(), cd_ref = character(), genus = character()))
  readr::read_csv(p, show_col_types = FALSE) %>%
    transmute(taxon_unit, cd_ref, genus, family, accepted_name) %>%
    mutate(across(c(taxon_unit, cd_ref, genus, family, accepted_name), as.character)) %>%
    distinct()
}

resolve_species_unit <- function(cd_ref_value, observed_audit, pool) {
  cd_ref_value <- clean_chr(cd_ref_value)
  if (is.na(cd_ref_value)) return(NA_character_)

  observed <- observed_audit %>%
    filter(.data$cd_ref == cd_ref_value) %>%
    pull(observed_taxon_unit) %>% unique()

  if (length(observed)) return(observed[[1]])

  candidates <- pool %>%
    filter(.data$cd_ref == cd_ref_value) %>%
    pull(taxon_unit) %>% unique()

  if (length(candidates)) return(candidates[[1]])
  NA_character_
}

project_expert_rules_one <- function(assemblage_row, taxon_row, settings) {
  assemblage_id_value <- assemblage_row$assemblage_id[[1]]
  taxon_key_value <- assemblage_row$taxon_key[[1]]
  rules <- read_expert_rules_for_taxon(taxon_row, settings)

  empty <- list(
    species_map = tibble(source_taxon_unit = character(), target_taxon_unit = character(), weight = numeric(), comment = character()),
    reporting_rules = tibble(source_taxon_unit = character(), target_taxon_unit = character(), comment = character()),
    rtu_map = tibble(source_taxon_unit = character(), target_taxon_unit = character(), weight = numeric(), comment = character()),
    audit = tibble(assemblage_id = character(), taxon_key = character(), rule_id = character(), rule_type = character(), status = character(), detail = character())
  )
  if (!nrow(rules)) return(empty)

  lookup_path <- file.path(settings$outputs_dir, "assemblages", paste0(assemblage_id_value, "__taxon_lookup.csv"))
  lookup <- readr::read_csv(lookup_path, show_col_types = FALSE)
  observed_audit <- read_assemblage_taxref_audit(assemblage_id_value, settings)
  pool <- read_assemblage_pool(assemblage_id_value, settings)

  rtu_units <- unique(lookup$taxon_unit_rtu)
  species_units <- unique(lookup$taxon_unit_species)
  rtu_units <- rtu_units[!is.na(rtu_units)]
  species_units <- species_units[!is.na(species_units)]

  resolved <- rules %>%
    rowwise() %>%
    mutate(
      source_resolved = case_when(
        rule_type %in% c("species_confusion", "reporting_rule") && !is.na(source_cd_ref) ~ resolve_species_unit(source_cd_ref, observed_audit, pool),
        !is.na(source_unit) ~ source_unit,
        TRUE ~ NA_character_
      ),
      target_resolved = case_when(
        rule_type == "species_confusion" && !is.na(target_cd_ref) ~ resolve_species_unit(target_cd_ref, observed_audit, pool),
        !is.na(target_unit) ~ target_unit,
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    mutate(
      source_present = case_when(
        rule_type == "species_confusion" ~ source_resolved %in% species_units,
        TRUE ~ source_resolved %in% rtu_units
      ),
      target_available = !is.na(target_resolved),
      status = case_when(
        !source_present ~ "not_applicable_source_absent",
        !target_available ~ "unresolved_target",
        source_resolved == target_resolved ~ "self_transition",
        TRUE ~ "projected"
      )
    )

  audit <- resolved %>%
    mutate(
      assemblage_id = assemblage_id_value,
      taxon_key = taxon_key_value
    ) %>%
    transmute(
      assemblage_id, taxon_key, rule_id, rule_type, source_cd_ref, source_unit,
      target_cd_ref, target_unit, source_resolved, target_resolved,
      status,
      detail = comment
    )

  projected <- resolved %>% filter(status == "projected")

  list(
    species_map = projected %>%
      filter(rule_type == "species_confusion") %>%
      transmute(source_taxon_unit = source_resolved, target_taxon_unit = target_resolved,
                weight, comment) %>% distinct(),
    reporting_rules = projected %>%
      filter(rule_type == "reporting_rule") %>%
      transmute(source_taxon_unit = source_resolved, target_taxon_unit = target_resolved,
                comment) %>% distinct(),
    rtu_map = projected %>%
      filter(rule_type == "rtu_confusion") %>%
      transmute(source_taxon_unit = source_resolved, target_taxon_unit = target_resolved,
                weight, comment) %>% distinct(),
    audit = audit
  )
}
