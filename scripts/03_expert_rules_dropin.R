# R/expert_rules.R — drop-in replacement
# Robust handling of manually curated expert rules.
#
# Fixes:
#   - reads both comma- and semicolon-delimited CSV files;
#   - accepts Latin1/UTF-8 encoded expert files;
#   - resolves source/target IDs against both observed_cd_nom and cd_ref;
#   - normalises RTU targets such as "Oritoniscus" into "genus:Oritoniscus";
#   - optionally filters a master expert file by taxon_key if such a column exists.

empty_expert_rules <- function() {
  tibble::tibble(
    rule_id = character(), rule_type = character(),
    source_cd_ref = character(), source_unit = character(),
    target_cd_ref = character(), target_unit = character(),
    weight = numeric(), enabled = logical(), degree_confusion = numeric(),
    comment = character(), taxon_key = character()
  )
}

detect_delim <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (length(first) == 0L) return(",")
  if (stringr::str_count(first, ";") > stringr::str_count(first, ",")) ";" else ","
}

read_expert_csv_flexible <- function(path) {
  delim <- detect_delim(path)

  # Try UTF-8 first, then Latin1. readr can read both delimiters.
  out <- try(
    readr::read_delim(
      path,
      delim = delim,
      locale = readr::locale(encoding = "UTF-8"),
      show_col_types = FALSE,
      name_repair = "unique"
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    out <- readr::read_delim(
      path,
      delim = delim,
      locale = readr::locale(encoding = "Latin1"),
      show_col_types = FALSE,
      name_repair = "unique"
    )
  }

  janitor::clean_names(out)
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
    read_expert_csv_flexible(path)
  }

  required <- c(
    "rule_id", "rule_type", "source_cd_ref", "source_unit",
    "target_cd_ref", "target_unit", "weight", "enabled"
  )
  for (nm in setdiff(required, names(x))) x[[nm]] <- NA
  if (!"degree_confusion" %in% names(x)) x$degree_confusion <- NA_real_
  if (!"comment" %in% names(x)) x$comment <- NA_character_
  if (!"taxon_key" %in% names(x)) x$taxon_key <- NA_character_

  out <- x %>%
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
      comment = clean_chr(comment),
      taxon_key = clean_chr(taxon_key)
    ) %>%
    filter(enabled, !is.na(rule_type)) %>%
    mutate(
      # Backward-compatible aliases.
      rule_type = dplyr::recode(
        rule_type,
        genus_confusion = "rtu_confusion",
        genus_to_genus = "rtu_confusion",
        species_to_species = "species_confusion",
        report_to_genus = "reporting_rule",
        report_to_family = "reporting_rule",
        .default = rule_type
      )
    ) %>%
    distinct()

  # If a master expert table includes an explicit taxon_key column, use it.
  # If it does not, all rules are read and later filtered by source presence.
  this_taxon <- taxon_row$taxon_key[[1]]
  if (any(!is.na(out$taxon_key))) {
    out <- out %>% filter(is.na(taxon_key) | taxon_key == this_taxon)
  }

  out
}

read_assemblage_taxref_audit <- function(assemblage_id_value, settings) {
  p <- file.path(
    settings$regional_pool_dir,
    paste0(assemblage_id_value, "__observed_taxref_match_audit.csv")
  )

  if (!file.exists(p)) {
    return(tibble::tibble(
      observed_taxon_unit = character(),
      observed_cd_nom = character(),
      cd_ref = character()
    ))
  }

  x <- readr::read_csv(p, show_col_types = FALSE) %>% janitor::clean_names()
  if (!"observed_cd_nom" %in% names(x)) x$observed_cd_nom <- NA_character_
  if (!"cd_ref" %in% names(x)) x$cd_ref <- NA_character_

  x %>%
    transmute(
      observed_taxon_unit = as.character(observed_taxon_unit),
      observed_cd_nom = clean_chr(observed_cd_nom),
      cd_ref = clean_chr(cd_ref)
    ) %>%
    filter(!is.na(observed_taxon_unit)) %>%
    distinct()
}

read_assemblage_pool <- function(assemblage_id_value, settings) {
  p <- file.path(
    settings$regional_pool_dir,
    paste0(assemblage_id_value, "__regional_pool.csv")
  )

  if (!file.exists(p)) {
    return(tibble::tibble(
      taxon_unit = character(),
      cd_ref = character(),
      genus = character(),
      family = character(),
      accepted_name = character()
    ))
  }

  readr::read_csv(p, show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    transmute(
      taxon_unit = as.character(taxon_unit),
      cd_ref = clean_chr(cd_ref),
      genus = clean_chr(genus),
      family = clean_chr(family),
      accepted_name = clean_chr(accepted_name)
    ) %>%
    distinct()
}

resolve_species_unit <- function(id_value, observed_audit, pool) {
  # Despite the historical column name source_cd_ref/target_cd_ref, manually
  # curated files often contain CD_NOM values. We therefore resolve against:
  #   1. observed CD_NOM;
  #   2. accepted CD_REF in the observed audit;
  #   3. CD_REF in the broader regional pool.
  id_value <- clean_chr(id_value)
  if (is.na(id_value)) return(NA_character_)

  observed_by_cd_nom <- observed_audit %>%
    filter(.data$observed_cd_nom == id_value) %>%
    pull(observed_taxon_unit) %>%
    unique()

  if (length(observed_by_cd_nom)) return(observed_by_cd_nom[[1]])

  observed_by_cd_ref <- observed_audit %>%
    filter(.data$cd_ref == id_value) %>%
    pull(observed_taxon_unit) %>%
    unique()

  if (length(observed_by_cd_ref)) return(observed_by_cd_ref[[1]])

  candidates <- pool %>%
    filter(.data$cd_ref == id_value) %>%
    pull(taxon_unit) %>%
    unique()

  if (length(candidates)) return(candidates[[1]])

  NA_character_
}

normalise_expert_unit <- function(x) {
  x <- clean_chr(x)
  if (is.na(x)) return(NA_character_)

  # Already in internal format.
  if (stringr::str_detect(x, "^(species|genus|family|order):")) return(x)

  # A bare numeric ID is not an RTU unit.
  if (stringr::str_detect(x, "^\\d+$")) return(NA_character_)

  # Remove simple authorship/details when a full name was entered, keeping the
  # first word as genus unless it looks like a family.
  candidate <- stringr::str_extract(x, "^[A-Z][A-Za-zÀ-ÖØ-öø-ÿ-]+")
  if (is.na(candidate)) return(NA_character_)

  if (stringr::str_detect(candidate, "idae$")) {
    paste0("family:", candidate)
  } else {
    paste0("genus:", candidate)
  }
}

project_expert_rules_one <- function(assemblage_row, taxon_row, settings) {
  assemblage_id_value <- assemblage_row$assemblage_id[[1]]
  taxon_key_value <- assemblage_row$taxon_key[[1]]
  rules <- read_expert_rules_for_taxon(taxon_row, settings)

  empty <- list(
    species_map = tibble::tibble(
      source_taxon_unit = character(),
      target_taxon_unit = character(),
      weight = numeric(),
      comment = character()
    ),
    reporting_rules = tibble::tibble(
      source_taxon_unit = character(),
      target_taxon_unit = character(),
      comment = character()
    ),
    rtu_map = tibble::tibble(
      source_taxon_unit = character(),
      target_taxon_unit = character(),
      weight = numeric(),
      comment = character()
    ),
    audit = tibble::tibble(
      assemblage_id = character(),
      taxon_key = character(),
      rule_id = character(),
      rule_type = character(),
      status = character(),
      detail = character()
    )
  )

  if (!nrow(rules)) return(empty)

  lookup_path <- file.path(
    settings$outputs_dir,
    "assemblages",
    paste0(assemblage_id_value, "__taxon_lookup.csv")
  )
  lookup <- readr::read_csv(lookup_path, show_col_types = FALSE) %>% janitor::clean_names()

  observed_audit <- read_assemblage_taxref_audit(assemblage_id_value, settings)
  pool <- read_assemblage_pool(assemblage_id_value, settings)

  rtu_units <- unique(as.character(lookup$taxon_unit_rtu))
  species_units <- unique(as.character(lookup$taxon_unit_species))
  rtu_units <- rtu_units[!is.na(rtu_units)]
  species_units <- species_units[!is.na(species_units)]

  resolved <- rules %>%
    rowwise() %>%
    mutate(
      source_resolved = case_when(
        rule_type %in% c("species_confusion", "reporting_rule") &&
          !is.na(source_cd_ref) ~ resolve_species_unit(source_cd_ref, observed_audit, pool),
        !is.na(source_unit) ~ normalise_expert_unit(source_unit),
        TRUE ~ NA_character_
      ),
      target_resolved = case_when(
        rule_type == "species_confusion" &&
          !is.na(target_cd_ref) ~ resolve_species_unit(target_cd_ref, observed_audit, pool),
        !is.na(target_unit) ~ normalise_expert_unit(target_unit),
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
        is.na(source_resolved) ~ "unresolved_source",
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
      assemblage_id,
      taxon_key,
      rule_id,
      rule_type,
      source_cd_ref,
      source_unit,
      target_cd_ref,
      target_unit,
      source_resolved,
      target_resolved,
      status,
      detail = comment
    )

  projected <- resolved %>% filter(status == "projected")

  list(
    species_map = projected %>%
      filter(rule_type == "species_confusion") %>%
      transmute(
        source_taxon_unit = source_resolved,
        target_taxon_unit = target_resolved,
        weight,
        comment
      ) %>%
      distinct(),
    reporting_rules = projected %>%
      filter(rule_type == "reporting_rule") %>%
      transmute(
        source_taxon_unit = source_resolved,
        target_taxon_unit = target_resolved,
        comment
      ) %>%
      distinct(),
    rtu_map = projected %>%
      filter(rule_type == "rtu_confusion") %>%
      transmute(
        source_taxon_unit = source_resolved,
        target_taxon_unit = target_resolved,
        weight,
        comment
      ) %>%
      distinct(),
    audit = audit
  )
}
