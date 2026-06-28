read_taxref_flexible <- function(path) {
  if (!file.exists(path)) stop("TAXREF file not found: ", path)
  first <- readLines(path, n = 1L, warn = FALSE, encoding = "UTF-8")
  delim <- if (stringr::str_count(first, ";") >= stringr::str_count(first, "\t")) ";" else "\t"
  readr::read_delim(path, delim = delim, locale = readr::locale(encoding = "Latin1"),
                    show_col_types = FALSE, name_repair = "unique") %>% janitor::clean_names()
}

first_taxref_col <- function(x, candidates, required = TRUE) {
  hit <- intersect(candidates, names(x))
  if (!length(hit)) {
    if (required) stop("Missing TAXREF column; tried: ", paste(candidates, collapse = ", "))
    return(rep(NA_character_, nrow(x)))
  }
  x[[hit[[1]]]]
}

make_taxref_objects <- function(raw, settings) {
  x <- tibble(
    cd_nom = clean_chr(first_taxref_col(raw, c("cd_nom"))),
    cd_ref = clean_chr(first_taxref_col(raw, c("cd_ref"))),
    rang = toupper(clean_chr(first_taxref_col(raw, c("rang", "rank")))),
    lb_nom = clean_chr(first_taxref_col(raw, c("lb_nom", "nom_complet"), FALSE)),
    nom_valide = clean_chr(first_taxref_col(raw, c("nom_valide"), FALSE)),
    classe = clean_chr(first_taxref_col(raw, c("classe", "class"), FALSE)),
    ordre = clean_chr(first_taxref_col(raw, c("ordre", "order"), FALSE)),
    famille = clean_chr(first_taxref_col(raw, c("famille", "family"), FALSE)),
    fr = toupper(clean_chr(first_taxref_col(raw, c("fr"), FALSE)))
  ) %>%
    mutate(
      cd_ref = coalesce(cd_ref, cd_nom),
      accepted_name = coalesce(extract_binomial(nom_valide), extract_binomial(lb_nom)),
      genus = extract_genus(accepted_name),
      is_accepted = cd_nom == cd_ref
    )

  accepted <- x %>%
    filter(
      rang %in% settings$taxref_species_ranks,
      fr %in% settings$taxref_fr_status_keep,
      !is.na(cd_ref), !is.na(accepted_name), !is.na(genus)
    ) %>%
    arrange(desc(is_accepted)) %>%
    distinct(cd_ref, .keep_all = TRUE) %>%
    select(cd_ref, accepted_name, genus, famille, classe, ordre, fr)

  aliases <- x %>%
    filter(!is.na(cd_nom), !is.na(cd_ref)) %>%
    inner_join(accepted %>% select(cd_ref), by = "cd_ref") %>%
    transmute(cd_nom, cd_ref, alias_name = coalesce(extract_binomial(lb_nom), extract_binomial(nom_valide))) %>%
    distinct()

  list(accepted = accepted, aliases = aliases)
}

read_observed_species <- function(assemblage_id, settings) {
  p <- file.path(settings$outputs_dir, "assemblages", paste0(assemblage_id, "__taxon_lookup.csv"))
  if (!file.exists(p)) stop("Missing lookup: ", p)
  readr::read_csv(p, show_col_types = FALSE) %>%
    transmute(
      observed_taxon_unit = as.character(taxon_unit_species),
      observed_cd_nom = clean_chr(observed_cd_nom),
      observed_label = coalesce(binomial, valid_name, name),
      observed_genus = genus,
      observed_family = family
    ) %>%
    filter(!is.na(observed_taxon_unit)) %>%
    distinct()
}

filter_taxref_scope <- function(accepted, filter_col, filter_value) {
  if (!filter_col %in% names(accepted)) stop("Unknown TAXREF filter column: ", filter_col)
  accepted %>% filter(normalise_taxon(.data[[filter_col]]) == normalise_taxon(filter_value))
}

build_regional_pool_one <- function(assemblage_row, taxref, settings) {
  assemblage_id <- assemblage_row$assemblage_id[[1]]
  observed <- read_observed_species(assemblage_id, settings)
  scope <- filter_taxref_scope(
    taxref$accepted,
    assemblage_row$taxref_filter_col[[1]],
    assemblage_row$taxref_filter_value[[1]]
  )

  observed_match <- observed %>%
    left_join(taxref$aliases, by = c("observed_cd_nom" = "cd_nom")) %>%
    mutate(name_key = normalise_taxon(extract_binomial(observed_label))) %>%
    left_join(
      scope %>% transmute(name_key = normalise_taxon(accepted_name), cd_ref_name = cd_ref,
                          genus_name = genus, family_name = famille, accepted_name_name = accepted_name),
      by = "name_key"
    ) %>%
    mutate(
      cd_ref = coalesce(cd_ref, cd_ref_name),
      genus = coalesce(genus_name, observed_genus),
      family = coalesce(family_name, observed_family),
      accepted_name = coalesce(accepted_name_name, observed_label),
      match_method = case_when(!is.na(cd_ref) ~ "CD_NOM_or_name", TRUE ~ "unmatched")
    ) %>%
    select(observed_taxon_unit, observed_cd_nom, observed_label, observed_genus,
           observed_family, cd_ref, genus, family, accepted_name, match_method)

  observed_by_ref <- observed_match %>%
    filter(!is.na(cd_ref)) %>%
    group_by(cd_ref) %>%
    summarise(observed_taxon_unit = first(observed_taxon_unit), .groups = "drop")

  pool <- scope %>%
    left_join(observed_by_ref, by = "cd_ref") %>%
    transmute(
      taxon_unit = coalesce(observed_taxon_unit, paste0("species:cdref:", cd_ref)),
      cd_ref, genus, family = famille, accepted_name,
      observed_in_assemblage = !is.na(observed_taxon_unit),
      candidate_origin = if_else(observed_in_assemblage, "observed_RMQS", "TAXREF_mainland")
    ) %>%
    distinct(taxon_unit, .keep_all = TRUE) %>%
    arrange(genus, accepted_name)

  audit <- observed_match %>%
    left_join(pool %>% count(genus, name = "n_regional_candidates"), by = "genus") %>%
    mutate(
      assemblage_id = assemblage_id,
      in_regional_pool = !is.na(cd_ref),
      n_alternative_candidates = pmax(coalesce(n_regional_candidates, 0L) - 1L, 0L)
    )

  list(pool = pool, audit = audit)
}
