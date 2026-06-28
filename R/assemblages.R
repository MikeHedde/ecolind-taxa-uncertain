rank_to_level <- function(rank) {
  dplyr::case_when(
    rank %in% c("S", "E", "SE", "sE") ~ "species",
    rank %in% c("G", "sG") ~ "genus",
    rank %in% c("F", "sF") ~ "family",
    rank %in% c("O", "sO") ~ "order",
    TRUE ~ "other"
  )
}

read_and_standardise_raw <- function(settings) {
  if (!file.exists(settings$raw_input_file)) stop("Raw input not found: ", settings$raw_input_file)

  raw <- readr::read_delim(
    settings$raw_input_file, delim = ";",
    locale = readr::locale(encoding = "Latin1"),
    show_col_types = FALSE, name_repair = "unique"
  ) %>% janitor::clean_names()

  needed <- c("code_ech", "project", "yr", "station", "method", "phylum", "class",
              "order", "family", "cd_nom", "name", "valid_name", "rank", "tot_abund")
  missing <- setdiff(needed, names(raw))
  if (length(missing)) stop("Raw input missing columns: ", paste(missing, collapse = ", "))

  optional_stage <- c("male", "female", "juv", "indet", "remarque_individu")
  for (nm in setdiff(optional_stage, names(raw))) raw[[nm]] <- NA

  raw %>%
    mutate(
      across(c(code_ech, project, station, method, phylum, subphylum, class, order,
               family, cd_nom, name, valid_name, rank, remarque_individu), clean_chr),
      yr = suppressWarnings(as.integer(yr)),
      across(any_of(c("male", "female", "juv", "indet", "tot_abund")), parse_num),
      abundance = coalesce(tot_abund, 0),
      taxo_level = rank_to_level(rank),
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
      taxon_unit_species = if_else(taxo_level == "species" & !is.na(species_key), paste0("species:", species_key), NA_character_),
      taxon_unit_genus = if_else(taxo_level %in% c("species", "genus") & !is.na(genus), paste0("genus:", genus), NA_character_),
      taxon_unit_family = if_else(taxo_level %in% c("species", "genus", "family") & !is.na(family), paste0("family:", family), NA_character_),
      adults_definite = coalesce(male, 0) + coalesce(female, 0),
      juveniles_explicit = coalesce(juv, 0),
      n_stage_fields = rowSums(!is.na(cbind(male, female, juv, indet))),
      stage_sum = rowSums(cbind(coalesce(male, 0), coalesce(female, 0), coalesce(juv, 0), coalesce(indet, 0))),
      abundance_not_stage_resolved = if_else(n_stage_fields == 0, abundance, pmax(abundance - stage_sum, 0)),
      stage_accounting_difference = abundance - stage_sum
    ) %>%
    filter(project == settings$target_project, yr == settings$target_year)
}

discover_assemblages <- function(records, registry) {
  out <- purrr::map_dfr(seq_len(nrow(registry)), function(i) {
    rule <- registry[i, ]
    col <- rule$selector_col[[1]]
    if (!col %in% names(records)) stop("Registry selector column absent from raw data: ", col)

    x <- records %>%
      filter(.data[[col]] == rule$selector_value[[1]], !is.na(method))

    methods <- sort(unique(x$method))
    allowed <- split_methods(rule$method_include[[1]])
    if (length(allowed)) methods <- intersect(methods, allowed)
    if (!length(methods)) return(tibble())

    tibble(
      taxon_key = rule$taxon_key[[1]],
      display_label = rule$display_label[[1]],
      selector_col = col,
      selector_value = rule$selector_value[[1]],
      taxref_filter_col = rule$taxref_filter_col[[1]],
      taxref_filter_value = rule$taxref_filter_value[[1]],
      expert_workbook = rule$expert_workbook[[1]],
      preferred_role = rule$preferred_role[[1]],
      stage_scenario_candidate = rule$stage_scenario_candidate[[1]],
      method = methods,
      assemblage_id = paste(rule$taxon_key[[1]], slugify(methods), sep = "__")
    )
  })

  if (!nrow(out)) stop("No assemblages discovered. Check registry selectors and raw-data values.")
  out %>% arrange(taxon_key, method)
}

make_station_matrix <- function(comm_long, station_frame) {
  stations <- as.character(station_frame$station)
  taxa <- sort(unique(comm_long$taxon_unit))
  if (!length(taxa)) return(matrix(0, nrow = length(stations), ncol = 0, dimnames = list(stations, character())))
  comm_long %>%
    filter(!is.na(station), !is.na(taxon_unit), abundance > 0) %>%
    group_by(station, taxon_unit) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    tidyr::complete(station = stations, taxon_unit = taxa, fill = list(abundance = 0)) %>%
    tidyr::pivot_wider(names_from = taxon_unit, values_from = abundance, values_fill = 0) %>%
    arrange(match(station, stations)) %>%
    tibble::column_to_rownames("station") %>%
    as.matrix()
}

analysis_flag <- function(n_stations, abundance, n_taxa, settings) {
  dplyr::case_when(
    n_stations >= settings$min_main_stations && abundance >= settings$min_main_abundance && n_taxa >= settings$min_main_taxa ~ "main_analysis_candidate",
    n_stations >= 20 && abundance >= 20 && n_taxa >= 8 ~ "secondary_or_sensitivity_candidate",
    TRUE ~ "audit_only_or_insufficient"
  )
}

export_assemblages <- function(records, assemblages, settings) {
  out_dir <- file.path(settings$outputs_dir, "assemblages")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  sampling_frame <- records %>%
    distinct(code_ech, station, method) %>%
    filter(!is.na(station), !is.na(method))

  manifest <- vector("list", nrow(assemblages))
  audit_tax <- vector("list", nrow(assemblages))
  audit_stage <- vector("list", nrow(assemblages))

  for (i in seq_len(nrow(assemblages))) {
    def <- assemblages[i, ]
    message("Preparing: ", def$assemblage_id)

    station_frame <- sampling_frame %>%
      filter(method == def$method) %>%
      distinct(station) %>%
      arrange(station)

    x <- records %>%
      filter(method == def$method, .data[[def$selector_col]] == def$selector_value, abundance > 0)

    usable <- x %>% filter(taxo_level %in% c("species", "genus", "family"), !is.na(taxon_unit_rtu))

    tax_audit <- x %>%
      mutate(resolution = case_when(
        taxo_level == "species" ~ "species",
        taxo_level == "genus" ~ "genus",
        taxo_level == "family" ~ "family",
        taxo_level == "order" ~ "order",
        TRUE ~ "other_or_unusable"
      )) %>%
      group_by(resolution) %>%
      summarise(
        n_records = n(), total_abundance = sum(abundance),
        n_reported_names = n_distinct(coalesce(valid_name, name), na.rm = TRUE),
        n_stations_positive = n_distinct(station),
        .groups = "drop"
      ) %>%
      mutate(
        assemblage_id = def$assemblage_id,
        taxon_key = def$taxon_key,
        group_label = def$display_label,
        method = def$method,
        abundance_share = total_abundance / sum(total_abundance),
        record_share = n_records / sum(n_records)
      )

    comm_rtu <- usable %>%
      transmute(station, taxon_unit = taxon_unit_rtu, abundance) %>%
      group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop")
    comm_species <- usable %>%
      filter(taxo_level == "species", !is.na(taxon_unit_species)) %>%
      transmute(station, taxon_unit = taxon_unit_species, abundance) %>%
      group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop")
    comm_genus <- usable %>%
      filter(!is.na(taxon_unit_genus)) %>%
      transmute(station, taxon_unit = taxon_unit_genus, abundance) %>%
      group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop")
    comm_family <- usable %>%
      filter(!is.na(taxon_unit_family)) %>%
      transmute(station, taxon_unit = taxon_unit_family, abundance) %>%
      group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop")

    lookup <- usable %>%
      transmute(
        taxon_unit_rtu, taxon_unit_species, taxon_unit_genus, taxon_unit_family,
        taxo_level, species_key, binomial, genus, family, order, class, valid_name, name,
        observed_cd_nom = str_match(taxon_unit_species, "^species:cdnom:(.+)$")[, 2]
      ) %>% distinct()

    stage_audit <- x %>%
      summarise(
        assemblage_id = def$assemblage_id,
        taxon_key = def$taxon_key,
        group_label = def$display_label,
        method = def$method,
        total_abundance = sum(abundance),
        male = sum(coalesce(male, 0)),
        female = sum(coalesce(female, 0)),
        juvenile_explicit = sum(coalesce(juv, 0)),
        indeterminate_stage_or_sex = sum(coalesce(indet, 0)),
        juvenile_share = juvenile_explicit / total_abundance,
        stage_scenario_candidate = def$stage_scenario_candidate
      )

    adult <- if (isTRUE(def$stage_scenario_candidate)) {
      usable %>%
        mutate(adult_abundance = adults_definite) %>%
        filter(adult_abundance > 0) %>%
        transmute(station, taxon_unit = taxon_unit_rtu, abundance = adult_abundance) %>%
        group_by(station, taxon_unit) %>% summarise(abundance = sum(abundance), .groups = "drop")
    } else tibble(station = character(), taxon_unit = character(), abundance = numeric())

    station_audit <- station_frame %>%
      left_join(
        comm_rtu %>% group_by(station) %>% summarise(total_abundance = sum(abundance), rtu_richness = n_distinct(taxon_unit), .groups = "drop"),
        by = "station"
      ) %>%
      mutate(
        assemblage_id = def$assemblage_id,
        total_abundance = coalesce(total_abundance, 0),
        rtu_richness = coalesce(rtu_richness, 0L),
        detected = total_abundance > 0
      )

    summary <- tibble(
      assemblage_id = def$assemblage_id,
      taxon_key = def$taxon_key,
      group_label = def$display_label,
      method = def$method,
      role = def$preferred_role,
      selector = paste0(def$selector_col, " = ", def$selector_value),
      taxref_filter_col = def$taxref_filter_col,
      taxref_filter_value = def$taxref_filter_value,
      expert_workbook = def$expert_workbook,
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
      family_coarsening_applicable = n_distinct(usable$family, na.rm = TRUE) >= 2,
      analysis_flag = analysis_flag(nrow(station_frame), sum(comm_rtu$abundance), n_distinct(comm_rtu$taxon_unit), settings)
    )

    prefix <- file.path(out_dir, def$assemblage_id)
    readr::write_csv(station_frame %>% mutate(assemblage_id = def$assemblage_id), paste0(prefix, "__station_frame.csv"))
    readr::write_csv(comm_rtu, paste0(prefix, "__rtu_long.csv"))
    readr::write_csv(comm_species, paste0(prefix, "__species_long.csv"))
    readr::write_csv(comm_genus, paste0(prefix, "__genus_long.csv"))
    readr::write_csv(comm_family, paste0(prefix, "__family_long.csv"))
    readr::write_csv(lookup, paste0(prefix, "__taxon_lookup.csv"))
    readr::write_csv(tax_audit, paste0(prefix, "__taxonomic_resolution_audit.csv"))
    readr::write_csv(stage_audit, paste0(prefix, "__stage_audit.csv"))
    readr::write_csv(station_audit, paste0(prefix, "__station_audit.csv"))
    if (nrow(adult)) readr::write_csv(adult, paste0(prefix, "__adult_only_rtu_long.csv"))
    saveRDS(make_station_matrix(comm_rtu, station_frame), paste0(prefix, "__rtu_matrix.rds"))

    manifest[[i]] <- summary
    audit_tax[[i]] <- tax_audit
    audit_stage[[i]] <- stage_audit
  }

  manifest <- bind_rows(manifest) %>% arrange(role, taxon_key, method)
  readr::write_csv(assemblages, file.path(settings$outputs_dir, "assemblage_definitions.csv"))
  readr::write_csv(manifest, file.path(settings$outputs_dir, "assemblage_manifest.csv"))
  readr::write_csv(bind_rows(audit_tax), file.path(settings$outputs_dir, "taxonomic_resolution_audit_all.csv"))
  readr::write_csv(bind_rows(audit_stage), file.path(settings$outputs_dir, "stage_audit_all.csv"))
  invisible(manifest)
}
