# Shared helpers for RMQS taxonomic-sensitivity figures (v6).
# Source R/utils.R first; it must provide read_settings() and read_scenario_catalog().

load_figure_packages <- function() {
  pkgs <- c("dplyr", "tidyr", "readr", "stringr", "purrr", "ggplot2", "scales")
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) stop("Install required figure packages: ", paste(missing, collapse = ", "))
  invisible(lapply(pkgs, library, character.only = TRUE))
}

read_csv_if_exists <- function(path, empty) {
  if (file.exists(path)) readr::read_csv(path, show_col_types = FALSE) else empty
}

clean_label <- function(x) {
  x <- as.character(x)
  stringr::str_replace_all(x, "\\n", "\n")
}

read_editorial_overrides <- function() {
  assemblage <- read_csv_if_exists(
    "config/figure_assemblage_overrides.csv",
    tibble::tibble(
      assemblage_id = character(), display_label = character(), show_main = logical(),
      main_order = numeric(), show_supplement = logical(), supp_order = numeric(), comment = character()
    )
  ) %>%
    mutate(
      display_label = clean_label(display_label),
      show_main = as.logical(show_main),
      show_supplement = as.logical(show_supplement)
    )

  scenario <- read_csv_if_exists(
    "config/figure_scenario_overrides.csv",
    tibble::tibble(
      scenario_id = character(), short_label = character(), scenario_class = character(),
      show_main_override = logical(), main_order = numeric(), comment = character()
    )
  ) %>%
    mutate(
      short_label = clean_label(short_label),
      show_main_override = as.logical(show_main_override)
    )

  list(assemblage = assemblage, scenario = scenario)
}

read_figure_inputs <- function(settings) {
  sensitivity_dir <- file.path(settings$outputs_dir, "sensitivity")
  paths <- list(
    results = file.path(sensitivity_dir, "results_summary.csv"),
    manifest = file.path(settings$outputs_dir, "assemblage_manifest.csv"),
    resolution = file.path(settings$outputs_dir, "taxonomic_resolution_audit_all.csv"),
    stage = file.path(settings$outputs_dir, "stage_audit_all.csv")
  )
  missing <- names(paths)[!file.exists(unlist(paths))]
  if (length(missing)) {
    stop("Missing figure input(s): ", paste(missing, collapse = ", "), ". Run scripts 01–05 first.")
  }

  catalog <- read_scenario_catalog("config/scenario_catalog.csv")
  editorial <- read_editorial_overrides()

  manifest <- readr::read_csv(paths$manifest, show_col_types = FALSE) %>%
    left_join(editorial$assemblage, by = "assemblage_id") %>%
    mutate(
      display_label = coalesce(display_label, stringr::str_replace_all(assemblage_id, "__", "\n")),
      show_main_figure = coalesce(show_main, FALSE),
      show_supplement_figure = coalesce(show_supplement, TRUE),
      main_order = coalesce(main_order, 9999),
      supp_order = coalesce(supp_order, 9999)
    )

  results <- readr::read_csv(paths$results, show_col_types = FALSE) %>%
    left_join(
      catalog %>% select(scenario_id, catalog_label = label, description, catalog_show_main = show_main, show_supplement),
      by = c("scenario" = "scenario_id")
    ) %>%
    left_join(
      manifest %>% transmute(
        assemblage_id,
        display_label,
        show_main_figure,
        show_supplement_figure,
        assembly_main_order = main_order,
        assembly_supp_order = supp_order
      ),
      by = "assemblage_id"
    ) %>%
    left_join(editorial$scenario, by = c("scenario" = "scenario_id")) %>%
    mutate(
      scenario_short = coalesce(short_label, catalog_label, scenario),
      scenario_class = coalesce(scenario_class, scenario_family, "Other"),
      show_main_figure = coalesce(show_main_override, catalog_show_main, FALSE),
      scenario_order = coalesce(main_order, 9999),
      display_label = clean_label(display_label),
      scenario_short = clean_label(scenario_short)
    )

  scenario_levels <- results %>%
    distinct(scenario, scenario_short, scenario_order) %>%
    arrange(scenario_order, scenario_short) %>%
    pull(scenario_short)

  results <- results %>%
    mutate(scenario_short = factor(scenario_short, levels = unique(scenario_levels)))

  list(
    results = results,
    catalog = catalog,
    manifest = manifest,
    resolution = readr::read_csv(paths$resolution, show_col_types = FALSE),
    stage = readr::read_csv(paths$stage, show_col_types = FALSE),
    sensitivity_dir = sensitivity_dir
  )
}

main_results <- function(results) {
  results %>%
    filter(show_main_figure, scenario != baseline_scenario) %>%
    arrange(assembly_main_order, scenario_order)
}

supplement_results <- function(results) {
  results %>%
    filter(show_supplement, scenario != baseline_scenario) %>%
    arrange(assembly_supp_order, scenario_order)
}

summary_to_long <- function(data, metrics, labels, id_cols) {
  stopifnot(is.character(metrics), is.character(id_cols))

  purrr::map_dfr(metrics, function(metric_id) {
    cols <- paste0(metric_id, c("_median", "_p10", "_p90"))
    if (!all(cols %in% names(data))) return(tibble::tibble())

    label_value <- unname(labels[[metric_id]])
    if (length(label_value) != 1L || is.na(label_value)) label_value <- metric_id

    data %>%
      select(dplyr::all_of(id_cols), dplyr::all_of(cols)) %>%
      transmute(
        across(dplyr::all_of(id_cols)),
        metric = metric_id,
        metric_label = label_value,
        median = .data[[cols[[1]]]],
        p10 = .data[[cols[[2]]]],
        p90 = .data[[cols[[3]]]]
      )
  }) %>%
    filter(is.finite(median) | is.finite(p10) | is.finite(p90))
}

save_figure <- function(plot, path, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = height, units = "mm", dpi = 450)
}

clear_generated_figures <- function(dir) {
  if (!dir.exists(dir)) return(invisible(NULL))
  files <- list.files(dir, pattern = "^Fig.*\\.(png|pdf)$", full.names = TRUE)
  if (length(files)) unlink(files)
  invisible(NULL)
}

reference_value <- function(metric) {
  dplyr::case_when(
    stringr::str_detect(metric, "stability") ~ 1,
    TRUE ~ 0
  )
}

error_workflow_divider <- function(data) {
  error_orders <- data %>%
    filter(scenario_class %in% c("Identification error", "Expert-informed", "Expert component")) %>%
    distinct(scenario_order) %>%
    pull(scenario_order)

  if (!length(error_orders)) return(NA_real_)
  max(error_orders, na.rm = TRUE) + 0.5
}
