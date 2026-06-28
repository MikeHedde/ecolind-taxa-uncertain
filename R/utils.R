required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble",
  "janitor", "vegan", "ggplot2", "sf", "gdm", "readxl"
)

load_pipeline_packages <- function(require_gdm = TRUE) {
  pkgs <- required_packages
  if (!require_gdm) pkgs <- setdiff(pkgs, c("sf", "gdm"))
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    stop("Install required packages: ", paste(missing, collapse = ", "))
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

clean_chr <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

parse_bool <- function(x, default = FALSE) {
  x <- tolower(clean_chr(x))
  out <- x %in% c("true", "t", "1", "yes", "y")
  out[is.na(x)] <- default
  out
}

parse_num <- function(x) suppressWarnings(readr::parse_number(as.character(x), na = c("", "NA", "NaN")))

slugify <- function(x) {
  x <- stringr::str_to_lower(clean_chr(x))
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- stringr::str_replace_all(x, "[^a-z0-9]+", "_")
  stringr::str_replace_all(x, "^_+|_+$", "")
}

extract_binomial <- function(x) {
  clean_chr(stringr::str_extract(clean_chr(x), "^[A-Z][A-Za-zÀ-ÖØ-öø-ÿ-]+\\s+[a-z][A-Za-zÀ-ÖØ-öø-ÿ-]+"))
}

extract_genus <- function(x) {
  clean_chr(stringr::str_extract(clean_chr(x), "^[A-Z][A-Za-zÀ-ÖØ-öø-ÿ-]+"))
}

normalise_taxon <- function(x) {
  x <- clean_chr(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  stringr::str_to_lower(x)
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) x[[1]] else NA
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || dplyr::n_distinct(x[keep]) < 2L || dplyr::n_distinct(y[keep]) < 2L) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

source_pipeline_files <- function(files = c(
  "R/utils.R", "R/assemblages.R", "R/taxref.R", "R/expert_rules.R",
  "R/metrics_gdm.R", "R/scenarios.R"
)) {
  for (path in files) {
    if (!file.exists(path)) stop("Missing pipeline file: ", path)
    source(path)
  }
}

read_settings <- function(path = "config/analysis_settings.R") {
  e <- new.env(parent = baseenv())
  sys.source(path, envir = e)
  
  if (!exists("settings", envir = e, inherits = FALSE)) {
    stop("Config file must define `settings`.")
  }
  
  e$settings
}

read_registry <- function(path = "config/taxon_registry.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  needed <- c(
    "taxon_key", "display_label", "selector_col", "selector_value",
    "taxref_filter_col", "taxref_filter_value", "expert_workbook",
    "include", "preferred_role", "stage_scenario_candidate"
  )
  missing <- setdiff(needed, names(x))
  if (length(missing)) stop("taxon_registry.csv missing: ", paste(missing, collapse = ", "))
  if (!"method_include" %in% names(x)) x$method_include <- NA_character_
  if (!"notes" %in% names(x)) x$notes <- NA_character_
  x %>%
    mutate(
      across(c(taxon_key, display_label, selector_col, selector_value, taxref_filter_col,
               taxref_filter_value, expert_workbook, preferred_role, method_include, notes), clean_chr),
      include = parse_bool(include),
      stage_scenario_candidate = parse_bool(stage_scenario_candidate)
    ) %>%
    filter(include) %>%
    distinct(taxon_key, .keep_all = TRUE)
}

read_scenario_catalog <- function(path = "config/scenario_catalog.csv") {
  x <- readr::read_csv(path, show_col_types = FALSE) %>% janitor::clean_names()
  needed <- c("scenario_id", "handler", "scenario_family", "baseline_id", "run", "show_main", "show_supplement", "label")
  missing <- setdiff(needed, names(x))
  if (length(missing)) stop("scenario_catalog.csv missing: ", paste(missing, collapse = ", "))
  if (!"error_rate" %in% names(x)) x$error_rate <- NA_real_
  if (!"description" %in% names(x)) x$description <- NA_character_
  x %>%
    mutate(
      across(c(scenario_id, handler, scenario_family, baseline_id, label, description), clean_chr),
      error_rate = suppressWarnings(as.numeric(error_rate)),
      run = parse_bool(run),
      show_main = parse_bool(show_main),
      show_supplement = parse_bool(show_supplement)
    )
}

split_methods <- function(x) {
  x <- clean_chr(x)
  if (is.na(x)) return(character())
  stringr::str_split(x, "\\|")[[1]] %>% clean_chr() %>% stats::na.omit()
}
