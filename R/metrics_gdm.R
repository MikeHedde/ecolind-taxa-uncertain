make_matrix <- function(comm_long, station_frame) {
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

alpha_metrics <- function(mat) {
  if (!nrow(mat) || !ncol(mat)) return(tibble(station = rownames(mat), total_abundance = 0, q0 = 0, q1 = 0, q2 = 0))
  total <- rowSums(mat)
  q0 <- rowSums(mat > 0)
  q1 <- exp(vegan::diversity(mat, "shannon"))
  q2 <- vegan::diversity(mat, "invsimpson")
  q1[!is.finite(q1)] <- 0
  q2[!is.finite(q2)] <- 0
  tibble(station = rownames(mat), total_abundance = total, q0 = q0, q1 = q1, q2 = q2)
}

community_gamma <- function(mat) if (!ncol(mat)) 0L else sum(colSums(mat) > 0)

prepare_beta_matrix <- function(mat, settings) {
  if (!nrow(mat) || !ncol(mat)) return(NULL)
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (nrow(mat) < settings$min_beta_sites || ncol(mat) < 2L) return(NULL)
  mat
}

compute_distance <- function(mat, metric, settings) {
  mat <- prepare_beta_matrix(mat, settings)
  if (is.null(mat)) return(NULL)
  if (metric == "bray") return(vegan::vegdist(mat, method = "bray"))
  if (metric == "sorensen") return(vegan::vegdist(mat > 0, method = "bray", binary = TRUE))
  vegan::vegdist(mat > 0, method = "jaccard", binary = TRUE)
}

procrustes_r2 <- function(d1, d2, settings) {
  if (is.null(d1) || is.null(d2) || length(d1) != length(d2)) return(NA_real_)
  points <- function(d) {
    x <- try(stats::cmdscale(d, k = 2, eig = TRUE, add = TRUE), silent = TRUE)
    if (inherits(x, "try-error")) return(NULL)
    x <- if (is.list(x) && !is.null(x$points)) x$points else x
    x <- as.matrix(x)
    if (nrow(x) < settings$min_beta_sites || ncol(x) < 2L || any(!is.finite(x[, 1:2]))) return(NULL)
    x[, 1:2, drop = FALSE]
  }
  p1 <- points(d1); p2 <- points(d2)
  if (is.null(p1) || is.null(p2)) return(NA_real_)
  common <- intersect(rownames(p1), rownames(p2))
  if (length(common) < settings$min_beta_sites) return(NA_real_)
  fit <- try(vegan::procrustes(p1[common,,drop=FALSE], p2[common,,drop=FALSE], symmetric = TRUE), silent = TRUE)
  if (inherits(fit, "try-error") || !is.finite(fit$ss)) return(NA_real_)
  max(0, 1 - fit$ss)
}

compare_beta <- function(base, scenario, settings) {
  sites <- intersect(rownames(base)[rowSums(base) > 0], rownames(scenario)[rowSums(scenario) > 0])
  if (length(sites) < settings$min_beta_sites) return(tibble(
    n_common_sites = length(sites), bray_stability = NA_real_, sorensen_stability = NA_real_,
    jaccard_stability = NA_real_, bray_mean_change_pct = NA_real_,
    sorensen_mean_change_pct = NA_real_, jaccard_mean_change_pct = NA_real_,
    ordination_procrustes_r2 = NA_real_
  ))
  b <- base[sites,,drop=FALSE]; s <- scenario[sites,,drop=FALSE]
  db <- compute_distance(b, "bray", settings); ds <- compute_distance(s, "bray", settings)
  sb <- compute_distance(b, "sorensen", settings); ss <- compute_distance(s, "sorensen", settings)
  jb <- compute_distance(b, "jaccard", settings); js <- compute_distance(s, "jaccard", settings)
  tibble(
    n_common_sites = length(sites),
    bray_stability = if (!is.null(db) && !is.null(ds)) safe_cor(as.vector(db), as.vector(ds)) else NA_real_,
    sorensen_stability = if (!is.null(sb) && !is.null(ss)) safe_cor(as.vector(sb), as.vector(ss)) else NA_real_,
    jaccard_stability = if (!is.null(jb) && !is.null(js)) safe_cor(as.vector(jb), as.vector(js)) else NA_real_,
    bray_mean_change_pct = if (!is.null(db) && !is.null(ds) && mean(db) > 0) 100*(mean(ds)/mean(db)-1) else NA_real_,
    sorensen_mean_change_pct = if (!is.null(sb) && !is.null(ss) && mean(sb) > 0) 100*(mean(ss)/mean(sb)-1) else NA_real_,
    jaccard_mean_change_pct = if (!is.null(jb) && !is.null(js) && mean(jb) > 0) 100*(mean(js)/mean(jb)-1) else NA_real_,
    ordination_procrustes_r2 = procrustes_r2(db, ds, settings)
  )
}

read_env_meta <- function(settings) {
  if (!settings$run_gdm) return(tibble())
  if (!file.exists(settings$env_input_file)) stop("Environmental input missing: ", settings$env_input_file)
  x <- readr::read_csv(settings$env_input_file, show_col_types = FALSE) %>% janitor::clean_names()
  needed <- c("station", "longitude", "latitude", settings$env_predictors)
  missing <- setdiff(needed, names(x))
  if (length(missing)) stop("Environmental input missing: ", paste(missing, collapse=", "))
  x <- x %>%
    mutate(station = as.character(station), across(all_of(c("longitude","latitude",settings$env_predictors)), as.numeric)) %>%
    group_by(station) %>%
    summarise(longitude = first_non_missing(longitude), latitude = first_non_missing(latitude),
              across(all_of(settings$env_predictors), first_non_missing), .groups="drop") %>%
    filter(complete.cases(.))
  sf_obj <- sf::st_as_sf(x, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>% sf::st_transform(2154)
  xy <- sf::st_coordinates(sf_obj)
  x %>% mutate(x_l93 = xy[,1], y_l93 = xy[,2])
}

empty_gdm <- function() tibble(
  gdm_n_sites = NA_real_, gdm_explained_baseline = NA_real_, gdm_explained_scenario = NA_real_,
  gdm_delta_explained = NA_real_, gdm_predicted_stability = NA_real_, gdm_term_rank_stability = NA_real_
)

fit_gdm <- function(mat, env, metric, settings) {
  if (!settings$run_gdm || !nrow(env)) return(NULL)
  common <- intersect(rownames(mat), env$station)
  mat <- mat[common,,drop=FALSE]
  meta <- env %>% filter(station %in% common) %>% arrange(match(station, common))
  mat <- mat[meta$station,,drop=FALSE]
  keep <- rowSums(mat)>0
  mat <- mat[keep,,drop=FALSE]; meta <- meta[keep,,drop=FALSE]
  mat <- mat[,colSums(mat)>0,drop=FALSE]
  if (nrow(mat) < settings$min_gdm_sites || ncol(mat)<2L) return(NULL)
  bio <- if (metric=="bray") mat else (mat>0)*1L
  bio_df <- data.frame(station=rownames(bio), x_l93=meta$x_l93, y_l93=meta$y_l93,
                       as.data.frame(bio, check.names=FALSE), check.names=FALSE)
  pred <- data.frame(meta[,c("station","x_l93","y_l93",settings$env_predictors),drop=FALSE], check.names=FALSE)
  sp <- try(gdm::formatsitepair(bioData=bio_df, bioFormat=1, dist="bray",
                                abundance=identical(metric,"bray"), siteColumn="station",
                                XColumn="x_l93",YColumn="y_l93",predData=pred,verbose=FALSE), silent=TRUE)
  if (inherits(sp,"try-error") || is.null(sp) || !nrow(sp)) return(NULL)
  mod <- try(gdm::gdm(sp, geo=TRUE), silent=TRUE)
  if (inherits(mod,"try-error") || is.null(mod$explained) || !is.finite(mod$explained)) return(NULL)
  list(model=mod,n_sites=nrow(mat))
}

compare_gdm <- function(base, scenario, env, metric, settings) {
  if (!settings$run_gdm) return(empty_gdm())
  common <- intersect(rownames(base)[rowSums(base)>0], rownames(scenario)[rowSums(scenario)>0])
  if (length(common)<settings$min_gdm_sites) return(empty_gdm())
  b <- fit_gdm(base[common,,drop=FALSE],env,metric,settings)
  s <- fit_gdm(scenario[common,,drop=FALSE],env,metric,settings)
  if (is.null(b)||is.null(s)) return(empty_gdm())
  tibble(
    gdm_n_sites = min(b$n_sites,s$n_sites),
    gdm_explained_baseline=b$model$explained,
    gdm_explained_scenario=s$model$explained,
    gdm_delta_explained=s$model$explained-b$model$explained,
    gdm_predicted_stability=safe_cor(b$model$predicted,s$model$predicted),
    gdm_term_rank_stability=NA_real_
  )
}
