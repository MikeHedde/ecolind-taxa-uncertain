read_assemblage <- function(assemblage_id, settings) {
  p <- file.path(settings$outputs_dir, "assemblages", assemblage_id)
  list(
    station_frame = readr::read_csv(paste0(p, "__station_frame.csv"), show_col_types = FALSE) %>% mutate(station=as.character(station)),
    rtu = readr::read_csv(paste0(p, "__rtu_long.csv"), show_col_types = FALSE) %>% mutate(station=as.character(station)),
    species = readr::read_csv(paste0(p, "__species_long.csv"), show_col_types = FALSE) %>% mutate(station=as.character(station)),
    lookup = readr::read_csv(paste0(p, "__taxon_lookup.csv"), show_col_types = FALSE),
    adult_path = paste0(p, "__adult_only_rtu_long.csv")
  )
}

lookup_maps <- function(lookup) {
  lookup %>%
    transmute(taxon_unit=taxon_unit_rtu, species_unit=taxon_unit_species,
              genus_unit=taxon_unit_genus, family_unit=taxon_unit_family,
              taxo_level, genus, family) %>%
    filter(!is.na(taxon_unit)) %>% distinct(taxon_unit,.keep_all=TRUE)
}

collapse_long <- function(x) x %>% group_by(station,taxon_unit) %>% summarise(abundance=sum(abundance),.groups="drop")

build_rare_to_genus <- function(rtu, maps, settings) {
  x <- rtu %>% left_join(maps,by="taxon_unit")
  rare <- x %>% filter(taxo_level=="species",!is.na(species_unit)) %>%
    group_by(species_unit) %>% summarise(total=sum(abundance),.groups="drop") %>%
    filter(total<=settings$rare_max_total_abundance) %>% pull(species_unit)
  x %>% mutate(taxon_unit=if_else(taxo_level=="species"&species_unit%in%rare&!is.na(genus_unit),genus_unit,taxon_unit)) %>%
    select(station,taxon_unit,abundance) %>% collapse_long()
}

build_rank_pair <- function(rtu,maps,rank=c("genus","family")) {
  rank <- match.arg(rank)
  target <- if(rank=="genus") "genus_unit" else "family_unit"
  x <- rtu %>% left_join(maps,by="taxon_unit") %>% filter(!is.na(.data[[target]]))
  list(
    baseline=x%>%select(station,taxon_unit,abundance)%>%collapse_long(),
    scenario=x%>%transmute(station,taxon_unit=.data[[target]],abundance)%>%collapse_long()
  )
}

apply_reporting_rules <- function(rtu, rules) {
  if (is.null(rules)||!nrow(rules)) return(rtu)
  rtu %>% left_join(rules %>% select(source_taxon_unit,target_taxon_unit),by=c("taxon_unit"="source_taxon_unit")) %>%
    mutate(taxon_unit=coalesce(target_taxon_unit,taxon_unit)) %>%
    select(station,taxon_unit,abundance) %>% collapse_long()
}

draw_map <- function(map, seed=NULL) {
  if (is.null(map)||!nrow(map)) return(tibble(source_taxon_unit=character(),target_taxon_unit=character()))
  if(!is.null(seed)) set.seed(seed)
  map %>% group_by(source_taxon_unit) %>%
    group_modify(~ slice_sample(.x,n=1,weight_by=weight)) %>%
    ungroup() %>% select(source_taxon_unit,target_taxon_unit)
}

simulate_reassignment <- function(comm, map, rate, seed=NULL, probabilities=NULL) {
  if(!is.null(seed)) set.seed(seed)
  if(is.null(map)||!nrow(map)) return(comm)
  x <- comm %>%
    left_join(map,by=c("taxon_unit"="source_taxon_unit")) %>%
    left_join(probabilities %||% tibble(taxon_unit=character(),p_error=numeric()),by="taxon_unit") %>%
    mutate(p_error=coalesce(p_error,rate), abundance=as.integer(round(abundance)))
  out <- purrr::pmap_dfr(x,function(station,taxon_unit,abundance,target_taxon_unit,p_error,...){
    if(is.na(target_taxon_unit)||abundance<=0||!is.finite(p_error)||p_error<=0)
      return(tibble(station=station,taxon_unit=taxon_unit,abundance=abundance))
    swapped <- stats::rbinom(1,size=abundance,prob=min(p_error,1))
    bind_rows(
      if(abundance-swapped>0) tibble(station=station,taxon_unit=taxon_unit,abundance=abundance-swapped),
      if(swapped>0) tibble(station=station,taxon_unit=target_taxon_unit,abundance=swapped)
    )
  })
  collapse_long(out)
}

observed_pool_map <- function(species_meta, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  purrr::map_dfr(seq_len(nrow(species_meta)),function(i){
    src<-species_meta[i,,drop=FALSE]
    candidates<-species_meta%>%filter(genus==src$genus,taxon_unit!=src$taxon_unit)
    tibble(source_taxon_unit=src$taxon_unit,target_taxon_unit=if(nrow(candidates)) sample(candidates$taxon_unit,1) else NA_character_)
  })
}

regional_pool_map <- function(species_meta, pool, source_audit, seed=NULL) {
  if(is.null(pool)||!nrow(pool)) return(tibble(source_taxon_unit=character(),target_taxon_unit=character()))
  if(!is.null(seed)) set.seed(seed)
  source_meta<-species_meta%>%left_join(source_audit%>%transmute(taxon_unit=observed_taxon_unit,source_cd_ref=cd_ref),by="taxon_unit")
  purrr::map_dfr(seq_len(nrow(source_meta)),function(i){
    src<-source_meta[i,,drop=FALSE]
    candidates<-pool%>%filter(genus==src$genus,taxon_unit!=src$taxon_unit,(is.na(src$source_cd_ref)|is.na(cd_ref)|cd_ref!=src$source_cd_ref))
    tibble(source_taxon_unit=src$taxon_unit,target_taxon_unit=if(nrow(candidates)) sample(candidates$taxon_unit,1) else NA_character_)
  })
}

rare_probabilities <- function(species, rate, cap) {
  x<-species%>%group_by(taxon_unit)%>%summarise(total=sum(abundance),.groups="drop")
  raw <- 1/sqrt(x$total)
  sf <- rate/stats::weighted.mean(raw,w=x$total)
  x%>%transmute(taxon_unit,p_error=pmin(raw*sf,cap))
}

scenario_long <- function(id, data, expert, regional_pool, source_audit, settings, rate, iter) {
  maps <- lookup_maps(data$lookup)
  species_meta <- maps %>% filter(!is.na(species_unit),!is.na(genus)) %>% transmute(taxon_unit=species_unit,genus)%>%distinct()
  if (id=="rtu_best_available") return(data$rtu)
  if (id=="rtu_drop_unresolved") return(data$species)
  if (id=="rtu_rare_to_genus") return(build_rare_to_genus(data$rtu,maps,settings))
  if (id=="adult_only") return(if(file.exists(data$adult_path)) readr::read_csv(data$adult_path,show_col_types=FALSE)%>%mutate(station=as.character(station)) else NULL)
  if (id %in% c("rtu_genus_resolvable","rtu_genus_level")) {
    x<-build_rank_pair(data$rtu,maps,"genus"); return(if(id=="rtu_genus_resolvable")x$baseline else x$scenario)
  }
  if (id %in% c("rtu_family_resolvable","rtu_family_level")) {
    x<-build_rank_pair(data$rtu,maps,"family"); return(if(id=="rtu_family_resolvable")x$baseline else x$scenario)
  }
  if (id=="species_baseline") return(data$species)

  requires_species_matrix <- id %in% c(
    "species_observed_pool_10",
    "species_regional_pool_10",
    "species_rare_weighted_10",
    "expert_species_10"
  )
  if (requires_species_matrix && (!nrow(data$species) || !nrow(species_meta))) return(NULL)

  if(id=="species_observed_pool_10") {
    return(simulate_reassignment(data$species,observed_pool_map(species_meta,10000+iter),rate,11000+iter))
  }
  if(id=="species_regional_pool_10") {
    return(simulate_reassignment(data$species,regional_pool_map(species_meta,regional_pool,source_audit,20000+iter),rate,21000+iter))
  }
  if(id=="species_rare_weighted_10") {
    return(simulate_reassignment(data$species,observed_pool_map(species_meta,30000+iter),rate,31000+iter,
                                 rare_probabilities(data$species,rate,settings$rare_error_cap)))
  }
  if(id=="expert_species_10") {
    if(is.null(expert$species_map) || !nrow(expert$species_map)) return(NULL)
    return(simulate_reassignment(data$species,draw_map(expert$species_map,40000+iter),rate,41000+iter))
  }
  if(id=="expert_reporting") {
    if(is.null(expert$reporting_rules) || !nrow(expert$reporting_rules)) return(NULL)
    return(apply_reporting_rules(data$rtu,expert$reporting_rules))
  }
  if(id=="expert_rtu_10") {
    if(is.null(expert$rtu_map) || !nrow(expert$rtu_map)) return(NULL)
    return(simulate_reassignment(data$rtu,draw_map(expert$rtu_map,50000+iter),rate,51000+iter))
  }
  if(id=="expert_integrated_10") {
    has_any <- (nrow(expert$species_map) + nrow(expert$reporting_rules) + nrow(expert$rtu_map)) > 0
    if(!has_any) return(NULL)
    x <- simulate_reassignment(data$rtu,draw_map(expert$species_map,60000+iter),rate,61000+iter)
    x <- apply_reporting_rules(x,expert$reporting_rules)
    x <- simulate_reassignment(x,draw_map(expert$rtu_map,62000+iter),rate,63000+iter)
    return(x)
  }
  stop("No scenario handler implemented for: ",id)
}
