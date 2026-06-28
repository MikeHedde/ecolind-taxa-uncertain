source("R/utils.R")
source("R/expert_rules.R")
settings <- read_settings()
load_pipeline_packages(require_gdm = FALSE)
dir.create(settings$expert_internal_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- readr::read_csv(file.path(settings$outputs_dir,"assemblage_manifest.csv"),show_col_types=FALSE)
registry <- read_registry()

all_audits <- list()
summary <- list()
k <- 1L

for(i in seq_len(nrow(manifest))) {
  a <- manifest[i,]
  taxon_row <- registry %>% filter(taxon_key == a$taxon_key)
  if(!nrow(taxon_row)) next
  out <- project_expert_rules_one(a,taxon_row,settings)
  prefix <- file.path(settings$expert_internal_dir,a$assemblage_id)
  readr::write_csv(out$species_map,paste0(prefix,"__species_confusions.csv"))
  readr::write_csv(out$reporting_rules,paste0(prefix,"__reporting_rules.csv"))
  readr::write_csv(out$rtu_map,paste0(prefix,"__rtu_confusions.csv"))
  readr::write_csv(out$audit,paste0(prefix,"__expert_projection_audit.csv"))
  all_audits[[k]] <- out$audit
  summary[[k]] <- tibble(
    assemblage_id=a$assemblage_id,taxon_key=a$taxon_key,
    n_species_rules=nrow(out$species_map),
    n_reporting_rules=nrow(out$reporting_rules),
    n_rtu_rules=nrow(out$rtu_map),
    n_unprojected=sum(out$audit$status!="projected")
  )
  k <- k+1L
}
readr::write_csv(bind_rows(all_audits),file.path(settings$expert_internal_dir,"expert_projection_audit_all.csv"))
readr::write_csv(bind_rows(summary),file.path(settings$expert_internal_dir,"expert_projection_summary.csv"))
message("Expert rules projected. Review: ",file.path(settings$expert_internal_dir,"expert_projection_summary.csv"))
