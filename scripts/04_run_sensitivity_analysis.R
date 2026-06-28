source("R/utils.R")
source_pipeline_files()
settings <- read_settings()
load_pipeline_packages(require_gdm = settings$run_gdm)

out_dir <- file.path(settings$outputs_dir,"sensitivity")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out_dir,"per_assemblage"),recursive=TRUE,showWarnings=FALSE)

manifest <- readr::read_csv(file.path(settings$outputs_dir,"assemblage_manifest.csv"),show_col_types=FALSE)
catalog <- read_scenario_catalog() %>% filter(run)

env <- if(settings$run_gdm) read_env_meta(settings) else tibble()

read_map <- function(path) if(file.exists(path)) readr::read_csv(path,show_col_types=FALSE) else tibble(
  source_taxon_unit=character(),target_taxon_unit=character(),weight=numeric(),comment=character()
)

all_results <- list(); all_state <- list(); all_alpha <- list(); all_index <- list(); z <- 1L

for(i in seq_len(nrow(manifest))) {
  a <- manifest[i,]
  message("\nAnalysing ",a$assemblage_id)
  data <- read_assemblage(a$assemblage_id,settings)
  expert_prefix <- file.path(settings$expert_internal_dir,a$assemblage_id)
  expert <- list(
    species_map=read_map(paste0(expert_prefix,"__species_confusions.csv")),
    reporting_rules=read_map(paste0(expert_prefix,"__reporting_rules.csv")),
    rtu_map=read_map(paste0(expert_prefix,"__rtu_confusions.csv"))
  )
  pool_path <- file.path(settings$regional_pool_dir,paste0(a$assemblage_id,"__regional_pool.csv"))
  audit_path <- file.path(settings$regional_pool_dir,paste0(a$assemblage_id,"__observed_taxref_match_audit.csv"))
  pool <- if(file.exists(pool_path)) readr::read_csv(pool_path,show_col_types=FALSE) else tibble(taxon_unit=character(),cd_ref=character(),genus=character())
  source_audit <- if(file.exists(audit_path)) readr::read_csv(audit_path,show_col_types=FALSE) else tibble(observed_taxon_unit=character(),cd_ref=character())

  # Build every enabled scenario once per requested iteration.
  scen_objects <- list()
  for(j in seq_len(nrow(catalog))) {
    sc <- catalog[j,]
    n_iter <- if(str_detect(sc$handler,"error|expert") && !sc$handler %in% c("expert_reporting")) settings$n_sim else 1L
    for(iter in seq_len(n_iter)) {
      long <- try(scenario_long(sc$scenario_id,data,expert,pool,source_audit,settings,sc$error_rate,iter),silent=TRUE)
      if(inherits(long,"try-error") || is.null(long) || !nrow(long)) next
      scen_objects[[length(scen_objects)+1]] <- list(
        scenario=sc$scenario_id,scenario_family=sc$scenario_family,
        baseline_id=sc$baseline_id,unit_type=sc$handler,iter=iter,comm_long=long
      )
    }
  }
  # Keep scenario matrices keyed by scenario+iteration.
  for(obj in scen_objects) {
    base_candidates <- Filter(function(x) x$scenario==obj$baseline_id,scen_objects)
    if(!length(base_candidates)) next
    base_obj <- base_candidates[[min(obj$iter,length(base_candidates))]]
    mat <- make_matrix(obj$comm_long,data$station_frame)
    base_mat <- make_matrix(base_obj$comm_long,data$station_frame)
    alpha <- alpha_metrics(mat); alpha_base <- alpha_metrics(base_mat)
    beta <- compare_beta(base_mat,mat,settings)
    gdm_bray <- if(settings$run_gdm && (obj$iter<=settings$gdm_max_stochastic_iters || obj$iter==1L)) compare_gdm(base_mat,mat,env,"bray",settings) else empty_gdm()
    gdm_sor <- if(settings$run_gdm && (obj$iter<=settings$gdm_max_stochastic_iters || obj$iter==1L)) compare_gdm(base_mat,mat,env,"sorensen",settings) else empty_gdm()

    all_alpha[[z]] <- alpha %>% mutate(assemblage_id=a$assemblage_id,scenario=obj$scenario,iter=obj$iter)
    all_state[[z]] <- tibble(
      assemblage_id=a$assemblage_id,taxon_key=a$taxon_key,group_label=a$group_label,method=a$method,
      scenario=obj$scenario,scenario_family=obj$scenario_family,baseline_scenario=obj$baseline_id,iter=obj$iter,
      mean_total_abundance=mean(alpha$total_abundance),mean_q0=mean(alpha$q0),mean_q1=mean(alpha$q1),mean_q2=mean(alpha$q2),
      gamma=community_gamma(mat)
    )
    all_results[[z]] <- bind_cols(
      tibble(
        assemblage_id=a$assemblage_id,taxon_key=a$taxon_key,group_label=a$group_label,method=a$method,role=a$role,
        scenario=obj$scenario,scenario_family=obj$scenario_family,baseline_scenario=obj$baseline_id,
        unit_type=obj$unit_type,iter=obj$iter,
        gamma_change_pct=100*(community_gamma(mat)/community_gamma(base_mat)-1),
        mean_q0_change_pct=100*(mean(alpha$q0)/mean(alpha_base$q0)-1),
        mean_q1_change_pct=100*(mean(alpha$q1)/mean(alpha_base$q1)-1),
        mean_q2_change_pct=100*(mean(alpha$q2)/mean(alpha_base$q2)-1),
        q0_stability=safe_cor(alpha_base$q0,alpha$q0),
        q1_stability=safe_cor(alpha_base$q1,alpha$q1),
        q2_stability=safe_cor(alpha_base$q2,alpha$q2)
      ),beta,
      gdm_bray %>% rename_with(~paste0("gdm_bray_",sub("^gdm_","",.x))),
      gdm_sor %>% rename_with(~paste0("gdm_sorensen_",sub("^gdm_","",.x)))
    )
    z <- z+1L
  }
}

results <- bind_rows(all_results); state <- bind_rows(all_state); alpha <- bind_rows(all_alpha)
readr::write_csv(results,file.path(out_dir,"results_by_iter_all.csv"))
readr::write_csv(state,file.path(out_dir,"state_by_iter_all.csv"))
readr::write_csv(alpha,file.path(out_dir,"alpha_by_site_all.csv"))

numeric_cols <- names(results)[vapply(results,is.numeric,logical(1))]
numeric_cols <- setdiff(numeric_cols,c("iter"))
summary <- results %>% group_by(assemblage_id,taxon_key,group_label,method,role,scenario,scenario_family,baseline_scenario,unit_type) %>%
  summarise(across(all_of(numeric_cols),list(median=~median(.x,na.rm=TRUE),p10=~quantile(.x,.1,na.rm=TRUE,names=FALSE),p90=~quantile(.x,.9,na.rm=TRUE,names=FALSE)),.names="{.col}_{.fn}"),n_iter=n(),.groups="drop") %>%
  mutate(across(ends_with("_median")|ends_with("_p10")|ends_with("_p90"),~ifelse(is.nan(.x),NA_real_,.x)))
readr::write_csv(summary,file.path(out_dir,"results_summary.csv"))
readr::write_csv(catalog,file.path(out_dir,"scenario_catalog_used.csv"))
message("\nFinished. Results: ",out_dir)
