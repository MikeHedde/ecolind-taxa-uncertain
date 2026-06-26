# =============================================================================
# 06_refresh_multitaxa_figures_v3_with_taxref.R
# Regenerate manuscript Figures 2–5 after running the coherent-confusion v3
# pipeline, including regional TAXREF and curated expert-confusion scenarios.
#
# Prerequisite:
#   01_run_multitaxa_uncertainty_2024_v3_taxref_coherent_confusion.R
#
# This script reads results only. It does not rerun simulations or GDMs.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(stringr)
})

# ---- paths -------------------------------------------------------------------
INPUT_DIR <- "outputs_multitaxa_2024"
RESULT_DIR <- file.path(INPUT_DIR, "uncertainty_results")
FIG_DIR <- file.path(RESULT_DIR, "publication_figures_v3")
MAIN_DIR <- file.path(FIG_DIR, "main_text")
SUPP_DIR <- file.path(FIG_DIR, "supplement")
dir.create(MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

MAIN_ASSEMBLAGES <- c(
  "collembola_soil_core", "araneae_pitfall", "carabidae_pitfall",
  "formicidae_pitfall", "isopoda_pitfall"
)

# Scenario labels and conceptual families.
scenario_info <- tibble::tribble(
  ~scenario, ~scenario_label, ~scenario_family_display,
  "species_observed_congeneric_10pct", "Observed-pool\nerror 10%", "Observed-pool identification error",
  "species_observed_rare_weighted_10pct", "Rare-weighted\nerror 10%", "Observed-pool identification error",
  "species_regional_congeneric_10pct", "Regional TAXREF-pool\nerror 10%", "Regional-pool identification error",
  "species_expert_confusion_10pct", "Expert\nconfusion 10%", "Expert identification error",
  "rtu_rare_species_to_genus", "Rare taxa\nreported at genus", "Reporting decision",
  "rtu_species_only_drop_unresolved", "Drop unresolved\nrecords", "Reporting decision",
  "rtu_genus_level", "Genus-level\nresolution", "Taxonomic coarsening",
  "rtu_family_level", "Family-level\nresolution", "Taxonomic coarsening"
)

workflow_scenarios <- c(
  "rtu_rare_species_to_genus",
  "rtu_species_only_drop_unresolved",
  "rtu_genus_level",
  "rtu_family_level"
)
observed_error_scenarios <- c(
  "species_observed_congeneric_10pct",
  "species_observed_rare_weighted_10pct"
)
regional_error_scenarios <- c(
  "species_regional_congeneric_10pct",
  "species_expert_confusion_10pct"
)

assemblage_labels <- c(
  collembola_soil_core = "Collembola\nsoil cores",
  araneae_pitfall = "Araneae\npitfall traps",
  carabidae_pitfall = "Carabidae\npitfall traps",
  formicidae_pitfall = "Formicidae\npitfall traps",
  isopoda_pitfall = "Isopoda\npitfall traps",
  isopoda_hand_sorting = "Isopoda\nhand sorting",
  diplopoda_pitfall = "Diplopoda\npitfall traps",
  diplopoda_hand_sorting = "Diplopoda\nhand sorting"
)

cols_id <- c(
  collembola_soil_core = "#009E73",
  araneae_pitfall = "#CC79A7",
  carabidae_pitfall = "#D55E00",
  formicidae_pitfall = "#0072B2",
  isopoda_pitfall = "#56B4E9",
  isopoda_hand_sorting = "#E69F00",
  diplopoda_pitfall = "#7A7A7A",
  diplopoda_hand_sorting = "#4D4D4D"
)
cols <- setNames(unname(cols_id), unname(assemblage_labels[names(cols_id)]))

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = rel(1.15)),
      plot.subtitle = element_text(hjust = 0, size = rel(.95)),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold", size = rel(.92)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(colour = "grey91", linewidth = .25),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(7, 8, 7, 8)
    )
}

save_main <- function(p, filename, width, height) {
  ggsave(file.path(MAIN_DIR, paste0(filename, ".pdf")), p,
         width = width, height = height, units = "mm", device = grDevices::pdf)
  ggsave(file.path(MAIN_DIR, paste0(filename, ".png")), p,
         width = width, height = height, units = "mm", dpi = 500)
}

summarise_iter <- function(data, value, groups) {
  data %>%
    filter(is.finite(.data[[value]])) %>%
    group_by(across(all_of(groups))) %>%
    summarise(
      median = median(.data[[value]], na.rm = TRUE),
      p10 = quantile(.data[[value]], .10, na.rm = TRUE, names = FALSE),
      p90 = quantile(.data[[value]], .90, na.rm = TRUE, names = FALSE),
      n_iter = n(),
      .groups = "drop"
    )
}

# ---- inputs ------------------------------------------------------------------
raw_path <- file.path(RESULT_DIR, "results_by_iter_all.csv")
manifest_path <- file.path(INPUT_DIR, "assemblage_manifest.csv")
if (!file.exists(raw_path)) stop("Missing v2 result table: ", raw_path)

raw <- readr::read_csv(raw_path, show_col_types = FALSE)
manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)

meta <- manifest %>%
  transmute(
    assemblage_id,
    assemblage_label = recode(assemblage_id, !!!assemblage_labels, .default = assemblage_id)
  )

available_scenarios <- sort(unique(raw$scenario))
scenario_info <- scenario_info %>% filter(scenario %in% available_scenarios)
plot_scenarios <- scenario_info$scenario

raw_contrasts <- raw %>%
  filter(!is.na(baseline_scenario), scenario != baseline_scenario) %>%
  filter(assemblage_id %in% MAIN_ASSEMBLAGES, scenario %in% plot_scenarios) %>%
  left_join(meta, by = "assemblage_id") %>%
  left_join(scenario_info, by = "scenario") %>%
  mutate(
    assemblage_label = factor(assemblage_label, levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])),
    scenario_label = factor(scenario_label, levels = scenario_info$scenario_label)
  )

# =============================================================================
# Fig 2: balanced robustness heatmap
# =============================================================================
metric_map <- tibble::tribble(
  ~metric_key, ~metric_label, ~domain,
  "q0_stability", "q0 richness", "Alpha diversity",
  "q1_stability", "Hill q1", "Alpha diversity",
  "q2_stability", "Hill q2", "Alpha diversity",
  "bray_stability", "Bray–Curtis", "Beta diversity",
  "sorensen_stability", "Sørensen", "Beta diversity",
  "ordination_procrustes_r2", "Ordination", "Beta diversity",
  "gdm_bray_predicted_stability", "GDM Bray", "Ecological inference",
  "gdm_sorensen_predicted_stability", "GDM Sørensen", "Ecological inference"
)

heat <- raw_contrasts %>%
  select(
    assemblage_id, assemblage_label, scenario, scenario_label,
    q0_stability, q1_stability, q2_stability,
    bray_stability, sorensen_stability, ordination_procrustes_r2,
    gdm_bray_predicted_stability, gdm_sorensen_predicted_stability
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, assemblage_label, scenario, scenario_label),
    names_to = "metric_key", values_to = "stability"
  ) %>%
  left_join(metric_map, by = "metric_key") %>%
  summarise_iter(
    value = "stability",
    groups = c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "metric_key", "metric_label", "domain")
  )

grid <- tidyr::expand_grid(
  assemblage_id = MAIN_ASSEMBLAGES,
  scenario = scenario_info$scenario,
  metric_key = metric_map$metric_key
) %>%
  left_join(meta, by = "assemblage_id") %>%
  left_join(scenario_info %>% select(scenario, scenario_label), by = "scenario") %>%
  left_join(metric_map, by = "metric_key") %>%
  left_join(heat %>% select(assemblage_id, scenario, metric_key, median),
            by = c("assemblage_id", "scenario", "metric_key")) %>%
  mutate(
    loss = pmax(0, pmin(1, 1 - median)),
    assemblage_label = factor(assemblage_label, levels = unname(assemblage_labels[MAIN_ASSEMBLAGES])),
    scenario_label = factor(scenario_label, levels = scenario_info$scenario_label),
    metric_label = factor(metric_label, levels = rev(metric_map$metric_label)),
    domain = factor(domain, levels = rev(c("Alpha diversity", "Beta diversity", "Ecological inference")))
  )

p2 <- ggplot(grid, aes(x = scenario_label, y = metric_label, fill = loss)) +
  geom_tile(colour = "white", linewidth = .35) +
  facet_grid(domain ~ assemblage_label, scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = c("#FFF7EC", "#FEC44F", "#F03B20", "#7F0000"),
    limits = c(0, 1), oob = scales::squish, na.value = "grey86",
    name = "Loss of\nstability"
  ) +
  labs(
    title = "Identification errors and taxonomic workflows affect ecological inference differently",
    subtitle = "Stability is evaluated relative to the appropriate scenario-specific baseline; grey cells are non-applicable or non-estimable.",
    x = NULL, y = NULL
  ) +
  theme_paper(8.1) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.spacing = grid::unit(1.4, "mm"))

save_main(p2, "Fig2_cross_taxon_robustness_v3", 210, 155)

# =============================================================================
# Fig 3: directional effects by process family
# =============================================================================
directional <- raw_contrasts %>%
  transmute(
    assemblage_id, assemblage_label, scenario, scenario_label, scenario_family_display,
    `Regional richness (γ)` = gamma_change_pct,
    `Mean Bray–Curtis` = bray_mean_change_pct,
    `Mean Sørensen` = sorensen_mean_change_pct,
    `Environmental turnover explained (GDM)` = gdm_bray_delta_explained
  ) %>%
  pivot_longer(
    cols = -c(assemblage_id, assemblage_label, scenario, scenario_label, scenario_family_display),
    names_to = "response", values_to = "change"
  ) %>%
  summarise_iter(
    value = "change",
    groups = c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "scenario_family_display", "response")
  ) %>%
  mutate(
    response = factor(
      response,
      levels = c("Regional richness (γ)", "Mean Bray–Curtis", "Mean Sørensen", "Environmental turnover explained (GDM)")
    ),
    scenario_label = factor(scenario_label, levels = scenario_info$scenario_label)
  )

p3 <- ggplot(
  directional,
  aes(x = scenario_label, y = median, colour = assemblage_label)
) +
  geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.35) +
  geom_vline(
    xintercept = c(2.5, 4.5),
    colour = "grey60", linewidth = 0.35, linetype = "dashed"
  ) +
  geom_errorbar(
    aes(ymin = p10, ymax = p90),
    width = 0, linewidth = 0.4,
    position = position_dodge(width = 0.45),
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.0,
    position = position_dodge(width = 0.45),
    na.rm = TRUE
  ) +
  facet_wrap(~ response, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = cols, name = "Assemblage") +
  labs(
    title = "The direction of taxonomic artefacts depends on the error or workflow mechanism",
    subtitle = paste(
      "Scenarios are ordered as observed-pool errors, regional/expert errors and workflow decisions.",
      "Points are medians; intervals show 10th–90th percentiles."
    ),
    x = NULL,
    y = "Change relative to the appropriate baseline (%)\n(GDM: percentage points)"
  ) +
  theme_paper(8.0) +
  theme(
    axis.text.x = element_text(angle = 38, hjust = 1),
    panel.spacing = grid::unit(5, "mm")
  )

save_main(p3, "Fig3_directional_process_effects_v3", 205, 160)

# =============================================================================
# Fig 4: GDM ecological-inference heatmaps
# =============================================================================
gdm <- raw_contrasts %>%
  transmute(
    assemblage_id, assemblage_label, scenario, scenario_label,
    `Bray–Curtis` = gdm_bray_delta_explained,
    `Sørensen` = gdm_sorensen_delta_explained,
    `Bray–Curtis__prediction_loss` = 1 - gdm_bray_predicted_stability,
    `Sørensen__prediction_loss` = 1 - gdm_sorensen_predicted_stability
  )

gdm_fit <- gdm %>%
  pivot_longer(cols = c(`Bray–Curtis`, `Sørensen`), names_to = "distance", values_to = "delta_explained") %>%
  summarise_iter("delta_explained", c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "distance")) %>%
  mutate(
    assemblage_label = factor(assemblage_label, levels = rev(unname(assemblage_labels[MAIN_ASSEMBLAGES]))),
    scenario_label = factor(scenario_label, levels = scenario_info$scenario_label)
  )

gdm_pred <- gdm %>%
  pivot_longer(cols = ends_with("prediction_loss"), names_to = "distance", values_to = "prediction_loss") %>%
  mutate(distance = str_remove(distance, "__prediction_loss")) %>%
  summarise_iter("prediction_loss", c("assemblage_id", "assemblage_label", "scenario", "scenario_label", "distance")) %>%
  mutate(
    assemblage_label = factor(assemblage_label, levels = rev(unname(assemblage_labels[MAIN_ASSEMBLAGES]))),
    scenario_label = factor(scenario_label, levels = scenario_info$scenario_label)
  )

p4a <- ggplot(gdm_fit, aes(x = scenario_label, y = assemblage_label, fill = median)) +
  geom_tile(colour = "white", linewidth = .38) +
  geom_text(aes(label = ifelse(is.na(median), "", sprintf("%+.1f", median))), size = 2.5) +
  facet_wrap(~ distance, ncol = 1) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, na.value = "grey86",
                       name = "Δ explained deviance\n(percentage points)") +
  labs(title = "GDM fit", x = NULL, y = NULL) +
  theme_paper(7.8) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

p4b <- ggplot(gdm_pred, aes(x = scenario_label, y = assemblage_label, fill = median)) +
  geom_tile(colour = "white", linewidth = .38) +
  geom_text(aes(label = ifelse(is.na(median), "", sprintf("%.2f", median))), size = 2.5) +
  facet_wrap(~ distance, ncol = 1) +
  scale_fill_gradientn(
    colours = c("#FFF7EC", "#FEC44F", "#F03B20", "#7F0000"),
    limits = c(0, 1), oob = scales::squish, na.value = "grey86",
    name = "Loss of predicted-\ndissimilarity stability"
  ) +
  labs(title = "GDM prediction stability", x = NULL, y = NULL) +
  theme_paper(7.8) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

p4 <- (p4a + p4b) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Taxonomic uncertainty can modify environmental turnover inference",
    subtitle = "Each scenario is compared with the GDM fitted to the same assemblage and sites.",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_main(p4, "Fig4_GDM_ecological_inference_v3", 210, 155)

# =============================================================================
# Fig 5: Blowes alpha-gamma-occupancy synthesis
# =============================================================================
summarise_blowes <- function(x) {
  x %>%
    filter(is.finite(mean_q0_change_pct), is.finite(gamma_change_pct)) %>%
    mutate(
      alpha_ratio = 1 + mean_q0_change_pct / 100,
      gamma_ratio = 1 + gamma_change_pct / 100
    ) %>%
    filter(alpha_ratio > 0, gamma_ratio > 0) %>%
    mutate(occupancy_ratio = alpha_ratio / gamma_ratio) %>%
    group_by(assemblage_id, assemblage_label, scenario, scenario_label) %>%
    summarise(
      alpha = 100 * (median(alpha_ratio) - 1),
      gamma = 100 * (median(gamma_ratio) - 1),
      alpha_p10 = 100 * (quantile(alpha_ratio, .10) - 1),
      alpha_p90 = 100 * (quantile(alpha_ratio, .90) - 1),
      gamma_p10 = 100 * (quantile(gamma_ratio, .10) - 1),
      gamma_p90 = 100 * (quantile(gamma_ratio, .90) - 1),
      occupancy_change = 100 * (median(occupancy_ratio) - 1),
      .groups = "drop"
    )
}

plot_blowes <- function(dat, title, subtitle, ncol = 2, minspan = 30) {
  vals <- c(dat$alpha_p10, dat$alpha_p90, dat$gamma_p10, dat$gamma_p90, 0)
  vals <- vals[is.finite(vals)]
  lo <- min(vals); hi <- max(vals)
  span <- max(hi - lo, minspan)
  lims <- c(floor((lo - .08 * span) / 10) * 10, ceiling((hi + .08 * span) / 10) * 10)
  
  ggplot(dat, aes(x = alpha, y = gamma, colour = assemblage_label)) +
    geom_hline(yintercept = 0, colour = "grey45", linewidth = .32) +
    geom_vline(xintercept = 0, colour = "grey45", linewidth = .32) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = .6) +
    geom_segment(aes(x = alpha_p10, xend = alpha_p90, y = gamma, yend = gamma), linewidth = .45) +
    geom_segment(aes(x = alpha, xend = alpha, y = gamma_p10, yend = gamma_p90), linewidth = .45) +
    geom_point(size = 2.8) +
    facet_wrap(~ scenario_label, ncol = ncol) +
    coord_equal(xlim = lims, ylim = lims, expand = FALSE) +
    scale_colour_manual(values = cols, name = "Assemblage") +
    labs(title = title, subtitle = subtitle,
         x = "Change in mean local richness, Δα (%)",
         y = "Change in regional richness, Δγ (%)") +
    theme_paper(8.2)
}

mk_blowes <- function(scenarios) {
  raw_contrasts %>%
    filter(scenario %in% scenarios) %>%
    summarise_blowes() %>%
    mutate(
      scenario_label = factor(scenario_label, levels = scenario_info$scenario_label[match(scenarios, scenario_info$scenario)]),
      assemblage_label = factor(assemblage_label, levels = unname(assemblage_labels[MAIN_ASSEMBLAGES]))
    )
}

b_workflow <- mk_blowes(intersect(workflow_scenarios, available_scenarios))
b_observed <- mk_blowes(intersect(observed_error_scenarios, available_scenarios))
b_regional <- mk_blowes(intersect(regional_error_scenarios, available_scenarios))

p5a <- plot_blowes(
  b_workflow,
  "Workflow decisions and taxonomic coarsening",
  "Below the dashed line: apparent homogenisation through a greater loss of γ than α.",
  ncol = 2, minspan = 40
)
p5b <- plot_blowes(
  b_observed,
  "Observed-pool identification error",
  "Errors are constrained to species already detected in the RMQS assemblage.",
  ncol = 2, minspan = 30
)

p5_parts <- list(p5a, p5b)
if (nrow(b_regional) > 0L) {
  p5c <- plot_blowes(
    b_regional,
    "Regional TAXREF-pool and expert identification error",
    "New regional candidate taxa can yield a contrasting alpha–gamma signature.",
    ncol = 2, minspan = 40
  )
  p5_parts <- append(p5_parts, list(p5c))
}

p5 <- wrap_plots(p5_parts, ncol = 1) +
  plot_annotation(
    title = "Taxonomic processes can shift apparent occupancy patterns in alpha–gamma space",
    subtitle = "The dashed line denotes unchanged mean occupancy; each point is a scenario median and crosshairs are 10th–90th percentiles.",
    tag_levels = "A"
  ) &
  theme(legend.position = "bottom")

save_main(p5, "Fig5_multitaxa_Blowes_synthesis_v3", 195, if (length(p5_parts) == 3) 300 else 225)

message("\nFinished v2 figures:")
message(MAIN_DIR)
