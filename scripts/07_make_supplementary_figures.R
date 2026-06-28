# Supplementary figures (v6).
# Results in outputs/sensitivity/figures/supplement/.

source("R/utils.R")
source("R/figure_helpers.R")
load_figure_packages()
settings <- read_settings()
inp <- read_figure_inputs(settings)

res <- inp$results
manifest <- inp$manifest
resolution <- inp$resolution
stage <- inp$stage

fig_dir <- file.path(inp$sensitivity_dir, "figures", "supplement")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
clear_generated_figures(fig_dir)

# -----------------------------------------------------------------------------
# Fig. S1 — coverage audit across all discovered assemblages.
# -----------------------------------------------------------------------------
audit <- manifest %>%
  mutate(display_label = factor(display_label, levels = display_label[order(supp_order, display_label)])) %>%
  select(display_label, n_sampled_stations, n_positive_stations, n_rtu, n_species) %>%
  pivot_longer(-display_label, names_to = "measure", values_to = "value") %>%
  mutate(measure = recode(
    measure,
    n_sampled_stations = "Sampled stations",
    n_positive_stations = "Positive stations",
    n_rtu = "Best-available units",
    n_species = "Species-level units"
  ))

p_s1 <- ggplot(audit, aes(reorder(display_label, value), value)) +
  geom_col() +
  facet_wrap(~ measure, scales = "free_x", ncol = 2) +
  coord_flip() +
  labs(title = "Assemblage audit: sampling coverage and taxonomic richness", x = NULL, y = NULL) +
  theme_classic(9) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA))
save_figure(p_s1, file.path(fig_dir, "FigS1_assemblage_audit.png"), 220, 180)

# -----------------------------------------------------------------------------
# Fig. S2 — resolution profiles across all assemblages.
# -----------------------------------------------------------------------------
resol_all <- resolution %>%
  filter(resolution %in% c("species", "genus", "family", "order", "other_or_unusable")) %>%
  left_join(manifest %>% select(assemblage_id, display_label, supp_order), by = "assemblage_id") %>%
  mutate(
    display_label = factor(display_label, levels = manifest %>% arrange(supp_order, display_label) %>% pull(display_label)),
    resolution = factor(
      resolution,
      levels = c("species", "genus", "family", "order", "other_or_unusable"),
      labels = c("Species", "Genus", "Family", "Order", "Other / unusable")
    )
  )

p_s2 <- ggplot(resol_all, aes(display_label, abundance_share, fill = resolution)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Full audit of taxonomic resolution", x = NULL, y = "Share of individuals", fill = "Reported rank") +
  theme_classic(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_figure(p_s2, file.path(fig_dir, "FigS2_taxonomic_resolution_all_assemblages.png"), 250, 145)

supp_res <- supplement_results(res)

# -----------------------------------------------------------------------------
# Fig. S3 — full heatmap.
# -----------------------------------------------------------------------------
stability_labels <- c(
  q0_stability = "Local richness (q0)",
  q1_stability = "Hill q1",
  q2_stability = "Hill q2",
  bray_stability = "Bray–Curtis",
  sorensen_stability = "Sørensen",
  jaccard_stability = "Jaccard",
  ordination_procrustes_r2 = "Ordination",
  gdm_bray_predicted_stability = "Bray–Curtis GDM fitted turnover",
  gdm_sorensen_predicted_stability = "Sørensen GDM fitted turnover"
)

heat <- summary_to_long(
  supp_res, names(stability_labels), stability_labels,
  c("assemblage_id", "display_label", "assembly_supp_order", "scenario", "scenario_short", "scenario_order", "scenario_class")
) %>%
  mutate(
    instability = pmax(0, pmin(1, 1 - median)),
    display_label = factor(display_label, levels = manifest %>% arrange(supp_order, display_label) %>% pull(display_label)),
    metric_label = factor(metric_label, levels = unname(stability_labels))
  )

if (nrow(heat)) {
  p_s3 <- ggplot(heat, aes(scenario_short, metric_label, fill = instability)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_grid(scenario_class ~ display_label, scales = "free_x", space = "free_x") +
    scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), na.value = "grey82", name = "Loss of\nstability") +
    labs(title = "Complete robustness heatmap across taxonomic scenarios", x = NULL, y = NULL) +
    theme_classic(7.5) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey92", colour = NA),
      panel.spacing = unit(2, "mm")
    )
  save_figure(p_s3, file.path(fig_dir, "FigS3_full_robustness_heatmap.png"), 310, 245)
}

# -----------------------------------------------------------------------------
# Fig. S4 — expert components, only when expert scenarios were actually run.
# -----------------------------------------------------------------------------
expert_ids <- c("expert_species_10", "expert_reporting", "expert_rtu_10", "expert_integrated_10")
expert <- res %>% filter(scenario %in% expert_ids, scenario != baseline_scenario)

expert_labels <- c(
  gamma_change_pct = "Regional richness (γ)",
  mean_q0_change_pct = "Mean local richness (q0)",
  bray_mean_change_pct = "Mean Bray–Curtis",
  sorensen_mean_change_pct = "Mean Sørensen",
  gdm_bray_delta_explained = "Bray–Curtis GDM explained deviance (pp)"
)

expert_long <- summary_to_long(
  expert, names(expert_labels), expert_labels,
  c("assemblage_id", "display_label", "scenario", "scenario_short", "scenario_order")
)

if (nrow(expert_long)) {
  p_s4 <- ggplot(expert_long, aes(scenario_short, median, colour = display_label)) +
    geom_hline(yintercept = 0, colour = "grey45") +
    geom_errorbar(aes(ymin = p10, ymax = p90), width = 0.12, na.rm = TRUE) +
    geom_point(size = 2.1) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    labs(
      title = "Decomposition of expert-informed taxonomic uncertainty",
      subtitle = "Expert species confusion, conservative reporting, RTU confusion and their integrated scenario.",
      x = NULL, y = "Change relative to the appropriate baseline", colour = "Assemblage"
    ) +
    theme_classic(9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
  save_figure(p_s4, file.path(fig_dir, "FigS4_expert_mechanism_decomposition.png"), 245, 180)
}

# -----------------------------------------------------------------------------
# Fig. S5 — 1–20% gradients, now including rare-weighted observed-pool error.
# -----------------------------------------------------------------------------
gradient_path <- file.path(settings$outputs_dir, "error_gradient", "error_gradient_summary.csv")

if (file.exists(gradient_path)) {
  gradient <- readr::read_csv(gradient_path, show_col_types = FALSE)

  gradient_labels <- c(
    gamma_change_pct = "Regional richness (γ)",
    mean_q0_change_pct = "Mean local richness (q0)",
    bray_stability = "Bray–Curtis stability",
    sorensen_stability = "Sørensen stability"
  )

  grad <- summary_to_long(
    gradient, names(gradient_labels), gradient_labels,
    c("assemblage_id", "taxon_key", "mechanism", "error_rate")
  ) %>%
    left_join(manifest %>% select(assemblage_id, display_label, supp_order), by = "assemblage_id") %>%
    mutate(
      display_label = factor(display_label, levels = manifest %>% arrange(supp_order, display_label) %>% pull(display_label)),
      mechanism = recode(
        mechanism,
        observed_pool = "Observed pool",
        rare_weighted_observed_pool = "Rare-weighted observed pool",
        regional_pool = "Regional TAXREF pool",
        .default = mechanism
      ),
      mechanism = factor(
        mechanism,
        levels = c("Observed pool", "Rare-weighted observed pool", "Regional TAXREF pool")
      ),
      reference = reference_value(metric)
    )

  make_gradient <- function(data, title, path) {
    if (!nrow(data)) return(invisible(NULL))

    p <- ggplot(data, aes(error_rate * 100, median, colour = mechanism, group = mechanism)) +
      geom_hline(
        data = data %>% distinct(metric, metric_label, reference),
        aes(yintercept = reference), inherit.aes = FALSE, colour = "grey55"
      ) +
      geom_ribbon(
        aes(ymin = p10, ymax = p90, fill = mechanism),
        alpha = 0.14, colour = NA, show.legend = FALSE
      ) +
      geom_line(linewidth = 0.55) +
      facet_grid(metric_label ~ display_label, scales = "free_y") +
      labs(
        title = title,
        x = "Nominal identification-error rate (%)",
        y = NULL,
        colour = "Error mechanism"
      ) +
      theme_classic(7.5) +
      theme(
        strip.background = element_rect(fill = "grey92", colour = NA),
        axis.text.x = element_text(angle = 0),
        legend.position = "bottom"
      )

    save_figure(p, path, 310, 180)
  }

  make_gradient(
    grad %>% filter(metric %in% c("gamma_change_pct", "mean_q0_change_pct")),
    "Sensitivity to increasing identification-error rates: richness metrics",
    file.path(fig_dir, "FigS5A_error_rate_gradients_richness.png")
  )

  make_gradient(
    grad %>% filter(metric %in% c("bray_stability", "sorensen_stability")),
    "Sensitivity to increasing identification-error rates: beta-diversity stability",
    file.path(fig_dir, "FigS5B_error_rate_gradients_beta_stability.png")
  )
}

# -----------------------------------------------------------------------------
# Fig. S6 — adults-only sensitivity.
# -----------------------------------------------------------------------------
adult_labels <- c(
  gamma_change_pct = "Regional richness (γ)",
  mean_q0_change_pct = "Mean local richness (q0)",
  bray_mean_change_pct = "Mean Bray–Curtis",
  sorensen_mean_change_pct = "Mean Sørensen",
  q0_stability = "Local richness stability",
  bray_stability = "Bray–Curtis stability"
)

adult <- res %>% filter(scenario == "adult_only", scenario != baseline_scenario)
adult_long <- summary_to_long(
  adult, names(adult_labels), adult_labels,
  c("assemblage_id", "display_label", "scenario_short")
) %>% mutate(reference = reference_value(metric))

if (nrow(adult_long)) {
  p_s6 <- ggplot(adult_long, aes(display_label, median)) +
    geom_hline(
      data = adult_long %>% distinct(metric, metric_label, reference),
      aes(yintercept = reference), inherit.aes = FALSE, colour = "grey45"
    ) +
    geom_point(size = 2.4) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    labs(title = "Sensitivity to explicit removal of non-adults", x = NULL, y = "Change or stability relative to best-available RTUs") +
    theme_classic(9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_figure(p_s6, file.path(fig_dir, "FigS6_adult_only_sensitivity.png"), 220, 165)
}

# -----------------------------------------------------------------------------
# Fig. S7 — Jaccard versus Sørensen.
# -----------------------------------------------------------------------------
jaccard_labels <- c(
  sorensen_stability = "Sørensen stability",
  jaccard_stability = "Jaccard stability",
  sorensen_mean_change_pct = "Mean Sørensen change",
  jaccard_mean_change_pct = "Mean Jaccard change"
)

jaccard <- summary_to_long(
  supp_res, names(jaccard_labels), jaccard_labels,
  c("assemblage_id", "display_label", "scenario", "scenario_short", "scenario_order")
) %>% mutate(reference = reference_value(metric))

if (nrow(jaccard)) {
  p_s7 <- ggplot(jaccard, aes(scenario_short, median, colour = display_label)) +
    geom_hline(
      data = jaccard %>% distinct(metric, metric_label, reference),
      aes(yintercept = reference), inherit.aes = FALSE, colour = "grey45"
    ) +
    geom_point(size = 1.85) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 2) +
    labs(title = "Binary-dissimilarity sensitivity: Sørensen versus Jaccard", x = NULL, y = NULL, colour = "Assemblage") +
    theme_classic(8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  save_figure(p_s7, file.path(fig_dir, "FigS7_jaccard_sensitivity.png"), 250, 185)
}

# -----------------------------------------------------------------------------
# Fig. S8A/S8B — GDM diagnostics for both Bray–Curtis and Sørensen.
# -----------------------------------------------------------------------------
gdm_delta_labels <- c(
  gdm_bray_delta_explained = "Bray–Curtis GDM",
  gdm_sorensen_delta_explained = "Sørensen GDM"
)

gdm_delta <- res %>%
  filter(scenario != baseline_scenario) %>%
  summary_to_long(
    names(gdm_delta_labels),
    gdm_delta_labels,
    c("assemblage_id", "display_label", "assembly_supp_order", "scenario", "scenario_short", "scenario_order", "scenario_class")
  ) %>%
  mutate(
    display_label = factor(display_label, levels = manifest %>% arrange(supp_order, display_label) %>% pull(display_label)),
    metric_label = factor(metric_label, levels = unname(gdm_delta_labels))
  )

if (nrow(gdm_delta)) {
  p_s8a <- ggplot(gdm_delta, aes(scenario_short, display_label, fill = median)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    facet_grid(metric_label ~ scenario_class, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(midpoint = 0, name = "Δ explained\ndeviance (pp)", na.value = "grey82") +
    labs(title = "GDM sensitivity: change in explained deviance", x = NULL, y = NULL) +
    theme_classic(7.5) +
    theme(axis.text.x = element_text(angle = 42, hjust = 1), strip.background = element_rect(fill = "grey92", colour = NA))
  save_figure(p_s8a, file.path(fig_dir, "FigS8A_GDM_delta_explained_heatmap.png"), 290, 220)
}

gdm_stability_labels <- c(
  gdm_bray_predicted_stability = "Bray–Curtis GDM",
  gdm_sorensen_predicted_stability = "Sørensen GDM"
)

gdm_stability <- res %>%
  filter(scenario != baseline_scenario) %>%
  summary_to_long(
    names(gdm_stability_labels),
    gdm_stability_labels,
    c("assemblage_id", "display_label", "assembly_supp_order", "scenario", "scenario_short", "scenario_order", "scenario_class")
  ) %>%
  mutate(
    instability = pmax(0, pmin(1, 1 - median)),
    display_label = factor(display_label, levels = manifest %>% arrange(supp_order, display_label) %>% pull(display_label)),
    metric_label = factor(metric_label, levels = unname(gdm_stability_labels))
  )

if (nrow(gdm_stability)) {
  p_s8b <- ggplot(gdm_stability, aes(scenario_short, display_label, fill = instability)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    facet_grid(metric_label ~ scenario_class, scales = "free_x", space = "free_x") +
    scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), name = "Loss of GDM\nstability", na.value = "grey82") +
    labs(title = "GDM sensitivity: stability of fitted environmental turnover", x = NULL, y = NULL) +
    theme_classic(7.5) +
    theme(axis.text.x = element_text(angle = 42, hjust = 1), strip.background = element_rect(fill = "grey92", colour = NA))
  save_figure(p_s8b, file.path(fig_dir, "FigS8B_GDM_fitted_turnover_stability_heatmap.png"), 290, 220)
}

message("Refined supplementary figures written to: ", fig_dir)
