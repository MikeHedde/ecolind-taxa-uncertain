# Main-text figure layer (v6).
# Results in outputs/sensitivity/figures/main_text/.

source("R/utils.R")
source("R/figure_helpers.R")
load_figure_packages()
settings <- read_settings()
inp <- read_figure_inputs(settings)

res <- inp$results
manifest <- inp$manifest
resolution <- inp$resolution
stage <- inp$stage

fig_dir <- file.path(inp$sensitivity_dir, "figures", "main_text")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
clear_generated_figures(fig_dir)

focal_ids <- manifest %>% filter(show_main_figure) %>% arrange(main_order) %>% pull(assemblage_id)
if (!length(focal_ids)) stop("No focal assemblages selected. Edit config/figure_assemblage_overrides.csv.")

focal_labels <- manifest %>%
  filter(assemblage_id %in% focal_ids) %>%
  arrange(main_order) %>%
  pull(display_label)

# -----------------------------------------------------------------------------
# Fig. 1A — reported taxonomic resolution.
# -----------------------------------------------------------------------------
resolution_plot <- resolution %>%
  filter(
    assemblage_id %in% focal_ids,
    resolution %in% c("species", "genus", "family", "order", "other_or_unusable")
  ) %>%
  left_join(manifest %>% select(assemblage_id, display_label, main_order), by = "assemblage_id") %>%
  mutate(
    display_label = factor(display_label, levels = focal_labels),
    resolution = factor(
      resolution,
      levels = c("species", "genus", "family", "order", "other_or_unusable"),
      labels = c("Species", "Genus", "Family", "Order", "Other / unusable")
    )
  )

p1a <- ggplot(resolution_plot, aes(display_label, abundance_share, fill = resolution)) +
  geom_col(width = 0.78) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Taxonomic resolution of the focal assemblages",
    subtitle = "The best-available dataset retains records at their reported resolution.",
    x = NULL, y = "Share of individuals", fill = "Reported rank"
  ) +
  theme_classic(10) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
save_figure(p1a, file.path(fig_dir, "Fig1A_taxonomic_resolution.png"), 185, 130)

# -----------------------------------------------------------------------------
# Fig. 1B — direct operational sources of uncertainty, restored from the
# pre-factorisation figure set.
# -----------------------------------------------------------------------------
unresolved <- resolution_plot %>%
  group_by(assemblage_id, display_label) %>%
  summarise(
    share = sum(abundance_share[resolution != "Species"], na.rm = TRUE),
    source = "Not resolved to species",
    .groups = "drop"
  )

stage_source <- tibble()
if ("juvenile_share" %in% names(stage)) {
  stage_source <- stage %>%
    filter(assemblage_id %in% focal_ids, is.finite(juvenile_share)) %>%
    left_join(manifest %>% select(assemblage_id, display_label, main_order), by = "assemblage_id") %>%
    group_by(assemblage_id, display_label) %>%
    summarise(share = dplyr::first(juvenile_share), .groups = "drop") %>%
    mutate(source = "Explicit non-adults")
}

uncertainty_sources <- bind_rows(unresolved, stage_source) %>%
  mutate(
    display_label = factor(display_label, levels = rev(focal_labels)),
    source = factor(source, levels = c("Not resolved to species", "Explicit non-adults"))
  )

if (nrow(uncertainty_sources)) {
  p1b <- ggplot(uncertainty_sources, aes(share, display_label)) +
    geom_vline(xintercept = 0, colour = "grey55") +
    geom_point(size = 2.8) +
    facet_wrap(~ source, ncol = 1, scales = "free_x") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Operational sources of taxonomic uncertainty",
      subtitle = "Non-specific reporting and explicit non-adults are distinct sources of uncertainty.",
      x = "Share of individuals", y = NULL
    ) +
    theme_classic(10) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA))

  save_figure(p1b, file.path(fig_dir, "Fig1B_sources_taxonomic_uncertainty.png"), 175, 130)
}

main_res <- main_results(res) %>% filter(assemblage_id %in% focal_ids)
divider <- error_workflow_divider(main_res)

# -----------------------------------------------------------------------------
# Fig. 2 — robustness heatmap.
# -----------------------------------------------------------------------------
stability_labels <- c(
  q0_stability = "Local richness (q0)",
  q1_stability = "Hill q1",
  q2_stability = "Hill q2",
  bray_stability = "Bray–Curtis",
  sorensen_stability = "Sørensen",
  ordination_procrustes_r2 = "Ordination",
  gdm_bray_predicted_stability = "Bray–Curtis GDM fitted turnover"
)

heat <- summary_to_long(
  main_res, names(stability_labels), stability_labels,
  c("assemblage_id", "display_label", "assembly_main_order", "scenario", "scenario_short", "scenario_order")
) %>%
  mutate(
    instability = pmax(0, pmin(1, 1 - median)),
    display_label = factor(display_label, levels = focal_labels),
    metric_label = factor(metric_label, levels = unname(stability_labels))
  )

if (nrow(heat)) {
  p2 <- ggplot(heat, aes(scenario_short, metric_label, fill = instability)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    { if (is.finite(divider)) geom_vline(xintercept = divider, linetype = "dashed", colour = "grey40") } +
    facet_wrap(~ display_label, nrow = 1) +
    scale_fill_gradient(low = "white", high = "black", limits = c(0, 1), na.value = "grey82", name = "Loss of\nstability") +
    labs(
      title = "Most inference is robust to moderate identification error, but not to broad taxonomic coarsening",
      subtitle = "Vertical divider: identification-error mechanisms (left) versus taxonomic workflow decisions (right). Grey cells indicate unavailable estimates.",
      x = NULL, y = NULL
    ) +
    theme_classic(9) +
    theme(
      axis.text.x = element_text(angle = 38, hjust = 1),
      strip.background = element_rect(fill = "grey92", colour = NA),
      panel.spacing.x = unit(3, "mm")
    )
  save_figure(p2, file.path(fig_dir, "Fig2_main_robustness_heatmap.png"), 255, 145)
}

# -----------------------------------------------------------------------------
# Fig. 3 — community-level directional effects.
# -----------------------------------------------------------------------------
effect_labels <- c(
  gamma_change_pct = "Regional richness (γ)",
  bray_mean_change_pct = "Mean Bray–Curtis",
  sorensen_mean_change_pct = "Mean Sørensen"
)

effects <- summary_to_long(
  main_res, names(effect_labels), effect_labels,
  c("assemblage_id", "display_label", "scenario", "scenario_short", "scenario_order")
) %>%
  mutate(metric_label = factor(metric_label, levels = unname(effect_labels)))

if (nrow(effects)) {
  p3 <- ggplot(effects, aes(scenario_short, median, colour = display_label)) +
    geom_hline(yintercept = 0, colour = "grey45") +
    { if (is.finite(divider)) geom_vline(xintercept = divider, linetype = "dashed", colour = "grey40") } +
    geom_errorbar(aes(ymin = p10, ymax = p90), width = 0.12, na.rm = TRUE) +
    geom_point(size = 2.25, na.rm = TRUE) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 3) +
    labs(
      title = "Workflow decisions shift richness and community dissimilarity more than moderate identification errors",
      subtitle = "Points are medians; intervals show 10th–90th percentiles for stochastic scenarios.",
      x = NULL, y = "Change relative to the matched baseline (%)", colour = "Assemblage"
    ) +
    theme_classic(9) +
    theme(
      axis.text.x = element_text(angle = 38, hjust = 1),
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.position = "bottom"
    )
  save_figure(p3, file.path(fig_dir, "Fig3_directional_process_effects.png"), 255, 125)
}

# -----------------------------------------------------------------------------
# Fig. 4 — Bray–Curtis GDM only. Sørensen GDM diagnostics remain supplementary.
# -----------------------------------------------------------------------------
gdm_labels <- c(
  gdm_bray_delta_explained = "Bray–Curtis GDM: change in explained deviance (pp)",
  gdm_bray_predicted_stability = "Bray–Curtis GDM: stability of fitted turnover"
)

gdm <- summary_to_long(
  main_res, names(gdm_labels), gdm_labels,
  c("assemblage_id", "display_label", "scenario", "scenario_short", "scenario_order")
) %>%
  mutate(
    reference = reference_value(metric),
    metric_label = factor(metric_label, levels = unname(gdm_labels))
  )

if (nrow(gdm)) {
  p4 <- ggplot(gdm, aes(scenario_short, median, colour = display_label)) +
    geom_hline(
      data = gdm %>% distinct(metric, metric_label, reference),
      aes(yintercept = reference), inherit.aes = FALSE, colour = "grey45"
    ) +
    { if (is.finite(divider)) geom_vline(xintercept = divider, linetype = "dashed", colour = "grey40") } +
    geom_errorbar(aes(ymin = p10, ymax = p90), width = 0.12, na.rm = TRUE) +
    geom_point(size = 2.25, na.rm = TRUE) +
    facet_wrap(~ metric_label, scales = "free_y", ncol = 1) +
    labs(
      title = "Bray–Curtis environmental-turnover inference is affected chiefly by taxonomic coarsening",
      subtitle = "Explained deviance is in percentage points; fitted turnover is compared by correlation.",
      x = NULL, y = NULL, colour = "Assemblage"
    ) +
    theme_classic(9) +
    theme(axis.text.x = element_text(angle = 38, hjust = 1), legend.position = "bottom")
  save_figure(p4, file.path(fig_dir, "Fig4_Bray_GDM_summary.png"), 220, 175)
}

# -----------------------------------------------------------------------------
# Fig. 5A and Fig. 5B — alpha–gamma synthesis split by mechanism family.
# -----------------------------------------------------------------------------
blowes <- main_res %>%
  filter(is.finite(mean_q0_change_pct_median), is.finite(gamma_change_pct_median)) %>%
  mutate(
    scenario_group = case_when(
      scenario_class %in% c("Identification error", "Expert-informed", "Expert component") ~ "Identification-error mechanisms",
      TRUE ~ "Workflow decisions"
    )
  )

make_blowes <- function(data, title, path, width) {
  if (!nrow(data)) return(invisible(NULL))
  p <- ggplot(data, aes(mean_q0_change_pct_median, gamma_change_pct_median, colour = display_label)) +
    geom_hline(yintercept = 0, colour = "grey45") +
    geom_vline(xintercept = 0, colour = "grey45") +
    geom_abline(slope = 1, linetype = "dashed") +
    geom_point(size = 2.7) +
    facet_wrap(~ scenario_short, scales = "fixed") +
    labs(
      title = title,
      subtitle = "Dashed diagonal: unchanged mean occupancy. Below: apparent homogenisation. Above: apparent differentiation.",
      x = "Change in mean local richness, Δα (%)",
      y = "Change in regional richness, Δγ (%)",
      colour = "Assemblage"
    ) +
    theme_classic(9) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = NA),
      legend.position = "bottom"
    )
  save_figure(p, path, width, 150)
}

make_blowes(
  blowes %>% filter(scenario_group == "Identification-error mechanisms"),
  "Alpha–gamma synthesis: identification-error mechanisms",
  file.path(fig_dir, "Fig5A_alpha_gamma_identification_errors.png"),
  210
)
make_blowes(
  blowes %>% filter(scenario_group == "Workflow decisions"),
  "Alpha–gamma synthesis: workflow decisions",
  file.path(fig_dir, "Fig5B_alpha_gamma_workflow_decisions.png"),
  245
)

message("Refined main-text figures written to: ", fig_dir)
