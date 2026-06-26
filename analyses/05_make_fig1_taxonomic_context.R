# =============================================================================
# 04_make_fig1_taxonomic_context.R
# Main-text Figure 1: taxonomic and ontogenetic context of the focal RMQS 2024
# assemblages. Run after 00_prepare_multitaxa_inputs_2024.R.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

INPUT_DIR <- "outputs_multitaxa_2024"
OUT_DIR <- file.path(INPUT_DIR, "uncertainty_results", "publication_figures_v3", "main_text")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MAIN_ASSEMBLAGES <- c(
  "collembola_soil_core",
  "araneae_pitfall",
  "carabidae_pitfall",
  "formicidae_pitfall",
  "isopoda_pitfall"
)

assemblage_labels <- c(
  collembola_soil_core = "Collembola\nsoil cores",
  araneae_pitfall = "Araneae\npitfall traps",
  carabidae_pitfall = "Carabidae\npitfall traps",
  formicidae_pitfall = "Formicidae\npitfall traps",
  isopoda_pitfall = "Isopoda\npitfall traps"
)

resolution_colours <- c(
  "Species" = "#0072B2",
  "Genus" = "#E69F00",
  "Family" = "#009E73",
  "Order" = "#CC79A7",
  "Other / unusable" = "#999999"
)

theme_paper <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = rel(1.15)),
      plot.subtitle = element_text(hjust = 0),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "grey20"),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}

manifest <- readr::read_csv(file.path(INPUT_DIR, "assemblage_manifest.csv"), show_col_types = FALSE)
audit <- readr::read_csv(file.path(INPUT_DIR, "taxonomic_resolution_audit_all.csv"), show_col_types = FALSE)

meta <- manifest %>%
  filter(assemblage_id %in% MAIN_ASSEMBLAGES) %>%
  transmute(
    assemblage_id,
    assemblage_label = factor(
      recode(assemblage_id, !!!assemblage_labels),
      levels = rev(unname(assemblage_labels))
    ),
    unresolved_pct = 100 * (share_genus_abundance + share_family_or_coarser_abundance),
    non_adult_pct = 100 * juvenile_share
  )

p_a <- audit %>%
  filter(
    assemblage_id %in% MAIN_ASSEMBLAGES,
    resolution %in% c("species", "genus", "family", "order", "other_or_unusable")
  ) %>%
  left_join(meta %>% select(assemblage_id, assemblage_label), by = "assemblage_id") %>%
  mutate(
    resolution = recode(
      resolution,
      species = "Species",
      genus = "Genus",
      family = "Family",
      order = "Order",
      other_or_unusable = "Other / unusable"
    ),
    resolution = factor(resolution, levels = names(resolution_colours))
  ) %>%
  ggplot(aes(x = 100 * abundance_share, y = assemblage_label, fill = resolution)) +
  geom_col(width = .74, colour = "white", linewidth = .28) +
  scale_fill_manual(values = resolution_colours, name = "Reported resolution") +
  scale_x_continuous(
    limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
    labels = \(x) paste0(x, "%"), expand = c(0, 0)
  ) +
  labs(title = "Reported taxonomic resolution", x = "Share of total abundance", y = NULL) +
  theme_paper()

p_b <- meta %>%
  select(assemblage_label, `Not resolved to species` = unresolved_pct, `Explicit non-adults` = non_adult_pct) %>%
  pivot_longer(-assemblage_label, names_to = "source", values_to = "share_pct") %>%
  mutate(source = factor(source, levels = c("Not resolved to species", "Explicit non-adults"))) %>%
  ggplot(aes(x = share_pct, y = assemblage_label)) +
  geom_segment(aes(x = 0, xend = share_pct, yend = assemblage_label), colour = "grey78", linewidth = .45) +
  geom_point(size = 2.6, colour = "#333333") +
  geom_text(aes(label = sprintf("%.0f%%", share_pct)), hjust = -.25, size = 2.7) +
  facet_wrap(~ source, ncol = 1, scales = "free_x") +
  scale_x_continuous(
    labels = \(x) paste0(x, "%"),
    expand = expansion(mult = c(0, .18))
  ) +
  labs(title = "Sources of taxonomic uncertainty", x = "Share of total abundance", y = NULL) +
  theme_paper() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p <- (p_a + p_b) +
  patchwork::plot_layout(widths = c(1.55, 1)) +
  patchwork::plot_annotation(
    title = "Taxonomic and ontogenetic context of the focal RMQS assemblages",
    subtitle = "Each group × protocol dataset is analysed independently.",
    tag_levels = "A"
  )

ggsave(file.path(OUT_DIR, "Fig1_taxonomic_context.pdf"), p,
       width = 180, height = 105, units = "mm", device = grDevices::pdf)
ggsave(file.path(OUT_DIR, "Fig1_taxonomic_context.png"), p,
       width = 180, height = 105, units = "mm", dpi = 500)
