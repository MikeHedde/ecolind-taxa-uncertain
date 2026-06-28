# RMQS taxonomic-sensitivity pipeline

A configuration-driven R workflow for quantifying how taxonomic uncertainty and
taxonomic data-processing decisions affect biodiversity indicators, community
structure, and environmental-turnover inference.

## Design principles

- An **assemblage** is generated automatically as `taxon × protocol`.
- Taxa are registered once in `config/taxon_registry.csv`; protocols are
  discovered from the raw dataset.
- One manually curated **expert workbook per taxon** is sufficient. It is
  projected automatically onto every matching assemblage.
- Expert inputs remain distinct internally as species confusions, reporting
  rules, and RTU confusions, but the manuscript can show one readable
  `expert_integrated_10pct` scenario while retaining components in supplements.
- Scenarios to run and scenarios to display are controlled in
  `config/scenario_catalog.csv`, not hard-coded in figures.

## First-time setup

1. Open the project from its `.Rproj` file or set the working directory to the
   project root.
2. Edit `config/analysis_settings.R`.
3. Edit `config/taxon_registry.csv`.
4. Add manually curated expert workbooks to `expert_input/`.
5. Run scripts in order:

```r
source("scripts/01_prepare_assemblages.R")
source("scripts/02_build_taxref_pools.R")
source("scripts/03_project_expert_rules.R")
source("scripts/04_run_sensitivity_analysis.R")
source("scripts/05_run_error_gradient_appendix.R")
source("scripts/06_make_figures.R")
```

## Adding a new taxon

Add **one row** to `config/taxon_registry.csv`. Do not add protocol-specific
rows. The preparation script will detect all protocols containing records that
match the selector.

For example, a Staphylinidae row can use:

```text
staphylinidae,Staphylinidae,family,Staphylinidae,famille,Staphylinidae,...
```

If Earthworms are represented by a different selector in a future raw dataset,
change only the selector columns in that same registry row.

## Expert workbook: one file per taxon

Store one workbook per taxon under `expert_input/`, for example:

```text
expert_input/isopoda_expert_confusions.xlsx
```

It contains a `rules` sheet with manually resolved TAXREF concepts and/or RTU
units. See `config/expert_rules_template.csv`.

The pipeline does **not** infer taxonomic matches from names. The fields
`source_cd_ref`, `target_cd_ref`, `source_unit`, and `target_unit` are expected
to have been curated manually as expertise accumulates.

Supported rule types:

- `species_confusion`:
  `source_cd_ref → target_cd_ref`
- `reporting_rule`:
  an observed source species or RTU is deliberately reported at a coarser
  target unit, such as `genus:Oritoniscus`.
- `rtu_confusion`:
  a genus/family RTU can be confused with another genus/family RTU.

The projection script filters rules to the taxa actually present in each
assemblage and produces auditable per-assemblage internal maps.

## Scenario presentation

The scenario catalogue separates computation from presentation:

- `run = TRUE`: scenario is calculated.
- `show_main = TRUE`: scenario is available to main-text plots.
- `show_supplement = TRUE`: scenario is available to supplementary plots.

The default main-text set is deliberately small:

1. observed-pool error at 10%;
2. regional TAXREF-pool error at 10%;
3. expert-integrated uncertainty at 10%;
4. rare taxa reported at genus;
5. unresolved records dropped;
6. genus-level aggregation;
7. family-level aggregation where applicable.

Detailed expert components, rare-weighted errors, life-stage sensitivity,
Jaccard sensitivity, and the 1–20% error gradient remain supplementary.

## Outputs

- `outputs/assemblages/`: one set of inputs per discovered taxon × protocol.
- `regional_pools/`: TAXREF candidate pools and matching audits.
- `expert_internal/`: projected expert maps and audits by assemblage.
- `outputs/sensitivity/`: scenario metrics, GDM outputs, and figures.
