# Taxonomic uncertainty in multi-taxa soil-biodiversity monitoring

## Overview

This repository contains the R workflow used to quantify how taxonomic uncertainty and taxonomic data-processing choices affect biodiversity inference in the 2024 French Soil Quality Monitoring Network (RMQS) dataset.

The analysis is conducted separately for each **taxonomic group × sampling protocol** assemblage. It evaluates consequences for:

* local diversity: Hill numbers q0, q1 and q2;
* regional richness (gamma diversity);
* beta diversity: Bray–Curtis, Sørensen and Jaccard dissimilarities;
* multivariate community structure and ordination stability;
* environmental turnover inference using Generalised Dissimilarity Models (GDMs);
* apparent changes in mean taxon occupancy using an alpha–gamma synthesis.

The central aim is to distinguish the effects of:

1. moderate identification errors among congeners already observed in the monitoring network;
2. errors directed towards a broader mainland-France TAXREF pool;
3. expert-validated source-to-target confusion pathways;
4. analytical choices used to handle uncertain records, including removal of unresolved records, conservative genus-level reporting, and taxonomic aggregation.

## Study design

Each assemblage is defined as a unique combination of taxonomic group and collection method. Data from different protocols are never merged into a single community matrix.

Main assemblages:

* Collembola — soil cores;
* Araneae — pitfall traps;
* Carabidae — pitfall traps;
* Formicidae — pitfall traps;
* Isopoda — pitfall traps.

Secondary assemblages are retained for sensitivity analyses:

* Isopoda — hand sorting;
* Diplopoda — pitfall traps;
* Diplopoda — hand sorting.

## Repository structure

```text
.
├── analyses/
│   └── 4.multitaxa ID precision/
│       ├── 00_prepare_multitaxa_inputs_2024.R
│       ├── 01_run_multitaxa_uncertainty_2024_v3_taxref_coherent_confusion.R
│       ├── 02_run_multitaxa_error_gradients_appendix_v3_1_fast.R
│       ├── 04_make_fig1_taxonomic_context.R
│       ├── 05_build_multitaxa_taxref_pools_v2.R
│       ├── 05b_check_taxref_pools_and_expert_maps.R
│       ├── 05c_prepare_expert_confusion_questionnaires_v3.R
│       └── 06_refresh_multitaxa_figures_v3_with_taxref.R
│
├── data/
│   └── README.md
│
├── regional_pools/
│   ├── observed_taxref_manual_overrides.csv
│   └── [generated TAXREF candidate pools]
│
├── expert_confusion_maps/
│   └── [expert-reviewed source-to-target confusion maps]
│
├── outputs_multitaxa_2024/
│   └── [generated analysis outputs]
│
├── README.md
└── .gitignore
```

Exact folder names can be adapted, but the paths defined at the beginning of the scripts must be updated accordingly.

## Data availability

This repository does **not** include:

* raw RMQS biological data;
* raw environmental data;
* full TAXREF exports;
* confidential metadata;
* large intermediate or figure-export files.

These files are excluded because of data-governance, licensing, or storage constraints.

The workflow expects the user to provide local input data in the structure required by `00_prepare_multitaxa_inputs_2024.R`.

The complete TAXREF export must be obtained separately from the official French taxonomic reference source. Its local path is defined in:

```r
TAXREF_FILE <- "TAXREFv17.txt"
```

within `05_build_multitaxa_taxref_pools_v2.R`.

The TAXREF file is intentionally ignored by Git.

## Main analytical scenarios

### Taxonomic workflow scenarios

The workflow compares several deterministic decisions for handling uncertain records:

* best-available mixed-resolution RTU matrix;
* species-only matrix after removing unresolved records;
* conservative reporting of rare species at genus level;
* genus-level aggregation;
* family-level aggregation, when meaningful;
* adult-only sensitivity scenario for assemblages with explicit non-adult records.

### Identification-error scenarios

Identification errors are simulated from species-level community matrices.

For each stochastic iteration, a source species is assigned a target species and this source-to-target association remains constant across all sites. The realised number of reassigned individuals varies among sites through binomial sampling.

The main scenarios include:

* observed-pool congeneric error at 10%;
* rare-weighted observed-pool error at 10%;
* regional TAXREF-pool congeneric error at 10%;
* expert-validated confusion error at 10%, when an expert map is available.

The supplementary gradient analysis evaluates error rates from 1% to 20%.

## Expert confusion maps

Expert confusion scenarios are optional and are only activated when a reviewed source-to-target map is available.

To prepare taxonomist-facing review sheets:

```r
source("analyses/4.multitaxa ID precision/05c_prepare_expert_confusion_questionnaires_v3.R")
```

This generates files in:

```text
expert_review_to_send/
```

Each specialist receives:

* a candidate source-to-target confusion table;
* a source-taxa overview;
* an additional-pairs table for cross-genus or missing candidates;
* instructions for reviewing plausible directional confusions.

Completed files can then be compiled into:

```text
expert_confusion_maps/<assemblage_id>__expert_confusions.csv
```

The final map must contain at least:

```text
source_taxon_unit,target_taxon_unit,weight,enabled,comment
```

Only expert-reviewed and biologically plausible directional confusion pathways should be enabled.

## Recommended workflow

Run the scripts in the following order.

### 1. Prepare assemblage-level inputs

```r
source("analyses/4.multitaxa ID precision/00_prepare_multitaxa_inputs_2024.R")
```

This creates the standardised assemblage matrices, taxonomic lookup tables and audit tables.

### 2. Build TAXREF candidate pools

```r
source("analyses/4.multitaxa ID precision/05_build_multitaxa_taxref_pools_v2.R")
```

This creates:

* mainland-France TAXREF candidate pools for each assemblage;
* TAXREF matching audits;
* templates for expert confusion maps;
* a regional-pool summary table.

### 3. Review TAXREF matching

```r
source("analyses/4.multitaxa ID precision/05b_check_taxref_pools_and_expert_maps.R")
```

Inspect especially:

```text
regional_pools/regional_pool_build_summary.csv
regional_pools/*__observed_taxref_match_audit.csv
regional_pools/observed_taxref_manual_overrides.csv
```

Manual overrides should only be added when an observed RMQS unit can be confidently matched to a TAXREF accepted concept.

### 4. Prepare and collect expert confusion maps

```r
source("analyses/4.multitaxa ID precision/05c_prepare_expert_confusion_questionnaires_v3.R")
```

After expert review, compile the returned maps according to the instructions in that script.

### 5. Run the main multi-taxa analysis

```r
source("analyses/4.multitaxa ID precision/01_run_multitaxa_uncertainty_2024_v3_taxref_coherent_confusion.R")
```

This generates alpha-diversity, beta-diversity, ordination and GDM outputs for all applicable scenarios.

### 6. Run the 1–20% error-gradient appendix

```r
source("analyses/4.multitaxa ID precision/02_run_multitaxa_error_gradients_appendix_v3_1_fast.R")
```

This script uses checkpoints and can resume after interruption.

Checkpoints are stored in:

```text
outputs_multitaxa_2024/uncertainty_results/appendix_error_gradients_v3/checkpoints/
```

### 7. Generate figures

```r
source("analyses/4.multitaxa ID precision/04_make_fig1_taxonomic_context.R")
source("analyses/4.multitaxa ID precision/06_refresh_multitaxa_figures_v3_with_taxref.R")
```

## Key outputs

The main analysis produces:

* taxonomic-resolution audits by assemblage;
* scenario-level alpha and gamma diversity metrics;
* Bray–Curtis, Sørensen and Jaccard beta-diversity metrics;
* ordination stability metrics;
* GDM explained deviance and fitted-turnover stability;
* alpha–gamma–occupancy summaries;
* scenario-specific source-to-target confusion maps;
* publication-ready figures and figure data tables.

The error-gradient appendix produces:

* responses of richness and beta diversity from 1% to 20% error;
* stability curves for alpha and beta metrics;
* eligibility diagnostics showing the share of individuals and species that could actually be reassigned under each mechanism.

## Reproducibility notes

* The full TAXREF export is not included in the repository.
* RMQS raw data are not included in the repository.
* Large generated outputs should not be committed unless they are required for a specific release or archive.
* Scripts should be run from the root R project directory.
* Paths and package versions should be recorded before final manuscript submission.

## Software

Analyses were developed in R using, among others:

* `tidyverse`
* `vegan`
* `gdm`
* `janitor`
* `Matrix`
* `patchwork`
* `ggplot2`

A `renv` environment is recommended to record package versions:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

## Citation

A manuscript describing this workflow is in preparation.

Please cite the RMQS programme, TAXREF, and the relevant taxonomic experts when using or adapting this workflow.
