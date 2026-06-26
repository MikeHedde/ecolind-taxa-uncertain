# RMQS 2024 — reprise du pipeline d’incertitude taxonomique multi-taxa

Ce dossier reprend le chantier exactement dans la direction arrêtée avant l’interruption :

- campagne **RMQS 2024 uniquement** ;
- assemblages analysés séparément comme **groupe × protocole** ;
- équilibre entre **alpha-diversité**, **bêta-diversité** et **inférence environnementale** ;
- GDM avec Bray–Curtis et Sørensen ; Jaccard seulement en sensibilité ;
- stades non adultes uniquement là où le codage le justifie ;
- agrégation familiale uniquement lorsque l’assemblage couvre réellement plusieurs familles ;
- confusion d’identification **cohérente source → cible** dans tous les sites d’une itération ;
- scénario principal à **10 %** pour chaque mécanisme ;
- gradients **1–20 %** séparés en annexe ;
- pool régional national TAXREF et cartes d’experts curées ;
- couche Blowes reconstruite seulement après les analyses principales.

## Scripts actifs, dans l’ordre

### 0. Préparation des matrices par assemblage × protocole

```r
source("00_prepare_multitaxa_inputs_2024.R")
```

Produit les matrices station × RTU / espèce / genre / famille, les audits de résolution et les jeux `adult_only` pour Araignées et Diplopodes.

### 1. Construction des pools TAXREF

Dans `05_build_multitaxa_taxref_pools_v2.R`, renseigner :

```r
TAXREF_FILE <- "CHEMIN/VERS/TAXREF.csv"
```

Puis :

```r
source("05_build_multitaxa_taxref_pools_v2.R")
```

Inspecter impérativement :

- `regional_pools/regional_pool_build_summary.csv`
- `regional_pools/*__observed_taxref_match_audit.csv`
- `regional_pools/observed_taxref_manual_overrides.csv`

Un taxon RMQS non apparié doit être contrôlé avant le run principal. Les overrides servent seulement à documenter et corriger des synonymies/orthographes non résolues automatiquement.

### 1b. Contrôle pré-vol des pools et cartes

```r
source("05b_check_taxref_pools_and_expert_maps.R")
```

Ne pas interpréter un scénario régional ou expert avant que le tableau
`regional_pools/regional_pool_preflight_check.csv` soit satisfaisant. Les statuts
`REVIEW_MATCH_RATE`, `INVALID_EXPERT_SOURCE` et `INVALID_EXPERT_TARGET` demandent
une correction explicite.

### 2. Cartes de confusion expertes

Le script TAXREF produit des templates :

```text
expert_confusion_maps/<assemblage_id>__expert_confusions_TEMPLATE.csv
```

Les templates de Collemboles sont restreints aux genres déjà signalés comme difficiles. Ils sont des candidats de revue, pas une carte prête à utiliser.

Pour activer une carte, copier le template sous :

```text
expert_confusion_maps/<assemblage_id>__expert_confusions.csv
```

Puis mettre `enabled = TRUE` uniquement pour des paires source → cible défendables. Les cartes expertes peuvent inclure des confusions inter-genres si elles sont explicites dans le fichier final : elles ne sont pas forcées à être intra-génériques.

### 3. Analyse principale

```r
source("01_run_multitaxa_uncertainty_2024_v3_taxref_coherent_confusion.R")
```

Le main text utilise **un seul niveau nominal à 10 %** :

- erreur congénérique dans le pool observé ;
- erreur pondérée vers les rares, dans le pool observé ;
- erreur congénérique vers le pool TAXREF régional ;
- erreur experte, seulement si une carte curée existe ;
- reporting prudent, suppression des non-résolus, passage au genre et à la famille selon applicabilité.

Chaque itération tire une carte source → cible ; elle est ensuite appliquée à tous les sites. Le nombre d’individus réassignés varie entre sites, mais pas la direction de la confusion.

Par défaut, le script exporte une carte auditable pour l’itération 1 de chaque scénario × assemblage dans `uncertainty_results/per_assemblage/`. Pour tous les maps, activer :

```r
WRITE_ALL_COHERENT_CONFUSION_MAPS <- TRUE
```

### 4. Annexe : gradient d’erreur 1–20 %

```r
source("02_run_multitaxa_error_gradients_appendix_v2.R")
```

Cette analyse ne recalcule pas les GDM. Elle produit les courbes alpha/bêta le long du gradient 1–20 % et conserve **la même carte source → cible à toutes les intensités** au sein d’une itération. Ainsi, les courbes isolent l’effet de l’intensité d’erreur.

Lire aussi :

```text
error_gradient_eligibility_audit.csv
```

Car un taux de 10 % est demandé pour les sources ayant une cible de confusion disponible ; la fraction réellement réassignable dépend donc de la part des individus appartenant à des espèces pour lesquelles un congénère candidat existe.

### 5. Figures principales après TAXREF

```r
source("04_make_fig1_taxonomic_context.R")
source("06_refresh_multitaxa_figures_v3_with_taxref.R")
```

Ces scripts régénèrent :

- Fig. 1 : contexte de résolution taxonomique et de stades ;
- Fig. 2 : robustesse équilibrée alpha / bêta / GDM ;
- Fig. 3 : effets directionnels ;
- Fig. 4 : inférence GDM ;
- Fig. 5 : espace alpha–gamma–occupation de Blowes.

## Important pour l’interprétation

Le scénario TAXREF est une **borne taxonomique nationale**, pas une reconstitution du pool local effectivement accessible. Il sert à quantifier ce qui se passe lorsqu’une erreur intra-générique ouvre la porte à des espèces françaises non observées dans le réseau ; ce scénario peut donc produire de la différenciation apparente et une inflation de gamma.

Les scénarios experts ont une ambition différente : reproduire uniquement des confusions plausibles et documentées. Ils peuvent être plus réalistes mais sont volontairement restreints aux groupes pour lesquels une carte est disponible.

