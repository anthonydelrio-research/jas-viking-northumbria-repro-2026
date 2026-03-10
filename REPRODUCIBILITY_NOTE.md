# Reproducibility Note

## Scope

This note documents exact reruns for the manuscript's AD `750-954` package-level analyses.

- Primary run output: `outputs_750_primary_pasall_v1/`
- Sensitivity run output: `outputs_750_sensitivity_pasall_v2/`
- Canonical submission bundle: `outputs_submission_2026-03-03_pasall/`
- Main script: `scripts/run_jas_analysis.py`

The `outputs_*` paths above are the original author-side run directories used to produce the manuscript package. This repository redistributes the corresponding derived result tables under `data/derived/`.

## Computational Setup

- Coordinate reference system: EPSG:27700
- Randomized procedures in the script are seeded internally for deterministic reruns
- Full iteration settings were used (no reduced bootstrap/permutation shortcuts)

## Exact Commands Used

From repository root:

```bash
JAS_OUTPUT_DIR=$PWD/outputs_750_primary_pasall_v1 \
JAS_ROOT=$PWD \
JAS_INBOX=$PWD/incoming_arcgis/01_inbox \
JAS_TARGET_YEAR_FROM=750 \
JAS_TARGET_YEAR_TO=954 \
JAS_PAS_TEMPORAL_MODE=overlap \
JAS_HOLDOUT_GPKG=$PWD/incoming_arcgis/01_inbox/validation_holdout_sites_vetted_primary_2026-03-03.gpkg \
JAS_HOLDOUT_LAYER=validation_holdout_sites_vetted_primary \
JAS_HOLDOUT_TEMPORAL_MODE=overlap \
JAS_USE_OBSERVED_OPPORTUNITY=1 \
JAS_OBSERVED_OPPORTUNITY_MODE=pas_all \
JAS_PAS_ALL_PATH=$PWD/incoming_arcgis/01_inbox/PAS_ALL/PAS_ALL_PERIODS_Northumbria.shp \
python3 scripts/run_jas_analysis.py
```

```bash
JAS_OUTPUT_DIR=$PWD/outputs_750_sensitivity_pasall_v2 \
JAS_ROOT=$PWD \
JAS_INBOX=$PWD/incoming_arcgis/01_inbox \
JAS_TARGET_YEAR_FROM=750 \
JAS_TARGET_YEAR_TO=954 \
JAS_PAS_TEMPORAL_MODE=overlap \
JAS_HOLDOUT_GPKG=$PWD/incoming_arcgis/01_inbox/validation_holdout_sites_vetted_2026-03-03.gpkg \
JAS_HOLDOUT_LAYER=validation_holdout_sites_vetted \
JAS_INCLUDE_SENSITIVITY_HOLDOUTS=1 \
JAS_HOLDOUT_TEMPORAL_MODE=overlap \
JAS_USE_OBSERVED_OPPORTUNITY=1 \
JAS_OBSERVED_OPPORTUNITY_MODE=pas_all \
JAS_PAS_ALL_PATH=$PWD/incoming_arcgis/01_inbox/PAS_ALL/PAS_ALL_PERIODS_Northumbria.shp \
python3 scripts/run_jas_analysis.py
```

## Expected Key Outputs

### Primary (`outputs_750_primary_pasall_v1`)
- `tables/validation_metrics.tsv`: `n_validation=25`, `n_hits_top15=6`, `hit_rate=0.24`
- `tables/auc_metrics.tsv`: `auc_mean=0.617026`, `auc_temporal_weighted=0.6336099221`
- `tables/auc_opportunity_null_metrics.tsv`: `auc_opportunity_null=0.531501`
- `tables/baseline_gain_metrics.tsv`: `gain_vs_csr=1.6`, `gain_vs_opportunity=2.0`
- `tables/temporal_holdout_scenarios.tsv`:
  - strict-inside: `n=21`, Hit@15% `0.2381`, AUC `0.6393`
  - overlap: `n=25`, Hit@15% `0.24`, AUC `0.6170`

### Sensitivity (`outputs_750_sensitivity_pasall_v2`)
- `tables/validation_metrics.tsv`: `n_validation=27`, `n_hits_top15=6`, `hit_rate=0.222222...`
- `tables/auc_metrics.tsv`: `auc_mean=0.5947518519`, `auc_temporal_weighted=0.6135320131`
- `tables/auc_opportunity_null_metrics.tsv`: `auc_opportunity_null=0.5056453704`

### Comparison targets
- Primary holdout package: `outputs_750_primary_pasall_v1`
- Sensitivity holdout package: `outputs_750_sensitivity_pasall_v2`

## Current Reproducibility Constraints

1. `SICI_validation_sample.tsv` distributed local rater columns are blank, so in-package kappa/AC1 are thesis-locked reference values; current-package SICI grading is author-coded with supervisory expert review.
2. Opportunity baseline uses PAS all-finds/all-periods density as an observed reporting-intensity proxy, not direct detector-count telemetry.
3. External archival DOI minting is pending and can be added at revision/acceptance stage.

## Author

- Anthony Brian Del Rio
- ORCID: `0009-0002-3799-8225`
- Contact: `anthony.delrio@arch.ox.ac.uk`; `contact@anthonydelrio.com`
