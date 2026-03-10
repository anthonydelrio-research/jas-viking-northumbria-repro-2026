# TSV Data Dictionary

This document summarizes the key TSV outputs redistributed under `data/derived/`.

## Directory layout

- `data/derived/primary/tables/`: primary holdout run outputs.
- `data/derived/sensitivity/tables/`: sensitivity holdout run outputs.
- `data/derived/holdout/HOLDOUT_BREADTH_STRESS_TEST.tsv`: holdout breadth stress-test summary.

Most table schemas are the same between `primary` and `sensitivity`.

## Core tables used for manuscript headline metrics

### `validation_metrics.tsv`

One-row summary of top-area capture performance and uncertainty.

Key fields:
- `n_validation`: number of validation sites.
- `n_hits_top15`: number of validation sites in top 15% ranked area.
- `hit_rate`: `n_hits_top15 / n_validation`.
- `ci_lower`, `ci_upper`: Wilson interval bounds for `hit_rate`.
- `binomial_p_one_sided_vs_csr`: one-sided exact binomial p-value vs CSR baseline.
- `posterior_prob_hit_rate_gt_csr`: posterior probability that hit rate exceeds CSR baseline.

### `auc_metrics.tsv`

One-row summary of ranking discrimination.

Key fields:
- `auc_mean`: AUC from unweighted holdout evaluation.
- `auc_temporal_weighted`: AUC from temporal-overlap-weighted evaluation.
- `ci_lower`, `ci_upper`: bootstrap interval bounds for `auc_mean`.
- `negative_samples`: number of sampled negatives used for AUC estimation.

### `auc_opportunity_null_metrics.tsv`

One-row comparison against opportunity-weighted null.

Key fields:
- `auc_opportunity_null`: AUC against opportunity-weighted negatives.
- `auc_opportunity_null_temporal_weighted`: temporal-overlap-weighted variant.
- `sampling_mode`: identifier of opportunity sampling mode.

### `baseline_gain_metrics.tsv`

One-row gain statistics at the top-area threshold.

Key fields:
- `gain_vs_csr`: hit-rate gain vs complete spatial randomness expectation.
- `gain_vs_opportunity`: hit-rate gain vs opportunity baseline.
- `opportunity_baseline_hit_rate`: baseline hit rate under opportunity sampling.

### `temporal_holdout_scenarios.tsv`

Multi-scenario comparison across strict/overlap/weighted temporal handling.

Key fields:
- `scenario`: scenario label.
- `n_validation`, `n_hits_top15`, `hit_rate`, `auc`: core metrics.
- `hit_rate_temporal_weighted`, `auc_temporal_weighted`: weighted variants where applicable.
- `n_burial`, `n_settlement`, `n_hoard`: holdout composition counts.

### `region_holdout_transportability.tsv`

Region-held-out transfer diagnostics.

Key fields:
- `scenario`: train/test region direction.
- `train_region`, `test_region`.
- `n_train_pas`, `n_test_holdouts`.
- `hit_rate`, `auc`.
- `gain_vs_csr`, `gain_vs_opportunity`.

### `decile_lift_table.tsv`

Capture-vs-area table used for ranking utility diagnostics.

Key fields:
- `top_area_pct`: area threshold.
- `capture_rate`: observed capture rate.
- `expected_random_rate`: CSR expectation.
- `opportunity_capture_rate`: opportunity-baseline expectation.
- `lift_vs_random`, `lift_vs_opportunity`.

### `decile_bin_reliability.tsv`

Score-bin reliability summary.

Key fields:
- `score_decile`, `score_lower`, `score_upper`.
- `n_in_bin`.
- `observed_rate`, `expected_rate`.
- `lift_vs_uniform`.

## Additional TSV outputs

Additional files in each `tables/` directory provide:
- cross-validation summaries (`crossvalidation_summary.tsv`, `fold_metrics.tsv`, `fold_assignments.tsv`);
- robustness checks (`bandwidth_sensitivity_results.tsv`, `kernel_comparison_results.tsv`, `coordinate_perturbation_results.tsv`, `grade_perturbation_results.tsv`, `jackknife_results.tsv`);
- environmental and proximity diagnostics (`environmental_correlation_matrix.tsv`, `soil_chisquare_stats.tsv`, `elevation_chisquare_stats.tsv`, `roman_roads_proximity.tsv`, `waterways_proximity.tsv`, `infrastructure_proximity_results.tsv`);
- multiple-testing and permutation outputs (`fdr_adjusted_pvalues.tsv`, `permutation_test_results.tsv`, `permutation_tests_extended.tsv`);
- traceability summaries (`dataset_summary_computed.tsv`, `exclusion_counts.tsv`, `holdout_influence_loo.tsv`, `holdout_uncertainty_diagnostics.tsv`).

## Holdout breadth stress-test table

### `HOLDOUT_BREADTH_STRESS_TEST.tsv`

Row-per-holdout-set comparison across:
- `primary`
- `primary_plus_sensitivity`
- `primary_plus_sensitivity_plus_context`

Key fields:
- `n_validation`, `n_hits_top15`, `hit_rate`, `auc_uniform`, `gain_vs_csr`.
- holdout composition counts: `n_burial`, `n_settlement`, `n_hoard`.
