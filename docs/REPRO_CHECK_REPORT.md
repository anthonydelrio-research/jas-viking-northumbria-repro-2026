# Reproducibility Check Report

Derived-table verification against manuscript headline metrics.

## Headline Metrics

| Metric | Observed | Expected | Match |
|---|---:|---:|:---:|
| primary_hit_rate | 0.24 | 0.24 | PASS |
| primary_auc | 0.617026 | 0.617026 | PASS |
| primary_auc_opportunity | 0.531501 | 0.531501 | PASS |
| primary_gain_csr | 1.6 | 1.6 | PASS |
| primary_gain_opportunity | 2 | 2 | PASS |
| sensitivity_hit_rate | 0.222222222222 | 0.222222222222 | PASS |
| sensitivity_auc | 0.594751851852 | 0.594751851852 | PASS |

## Region-Held-Out Transportability

| Scenario | Hit Rate | AUC |
|---|---:|---:|
| train_yorkshire_test_northwest | 0.000000 | 0.168383 |
| train_northwest_test_yorkshire | 0.000000 | 0.374552 |
| train_yorkshire_test_yorkshire | 0.416667 | 0.713263 |
| train_northwest_test_northwest | 0.769231 | 0.898079 |

## Verdict

All headline metric checks PASS.
