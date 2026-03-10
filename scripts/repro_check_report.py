#!/usr/bin/env python3
"""
Generate a compact reproducibility report from included derived tables.

This script is designed for JAS reproducibility screening when source datasets
have third-party licensing constraints. It verifies that the packaged derived
outputs reproduce manuscript headline metrics.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "data" / "derived" / "primary" / "tables"
SENS = ROOT / "data" / "derived" / "sensitivity" / "tables"
OUT = ROOT / "docs" / "REPRO_CHECK_REPORT.md"


def read_first_row_tsv(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            return {k: (v or "") for k, v in row.items()}
    raise RuntimeError(f"No rows in {path}")


def read_rows_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        return [{k: (v or "") for k, v in r.items()} for r in rdr]


def fnum(v: str) -> float:
    return float(v)


def close(a: float, b: float, tol: float = 1e-9) -> bool:
    return math.isclose(a, b, rel_tol=0.0, abs_tol=tol)


def main() -> None:
    primary_val = read_first_row_tsv(PRIMARY / "validation_metrics.tsv")
    primary_auc = read_first_row_tsv(PRIMARY / "auc_metrics.tsv")
    primary_opp = read_first_row_tsv(PRIMARY / "auc_opportunity_null_metrics.tsv")
    primary_gain = read_first_row_tsv(PRIMARY / "baseline_gain_metrics.tsv")
    sens_val = read_first_row_tsv(SENS / "validation_metrics.tsv")
    sens_auc = read_first_row_tsv(SENS / "auc_metrics.tsv")
    region = read_rows_tsv(PRIMARY / "region_holdout_transportability.tsv")

    expected = {
        "primary_hit_rate": 0.24,
        "primary_auc": 0.617026,
        "primary_auc_opportunity": 0.531501,
        "primary_gain_csr": 1.6,
        "primary_gain_opportunity": 2.0,
        "sensitivity_hit_rate": 0.2222222222222222,
        "sensitivity_auc": 0.5947518518518519,
    }

    observed = {
        "primary_hit_rate": fnum(primary_val["hit_rate"]),
        "primary_auc": fnum(primary_auc["auc_mean"]),
        "primary_auc_opportunity": fnum(primary_opp["auc_opportunity_null"]),
        "primary_gain_csr": fnum(primary_gain["gain_vs_csr"]),
        "primary_gain_opportunity": fnum(primary_gain["gain_vs_opportunity"]),
        "sensitivity_hit_rate": fnum(sens_val["hit_rate"]),
        "sensitivity_auc": fnum(sens_auc["auc_mean"]),
    }

    checks = []
    for key, exp in expected.items():
        obs = observed[key]
        ok = close(obs, exp)
        checks.append((key, obs, exp, ok))

    with OUT.open("w", encoding="utf-8") as f:
        f.write("# Reproducibility Check Report\n\n")
        f.write("Derived-table verification against manuscript headline metrics.\n\n")
        f.write("## Headline Metrics\n\n")
        f.write("| Metric | Observed | Expected | Match |\n")
        f.write("|---|---:|---:|:---:|\n")
        for key, obs, exp, ok in checks:
            f.write(f"| {key} | {obs:.12g} | {exp:.12g} | {'PASS' if ok else 'FAIL'} |\n")

        f.write("\n## Region-Held-Out Transportability\n\n")
        f.write("| Scenario | Hit Rate | AUC |\n")
        f.write("|---|---:|---:|\n")
        for r in region:
            f.write(
                f"| {r['scenario']} | {float(r['hit_rate']):.6f} | {float(r['auc']):.6f} |\n"
            )

        n_fail = sum(0 if ok else 1 for _, _, _, ok in checks)
        f.write("\n## Verdict\n\n")
        if n_fail == 0:
            f.write("All headline metric checks PASS.\n")
        else:
            f.write(f"{n_fail} headline metric checks FAIL.\n")

    print(str(OUT))


if __name__ == "__main__":
    main()
