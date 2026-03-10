#!/usr/bin/env python3
"""
Targeted stress tests for unresolved submission blockers:
  - Opportunity baseline semantics (alternative observed PAS-all proxies)
  - Holdout sample breadth (primary vs expanded vetted roles)

This script reuses the same analytical primitives as run_jas_analysis.py
but evaluates only the validation/opportunity metrics required for the
blocker decision.
"""

from __future__ import annotations

import csv
import importlib.util
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(os.environ.get("JAS_ROOT", str(Path(__file__).resolve().parents[1]))).resolve()
SCRIPT_MAIN = ROOT / "scripts" / "run_jas_analysis.py"
INBOX = ROOT / "incoming_arcgis" / "01_inbox"
VALIDATED = ROOT / "incoming_arcgis" / "03_validated" / "2026-03-05_opportunity_baselines"
OUT = ROOT / "outputs_submission_2026-03-05_blocker_mitigation"
OUT_TABLES = OUT / "tables"

TARGET_YEAR_FROM = 750
TARGET_YEAR_TO = 954
TOP_AREA_PCT = 15.0
N_NEG = int(os.environ.get("JAS_BLOCKER_N_NEG", "10000"))
N_BOOT = int(os.environ.get("JAS_BLOCKER_N_BOOT", "2000"))
SEED = 20260305


def load_main_module():
    spec = importlib.util.spec_from_file_location("run_jas_analysis_mod", SCRIPT_MAIN)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module from {SCRIPT_MAIN}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)  # type: ignore[attr-defined]
    return mod


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    OUT_TABLES.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path, rows: List[Dict[str, object]], header: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def load_holdouts(mod, include_sensitivity: bool, include_context: bool) -> List[Dict[str, object]]:
    gpkg = INBOX / "validation_holdout_sites_vetted_2026-03-03.gpkg"
    rows = mod.read_holdouts_from_gpkg(
        gpkg,
        layer_name="validation_holdout_sites_vetted",
        include_sensitivity=include_sensitivity,
        include_context=include_context,
    )
    out: List[Dict[str, object]] = []
    for r in rows:
        a0 = mod.parse_year(r.get("date_from_ad"), TARGET_YEAR_FROM)
        a1 = mod.parse_year(r.get("date_to_ad"), TARGET_YEAR_TO)
        if not mod.temporal_pass(a0, a1, TARGET_YEAR_FROM, TARGET_YEAR_TO, "overlap"):
            continue
        rr = dict(r)
        rr["date_from_ad"] = a0
        rr["date_to_ad"] = a1
        rr["temporal_overlap_fraction"] = mod.temporal_overlap_fraction(a0, a1, TARGET_YEAR_FROM, TARGET_YEAR_TO)
        out.append(rr)
    return out


def holdout_type_counts(holdouts: List[Dict[str, object]]) -> Dict[str, int]:
    c: Dict[str, int] = {}
    for r in holdouts:
        k = str(r.get("site_type", "unknown")).strip().lower() or "unknown"
        c[k] = c.get(k, 0) + 1
    return c


def evaluate_scenario(
    mod,
    grid,
    surf_main: np.ndarray,
    surf_full: np.ndarray,
    threshold: float,
    holdouts: List[Dict[str, object]],
    opp_points_path: Path,
    opp_layer: str | None,
    scenario_label: str,
    rng_seed: int,
) -> Dict[str, object]:
    if not holdouts:
        raise RuntimeError(f"No holdouts available for scenario {scenario_label}")
    hold_x = np.array([float(r["x"]) for r in holdouts], dtype=np.float64)
    hold_y = np.array([float(r["y"]) for r in holdouts], dtype=np.float64)
    hold_w = np.array([float(r.get("temporal_overlap_fraction", 1.0)) for r in holdouts], dtype=np.float64)

    hold_scores = mod.sample_surface_at_points(surf_full, grid.gt, hold_x, hold_y)
    valid = np.isfinite(hold_scores)
    hits = (hold_scores >= threshold) & valid
    n_val = int(np.sum(valid))
    n_hits = int(np.sum(hits))
    hit_rate = (n_hits / n_val) if n_val > 0 else float("nan")
    ci_l, ci_u = mod.wilson_ci(n_hits, n_val)

    rng = np.random.default_rng(rng_seed)
    mask_n = surf_main.size
    neg_idx = rng.choice(mask_n, size=min(N_NEG, mask_n), replace=False)
    neg_scores = surf_main[neg_idx]
    auc = mod.auc_from_scores(hold_scores, neg_scores)
    if N_BOOT > 0:
        auc_ci_l, auc_ci_u = mod.bootstrap_auc_ci(hold_scores, neg_scores, N_BOOT, rng)
    else:
        auc_ci_l, auc_ci_u = float("nan"), float("nan")

    opp_x, opp_y = mod.read_vector_points(opp_points_path, opp_layer)
    opp_weight, opp_diag = mod.opportunity_surface_observed_pas_all(
        grid, opp_x, opp_y, mod.KDE_BANDWIDTH_M
    )
    opp_thr = float(np.percentile(opp_weight, 100.0 - TOP_AREA_PCT))
    opp_full = mod.array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, opp_weight, fill=np.nan)
    hold_opp = mod.sample_surface_at_points(opp_full, grid.gt, hold_x, hold_y)
    opp_hit_rate = float(np.mean(hold_opp >= opp_thr))
    gain_vs_csr = hit_rate / (TOP_AREA_PCT / 100.0)
    gain_vs_opp = (hit_rate / opp_hit_rate) if opp_hit_rate > 0 else float("nan")

    opp_prob = opp_weight.copy()
    opp_prob[~np.isfinite(opp_prob)] = 0.0
    ps = float(np.sum(opp_prob))
    if ps <= 0:
        neg_idx_opp = rng.choice(mask_n, size=min(N_NEG, mask_n), replace=False)
        opp_mode = "fallback_uniform_due_zero_opportunity_mass"
    else:
        opp_prob = opp_prob / ps
        neg_idx_opp = rng.choice(mask_n, size=min(N_NEG, mask_n), replace=False, p=opp_prob)
        opp_mode = "observed_detector_effort_pas_all_kde"
    neg_scores_opp = surf_main[neg_idx_opp]
    auc_opp = mod.auc_from_scores(hold_scores, neg_scores_opp)
    if N_BOOT > 0:
        auc_opp_ci_l, auc_opp_ci_u = mod.bootstrap_auc_ci(hold_scores, neg_scores_opp, N_BOOT, rng)
    else:
        auc_opp_ci_l, auc_opp_ci_u = float("nan"), float("nan")

    hw = np.where(np.isfinite(hold_w), hold_w, 0.0)
    den = float(np.sum(hw))
    num = float(np.sum(hw * hits.astype(np.float64)))
    hit_rate_temporal_weighted = (num / den) if den > 0 else float("nan")
    auc_temporal_weighted = mod.weighted_auc_from_scores(hold_scores, neg_scores, hw)
    auc_opp_temporal_weighted = mod.weighted_auc_from_scores(hold_scores, neg_scores_opp, hw)

    counts = holdout_type_counts(holdouts)
    return {
        "scenario": scenario_label,
        "n_validation": n_val,
        "n_hits_top15": n_hits,
        "hit_rate": hit_rate,
        "hit_rate_temporal_weighted": hit_rate_temporal_weighted,
        "hit_ci_lower": ci_l,
        "hit_ci_upper": ci_u,
        "auc_uniform": auc,
        "auc_uniform_temporal_weighted": auc_temporal_weighted,
        "auc_uniform_ci_lower": auc_ci_l,
        "auc_uniform_ci_upper": auc_ci_u,
        "auc_opportunity_null": auc_opp,
        "auc_opportunity_null_temporal_weighted": auc_opp_temporal_weighted,
        "auc_opp_ci_lower": auc_opp_ci_l,
        "auc_opp_ci_upper": auc_opp_ci_u,
        "gain_vs_csr": gain_vs_csr,
        "gain_vs_opportunity": gain_vs_opp,
        "opportunity_baseline_hit_rate": opp_hit_rate,
        "opportunity_mode": opp_mode,
        "opportunity_points": int(opp_diag.get("pas_all_points", np.nan))
        if np.isfinite(opp_diag.get("pas_all_points", np.nan))
        else int(opp_x.size),
        "n_burial": counts.get("burial", 0),
        "n_settlement": counts.get("settlement", 0),
        "n_hoard": counts.get("hoard", 0),
    }


def main() -> None:
    ensure_dirs()
    mod = load_main_module()

    boundary_path = INBOX / "northumbria_boundary.gpkg"
    grid = mod.build_analysis_grid(boundary_path, "northumbria_boundary", mod.CELL_SIZE_M)

    score_raster = ROOT / "outputs_submission_2026-03-03_pasall" / "primary" / "rasters" / "kde_main_1km.tif"
    surf_main = mod.sample_raster_to_mask_cells(score_raster, grid)
    surf_full = mod.array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, surf_main, fill=np.nan)
    threshold = float(np.percentile(surf_main[np.isfinite(surf_main)], 100.0 - TOP_AREA_PCT))

    hold_primary = load_holdouts(mod, include_sensitivity=False, include_context=False)
    hold_expanded = load_holdouts(mod, include_sensitivity=True, include_context=True)

    baselines: List[Tuple[str, Path, str | None]] = [
        (
            "opp_pas_all_full",
            INBOX / "PAS_ALL" / "PAS_ALL_PERIODS_Northumbria.shp",
            None,
        ),
        (
            "opp_pas_all_non_target_strict",
            VALIDATED / "PAS_ALL_non_target_strict.gpkg",
            "pas_all_non_target_strict",
        ),
        (
            "opp_pas_all_non_early_medieval",
            VALIDATED / "PAS_ALL_non_early_medieval.gpkg",
            "pas_all_non_early_medieval",
        ),
    ]

    rows: List[Dict[str, object]] = []
    seed_i = 0
    for b_name, b_path, b_layer in baselines:
        rows.append(
            evaluate_scenario(
                mod,
                grid,
                surf_main,
                surf_full,
                threshold,
                hold_primary,
                b_path,
                b_layer,
                f"{b_name}__holdouts_primary",
                SEED + seed_i,
            )
        )
        seed_i += 1
        rows.append(
            evaluate_scenario(
                mod,
                grid,
                surf_main,
                surf_full,
                threshold,
                hold_expanded,
                b_path,
                b_layer,
                f"{b_name}__holdouts_expanded",
                SEED + seed_i,
            )
        )
        seed_i += 1

    header = [
        "scenario",
        "n_validation",
        "n_hits_top15",
        "hit_rate",
        "hit_rate_temporal_weighted",
        "hit_ci_lower",
        "hit_ci_upper",
        "auc_uniform",
        "auc_uniform_temporal_weighted",
        "auc_uniform_ci_lower",
        "auc_uniform_ci_upper",
        "auc_opportunity_null",
        "auc_opportunity_null_temporal_weighted",
        "auc_opp_ci_lower",
        "auc_opp_ci_upper",
        "gain_vs_csr",
        "gain_vs_opportunity",
        "opportunity_baseline_hit_rate",
        "opportunity_mode",
        "opportunity_points",
        "n_burial",
        "n_settlement",
        "n_hoard",
    ]
    out_tsv = OUT_TABLES / "BLOCKER_MITIGATION_SCENARIOS.tsv"
    write_tsv(out_tsv, rows, header)

    summary_md = OUT / "BLOCKER_MITIGATION_SUMMARY.md"
    with summary_md.open("w", encoding="utf-8") as f:
        f.write("# Blocker Mitigation Stress Tests (2026-03-05)\n\n")
        f.write("- Score surface: `outputs_submission_2026-03-03_pasall/primary/rasters/kde_main_1km.tif`\n")
        f.write("- Holdout temporal mode: overlap with AD 750-954\n")
        f.write("- Top area threshold: 15%\n")
        f.write(f"- Negatives: {N_NEG}; Bootstrap: {N_BOOT}\n\n")
        f.write("| scenario | n_validation | n_hits_top15 | hit_rate | auc_uniform | auc_opportunity_null | gain_vs_csr | gain_vs_opportunity |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for r in rows:
            f.write(
                f"| {r['scenario']} | {r['n_validation']} | {r['n_hits_top15']} | "
                f"{float(r['hit_rate']):.3f} | {float(r['auc_uniform']):.3f} | "
                f"{float(r['auc_opportunity_null']):.3f} | {float(r['gain_vs_csr']):.3f} | "
                f"{float(r['gain_vs_opportunity']):.3f} |\n"
            )

    print(f"Wrote: {out_tsv}")
    print(f"Wrote: {summary_md}")


if __name__ == "__main__":
    main()
