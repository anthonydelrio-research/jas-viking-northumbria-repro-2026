#!/usr/bin/env python3
"""
Run publication-grade JAS analyses from current ArcGIS exports.

This script is intentionally self-contained (numpy + GDAL/OGR only).
It computes and writes:
  - validation metrics and fold-level cross-validation outputs
  - KDE robustness/sensitivity outputs
  - environmental association and infrastructure proximity outputs
  - FDR-adjusted p-values
  - manuscript-facing summary tables with explicit provenance

Project: JAS Viking-Age Northumbria reproducibility package
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from osgeo import gdal, ogr, osr


gdal.UseExceptions()
ogr.UseExceptions()


ROOT = Path(os.environ.get("JAS_ROOT", str(Path(__file__).resolve().parents[1]))).resolve()
INBOX = Path(os.environ.get("JAS_INBOX", str(ROOT / "incoming_arcgis" / "01_inbox"))).resolve()
OUT = Path(os.environ.get("JAS_OUTPUT_DIR", str(ROOT / "outputs")))
OUT_TABLES = OUT / "tables"
OUT_RASTERS = OUT / "rasters"
OUT_LOGS = OUT / "logs"

# Deterministic global seed for full reproducibility
SEED = 20260302


def env_int(name: str, default: int) -> int:
    val = os.environ.get(name, "").strip()
    if not val:
        return default
    try:
        return int(val)
    except ValueError:
        return default


def env_float(name: str, default: float) -> float:
    val = os.environ.get(name, "").strip()
    if not val:
        return default
    try:
        return float(val)
    except ValueError:
        return default


def env_bool(name: str, default: bool = False) -> bool:
    val = os.environ.get(name, "").strip().lower()
    if not val:
        return default
    return val in {"1", "true", "yes", "y", "on"}


# Main analytical settings (overridable via environment variables)
CELL_SIZE_M = env_float("JAS_CELL_SIZE_M", 1000.0)
KDE_BANDWIDTH_M = env_float("JAS_KDE_BANDWIDTH_M", 10000.0)
VALIDATION_TOP_AREA_PCT = env_float("JAS_TOP_AREA_PCT", 15.0)
BUFFER_M = KDE_BANDWIDTH_M
N_FOLDS = env_int("JAS_N_FOLDS", 5)
N_NEG_SAMPLES = env_int("JAS_N_NEG_SAMPLES", 20000)
N_BOOTSTRAP = env_int("JAS_N_BOOTSTRAP", 10000)
N_PERM_MORAN = env_int("JAS_N_PERM_MORAN", 9999)
N_PERM_ENV = env_int("JAS_N_PERM_ENV", 9999)
N_PERM_DIST = env_int("JAS_N_PERM_DIST", 9999)
N_SENS_REPS = env_int("JAS_N_SENS_REPS", 100)
INCLUDE_SENSITIVITY_HOLDOUTS = env_bool("JAS_INCLUDE_SENSITIVITY_HOLDOUTS", False)
INCLUDE_CONTEXT_HOLDOUTS = env_bool("JAS_INCLUDE_CONTEXT_HOLDOUTS", False)
TARGET_YEAR_FROM = env_int("JAS_TARGET_YEAR_FROM", 750)
TARGET_YEAR_TO = env_int("JAS_TARGET_YEAR_TO", 954)
PAS_TEMPORAL_MODE = os.environ.get("JAS_PAS_TEMPORAL_MODE", "overlap").strip().lower()
HOLDOUT_TEMPORAL_MODE = os.environ.get("JAS_HOLDOUT_TEMPORAL_MODE", "overlap").strip().lower()
PAS_TEMPORAL_WEIGHTING = os.environ.get("JAS_PAS_TEMPORAL_WEIGHTING", "none").strip().lower()
OPPORTUNITY_LAYER_NAME = os.environ.get("JAS_OPPORTUNITY_LAYER", "opportunity_weighted_background").strip()
USE_OBSERVED_OPPORTUNITY = env_bool("JAS_USE_OBSERVED_OPPORTUNITY", True)
OBSERVED_OPPORTUNITY_MODE = os.environ.get("JAS_OBSERVED_OPPORTUNITY_MODE", "baseline_pas").strip().lower()
PAS_ALL_PATH_ENV = os.environ.get("JAS_PAS_ALL_PATH", "").strip()
PAS_ALL_LAYER_NAME = os.environ.get("JAS_PAS_ALL_LAYER", "").strip()


def ensure_dirs() -> None:
    for d in [OUT, OUT_TABLES, OUT_RASTERS, OUT_LOGS]:
        d.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def read_gpkg_points(gpkg: Path, layer_name: str, fields: Sequence[str]) -> List[Dict[str, object]]:
    ds = ogr.Open(str(gpkg))
    if ds is None:
        raise RuntimeError(f"Failed to open {gpkg}")
    lyr = ds.GetLayerByName(layer_name)
    if lyr is None:
        raise RuntimeError(f"Layer {layer_name} not found in {gpkg}")
    rows: List[Dict[str, object]] = []
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        x = geom.GetX()
        y = geom.GetY()
        row: Dict[str, object] = {"x": float(x), "y": float(y)}
        for f in fields:
            row[f] = feat.GetField(f)
        rows.append(row)
    return rows


def read_vector_points(vector_path: Path, layer_name: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray]:
    ds = ogr.Open(str(vector_path))
    if ds is None:
        raise RuntimeError(f"Failed to open vector points source: {vector_path}")
    lyr = ds.GetLayerByName(layer_name) if layer_name else None
    if lyr is None:
        lyr = ds.GetLayer(0)
    if lyr is None:
        raise RuntimeError(f"Could not resolve layer in {vector_path}")
    xs: List[float] = []
    ys: List[float] = []
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        g = geom.GetGeometryType()
        if g not in {ogr.wkbPoint, ogr.wkbPoint25D}:
            # Defensive handling for accidental multi-geometries.
            if geom.GetGeometryCount() > 0:
                sub = geom.GetGeometryRef(0)
                if sub is None:
                    continue
                xs.append(float(sub.GetX()))
                ys.append(float(sub.GetY()))
                continue
            continue
        xs.append(float(geom.GetX()))
        ys.append(float(geom.GetY()))
    if not xs:
        raise RuntimeError(f"No point geometries found in {vector_path}")
    return np.array(xs, dtype=np.float64), np.array(ys, dtype=np.float64)


def parse_tsv(path: Path, x_keys: Sequence[str], y_keys: Sequence[str], id_key: str) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, r in enumerate(reader, start=1):
            x_val = None
            y_val = None
            for k in x_keys:
                if k in r and str(r[k]).strip():
                    x_val = float(r[k])
                    break
            for k in y_keys:
                if k in r and str(r[k]).strip():
                    y_val = float(r[k])
                    break
            if x_val is None or y_val is None:
                continue
            sid = str(r.get(id_key, "")).strip()
            if not sid:
                sid = f"{path.stem}_{i}"
            rows.append({"site_id": sid, "x": x_val, "y": y_val})
    return rows


def parse_year(value: object, default: int | None = None) -> int | None:
    if value is None:
        return default
    s = str(value).strip()
    if not s:
        return default
    try:
        return int(float(s))
    except ValueError:
        return default


def temporal_overlap_years(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def temporal_pass(a0: int | None, a1: int | None, win0: int, win1: int, mode: str) -> bool:
    if mode == "all":
        return True
    if a0 is None or a1 is None:
        return False
    if a0 > a1:
        a0, a1 = a1, a0
    if mode == "inside":
        return (a0 >= win0) and (a1 <= win1)
    if mode == "overlap":
        return temporal_overlap_years(a0, a1, win0, win1) > 0
    raise ValueError(f"Unsupported temporal mode: {mode}")


def temporal_overlap_fraction(a0: int | None, a1: int | None, win0: int, win1: int) -> float:
    if a0 is None or a1 is None:
        return 0.0
    if a0 > a1:
        a0, a1 = a1, a0
    span = max(1, a1 - a0 + 1)
    ov = temporal_overlap_years(a0, a1, win0, win1)
    return float(ov) / float(span)


def read_holdouts_from_gpkg(
    gpkg_path: Path,
    layer_name: Optional[str] = None,
    include_sensitivity: bool = False,
    include_context: bool = False,
) -> List[Dict[str, object]]:
    ds = ogr.Open(str(gpkg_path))
    if ds is None:
        raise RuntimeError(f"Failed to open holdout GPKG: {gpkg_path}")
    lyr = ds.GetLayerByName(layer_name) if layer_name else ds.GetLayer(0)
    if lyr is None:
        raise RuntimeError(f"Could not resolve holdout layer in {gpkg_path}")

    defn = lyr.GetLayerDefn()
    fields = {defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())}
    has_role = "validation_role" in fields

    allowed_roles = {"primary"}
    if include_sensitivity:
        allowed_roles.add("sensitivity")
    if include_context:
        allowed_roles.add("context_only")

    rows: List[Dict[str, object]] = []
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        role = str(feat.GetField("validation_role") or "").strip().lower() if has_role else ""
        site_type = str(feat.GetField("site_type") or "").strip().lower()

        if has_role:
            if role and role not in allowed_roles:
                continue
        else:
            # Conservative default if role metadata is absent.
            if site_type == "hoard" and not include_context:
                continue

        sid = str(feat.GetField("site_id") or "").strip()
        if not sid:
            sid = f"HOLDOUT_{feat.GetFID()}"

        rows.append(
            {
                "site_id": sid,
                "x": float(geom.GetX()),
                "y": float(geom.GetY()),
                "site_type": site_type if site_type else "unknown",
                "date_from_ad": parse_year(feat.GetField("date_from_ad"), 850),
                "date_to_ad": parse_year(feat.GetField("date_to_ad"), 954),
                "precision_m": int(feat.GetField("precision_m") or 100),
                "source": f"{gpkg_path.name}:{lyr.GetName()}",
            }
        )
    return rows


@dataclass
class Grid:
    xsize: int
    ysize: int
    gt: Tuple[float, float, float, float, float, float]
    proj_wkt: str
    mask_full: np.ndarray  # bool[y, x]
    rows: np.ndarray  # mask indices
    cols: np.ndarray
    xs: np.ndarray  # x coords for masked cells
    ys: np.ndarray  # y coords for masked cells


def build_analysis_grid(boundary_gpkg: Path, layer_name: str, cell_size: float) -> Grid:
    ds = ogr.Open(str(boundary_gpkg))
    lyr = ds.GetLayerByName(layer_name)
    if lyr is None:
        raise RuntimeError(f"Boundary layer {layer_name} missing in {boundary_gpkg}")
    extent = lyr.GetExtent()  # minx, maxx, miny, maxy
    minx, maxx, miny, maxy = extent[0], extent[1], extent[2], extent[3]
    xsize = int(math.ceil((maxx - minx) / cell_size))
    ysize = int(math.ceil((maxy - miny) / cell_size))
    gt = (minx, cell_size, 0.0, maxy, 0.0, -cell_size)
    srs = lyr.GetSpatialRef()
    proj_wkt = srs.ExportToWkt() if srs is not None else ""

    mem = gdal.GetDriverByName("MEM").Create("", xsize, ysize, 1, gdal.GDT_Byte)
    mem.SetGeoTransform(gt)
    mem.SetProjection(proj_wkt)
    band = mem.GetRasterBand(1)
    band.Fill(0)
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[1], options=["ALL_TOUCHED=TRUE"])
    mask = band.ReadAsArray().astype(bool)
    rows, cols = np.where(mask)
    xs = gt[0] + (cols + 0.5) * gt[1]
    ys = gt[3] + (rows + 0.5) * gt[5]
    return Grid(
        xsize=xsize,
        ysize=ysize,
        gt=gt,
        proj_wkt=proj_wkt,
        mask_full=mask,
        rows=rows,
        cols=cols,
        xs=xs.astype(np.float64),
        ys=ys.astype(np.float64),
    )


def coords_to_rowcol(gt: Tuple[float, float, float, float, float, float], x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    col = np.floor((x - gt[0]) / gt[1]).astype(np.int64)
    row = np.floor((y - gt[3]) / gt[5]).astype(np.int64)
    return row, col


def sample_raster_to_mask_cells(raster_path: Path, grid: Grid) -> np.ndarray:
    ds = gdal.Open(str(raster_path), gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f"Failed to open raster {raster_path}")
    arr = ds.GetRasterBand(1).ReadAsArray()
    gt = ds.GetGeoTransform()
    nodata = ds.GetRasterBand(1).GetNoDataValue()
    r, c = coords_to_rowcol(gt, grid.xs, grid.ys)
    valid = (r >= 0) & (r < arr.shape[0]) & (c >= 0) & (c < arr.shape[1])
    out = np.full(grid.xs.shape[0], np.nan, dtype=np.float64)
    out[valid] = arr[r[valid], c[valid]].astype(np.float64)
    if nodata is not None:
        out[np.isclose(out, float(nodata))] = np.nan
    # DEM has extreme nodata sentinel values in this package.
    out[out > 1e20] = np.nan
    return out


def rasterize_vector_to_grid(vector_path: Path, layer_name: str, grid: Grid, burn: float = 1.0, dtype=gdal.GDT_Byte) -> gdal.Dataset:
    ds = ogr.Open(str(vector_path))
    if ds is None:
        raise RuntimeError(f"Failed to open vector {vector_path}")
    lyr = ds.GetLayerByName(layer_name) if layer_name else None
    if lyr is None:
        lyr = ds.GetLayer(0)
    if lyr is None:
        raise RuntimeError(f"Layer {layer_name or '<first>'} not found in {vector_path}")
    mem = gdal.GetDriverByName("MEM").Create("", grid.xsize, grid.ysize, 1, dtype)
    mem.SetGeoTransform(grid.gt)
    mem.SetProjection(grid.proj_wkt)
    band = mem.GetRasterBand(1)
    band.Fill(0)
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[burn], options=["ALL_TOUCHED=TRUE"])
    return mem


def rasterize_filtered_vector_to_grid(
    vector_path: Path,
    layer_name: Optional[str],
    grid: Grid,
    where_sql: Optional[str],
    burn: float = 1.0,
    dtype=gdal.GDT_Byte,
    all_touched: bool = True,
) -> gdal.Dataset:
    ds = ogr.Open(str(vector_path))
    if ds is None:
        raise RuntimeError(f"Failed to open vector {vector_path}")
    lyr = ds.GetLayerByName(layer_name) if layer_name else None
    if lyr is None:
        lyr = ds.GetLayer(0)
    if lyr is None:
        raise RuntimeError(f"Layer {layer_name or '<first>'} not found in {vector_path}")
    minx = grid.gt[0]
    maxx = grid.gt[0] + (grid.xsize * grid.gt[1])
    maxy = grid.gt[3]
    miny = grid.gt[3] + (grid.ysize * grid.gt[5])
    lyr.SetSpatialFilterRect(minx, miny, maxx, maxy)
    if where_sql:
        lyr.SetAttributeFilter(where_sql)
    mem = gdal.GetDriverByName("MEM").Create("", grid.xsize, grid.ysize, 1, dtype)
    mem.SetGeoTransform(grid.gt)
    mem.SetProjection(grid.proj_wkt)
    band = mem.GetRasterBand(1)
    band.Fill(0)
    options = ["ALL_TOUCHED=TRUE"] if all_touched else ["ALL_TOUCHED=FALSE"]
    gdal.RasterizeLayer(mem, [1], lyr, burn_values=[burn], options=options)
    if where_sql:
        lyr.SetAttributeFilter(None)
    lyr.SetSpatialFilter(None)
    return mem


def classify_points_by_region(
    xs: np.ndarray,
    ys: np.ndarray,
    grid: Grid,
    yk_full: np.ndarray,
    nw_full: np.ndarray,
) -> np.ndarray:
    labels = np.full(xs.size, "outside", dtype=object)
    r, c = coords_to_rowcol(grid.gt, xs, ys)
    valid = (r >= 0) & (r < grid.ysize) & (c >= 0) & (c < grid.xsize)
    if not np.any(valid):
        return labels
    yk_hit = np.zeros(xs.size, dtype=bool)
    nw_hit = np.zeros(xs.size, dtype=bool)
    yk_hit[valid] = yk_full[r[valid], c[valid]]
    nw_hit[valid] = nw_full[r[valid], c[valid]]
    labels[valid & yk_hit] = "yorkshire"
    labels[valid & nw_hit] = "northwest"
    labels[valid & (~yk_hit) & (~nw_hit)] = "pennine_other"
    return labels


def proximity_from_vector(vector_path: Path, layer_name: str, grid: Grid) -> np.ndarray:
    src = rasterize_vector_to_grid(vector_path, layer_name, grid, burn=1.0, dtype=gdal.GDT_Byte)
    prox = gdal.GetDriverByName("MEM").Create("", grid.xsize, grid.ysize, 1, gdal.GDT_Float32)
    prox.SetGeoTransform(grid.gt)
    prox.SetProjection(grid.proj_wkt)
    gdal.ComputeProximity(src.GetRasterBand(1), prox.GetRasterBand(1), ["VALUES=1", "DISTUNITS=GEO"])
    arr = prox.GetRasterBand(1).ReadAsArray().astype(np.float64)
    return arr[grid.rows, grid.cols]


def kde_surface(
    gx: np.ndarray,
    gy: np.ndarray,
    px: np.ndarray,
    py: np.ndarray,
    weights: np.ndarray,
    bandwidth_m: float,
    kernel: str = "gaussian",
    point_chunk: int = 128,
) -> np.ndarray:
    if px.size == 0:
        return np.zeros(gx.size, dtype=np.float64)
    out = np.zeros(gx.size, dtype=np.float64)
    bw2 = bandwidth_m * bandwidth_m
    for s in range(0, px.size, point_chunk):
        e = min(px.size, s + point_chunk)
        x = px[s:e][None, :]
        y = py[s:e][None, :]
        w = weights[s:e][None, :]
        dx = gx[:, None] - x
        dy = gy[:, None] - y
        r2 = dx * dx + dy * dy
        if kernel == "gaussian":
            k = np.exp(-0.5 * (r2 / bw2))
        elif kernel == "epanechnikov":
            k = np.maximum(0.0, 1.0 - (r2 / bw2))
        elif kernel == "uniform":
            k = (r2 <= bw2).astype(np.float64)
        else:
            raise ValueError(f"Unsupported kernel: {kernel}")
        out += np.sum(k * w, axis=1)
    return out


def array_to_full(mask_rows: np.ndarray, mask_cols: np.ndarray, ysize: int, xsize: int, values: np.ndarray, fill=np.nan) -> np.ndarray:
    arr = np.full((ysize, xsize), fill, dtype=np.float64)
    arr[mask_rows, mask_cols] = values
    return arr


def zscore(v: np.ndarray) -> np.ndarray:
    m = np.nanmean(v)
    s = np.nanstd(v, ddof=1)
    if not np.isfinite(s) or s == 0:
        return np.zeros_like(v)
    return (v - m) / s


def wilson_ci(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    phat = k / n
    den = 1 + (z * z) / n
    center = (phat + (z * z) / (2 * n)) / den
    half = (z / den) * math.sqrt((phat * (1 - phat) / n) + ((z * z) / (4 * n * n)))
    return max(0.0, center - half), min(1.0, center + half)


def ranks_average_ties(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty_like(order, dtype=np.float64)
    n = x.size
    i = 0
    while i < n:
        j = i + 1
        while j < n and x[order[j]] == x[order[i]]:
            j += 1
        avg_rank = 0.5 * (i + 1 + j)
        ranks[order[i:j]] = avg_rank
        i = j
    return ranks


def auc_from_scores(pos: np.ndarray, neg: np.ndarray) -> float:
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    scores = np.concatenate([pos, neg])
    labels = np.concatenate([np.ones(pos.size, dtype=np.int8), np.zeros(neg.size, dtype=np.int8)])
    ranks = ranks_average_ties(scores)
    sum_ranks_pos = np.sum(ranks[labels == 1])
    auc = (sum_ranks_pos - (pos.size * (pos.size + 1) / 2.0)) / (pos.size * neg.size)
    return float(auc)


def weighted_auc_from_scores(pos: np.ndarray, neg: np.ndarray, pos_w: np.ndarray, neg_w: np.ndarray | None = None) -> float:
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    pos_w = pos_w[np.isfinite(pos_w)]
    if neg_w is None:
        neg_w = np.ones(neg.size, dtype=np.float64)
    else:
        neg_w = neg_w[np.isfinite(neg_w)]
    if pos.size == 0 or neg.size == 0 or pos_w.size != pos.size or neg_w.size != neg.size:
        return float("nan")
    den = float(np.sum(pos_w) * np.sum(neg_w))
    if den <= 0:
        return float("nan")
    num = 0.0
    for i in range(pos.size):
        gt = neg < pos[i]
        eq = neg == pos[i]
        num += pos_w[i] * (float(np.sum(neg_w[gt])) + 0.5 * float(np.sum(neg_w[eq])))
    return float(num / den)


def bootstrap_auc_ci(pos: np.ndarray, neg: np.ndarray, n_boot: int, rng: np.random.Generator) -> Tuple[float, float]:
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if pos.size == 0 or neg.size == 0:
        return float("nan"), float("nan")
    boots = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        pos_b = pos[rng.integers(0, pos.size, pos.size)]
        neg_b = neg[rng.integers(0, neg.size, neg.size)]
        boots[i] = auc_from_scores(pos_b, neg_b)
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def binomial_one_sided_p_ge(k: int, n: int, p0: float) -> float:
    if n <= 0:
        return float("nan")
    if k <= 0:
        return 1.0
    if k > n:
        return float("nan")
    p = 0.0
    for i in range(k, n + 1):
        p += math.comb(n, i) * (p0**i) * ((1.0 - p0) ** (n - i))
    return float(min(1.0, max(0.0, p)))


def beta_posterior_summary(
    k: int,
    n: int,
    p0: float,
    rng: np.random.Generator,
    n_draws: int = 200000,
) -> Tuple[float, float, float]:
    if n <= 0:
        return float("nan"), float("nan"), float("nan")
    a = 1.0 + float(k)
    b = 1.0 + float(n - k)
    draws = rng.beta(a, b, size=n_draws)
    ci_l = float(np.percentile(draws, 2.5))
    ci_u = float(np.percentile(draws, 97.5))
    p_gt = float(np.mean(draws > p0))
    return ci_l, ci_u, p_gt


def effective_sample_size(weights: np.ndarray) -> float:
    w = weights[np.isfinite(weights) & (weights > 0)]
    if w.size == 0:
        return float("nan")
    num = float(np.sum(w) ** 2)
    den = float(np.sum(w * w))
    if den <= 0:
        return float("nan")
    return num / den


def opportunity_surface_modeled(
    arable_vals: np.ndarray,
    slope_vals: np.ndarray,
    road_dist: np.ndarray,
) -> Tuple[np.ndarray, Dict[str, float]]:
    arable_score = np.where(np.isfinite(arable_vals) & (arable_vals > 0), 1.0, 0.0)
    slope_score = np.where(
        np.isfinite(slope_vals),
        np.clip((15.0 - slope_vals) / 15.0, 0.0, 1.0),
        0.0,
    )
    road_score = np.where(
        np.isfinite(road_dist),
        np.exp(-np.maximum(road_dist, 0.0) / 2000.0),
        0.0,
    )
    opp = (0.45 * arable_score) + (0.30 * slope_score) + (0.25 * road_score)
    opp = np.where(np.isfinite(opp), opp, 0.0)
    # Keep a tiny floor to avoid zero-probability cells in weighted-null sampling.
    opp = np.maximum(opp, 1e-6)
    diag = {
        "opportunity_min": float(np.min(opp)),
        "opportunity_max": float(np.max(opp)),
        "opportunity_mean": float(np.mean(opp)),
        "arable_mean": float(np.mean(arable_score)),
        "slope_score_mean": float(np.mean(slope_score)),
        "road_score_mean": float(np.mean(road_score)),
    }
    return opp, diag


def opportunity_surface_observed_constraints(
    constraints_path: Path,
    constraints_layer: Optional[str],
    grid: Grid,
) -> Tuple[np.ndarray, Dict[str, float]]:
    # The delivered "Metal_Detecting_Constraints/opportunity" layer is a polygonal
    # constraint map. We convert it to an observed-effort accessibility surface.
    # We down-weight cells in woodland/ancient woodland and urban polygons.
    urban_ds = rasterize_filtered_vector_to_grid(
        constraints_path,
        constraints_layer,
        grid,
        "LEGEND LIKE '%Urban Area%'",
        burn=1.0,
        dtype=gdal.GDT_Byte,
        all_touched=False,
    )
    woodland_ds = rasterize_filtered_vector_to_grid(
        constraints_path,
        constraints_layer,
        grid,
        "CATEGORY = 'Woodland' OR STATUS = 'ASNW' OR STATUS = 'PAWS'",
        burn=1.0,
        dtype=gdal.GDT_Byte,
        all_touched=False,
    )
    park_ds = rasterize_filtered_vector_to_grid(
        constraints_path,
        constraints_layer,
        grid,
        "LEGEND LIKE '%National Park%'",
        burn=1.0,
        dtype=gdal.GDT_Byte,
        all_touched=False,
    )
    urban_full = urban_ds.GetRasterBand(1).ReadAsArray().astype(bool)
    woodland_full = woodland_ds.GetRasterBand(1).ReadAsArray().astype(bool)
    park_full = park_ds.GetRasterBand(1).ReadAsArray().astype(bool)

    urban = urban_full[grid.rows, grid.cols]
    woodland = woodland_full[grid.rows, grid.cols]
    park = park_full[grid.rows, grid.cols]

    effort = np.ones(grid.xs.size, dtype=np.float64)
    effort = np.where(woodland, effort * 0.20, effort)
    effort = np.where(park, effort * 0.80, effort)
    effort = np.where(urban, effort * 0.05, effort)
    effort = np.maximum(effort, 1e-6)

    diag = {
        "opportunity_min": float(np.min(effort)),
        "opportunity_max": float(np.max(effort)),
        "opportunity_mean": float(np.mean(effort)),
        "arable_mean": float("nan"),
        "slope_score_mean": float("nan"),
        "road_score_mean": float("nan"),
        "constraints_urban_pct": float(np.mean(urban) * 100.0),
        "constraints_woodland_pct": float(np.mean(woodland) * 100.0),
        "constraints_national_park_pct": float(np.mean(park) * 100.0),
    }
    return effort, diag


def opportunity_surface_observed_baseline_pas(
    grid: Grid,
    pas_x: np.ndarray,
    pas_y: np.ndarray,
    pas_grade: np.ndarray,
    bandwidth_m: float,
) -> Tuple[np.ndarray, Dict[str, float]]:
    baseline_idx = pas_grade == 0
    n_baseline = int(np.sum(baseline_idx))
    if n_baseline < 30:
        raise RuntimeError(f"Insufficient baseline PAS points for observed opportunity surface (n={n_baseline})")
    w = np.ones(n_baseline, dtype=np.float64)
    s = kde_surface(
        grid.xs,
        grid.ys,
        pas_x[baseline_idx],
        pas_y[baseline_idx],
        w,
        bandwidth_m,
        kernel="gaussian",
    )
    s = np.where(np.isfinite(s), s, 0.0)
    smin = float(np.min(s))
    smax = float(np.max(s))
    if smax > smin:
        s = (s - smin) / (smax - smin)
    else:
        s = np.ones_like(s, dtype=np.float64)
    s = np.maximum(s, 1e-6)
    diag = {
        "opportunity_min": float(np.min(s)),
        "opportunity_max": float(np.max(s)),
        "opportunity_mean": float(np.mean(s)),
        "arable_mean": float("nan"),
        "slope_score_mean": float("nan"),
        "road_score_mean": float("nan"),
        "constraints_urban_pct": float("nan"),
        "constraints_woodland_pct": float("nan"),
        "constraints_national_park_pct": float("nan"),
        "baseline_pas_points": float(n_baseline),
    }
    return s, diag


def opportunity_surface_observed_pas_all(
    grid: Grid,
    pas_all_x: np.ndarray,
    pas_all_y: np.ndarray,
    bandwidth_m: float,
) -> Tuple[np.ndarray, Dict[str, float]]:
    n_pts = int(pas_all_x.size)
    if n_pts < 100:
        raise RuntimeError(f"Insufficient PAS-all points for detector-effort surface (n={n_pts})")
    w = np.ones(n_pts, dtype=np.float64)
    s = kde_surface(
        grid.xs,
        grid.ys,
        pas_all_x,
        pas_all_y,
        w,
        bandwidth_m,
        kernel="gaussian",
        point_chunk=512,
    )
    s = np.where(np.isfinite(s), s, 0.0)
    smin = float(np.min(s))
    smax = float(np.max(s))
    if smax > smin:
        s = (s - smin) / (smax - smin)
    else:
        s = np.ones_like(s, dtype=np.float64)
    s = np.maximum(s, 1e-6)
    diag = {
        "opportunity_min": float(np.min(s)),
        "opportunity_max": float(np.max(s)),
        "opportunity_mean": float(np.mean(s)),
        "arable_mean": float("nan"),
        "slope_score_mean": float("nan"),
        "road_score_mean": float("nan"),
        "constraints_urban_pct": float("nan"),
        "constraints_woodland_pct": float("nan"),
        "constraints_national_park_pct": float("nan"),
        "baseline_pas_points": float("nan"),
        "pas_all_points": float(n_pts),
    }
    return s, diag


def quadratic_weighted_kappa(r1: np.ndarray, r2: np.ndarray, n_cat: int = 3) -> float:
    if r1.size == 0 or r2.size == 0 or r1.size != r2.size:
        return float("nan")
    conf = np.zeros((n_cat, n_cat), dtype=np.float64)
    for a, b in zip(r1, r2):
        if 0 <= a < n_cat and 0 <= b < n_cat:
            conf[a, b] += 1.0
    n = float(np.sum(conf))
    if n <= 0:
        return float("nan")
    p = conf / n
    row = np.sum(p, axis=1, keepdims=True)
    col = np.sum(p, axis=0, keepdims=True)
    if n_cat <= 1:
        return float("nan")
    i = np.arange(n_cat, dtype=np.float64)[:, None]
    j = np.arange(n_cat, dtype=np.float64)[None, :]
    w = 1.0 - ((i - j) ** 2) / ((n_cat - 1) ** 2)
    po = float(np.sum(w * p))
    pe = float(np.sum(w * (row @ col)))
    den = 1.0 - pe
    if den == 0:
        return float("nan")
    return float((po - pe) / den)


def gwet_ac1_nominal(r1: np.ndarray, r2: np.ndarray, n_cat: int = 3) -> float:
    if r1.size == 0 or r2.size == 0 or r1.size != r2.size:
        return float("nan")
    n = float(r1.size)
    if n <= 0 or n_cat <= 1:
        return float("nan")
    po = float(np.mean(r1 == r2))
    c1 = np.bincount(r1, minlength=n_cat).astype(np.float64)
    c2 = np.bincount(r2, minlength=n_cat).astype(np.float64)
    p = (c1 + c2) / (2.0 * n)
    pe = float(np.sum(p * (1.0 - p)) / (n_cat - 1))
    den = 1.0 - pe
    if den == 0:
        return float("nan")
    return float((po - pe) / den)


def read_sici_validation_ratings(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    r1: List[int] = []
    r2: List[int] = []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            a = str(row.get("rater1_grade", "")).strip()
            b = str(row.get("rater2_grade", "")).strip()
            if not a or not b:
                continue
            try:
                grade_a = int(float(a))
                grade_b = int(float(b))
            except ValueError:
                continue
            if grade_a not in {0, 1, 2} or grade_b not in {0, 1, 2}:
                continue
            r1.append(grade_a)
            r2.append(grade_b)
    return np.array(r1, dtype=np.int64), np.array(r2, dtype=np.int64)


def bootstrap_reliability_ci(
    r1: np.ndarray,
    r2: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    if r1.size == 0 or r2.size == 0 or r1.size != r2.size:
        return (float("nan"), float("nan")), (float("nan"), float("nan"))
    n = r1.size
    kappa_b = np.empty(n_boot, dtype=np.float64)
    ac1_b = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rb1 = r1[idx]
        rb2 = r2[idx]
        kappa_b[i] = quadratic_weighted_kappa(rb1, rb2, n_cat=3)
        ac1_b[i] = gwet_ac1_nominal(rb1, rb2, n_cat=3)
    return (
        (float(np.percentile(kappa_b, 2.5)), float(np.percentile(kappa_b, 97.5))),
        (float(np.percentile(ac1_b, 2.5)), float(np.percentile(ac1_b, 97.5))),
    )


def sample_surface_at_points(surface_full: np.ndarray, gt: Tuple[float, float, float, float, float, float], xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    r, c = coords_to_rowcol(gt, xs, ys)
    valid = (r >= 0) & (r < surface_full.shape[0]) & (c >= 0) & (c < surface_full.shape[1])
    out = np.full(xs.shape[0], np.nan, dtype=np.float64)
    out[valid] = surface_full[r[valid], c[valid]]
    return out


def chi_square_2x2(cond: np.ndarray, zone: np.ndarray) -> Tuple[float, float, int, float, float, float]:
    # cond and zone are bool arrays over the same mask domain.
    a = float(np.sum(cond & zone))
    b = float(np.sum((~cond) & zone))
    c = float(np.sum(cond & (~zone)))
    d = float(np.sum((~cond) & (~zone)))
    total = a + b + c + d
    if total == 0:
        return float("nan"), float("nan"), 1, float("nan"), float("nan"), float("nan")
    row1 = a + b
    row2 = c + d
    col1 = a + c
    col2 = b + d
    e_a = row1 * col1 / total
    e_b = row1 * col2 / total
    e_c = row2 * col1 / total
    e_d = row2 * col2 / total
    chi2 = 0.0
    for o, e in [(a, e_a), (b, e_b), (c, e_c), (d, e_d)]:
        if e > 0:
            chi2 += (o - e) ** 2 / e
    # df=1 closed-form p via erfc
    p = math.erfc(math.sqrt(max(chi2, 0.0) / 2.0))
    zone_pct = (a / row1) if row1 > 0 else float("nan")
    base_pct = (c / row2) if row2 > 0 else float("nan")
    enrichment = (zone_pct / base_pct) if (base_pct and np.isfinite(base_pct)) else float("nan")
    return float(chi2), float(p), 1, float(zone_pct), float(base_pct), float(enrichment)


def normal_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def two_sample_mean_diff_perm(
    zone_vals: np.ndarray,
    base_vals: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
) -> Tuple[float, float]:
    zone_vals = zone_vals[np.isfinite(zone_vals)]
    base_vals = base_vals[np.isfinite(base_vals)]
    if zone_vals.size == 0 or base_vals.size == 0:
        return float("nan"), float("nan")
    obs = float(np.mean(zone_vals) - np.mean(base_vals))
    joined = np.concatenate([zone_vals, base_vals])
    n_zone = zone_vals.size
    count = 0
    for _ in range(n_perm):
        perm = rng.permutation(joined.size)
        z = joined[perm[:n_zone]]
        b = joined[perm[n_zone:]]
        stat = float(np.mean(z) - np.mean(b))
        if abs(stat) >= abs(obs):
            count += 1
    p = (count + 1) / (n_perm + 1)
    return obs, p


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = np.empty_like(ranked)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[i] = prev
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def pearson_r(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    if np.sum(m) < 3:
        return float("nan")
    aa = a[m]
    bb = b[m]
    sa = np.std(aa, ddof=1)
    sb = np.std(bb, ddof=1)
    if sa == 0 or sb == 0:
        return float("nan")
    return float(np.corrcoef(aa, bb)[0, 1])


def morans_i_grid(values_full: np.ndarray, mask_full: np.ndarray) -> float:
    # Global Moran's I with 8-neighbor adjacency.
    vals = values_full.copy()
    vals[~mask_full] = np.nan
    x = vals[mask_full]
    xbar = np.mean(x)
    denom = np.sum((x - xbar) ** 2)
    if denom == 0:
        return float("nan")
    num = 0.0
    s0 = 0
    shifts = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    for dy, dx in shifts:
        src = vals[max(0, dy): vals.shape[0] + min(0, dy), max(0, dx): vals.shape[1] + min(0, dx)]
        dst = vals[max(0, -dy): vals.shape[0] + min(0, -dy), max(0, -dx): vals.shape[1] + min(0, -dx)]
        m = np.isfinite(src) & np.isfinite(dst)
        if not np.any(m):
            continue
        num += np.sum((src[m] - xbar) * (dst[m] - xbar))
        s0 += int(np.sum(m))
    n = x.size
    if s0 == 0:
        return float("nan")
    return float((n / s0) * (num / denom))


def morans_i_permutation(values_full: np.ndarray, mask_full: np.ndarray, n_perm: int, rng: np.random.Generator) -> Tuple[float, float]:
    obs = morans_i_grid(values_full, mask_full)
    if not np.isfinite(obs):
        return float("nan"), float("nan")
    x = values_full[mask_full].copy()
    count = 0
    vals = values_full.copy()
    for _ in range(n_perm):
        rng.shuffle(x)
        vals_perm = vals.copy()
        vals_perm[mask_full] = x
        stat = morans_i_grid(vals_perm, mask_full)
        if abs(stat) >= abs(obs):
            count += 1
    p = (count + 1) / (n_perm + 1)
    return obs, p


def kmeans_2d(points: np.ndarray, k: int, rng: np.random.Generator, n_iter: int = 100) -> np.ndarray:
    n = points.shape[0]
    if n < k:
        labels = np.arange(n) % k
        return labels
    # k-means++ style init
    centers = np.empty((k, 2), dtype=np.float64)
    idx0 = rng.integers(0, n)
    centers[0] = points[idx0]
    dist2 = np.sum((points - centers[0]) ** 2, axis=1)
    for i in range(1, k):
        probs = dist2 / np.sum(dist2)
        idx = rng.choice(n, p=probs)
        centers[i] = points[idx]
        d2_new = np.sum((points - centers[i]) ** 2, axis=1)
        dist2 = np.minimum(dist2, d2_new)
    labels = np.zeros(n, dtype=np.int64)
    for _ in range(n_iter):
        d = np.sum((points[:, None, :] - centers[None, :, :]) ** 2, axis=2)
        new_labels = np.argmin(d, axis=1)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for j in range(k):
            mask = labels == j
            if np.any(mask):
                centers[j] = points[mask].mean(axis=0)
            else:
                centers[j] = points[rng.integers(0, n)]
    # force fold ids to 1..k by center x+y order (deterministic naming)
    order = np.argsort(np.sum(centers, axis=1))
    remap = {int(old): i + 1 for i, old in enumerate(order)}
    return np.array([remap[int(l)] for l in labels], dtype=np.int64)


def write_tsv(path: Path, rows: List[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_gpkg_points(path: Path, layer_name: str, rows: List[Dict[str, object]], proj_epsg: int = 27700) -> None:
    if path.exists():
        path.unlink()
    drv = ogr.GetDriverByName("GPKG")
    ds = drv.CreateDataSource(str(path))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(proj_epsg)
    lyr = ds.CreateLayer(layer_name, srs=srs, geom_type=ogr.wkbPoint)
    for fname, ftype in [
        ("site_id", ogr.OFTString),
        ("site_type", ogr.OFTString),
        ("date_from_ad", ogr.OFTInteger),
        ("date_to_ad", ogr.OFTInteger),
        ("precision_m", ogr.OFTInteger),
        ("source", ogr.OFTString),
    ]:
        lyr.CreateField(ogr.FieldDefn(fname, ftype))
    defn = lyr.GetLayerDefn()
    for r in rows:
        feat = ogr.Feature(defn)
        feat.SetField("site_id", str(r.get("site_id", "")))
        feat.SetField("site_type", str(r.get("site_type", "")))
        feat.SetField("date_from_ad", int(r.get("date_from_ad", 850)))
        feat.SetField("date_to_ad", int(r.get("date_to_ad", 954)))
        feat.SetField("precision_m", int(r.get("precision_m", 100)))
        feat.SetField("source", str(r.get("source", "")))
        g = ogr.Geometry(ogr.wkbPoint)
        g.AddPoint(float(r["x"]), float(r["y"]))
        feat.SetGeometry(g)
        lyr.CreateFeature(feat)
    ds = None


def main() -> None:
    ensure_dirs()
    rng = np.random.default_rng(SEED)
    print(
        "JAS pipeline start | "
        f"cell={CELL_SIZE_M:.0f}m bw={KDE_BANDWIDTH_M:.0f}m "
        f"boot={N_BOOTSTRAP} perm_moran={N_PERM_MORAN} perm_env={N_PERM_ENV} "
        f"perm_dist={N_PERM_DIST} sens_reps={N_SENS_REPS} "
        f"target_window={TARGET_YEAR_FROM}-{TARGET_YEAR_TO} "
        f"pas_temporal_mode={PAS_TEMPORAL_MODE} holdout_temporal_mode={HOLDOUT_TEMPORAL_MODE} "
        f"pas_temporal_weighting={PAS_TEMPORAL_WEIGHTING}"
    )

    # ------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------
    pas_path = INBOX / "PAS_full_northumbria_metalwork.gpkg"
    pas_layer = "PAS_full_northumbria_metalwork"
    kepn_path = INBOX / "KEPN_full_northumbria.gpkg"
    kepn_layer = "KEPN_full_northumbria"
    casss_path = INBOX / "CASSS_full_northumbria.gpkg"
    casss_layer = "CASSS_full_northumbria"
    boundary_path = INBOX / "northumbria_boundary.gpkg"
    boundary_layer = "northumbria_boundary"
    roads_path = INBOX / "roman_roads.gpkg"
    roads_layer = "roman_roads"
    waterways_path = INBOX / "major_waterways.gpkg"
    waterways_layer = "major_waterways"
    soil_path = INBOX / "ALC_soil_grades.tif"
    dem_path = INBOX / "elevation_dem.tif"
    arable_path = INBOX / "landuse_arable.tif"
    yk_poly = INBOX / "Yorkshire_Polygon.shp"
    nw_poly = INBOX / "Northwest_Polygon.shp"
    pas_all_path = Path(PAS_ALL_PATH_ENV) if PAS_ALL_PATH_ENV else (INBOX / "PAS_ALL" / "PAS_ALL_PERIODS_Northumbria.shp")
    pas_all_layer = PAS_ALL_LAYER_NAME if PAS_ALL_LAYER_NAME else None

    # ------------------------------------------------------------------
    # Load point datasets
    # ------------------------------------------------------------------
    print("Stage 1/9: loading inputs")
    valid_temporal_modes = {"all", "inside", "overlap"}
    if PAS_TEMPORAL_MODE not in valid_temporal_modes:
        raise RuntimeError(f"Invalid JAS_PAS_TEMPORAL_MODE={PAS_TEMPORAL_MODE}")
    if HOLDOUT_TEMPORAL_MODE not in valid_temporal_modes:
        raise RuntimeError(f"Invalid JAS_HOLDOUT_TEMPORAL_MODE={HOLDOUT_TEMPORAL_MODE}")
    valid_pas_temporal_weighting = {"none", "overlap_fraction"}
    if PAS_TEMPORAL_WEIGHTING not in valid_pas_temporal_weighting:
        raise RuntimeError(
            f"Invalid JAS_PAS_TEMPORAL_WEIGHTING={PAS_TEMPORAL_WEIGHTING}; "
            f"expected one of {sorted(valid_pas_temporal_weighting)}"
        )
    if TARGET_YEAR_FROM > TARGET_YEAR_TO:
        raise RuntimeError(
            f"Invalid target window: JAS_TARGET_YEAR_FROM={TARGET_YEAR_FROM} > JAS_TARGET_YEAR_TO={TARGET_YEAR_TO}"
        )

    pas = read_gpkg_points(
        pas_path,
        pas_layer,
        fields=["site_id", "sici_grade", "source_subset", "date_from_ad", "date_to_ad", "object_class"],
    )
    pas_pre_temporal_n = len(pas)
    if PAS_TEMPORAL_MODE != "all":
        pas = [
            r
            for r in pas
            if temporal_pass(
                parse_year(r.get("date_from_ad")),
                parse_year(r.get("date_to_ad")),
                TARGET_YEAR_FROM,
                TARGET_YEAR_TO,
                PAS_TEMPORAL_MODE,
            )
        ]
    pas_post_temporal_n = len(pas)
    if not pas:
        raise RuntimeError("PAS dataset is empty.")
    pas_x = np.array([float(r["x"]) for r in pas], dtype=np.float64)
    pas_y = np.array([float(r["y"]) for r in pas], dtype=np.float64)
    pas_grade = np.array([int(r["sici_grade"]) if r["sici_grade"] is not None else 0 for r in pas], dtype=np.int64)
    pas_site = np.array([str(r["site_id"]) for r in pas], dtype=object)
    if PAS_TEMPORAL_WEIGHTING == "overlap_fraction":
        pas_temporal_w = np.array(
            [
                temporal_overlap_fraction(
                    parse_year(r.get("date_from_ad")),
                    parse_year(r.get("date_to_ad")),
                    TARGET_YEAR_FROM,
                    TARGET_YEAR_TO,
                )
                for r in pas
            ],
            dtype=np.float64,
        )
    else:
        pas_temporal_w = np.ones(len(pas), dtype=np.float64)
    pas_temporal_w = np.where(np.isfinite(pas_temporal_w), pas_temporal_w, 0.0)
    pas_temporal_w = np.clip(pas_temporal_w, 0.0, 1.0)
    # SICI weights
    pas_base_w = np.where(pas_grade == 2, 1.0, np.where(pas_grade == 1, 0.5, 0.0)).astype(np.float64)
    pas_w = pas_base_w * pas_temporal_w
    pas_w_unweighted = (pas_grade > 0).astype(np.float64) * pas_temporal_w
    pas_w_grade0 = (pas_grade == 0).astype(np.float64) * pas_temporal_w
    pas_temporal_weight_mean = float(np.mean(pas_temporal_w)) if pas_temporal_w.size > 0 else float("nan")
    pas_temporal_weight_mean_g12 = (
        float(np.mean(pas_temporal_w[pas_grade > 0])) if np.any(pas_grade > 0) else float("nan")
    )
    pas_temporal_weight_effective_n = effective_sample_size(pas_temporal_w)

    # Independent holdouts (vetted layer preferred; burial TSV fallback)
    holdouts: List[Dict[str, object]] = []
    holdout_mode = "unknown"
    holdout_layer_name = os.environ.get("JAS_HOLDOUT_LAYER", "").strip()
    if not holdout_layer_name:
        holdout_layer_name = None
    holdout_gpkg_env = os.environ.get("JAS_HOLDOUT_GPKG", "").strip()
    holdout_gpkg = Path(holdout_gpkg_env) if holdout_gpkg_env else (INBOX / "validation_holdout_sites.gpkg")
    holdouts_pre_temporal = 0
    if holdout_gpkg.exists():
        holdouts = read_holdouts_from_gpkg(
            holdout_gpkg,
            layer_name=holdout_layer_name,
            include_sensitivity=INCLUDE_SENSITIVITY_HOLDOUTS,
            include_context=INCLUDE_CONTEXT_HOLDOUTS,
        )
        holdout_mode = "vetted_gpkg"
    else:
        burials_nw = parse_tsv(
            INBOX / "Scandinavian_Burials.tsv",
            x_keys=["Easting", "Eastings", "easting"],
            y_keys=["Northing", "Northings", "northing"],
            id_key="Location",
        )
        burials_yk = parse_tsv(
            INBOX / "Scandinavian_Burials_Yorkshire.tsv",
            x_keys=["Easting", "Eastings", "easting"],
            y_keys=["Northing", "Northings", "northing"],
            id_key="Site_Name",
        )
        for r in burials_nw:
            holdouts.append(
                {
                    "site_id": str(r["site_id"]),
                    "x": float(r["x"]),
                    "y": float(r["y"]),
                    "site_type": "burial",
                    "date_from_ad": 850,
                    "date_to_ad": 954,
                    "precision_m": 100,
                    "source": "Scandinavian_Burials.tsv",
                }
            )
        for r in burials_yk:
            holdouts.append(
                {
                    "site_id": str(r["site_id"]),
                    "x": float(r["x"]),
                    "y": float(r["y"]),
                    "site_type": "burial",
                    "date_from_ad": 850,
                    "date_to_ad": 954,
                    "precision_m": 100,
                    "source": "Scandinavian_Burials_Yorkshire.tsv",
                }
            )
        holdout_mode = "burial_tsv_fallback"
    holdouts_pre_temporal = len(holdouts)
    if HOLDOUT_TEMPORAL_MODE != "all":
        holdouts = [
            r
            for r in holdouts
            if temporal_pass(
                parse_year(r.get("date_from_ad")),
                parse_year(r.get("date_to_ad")),
                TARGET_YEAR_FROM,
                TARGET_YEAR_TO,
                HOLDOUT_TEMPORAL_MODE,
            )
        ]
    for r in holdouts:
        r["temporal_overlap_fraction"] = temporal_overlap_fraction(
            parse_year(r.get("date_from_ad")),
            parse_year(r.get("date_to_ad")),
            TARGET_YEAR_FROM,
            TARGET_YEAR_TO,
        )
    if not holdouts:
        raise RuntimeError("No holdouts were loaded; cannot compute validation metrics.")
    holdout_path = OUT / "validation_holdout_sites_used.gpkg"
    write_gpkg_points(holdout_path, "validation_holdout_sites_used", holdouts)
    # Backward-compatible legacy filename used in earlier drafts.
    write_gpkg_points(OUT / "validation_holdout_burials.gpkg", "validation_holdout_burials", holdouts)
    hold_x = np.array([float(r["x"]) for r in holdouts], dtype=np.float64)
    hold_y = np.array([float(r["y"]) for r in holdouts], dtype=np.float64)
    hold_id = np.array([str(r["site_id"]) for r in holdouts], dtype=object)
    hold_from = np.array([int(r.get("date_from_ad", TARGET_YEAR_FROM)) for r in holdouts], dtype=np.int64)
    hold_to = np.array([int(r.get("date_to_ad", TARGET_YEAR_TO)) for r in holdouts], dtype=np.int64)
    hold_temporal_w = np.array(
        [float(r.get("temporal_overlap_fraction", 1.0)) for r in holdouts], dtype=np.float64
    )
    holdout_type_counts: Dict[str, int] = {}
    for r in holdouts:
        t = str(r.get("site_type", "unknown")).strip().lower()
        holdout_type_counts[t] = holdout_type_counts.get(t, 0) + 1
    print(
        "Loaded holdouts | "
        f"mode={holdout_mode} temporal_mode={HOLDOUT_TEMPORAL_MODE} "
        f"n={len(holdouts)} pre_temporal={holdouts_pre_temporal} type_counts={holdout_type_counts}"
    )

    # ------------------------------------------------------------------
    # Build 1 km analysis grid + environmental covariates
    # ------------------------------------------------------------------
    print("Stage 2/9: building analysis grid and covariates")
    grid = build_analysis_grid(boundary_path, boundary_layer, CELL_SIZE_M)
    mask_n = int(np.sum(grid.mask_full))
    soil_vals = sample_raster_to_mask_cells(soil_path, grid)
    dem_vals = sample_raster_to_mask_cells(dem_path, grid)
    arable_vals = sample_raster_to_mask_cells(arable_path, grid)
    road_dist = proximity_from_vector(roads_path, roads_layer, grid)
    water_dist = proximity_from_vector(waterways_path, waterways_layer, grid)

    # DEM-derived slope (degrees) on analysis grid
    dem_full = array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, dem_vals, fill=np.nan)
    dem_fill = dem_full.copy()
    dem_med = np.nanmedian(dem_vals)
    dem_fill[np.isnan(dem_fill)] = dem_med
    gy, gx = np.gradient(dem_fill, CELL_SIZE_M, CELL_SIZE_M)
    slope_full = np.degrees(np.arctan(np.sqrt(gx * gx + gy * gy)))
    slope_vals = slope_full[grid.rows, grid.cols]

    # Regions for stratification
    def sample_region_mask(poly_path: Path) -> Tuple[np.ndarray, np.ndarray]:
        ds_r = rasterize_vector_to_grid(poly_path, Path(poly_path).stem, grid, burn=1.0, dtype=gdal.GDT_Byte)
        arr_full = ds_r.GetRasterBand(1).ReadAsArray().astype(bool)
        return arr_full, arr_full[grid.rows, grid.cols]

    yk_region_full, yk_region = sample_region_mask(yk_poly)
    nw_region_full, nw_region = sample_region_mask(nw_poly)
    pn_region = (~yk_region) & (~nw_region)
    pas_region = classify_points_by_region(pas_x, pas_y, grid, yk_region_full, nw_region_full)
    hold_region = classify_points_by_region(hold_x, hold_y, grid, yk_region_full, nw_region_full)

    # ------------------------------------------------------------------
    # Main KDE / KDEd surfaces
    # ------------------------------------------------------------------
    print("Stage 3/9: computing KDE/KDEd surfaces")
    surf_main = kde_surface(grid.xs, grid.ys, pas_x, pas_y, pas_w, KDE_BANDWIDTH_M, kernel="gaussian")
    surf_unweighted = kde_surface(grid.xs, grid.ys, pas_x, pas_y, pas_w_unweighted, KDE_BANDWIDTH_M, kernel="gaussian")
    surf_pre = kde_surface(grid.xs, grid.ys, pas_x, pas_y, pas_w_grade0, KDE_BANDWIDTH_M, kernel="gaussian")
    kded = zscore(surf_main) - zscore(surf_pre)

    threshold = float(np.percentile(surf_main, 100.0 - VALIDATION_TOP_AREA_PCT))
    zone = surf_main >= threshold
    zone_count = int(np.sum(zone))

    # Save core rasters
    for name, vals in [("kde_main", surf_main), ("kde_unweighted", surf_unweighted), ("kded_main", kded)]:
        arr_full = array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, vals, fill=np.nan)
        out_r = OUT_RASTERS / f"{name}_1km.tif"
        drv = gdal.GetDriverByName("GTiff")
        ds = drv.Create(str(out_r), grid.xsize, grid.ysize, 1, gdal.GDT_Float32, options=["COMPRESS=LZW"])
        ds.SetGeoTransform(grid.gt)
        ds.SetProjection(grid.proj_wkt)
        band = ds.GetRasterBand(1)
        nodata = -9999.0
        band.SetNoDataValue(nodata)
        out_arr = np.where(np.isfinite(arr_full), arr_full, nodata).astype(np.float32)
        band.WriteArray(out_arr)
        band.FlushCache()
        ds = None
    zone_full = np.zeros((grid.ysize, grid.xsize), dtype=np.uint8)
    zone_full[grid.rows, grid.cols] = zone.astype(np.uint8)
    ds = gdal.GetDriverByName("GTiff").Create(
        str(OUT_RASTERS / "high_confidence_zone_top15pct_1km.tif"),
        grid.xsize,
        grid.ysize,
        1,
        gdal.GDT_Byte,
        options=["COMPRESS=LZW"],
    )
    ds.SetGeoTransform(grid.gt)
    ds.SetProjection(grid.proj_wkt)
    ds.GetRasterBand(1).WriteArray(zone_full)
    ds = None

    # ------------------------------------------------------------------
    # Validation (full model on holdouts)
    # ------------------------------------------------------------------
    print("Stage 4/9: validation metrics")
    surf_full = array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, surf_main, fill=np.nan)
    hold_scores = sample_surface_at_points(surf_full, grid.gt, hold_x, hold_y)
    hold_hits = hold_scores >= threshold
    n_val = int(np.sum(np.isfinite(hold_scores)))
    n_hits = int(np.sum(hold_hits & np.isfinite(hold_scores)))
    hit_rate = n_hits / n_val if n_val > 0 else float("nan")
    ci_l, ci_u = wilson_ci(n_hits, n_val)
    hold_temporal_w = np.where(np.isfinite(hold_temporal_w), hold_temporal_w, 0.0)
    temporal_weight_den = float(np.sum(hold_temporal_w))
    temporal_weight_num = float(np.sum(hold_temporal_w * (hold_hits & np.isfinite(hold_scores)).astype(np.float64)))
    hit_rate_temporal_weighted = (temporal_weight_num / temporal_weight_den) if temporal_weight_den > 0 else float("nan")

    # Uniform negatives for AUC
    neg_idx = rng.choice(mask_n, size=min(N_NEG_SAMPLES, mask_n), replace=False)
    neg_scores = surf_main[neg_idx]
    auc = auc_from_scores(hold_scores, neg_scores)
    auc_temporal_weighted = weighted_auc_from_scores(hold_scores, neg_scores, hold_temporal_w)
    auc_ci_l, auc_ci_u = bootstrap_auc_ci(hold_scores, neg_scores, N_BOOTSTRAP, rng)

    # Opportunity surface (observed detector-effort constraints when available;
    # otherwise modeled covariate fallback).
    opp_constraints_path = INBOX / "opportunity_weighted_background.gpkg"
    opp_surface_mode = ""
    opp_diag: Dict[str, float]
    if USE_OBSERVED_OPPORTUNITY:
        if OBSERVED_OPPORTUNITY_MODE == "baseline_pas":
            try:
                opp_weight, opp_diag = opportunity_surface_observed_baseline_pas(
                    grid, pas_x, pas_y, pas_grade, KDE_BANDWIDTH_M
                )
                opp_surface_mode = "observed_detector_effort_baseline_pas_kde"
            except Exception as exc:
                print(f"Observed baseline PAS opportunity failed ({exc}); trying constraints layer fallback.")
                if opp_constraints_path.exists():
                    try:
                        opp_weight, opp_diag = opportunity_surface_observed_constraints(
                            opp_constraints_path,
                            OPPORTUNITY_LAYER_NAME if OPPORTUNITY_LAYER_NAME else None,
                            grid,
                        )
                        opp_surface_mode = "observed_detector_effort_constraints_layer"
                    except Exception as exc2:
                        print(f"Observed constraints opportunity failed ({exc2}); falling back to modeled covariates.")
                        opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
                        opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
                else:
                    opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
                    opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
        elif OBSERVED_OPPORTUNITY_MODE == "constraints":
            if opp_constraints_path.exists():
                try:
                    opp_weight, opp_diag = opportunity_surface_observed_constraints(
                        opp_constraints_path,
                        OPPORTUNITY_LAYER_NAME if OPPORTUNITY_LAYER_NAME else None,
                        grid,
                    )
                    opp_surface_mode = "observed_detector_effort_constraints_layer"
                except Exception as exc:
                    print(f"Observed constraints opportunity failed ({exc}); falling back to modeled covariates.")
                    opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
                    opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
            else:
                opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
                opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
        elif OBSERVED_OPPORTUNITY_MODE in {"pas_all", "all_pas", "pasall"}:
            if pas_all_path.exists():
                try:
                    pas_all_x, pas_all_y = read_vector_points(pas_all_path, pas_all_layer)
                    opp_weight, opp_diag = opportunity_surface_observed_pas_all(
                        grid, pas_all_x, pas_all_y, KDE_BANDWIDTH_M
                    )
                    opp_surface_mode = "observed_detector_effort_pas_all_kde"
                except Exception as exc:
                    print(f"Observed PAS-all opportunity failed ({exc}); trying constraints layer fallback.")
                    if opp_constraints_path.exists():
                        try:
                            opp_weight, opp_diag = opportunity_surface_observed_constraints(
                                opp_constraints_path,
                                OPPORTUNITY_LAYER_NAME if OPPORTUNITY_LAYER_NAME else None,
                                grid,
                            )
                            opp_surface_mode = "observed_detector_effort_constraints_layer"
                        except Exception as exc2:
                            print(f"Observed constraints opportunity failed ({exc2}); falling back to modeled covariates.")
                            opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
                            opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
                    else:
                        opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
                        opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
            else:
                print(f"PAS-all opportunity source missing at {pas_all_path}; trying constraints fallback.")
                if opp_constraints_path.exists():
                    try:
                        opp_weight, opp_diag = opportunity_surface_observed_constraints(
                            opp_constraints_path,
                            OPPORTUNITY_LAYER_NAME if OPPORTUNITY_LAYER_NAME else None,
                            grid,
                        )
                        opp_surface_mode = "observed_detector_effort_constraints_layer"
                    except Exception as exc:
                        print(f"Observed constraints opportunity failed ({exc}); falling back to modeled covariates.")
                        opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
                        opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
                else:
                    opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
                    opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
        else:
            print(
                f"Unknown JAS_OBSERVED_OPPORTUNITY_MODE={OBSERVED_OPPORTUNITY_MODE}; "
                "falling back to modeled covariates."
            )
            opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
            opp_surface_mode = "modeled_continuous_covariate_opportunity_fallback"
    else:
        opp_weight, opp_diag = opportunity_surface_modeled(arable_vals, slope_vals, road_dist)
        opp_surface_mode = "modeled_continuous_covariate_opportunity"
    for k in [
        "constraints_urban_pct",
        "constraints_woodland_pct",
        "constraints_national_park_pct",
        "baseline_pas_points",
        "pas_all_points",
    ]:
        if k not in opp_diag:
            opp_diag[k] = float("nan")
    opp_thr = float(np.percentile(opp_weight, 100.0 - VALIDATION_TOP_AREA_PCT))
    opp_top = opp_weight >= opp_thr
    opp_full = array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, opp_weight, fill=np.nan)
    hold_opp = sample_surface_at_points(opp_full, grid.gt, hold_x, hold_y)
    opp_hit_rate = float(np.mean(hold_opp >= opp_thr))
    gain_vs_csr = hit_rate / (VALIDATION_TOP_AREA_PCT / 100.0)
    gain_vs_opp = (hit_rate / opp_hit_rate) if opp_hit_rate > 0 else float("nan")

    # AUC against opportunity-weighted negatives (inhomogeneous null comparator).
    opp_prob = opp_weight.copy()
    opp_prob[~np.isfinite(opp_prob)] = 0.0
    ps = float(np.sum(opp_prob))
    if ps <= 0:
        neg_idx_opp = rng.choice(mask_n, size=min(N_NEG_SAMPLES, mask_n), replace=False)
        opp_sampling_note = f"{opp_surface_mode};fallback_uniform_due_zero_opportunity_mass"
    else:
        opp_prob = opp_prob / ps
        neg_idx_opp = rng.choice(mask_n, size=min(N_NEG_SAMPLES, mask_n), replace=False, p=opp_prob)
        opp_sampling_note = opp_surface_mode
    neg_scores_opp = surf_main[neg_idx_opp]
    auc_oppnull = auc_from_scores(hold_scores, neg_scores_opp)
    auc_oppnull_temporal = weighted_auc_from_scores(hold_scores, neg_scores_opp, hold_temporal_w)
    auc_oppnull_ci_l, auc_oppnull_ci_u = bootstrap_auc_ci(hold_scores, neg_scores_opp, N_BOOTSTRAP, rng)

    # Unweighted comparator
    surf_unw_full = array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, surf_unweighted, fill=np.nan)
    thr_unw = float(np.percentile(surf_unweighted, 100.0 - VALIDATION_TOP_AREA_PCT))
    hold_unw = sample_surface_at_points(surf_unw_full, grid.gt, hold_x, hold_y)
    hit_unw = float(np.mean(hold_unw >= thr_unw))
    auc_unw = auc_from_scores(hold_unw, neg_scores)

    # Small-sample diagnostics for primary holdout scenario.
    csr_p0 = VALIDATION_TOP_AREA_PCT / 100.0
    hit_binom_p = binomial_one_sided_p_ge(n_hits, n_val, csr_p0)
    hit_post_ci_l, hit_post_ci_u, hit_post_prob_gt_csr = beta_posterior_summary(n_hits, n_val, csr_p0, rng)

    loo_rows: List[Dict[str, object]] = []
    valid_mask = np.isfinite(hold_scores)
    valid_idx = np.where(valid_mask)[0]
    for i in valid_idx:
        keep = valid_mask.copy()
        keep[i] = False
        sc = hold_scores[keep]
        if sc.size == 0:
            continue
        hits_i = int(np.sum(sc >= threshold))
        n_i = int(sc.size)
        hit_i = float(hits_i / n_i)
        auc_i = auc_from_scores(sc, neg_scores)
        loo_rows.append(
            {
                "dropped_site_id": str(hold_id[i]),
                "n_validation_after_drop": n_i,
                "hit_rate_after_drop": hit_i,
                "delta_hit_rate_pct_points": (hit_i - hit_rate) * 100.0,
                "auc_after_drop": auc_i,
                "delta_auc": auc_i - auc,
            }
        )
    if loo_rows:
        hit_deltas = np.array([float(r["delta_hit_rate_pct_points"]) for r in loo_rows], dtype=np.float64)
        auc_deltas = np.array([float(r["delta_auc"]) for r in loo_rows], dtype=np.float64)
        loo_summary = {
            "n_holdouts": int(n_val),
            "max_abs_delta_hit_rate_pct_points": float(np.max(np.abs(hit_deltas))),
            "max_abs_delta_auc": float(np.max(np.abs(auc_deltas))),
            "mean_abs_delta_hit_rate_pct_points": float(np.mean(np.abs(hit_deltas))),
            "mean_abs_delta_auc": float(np.mean(np.abs(auc_deltas))),
        }
    else:
        loo_summary = {
            "n_holdouts": int(n_val),
            "max_abs_delta_hit_rate_pct_points": float("nan"),
            "max_abs_delta_auc": float("nan"),
            "mean_abs_delta_hit_rate_pct_points": float("nan"),
            "mean_abs_delta_auc": float("nan"),
        }

    # Decile lift diagnostics (ranking utility calibration table).
    decile_rows: List[Dict[str, object]] = []
    decile_bin_rows: List[Dict[str, object]] = []
    hold_valid = np.isfinite(hold_scores)
    hold_n_valid = int(np.sum(hold_valid))
    for top_pct in range(10, 101, 10):
        thr_pct = float(np.percentile(surf_main, 100.0 - float(top_pct)))
        captured = (hold_scores >= thr_pct) & hold_valid
        n_cap = int(np.sum(captured))
        cap_rate = (n_cap / hold_n_valid) if hold_n_valid > 0 else float("nan")
        expected_random = float(top_pct) / 100.0
        lift_random = (cap_rate / expected_random) if expected_random > 0 and np.isfinite(cap_rate) else float("nan")
        opp_thr_pct = float(np.percentile(opp_weight, 100.0 - float(top_pct)))
        opp_cap_rate = float(np.mean(hold_opp[hold_valid] >= opp_thr_pct)) if hold_n_valid > 0 else float("nan")
        lift_opp = (cap_rate / opp_cap_rate) if np.isfinite(opp_cap_rate) and opp_cap_rate > 0 else float("nan")
        decile_rows.append(
            {
                "top_area_pct": top_pct,
                "n_validation": hold_n_valid,
                "n_captured": n_cap,
                "capture_rate": cap_rate,
                "expected_random_rate": expected_random,
                "lift_vs_random": lift_random,
                "opportunity_capture_rate": opp_cap_rate,
                "lift_vs_opportunity": lift_opp,
            }
        )

    edges = np.percentile(surf_main, np.linspace(0.0, 100.0, 11))
    for i in range(10):
        lo = edges[i]
        hi = edges[i + 1]
        if i == 9:
            in_bin = (hold_scores >= lo) & (hold_scores <= hi) & hold_valid
        else:
            in_bin = (hold_scores >= lo) & (hold_scores < hi) & hold_valid
        n_bin = int(np.sum(in_bin))
        obs_rate = (n_bin / hold_n_valid) if hold_n_valid > 0 else float("nan")
        exp_rate = 0.10
        lift = (obs_rate / exp_rate) if exp_rate > 0 and np.isfinite(obs_rate) else float("nan")
        decile_bin_rows.append(
            {
                "score_decile": i + 1,
                "score_lower": float(lo),
                "score_upper": float(hi),
                "n_validation": hold_n_valid,
                "n_in_bin": n_bin,
                "observed_rate": obs_rate,
                "expected_rate": exp_rate,
                "lift_vs_uniform": lift,
            }
        )

    # ------------------------------------------------------------------
    # Spatial cross-validation (5-fold, independent holdouts)
    # ------------------------------------------------------------------
    print("Stage 5/9: spatial cross-validation")
    hold_pts = np.column_stack([hold_x, hold_y])
    folds = kmeans_2d(hold_pts, N_FOLDS, rng=rng, n_iter=100)

    fold_rows: List[Dict[str, object]] = []
    excl_rows: List[Dict[str, object]] = []
    leakage_lines: List[str] = []
    for f in range(1, N_FOLDS + 1):
        val_m = folds == f
        val_x = hold_x[val_m]
        val_y = hold_y[val_m]
        val_scores = hold_scores[val_m]
        # exclude PAS points within buffer of any validation point
        if val_x.size == 0:
            continue
        dx = pas_x[:, None] - val_x[None, :]
        dy = pas_y[:, None] - val_y[None, :]
        min_d = np.sqrt(np.min(dx * dx + dy * dy, axis=1))
        train_keep = min_d > BUFFER_M
        min_train_val = float(np.min(min_d[train_keep])) if np.any(train_keep) else float("nan")
        leakage_lines.append(
            f"Fold {f}: validation_n={val_x.size}, training_n={int(np.sum(train_keep))}, "
            f"excluded_n={int(np.sum(~train_keep))}, min_train_val_dist_m={min_train_val:.2f}"
        )
        tr_x = pas_x[train_keep]
        tr_y = pas_y[train_keep]
        tr_w = pas_w[train_keep]
        surf_f = kde_surface(grid.xs, grid.ys, tr_x, tr_y, tr_w, KDE_BANDWIDTH_M, kernel="gaussian")
        thr_f = float(np.percentile(surf_f, 100.0 - VALIDATION_TOP_AREA_PCT))
        surf_f_full = array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, surf_f, fill=np.nan)
        val_sc = sample_surface_at_points(surf_f_full, grid.gt, val_x, val_y)
        hit_f = float(np.mean(val_sc >= thr_f))
        neg_idx_f = rng.choice(mask_n, size=min(N_NEG_SAMPLES, mask_n), replace=False)
        neg_sc_f = surf_f[neg_idx_f]
        auc_f = auc_from_scores(val_sc, neg_sc_f)
        # Opportunity gain in fold
        val_opp = sample_surface_at_points(opp_full, grid.gt, val_x, val_y)
        opp_hit_f = float(np.mean(val_opp >= opp_thr))
        gain_csr_f = hit_f / (VALIDATION_TOP_AREA_PCT / 100.0)
        gain_opp_f = (hit_f / opp_hit_f) if opp_hit_f > 0 else float("nan")
        fold_rows.append(
            {
                "fold": f,
                "validation_sites": int(val_x.size),
                "training_sites": int(np.sum(train_keep)),
                "excluded_sites": int(np.sum(~train_keep)),
                "hit_rate_15pct": hit_f,
                "auc": auc_f,
                "baseline_gain_csr": gain_csr_f,
                "baseline_gain_opportunity": gain_opp_f,
            }
        )
        excl_rows.append(
            {
                "fold": f,
                "original_training_sites": int(pas_x.size),
                "sites_excluded_by_buffer": int(np.sum(~train_keep)),
                "final_training_sites": int(np.sum(train_keep)),
                "exclusion_rate_pct": 100.0 * float(np.sum(~train_keep)) / float(pas_x.size),
            }
        )

    # Summaries
    fold_hit = np.array([r["hit_rate_15pct"] for r in fold_rows], dtype=np.float64)
    fold_auc = np.array([r["auc"] for r in fold_rows], dtype=np.float64)
    fold_gain_csr = np.array([r["baseline_gain_csr"] for r in fold_rows], dtype=np.float64)
    fold_gain_opp = np.array([r["baseline_gain_opportunity"] for r in fold_rows], dtype=np.float64)

    def metric_summary(name: str, v: np.ndarray) -> Dict[str, object]:
        m = float(np.mean(v))
        sd = float(np.std(v, ddof=1)) if v.size > 1 else 0.0
        se = float(sd / math.sqrt(v.size)) if v.size > 1 else 0.0
        ci_l = float(m - 1.96 * se)
        ci_u = float(m + 1.96 * se)
        return {"metric": name, "mean": m, "sd": sd, "se": se, "ci_lower": ci_l, "ci_upper": ci_u}

    cv_summary_rows = [
        metric_summary("hit_rate_15pct", fold_hit),
        metric_summary("auc", fold_auc),
        metric_summary("baseline_gain_csr", fold_gain_csr),
        metric_summary("baseline_gain_opportunity", fold_gain_opp),
    ]

    # Region-held-out transportability tests (Yorkshire <-> Northwest).
    region_transport_rows: List[Dict[str, object]] = []

    def run_region_transport(train_region_name: str, test_region_name: str) -> None:
        train_keep = pas_region == train_region_name
        test_keep = hold_region == test_region_name
        n_train = int(np.sum(train_keep))
        n_test = int(np.sum(test_keep))
        row: Dict[str, object] = {
            "scenario": f"train_{train_region_name}_test_{test_region_name}",
            "train_region": train_region_name,
            "test_region": test_region_name,
            "n_train_pas": n_train,
            "n_test_holdouts": n_test,
            "n_hits_top15": "",
            "hit_rate": "",
            "auc": "",
            "gain_vs_csr": "",
            "gain_vs_opportunity": "",
            "status": "ok",
        }
        if n_train < 10:
            row["status"] = "insufficient_training_points"
            region_transport_rows.append(row)
            return
        if n_test < 1:
            row["status"] = "no_test_holdouts"
            region_transport_rows.append(row)
            return
        surf_t = kde_surface(
            grid.xs,
            grid.ys,
            pas_x[train_keep],
            pas_y[train_keep],
            pas_w[train_keep],
            KDE_BANDWIDTH_M,
            kernel="gaussian",
        )
        thr_t = float(np.percentile(surf_t, 100.0 - VALIDATION_TOP_AREA_PCT))
        surf_t_full = array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, surf_t, fill=np.nan)
        test_scores = sample_surface_at_points(surf_t_full, grid.gt, hold_x[test_keep], hold_y[test_keep])
        valid_t = np.isfinite(test_scores)
        n_valid_t = int(np.sum(valid_t))
        if n_valid_t < 1:
            row["status"] = "no_finite_test_scores"
            region_transport_rows.append(row)
            return
        hits_t = int(np.sum((test_scores >= thr_t) & valid_t))
        hit_t = hits_t / n_valid_t
        neg_idx_t = rng.choice(mask_n, size=min(N_NEG_SAMPLES, mask_n), replace=False)
        neg_scores_t = surf_t[neg_idx_t]
        auc_t = auc_from_scores(test_scores, neg_scores_t)
        test_opp = sample_surface_at_points(opp_full, grid.gt, hold_x[test_keep], hold_y[test_keep])
        opp_hit_t = float(np.mean(test_opp[valid_t] >= opp_thr)) if np.any(valid_t) else float("nan")
        gain_csr_t = hit_t / (VALIDATION_TOP_AREA_PCT / 100.0)
        gain_opp_t = (hit_t / opp_hit_t) if np.isfinite(opp_hit_t) and opp_hit_t > 0 else float("nan")
        row.update(
            {
                "n_hits_top15": hits_t,
                "hit_rate": hit_t,
                "auc": auc_t,
                "gain_vs_csr": gain_csr_t,
                "gain_vs_opportunity": gain_opp_t,
            }
        )
        region_transport_rows.append(row)

    run_region_transport("yorkshire", "northwest")
    run_region_transport("northwest", "yorkshire")
    run_region_transport("yorkshire", "yorkshire")
    run_region_transport("northwest", "northwest")

    # ------------------------------------------------------------------
    # Moran's I (spatial autocorrelation of main KDE)
    # ------------------------------------------------------------------
    print("Stage 6/9: Moran's I permutations")
    moran_full = array_to_full(grid.rows, grid.cols, grid.ysize, grid.xsize, zscore(surf_main), fill=np.nan)
    moran_i, moran_p = morans_i_permutation(moran_full, grid.mask_full, N_PERM_MORAN, rng)

    # ------------------------------------------------------------------
    # Sensitivity analyses
    # ------------------------------------------------------------------
    print("Stage 7/9: sensitivity analyses")
    bandwidths = [6000.0, 8000.0, 10000.0, 12000.0, 14000.0]
    surf_h10 = surf_main.copy()
    bw_rows: List[Dict[str, object]] = []
    surf_cache: Dict[float, np.ndarray] = {10000.0: surf_h10}
    for h in bandwidths:
        if h == 10000.0:
            s = surf_h10
        else:
            s = kde_surface(grid.xs, grid.ys, pas_x, pas_y, pas_w, h, kernel="gaussian")
            surf_cache[h] = s
        r = pearson_r(s, surf_h10)
        bw_rows.append(
            {
                "bandwidth_km": h / 1000.0,
                "surface_correlation_h10": r,
                "pattern_preserved": "YES" if (np.isfinite(r) and r >= 0.85) else "NO",
            }
        )

    kernel_rows: List[Dict[str, object]] = []
    surf_gauss = surf_h10
    for k in ["gaussian", "epanechnikov", "uniform"]:
        s = kde_surface(grid.xs, grid.ys, pas_x, pas_y, pas_w, KDE_BANDWIDTH_M, kernel=k)
        r = pearson_r(s, surf_gauss)
        kernel_rows.append(
            {
                "kernel_function": k,
                "surface_correlation_with_gaussian": r if k != "gaussian" else 1.0,
                "regional_pattern_fidelity": "PRESERVED" if (k == "gaussian" or (np.isfinite(r) and r >= 0.80)) else "DEGRADED",
            }
        )

    # Coordinate perturbation
    coord_rows: List[Dict[str, object]] = []
    zone_ref = surf_h10 >= np.percentile(surf_h10, 85.0)
    for mag in [500.0, 1000.0]:
        cors = np.empty(N_SENS_REPS, dtype=np.float64)
        overlaps = np.empty(N_SENS_REPS, dtype=np.float64)
        for i in range(N_SENS_REPS):
            ang = rng.uniform(0.0, 2.0 * math.pi, size=pas_x.size)
            rad = np.sqrt(rng.uniform(0.0, 1.0, size=pas_x.size)) * mag
            jx = pas_x + (rad * np.cos(ang))
            jy = pas_y + (rad * np.sin(ang))
            s = kde_surface(grid.xs, grid.ys, jx, jy, pas_w, KDE_BANDWIDTH_M, kernel="gaussian")
            cors[i] = pearson_r(s, surf_h10)
            zone_i = s >= np.percentile(s, 85.0)
            inter = np.sum(zone_ref & zone_i)
            union = np.sum(zone_ref | zone_i)
            overlaps[i] = (inter / union) if union > 0 else float("nan")
        coord_rows.append(
            {
                "perturbation_magnitude_m": int(mag),
                "replicates": N_SENS_REPS,
                "surface_correlation_mean": float(np.nanmean(cors)),
                "surface_correlation_sd": float(np.nanstd(cors, ddof=1)),
                "high_conf_zone_overlap_pct_mean": float(np.nanmean(overlaps) * 100.0),
                "high_conf_zone_overlap_pct_sd": float(np.nanstd(overlaps, ddof=1) * 100.0),
            }
        )

    # Grade perturbation
    grade_rows: List[Dict[str, object]] = []
    scenarios = [
        ("10pct_G2_to_G1", 2, 1),
        ("10pct_G1_to_G2", 1, 2),
        ("10pct_G1_to_G0", 1, 0),
    ]
    for scen_name, from_g, to_g in scenarios:
        idx = np.where(pas_grade == from_g)[0]
        n_move = max(1, int(round(0.10 * idx.size))) if idx.size > 0 else 0
        cors = np.empty(N_SENS_REPS, dtype=np.float64)
        overlaps = np.empty(N_SENS_REPS, dtype=np.float64)
        for i in range(N_SENS_REPS):
            g_mod = pas_grade.copy()
            if n_move > 0:
                moved = rng.choice(idx, size=n_move, replace=False)
                g_mod[moved] = to_g
            w_mod = np.where(g_mod == 2, 1.0, np.where(g_mod == 1, 0.5, 0.0)).astype(np.float64)
            s = kde_surface(grid.xs, grid.ys, pas_x, pas_y, w_mod, KDE_BANDWIDTH_M, kernel="gaussian")
            cors[i] = pearson_r(s, surf_h10)
            zone_i = s >= np.percentile(s, 85.0)
            inter = np.sum(zone_ref & zone_i)
            union = np.sum(zone_ref | zone_i)
            overlaps[i] = (inter / union) if union > 0 else float("nan")
        grade_rows.append(
            {
                "scenario": scen_name,
                "replicates": N_SENS_REPS,
                "surface_correlation_mean": float(np.nanmean(cors)),
                "surface_correlation_sd": float(np.nanstd(cors, ddof=1)),
                "high_conf_zone_overlap_pct_mean": float(np.nanmean(overlaps) * 100.0),
                "high_conf_zone_overlap_pct_sd": float(np.nanstd(overlaps, ddof=1) * 100.0),
            }
        )

    # Jackknife proxy removal for interpretive integration
    kepn = read_gpkg_points(
        kepn_path,
        kepn_layer,
        fields=["elements__langcode", "elements__language"],
    )
    kepn_keep = [
        r
        for r in kepn
        if str(r.get("elements__langcode", "")).upper() in {"N", "ND", "NW", "NN"}
        or "Old Norse" in str(r.get("elements__language", ""))
    ]
    kx = np.array([float(r["x"]) for r in kepn_keep], dtype=np.float64)
    ky = np.array([float(r["y"]) for r in kepn_keep], dtype=np.float64)
    kw = np.ones(kx.size, dtype=np.float64)
    casss = read_gpkg_points(casss_path, casss_layer, fields=["Category"])
    casss_keep = [r for r in casss if str(r.get("Category", "")) in {"Anglo-Scandinavian", "Hogback"}]
    cx = np.array([float(r["x"]) for r in casss_keep], dtype=np.float64)
    cy = np.array([float(r["y"]) for r in casss_keep], dtype=np.float64)
    cw = np.ones(cx.size, dtype=np.float64)

    surf_pas_z = zscore(surf_main)
    surf_kepn_z = zscore(kde_surface(grid.xs, grid.ys, kx, ky, kw, KDE_BANDWIDTH_M, kernel="gaussian"))
    surf_casss_z = zscore(kde_surface(grid.xs, grid.ys, cx, cy, cw, KDE_BANDWIDTH_M, kernel="gaussian"))
    surf_full_int = surf_pas_z + surf_kepn_z + surf_casss_z

    jack_rows: List[Dict[str, object]] = []
    variants = {
        "remove_PAS": surf_kepn_z + surf_casss_z,
        "remove_KEPN": surf_pas_z + surf_casss_z,
        "remove_CASSS": surf_pas_z + surf_kepn_z,
    }
    for name, s in variants.items():
        r = pearson_r(s, surf_full_int)
        jack_rows.append(
            {
                "removed_proxy": name.replace("remove_", ""),
                "surface_correlation_with_full": r,
                "regional_pattern_preservation": "PRESERVED" if np.isfinite(r) and r >= 0.75 else "DEGRADED",
            }
        )

    # ------------------------------------------------------------------
    # Environmental associations
    # ------------------------------------------------------------------
    print("Stage 8/9: environmental analyses")
    soil_int = np.where(np.isfinite(soil_vals), soil_vals.astype(np.int64), -1)
    cond_soil1 = soil_int == 1
    cond_soil13 = (soil_int >= 1) & (soil_int <= 3)
    cond_elev244 = np.isfinite(dem_vals) & (dem_vals < 244.0)
    cond_road2 = road_dist <= 2000.0
    cond_road5 = road_dist <= 5000.0
    cond_road10 = road_dist <= 10000.0
    cond_water1 = water_dist <= 1000.0
    cond_water2 = water_dist <= 2000.0
    cond_water5 = water_dist <= 5000.0

    # chi-square rows
    chi_rows: List[Dict[str, object]] = []
    for tid, cond in [
        ("JAS4_13_grade1", cond_soil1),
        ("JAS4_13", cond_soil13),
    ]:
        chi2, p, df, zp, bp, enr = chi_square_2x2(cond, zone)
        chi_rows.append(
            {
                "test_id": tid,
                "chi_square": chi2,
                "df": df,
                "p_value": p,
                "zone_proportion": zp,
                "baseline_proportion": bp,
                "enrichment_ratio": enr,
            }
        )
    elev_chi2, elev_p, elev_df, elev_zp, elev_bp, elev_enr = chi_square_2x2(cond_elev244, zone)
    elev_rows = [
        {
            "test_id": "JAS4_15",
            "chi_square": elev_chi2,
            "df": elev_df,
            "p_value": elev_p,
            "zone_proportion": elev_zp,
            "baseline_proportion": elev_bp,
            "enrichment_ratio": elev_enr,
        }
    ]

    # Distance mean-difference tests
    zone_road = road_dist[zone]
    base_road = road_dist[~zone]
    md_road_m, p_road = two_sample_mean_diff_perm(zone_road, base_road, N_PERM_DIST, rng)
    zone_water = water_dist[zone]
    base_water = water_dist[~zone]
    md_water_m, p_water = two_sample_mean_diff_perm(zone_water, base_water, N_PERM_DIST, rng)

    roman_rows = [
        {
            "mean_zone_km": float(np.nanmean(zone_road) / 1000.0),
            "mean_baseline_km": float(np.nanmean(base_road) / 1000.0),
            "mean_diff": float(md_road_m / 1000.0),
            "p_value": p_road,
        }
    ]
    water_rows = [
        {
            "mean_zone_km": float(np.nanmean(zone_water) / 1000.0),
            "mean_baseline_km": float(np.nanmean(base_water) / 1000.0),
            "mean_diff": float(md_water_m / 1000.0),
            "p_value": p_water,
        }
    ]

    # Infrastructure threshold tables
    infra_rows: List[Dict[str, object]] = []
    for infra_name, conds in [
        ("Roman Roads", [("lt_2km", cond_road2), ("lt_5km", cond_road5), ("lt_10km", cond_road10)]),
        ("Major Waterways", [("lt_1km", cond_water1), ("lt_2km", cond_water2), ("lt_5km", cond_water5)]),
    ]:
        for label, cond in conds:
            chi2, p, df, zp, bp, enr = chi_square_2x2(cond, zone)
            infra_rows.append(
                {
                    "infrastructure_type": infra_name,
                    "distance_threshold": label,
                    "zone_percent_within": zp * 100.0,
                    "landscape_percent_within": bp * 100.0,
                    "enrichment_ratio": enr,
                    "chi2_statistic": chi2,
                    "p_value": p,
                    "df": df,
                }
            )

    # Spatial permutation (toroidal shifts of zone mask)
    zone_full_bool = np.zeros((grid.ysize, grid.xsize), dtype=bool)
    zone_full_bool[grid.rows, grid.cols] = zone
    obs_zone_n = int(np.sum(zone))

    cond_map = {
        "Soil Grade 1 enrichment": cond_soil1,
        "Soil Grades 1-3 enrichment": cond_soil13,
        "Elevation <244 m enrichment": cond_elev244,
        "Roman Road <5 km enrichment": cond_road5,
        "Waterway <2 km enrichment": cond_water2,
    }
    cond_full_map: Dict[str, np.ndarray] = {}
    for k, cond in cond_map.items():
        cf = np.zeros((grid.ysize, grid.xsize), dtype=bool)
        cf[grid.rows, grid.cols] = cond
        cond_full_map[k] = cf

    perm_rows: List[Dict[str, object]] = []
    for name, cond in cond_map.items():
        obs_stat = float(np.mean(cond[zone]))
        null = np.empty(N_PERM_ENV, dtype=np.float64)
        cond_full = cond_full_map[name]
        for i in range(N_PERM_ENV):
            sx = int(rng.integers(0, grid.xsize))
            sy = int(rng.integers(0, grid.ysize))
            shifted = np.roll(np.roll(zone_full_bool, sy, axis=0), sx, axis=1) & grid.mask_full
            cur_n = int(np.sum(shifted))
            if cur_n > obs_zone_n:
                idx = np.where(shifted.ravel())[0]
                drop = rng.choice(idx, size=(cur_n - obs_zone_n), replace=False)
                flat = shifted.ravel()
                flat[drop] = False
                shifted = flat.reshape(shifted.shape)
            elif cur_n < obs_zone_n:
                idx = np.where((~shifted & grid.mask_full).ravel())[0]
                add = rng.choice(idx, size=(obs_zone_n - cur_n), replace=False)
                flat = shifted.ravel()
                flat[add] = True
                shifted = flat.reshape(shifted.shape)
            m = shifted[grid.rows, grid.cols]
            null[i] = float(np.mean(cond[m]))
        p_perm = (float(np.sum(null >= obs_stat)) + 1.0) / (N_PERM_ENV + 1.0)
        perm_rows.append(
            {
                "variable": name,
                "observed_statistic": obs_stat,
                "null_mean": float(np.mean(null)),
                "null_sd": float(np.std(null, ddof=1)),
                "p_perm": p_perm,
                "replicates": N_PERM_ENV,
                "interpretation": "SIGNIFICANT" if p_perm < 0.05 else "NOT_SIGNIFICANT",
            }
        )

    # Environmental correlation matrix + FDR
    env_tests = []
    env_vars = [
        ("Soil Grade 1", cond_soil1),
        ("Soil Grades 1-3", cond_soil13),
        ("Elevation <244 m", cond_elev244),
        ("Roman Road <5 km", cond_road5),
        ("Roman Road <10 km", cond_road10),
        ("Waterway <2 km", cond_water2),
        ("Waterway <5 km", cond_water5),
        ("Estuarine <5 km", cond_water5),  # proxy: waterway proximity (no explicit tidal limit layer in package)
    ]
    env_rows: List[Dict[str, object]] = []
    for name, cond in env_vars:
        chi2, p, df, zp, bp, enr = chi_square_2x2(cond, zone)
        # Pearson r on binary condition vs standardized density
        r = pearson_r(cond.astype(np.float64), zscore(surf_main))
        env_rows.append(
            {
                "variable": name,
                "zone_proportion_pct": zp * 100.0,
                "baseline_proportion_pct": bp * 100.0,
                "enrichment_ratio": enr,
                "pearson_r": r,
                "chi_square": chi2,
                "p_unadjusted": p,
            }
        )
        env_tests.append(p)
    env_p = np.array(env_tests, dtype=np.float64)
    env_q = bh_fdr(env_p)
    for i in range(len(env_rows)):
        env_rows[i]["p_fdr"] = float(env_q[i])
        env_rows[i]["fdr_significant"] = "YES" if env_q[i] < 0.05 else "NO"

    fdr_rows: List[Dict[str, object]] = []
    for r in env_rows:
        fdr_rows.append(
            {
                "test_id": r["variable"],
                "p_unadjusted": r["p_unadjusted"],
                "p_fdr": r["p_fdr"],
                "significant_after_fdr": r["fdr_significant"],
            }
        )
    n_sig = int(np.sum(env_q < 0.05))
    fdr_rows.append(
        {
            "test_id": "JAS4_19",
            "p_unadjusted": "",
            "p_fdr": "",
            "significant_after_fdr": "ALL_SIGNIFICANT" if n_sig == len(env_rows) else f"{n_sig}/{len(env_rows)}_SIGNIFICANT",
        }
    )

    # Regional stratification
    reg_rows: List[Dict[str, object]] = []
    regions = [("Yorkshire Lowlands", yk_region), ("Pennine Uplands", pn_region), ("Northwest Coastal/Estuarine", nw_region)]
    for reg_name, reg_mask in regions:
        m = reg_mask & np.isfinite(surf_main)
        if not np.any(m):
            continue
        zone_reg = zone & m
        denom = np.sum(m)
        reg_rows.append(
            {
                "region": reg_name,
                "soil_grade1_enrichment": (
                    (np.mean(cond_soil1[zone_reg]) / np.mean(cond_soil1[m]))
                    if np.sum(zone_reg) > 0 and np.mean(cond_soil1[m]) > 0
                    else float("nan")
                ),
                "elevation_lt244_pct": float(np.mean(cond_elev244[zone_reg]) * 100.0) if np.sum(zone_reg) > 0 else float("nan"),
                "roman_road_lt5km_pct": float(np.mean(cond_road5[zone_reg]) * 100.0) if np.sum(zone_reg) > 0 else float("nan"),
                "waterway_lt2km_pct": float(np.mean(cond_water2[zone_reg]) * 100.0) if np.sum(zone_reg) > 0 else float("nan"),
                "zone_cells_in_region": int(np.sum(zone_reg)),
                "total_cells_in_region": int(denom),
            }
        )

    # ------------------------------------------------------------------
    # Inter-rater agreement (compute from local ratings if available)
    # ------------------------------------------------------------------
    print("Stage 9/9: writing outputs")
    sici_validation_path = INBOX / "SICI_validation_sample.tsv"
    inter_rows: List[Dict[str, object]] = []
    conf_rows: List[Dict[str, object]] = []
    inter_source_note = "thesis_locked_value_due_empty_local_rater_columns"
    if sici_validation_path.exists():
        r1, r2 = read_sici_validation_ratings(sici_validation_path)
        if r1.size >= 30 and r2.size == r1.size:
            kappa_val = quadratic_weighted_kappa(r1, r2, n_cat=3)
            ac1_val = gwet_ac1_nominal(r1, r2, n_cat=3)
            (k_l, k_u), (a_l, a_u) = bootstrap_reliability_ci(r1, r2, N_BOOTSTRAP, rng)
            inter_rows = [
                {
                    "metric": "quadratic_weighted_cohens_kappa",
                    "value": kappa_val,
                    "ci_lower": k_l,
                    "ci_upper": k_u,
                    "n_validation_sample": int(r1.size),
                    "source": "computed_from_local_rater_columns",
                },
                {
                    "metric": "gwets_ac1",
                    "value": ac1_val,
                    "ci_lower": a_l,
                    "ci_upper": a_u,
                    "n_validation_sample": int(r1.size),
                    "source": "computed_from_local_rater_columns",
                },
            ]
            inter_source_note = "computed_from_local_rater_columns"
            conf = np.zeros((3, 3), dtype=np.int64)
            for a, b in zip(r1, r2):
                if 0 <= a < 3 and 0 <= b < 3:
                    conf[a, b] += 1
            for i in range(3):
                for j in range(3):
                    conf_rows.append(
                        {
                            "rater1_grade": i,
                            "rater2_grade": j,
                            "count": int(conf[i, j]),
                        }
                    )

    if not inter_rows:
        inter_rows = [
            {
                "metric": "quadratic_weighted_cohens_kappa",
                "value": 0.78,
                "ci_lower": 0.71,
                "ci_upper": 0.85,
                "n_validation_sample": 321,
                "source": "thesis_locked_value_due_empty_local_rater_columns",
            },
            {
                "metric": "gwets_ac1",
                "value": 0.81,
                "ci_lower": "",
                "ci_upper": "",
                "n_validation_sample": 321,
                "source": "thesis_locked_value_due_empty_local_rater_columns",
            },
        ]

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    write_tsv(
        OUT_TABLES / "validation_metrics.tsv",
        [
            {
                "n_validation": n_val,
                "n_hits_top15": n_hits,
                "hit_rate": hit_rate,
                "hit_rate_pct": hit_rate * 100.0,
                "ci_lower": ci_l,
                "ci_upper": ci_u,
                "ci_lower_pct": ci_l * 100.0,
                "ci_upper_pct": ci_u * 100.0,
                "top_area_pct": VALIDATION_TOP_AREA_PCT,
                "target_year_from": TARGET_YEAR_FROM,
                "target_year_to": TARGET_YEAR_TO,
                "holdout_temporal_mode": HOLDOUT_TEMPORAL_MODE,
                "hit_rate_temporal_weighted": hit_rate_temporal_weighted,
                "n_validation_temporal_weighted": temporal_weight_den,
                "binomial_p_one_sided_vs_csr": hit_binom_p,
                "posterior_hit_rate_ci_lower": hit_post_ci_l,
                "posterior_hit_rate_ci_upper": hit_post_ci_u,
                "posterior_prob_hit_rate_gt_csr": hit_post_prob_gt_csr,
            }
        ],
        [
            "n_validation",
            "n_hits_top15",
            "hit_rate",
            "hit_rate_pct",
            "ci_lower",
            "ci_upper",
            "ci_lower_pct",
            "ci_upper_pct",
            "top_area_pct",
            "target_year_from",
            "target_year_to",
            "holdout_temporal_mode",
            "hit_rate_temporal_weighted",
            "n_validation_temporal_weighted",
            "binomial_p_one_sided_vs_csr",
            "posterior_hit_rate_ci_lower",
            "posterior_hit_rate_ci_upper",
            "posterior_prob_hit_rate_gt_csr",
        ],
    )
    write_tsv(
        OUT_TABLES / "auc_metrics.tsv",
        [
            {
                "auc_mean": auc,
                "auc_temporal_weighted": auc_temporal_weighted,
                "auc_sd": float(np.std(fold_auc, ddof=1)) if fold_auc.size > 1 else 0.0,
                "ci_lower": auc_ci_l,
                "ci_upper": auc_ci_u,
                "negative_samples": int(min(N_NEG_SAMPLES, mask_n)),
            }
        ],
        ["auc_mean", "auc_temporal_weighted", "auc_sd", "ci_lower", "ci_upper", "negative_samples"],
    )
    write_tsv(
        OUT_TABLES / "auc_opportunity_null_metrics.tsv",
        [
            {
                "auc_opportunity_null": auc_oppnull,
                "auc_opportunity_null_temporal_weighted": auc_oppnull_temporal,
                "ci_lower": auc_oppnull_ci_l,
                "ci_upper": auc_oppnull_ci_u,
                "negative_samples": int(min(N_NEG_SAMPLES, mask_n)),
                "sampling_mode": opp_sampling_note,
            }
        ],
        [
            "auc_opportunity_null",
            "auc_opportunity_null_temporal_weighted",
            "ci_lower",
            "ci_upper",
            "negative_samples",
            "sampling_mode",
        ],
    )
    write_tsv(
        OUT_TABLES / "baseline_gain_metrics.tsv",
        [
            {
                "gain_vs_csr": gain_vs_csr,
                "gain_vs_opportunity": gain_vs_opp,
                "opportunity_baseline_hit_rate": opp_hit_rate,
                "opportunity_note": (
                    "Observed detector-effort baseline from independent baseline PAS (grade 0 / mid-Saxon) KDE surface."
                    if opp_surface_mode == "observed_detector_effort_baseline_pas_kde"
                    else
                    "Observed detector-effort baseline from full PAS all-finds/all-periods KDE surface."
                    if opp_surface_mode == "observed_detector_effort_pas_all_kde"
                    else
                    "Observed detector-effort constraints surface from delivered opportunity layer "
                    "(woodland/urban attenuation)."
                    if opp_surface_mode == "observed_detector_effort_constraints_layer"
                    else "Modeled continuous opportunity surface from arable+slope+road covariates "
                    "(fallback when observed constraints surface unavailable/unusable)."
                ),
            }
        ],
        ["gain_vs_csr", "gain_vs_opportunity", "opportunity_baseline_hit_rate", "opportunity_note"],
    )
    write_tsv(
        OUT_TABLES / "opportunity_model_diagnostics.tsv",
        [
            {
                "opportunity_threshold_top_area": opp_thr,
                "opportunity_top_area_pct": VALIDATION_TOP_AREA_PCT,
                "opportunity_sampling_mode": opp_sampling_note,
                **opp_diag,
            }
        ],
        [
            "opportunity_threshold_top_area",
            "opportunity_top_area_pct",
            "opportunity_sampling_mode",
            "opportunity_min",
            "opportunity_max",
            "opportunity_mean",
            "arable_mean",
            "slope_score_mean",
            "road_score_mean",
            "constraints_urban_pct",
            "constraints_woodland_pct",
            "constraints_national_park_pct",
            "baseline_pas_points",
            "pas_all_points",
        ],
    )
    write_tsv(
        OUT_TABLES / "decile_lift_table.tsv",
        decile_rows,
        [
            "top_area_pct",
            "n_validation",
            "n_captured",
            "capture_rate",
            "expected_random_rate",
            "lift_vs_random",
            "opportunity_capture_rate",
            "lift_vs_opportunity",
        ],
    )
    write_tsv(
        OUT_TABLES / "decile_bin_reliability.tsv",
        decile_bin_rows,
        [
            "score_decile",
            "score_lower",
            "score_upper",
            "n_validation",
            "n_in_bin",
            "observed_rate",
            "expected_rate",
            "lift_vs_uniform",
        ],
    )
    write_tsv(
        OUT_TABLES / "holdout_uncertainty_diagnostics.tsv",
        [
            {
                "n_validation_primary": n_val,
                "n_hits_top15_primary": n_hits,
                "hit_rate_primary": hit_rate,
                "binomial_p_one_sided_vs_csr": hit_binom_p,
                "posterior_hit_rate_ci_lower": hit_post_ci_l,
                "posterior_hit_rate_ci_upper": hit_post_ci_u,
                "posterior_prob_hit_rate_gt_csr": hit_post_prob_gt_csr,
                "max_abs_delta_hit_rate_pct_points_loo": loo_summary["max_abs_delta_hit_rate_pct_points"],
                "max_abs_delta_auc_loo": loo_summary["max_abs_delta_auc"],
                "mean_abs_delta_hit_rate_pct_points_loo": loo_summary["mean_abs_delta_hit_rate_pct_points"],
                "mean_abs_delta_auc_loo": loo_summary["mean_abs_delta_auc"],
            }
        ],
        [
            "n_validation_primary",
            "n_hits_top15_primary",
            "hit_rate_primary",
            "binomial_p_one_sided_vs_csr",
            "posterior_hit_rate_ci_lower",
            "posterior_hit_rate_ci_upper",
            "posterior_prob_hit_rate_gt_csr",
            "max_abs_delta_hit_rate_pct_points_loo",
            "max_abs_delta_auc_loo",
            "mean_abs_delta_hit_rate_pct_points_loo",
            "mean_abs_delta_auc_loo",
        ],
    )
    write_tsv(
        OUT_TABLES / "holdout_influence_loo.tsv",
        loo_rows,
        [
            "dropped_site_id",
            "n_validation_after_drop",
            "hit_rate_after_drop",
            "delta_hit_rate_pct_points",
            "auc_after_drop",
            "delta_auc",
        ],
    )
    write_tsv(
        OUT_TABLES / "unweighted_comparison_metrics.tsv",
        [
            {
                "hit_rate_unweighted": hit_unw,
                "auc_unweighted": auc_unw,
                "hit_rate_weighted": hit_rate,
                "auc_weighted": auc,
                "delta_hit_rate_pct_points": (hit_rate - hit_unw) * 100.0,
                "delta_auc": auc - auc_unw,
            }
        ],
        [
            "hit_rate_unweighted",
            "auc_unweighted",
            "hit_rate_weighted",
            "auc_weighted",
            "delta_hit_rate_pct_points",
            "delta_auc",
        ],
    )
    write_tsv(
        OUT_TABLES / "fold_assignments.tsv",
        [
            {
                "site_id": hold_id[i],
                "fold": int(folds[i]),
                "x": float(hold_x[i]),
                "y": float(hold_y[i]),
            }
            for i in range(hold_id.size)
        ],
        ["site_id", "fold", "x", "y"],
    )
    write_tsv(
        OUT_TABLES / "fold_metrics.tsv",
        fold_rows,
        [
            "fold",
            "validation_sites",
            "training_sites",
            "excluded_sites",
            "hit_rate_15pct",
            "auc",
            "baseline_gain_csr",
            "baseline_gain_opportunity",
        ],
    )
    write_tsv(
        OUT_TABLES / "exclusion_counts.tsv",
        excl_rows,
        ["fold", "original_training_sites", "sites_excluded_by_buffer", "final_training_sites", "exclusion_rate_pct"],
    )
    # Temporal holdout scenario audit for the target window.
    temporal_rows: List[Dict[str, object]] = []

    def add_temporal_row(label: str, mask: np.ndarray, weighted: bool = False) -> None:
        if mask.size == 0 or not np.any(mask):
            return
        sc = hold_scores[mask]
        hs = (sc >= threshold) & np.isfinite(sc)
        n_local = int(np.sum(np.isfinite(sc)))
        n_hits_local = int(np.sum(hs))
        hit_local = (n_hits_local / n_local) if n_local > 0 else float("nan")
        ci_l_local, ci_u_local = wilson_ci(n_hits_local, n_local)
        auc_local = auc_from_scores(sc, neg_scores)

        if weighted:
            w_local = hold_temporal_w[mask]
            denom_local = float(np.sum(w_local))
            hit_w = float(np.sum(w_local * hs.astype(np.float64)) / denom_local) if denom_local > 0 else float("nan")
            auc_w = weighted_auc_from_scores(sc, neg_scores, w_local)
            n_eff = effective_sample_size(w_local)
            binom_p = float("nan")
            post_ci_l = float("nan")
            post_ci_u = float("nan")
            post_prob = float("nan")
        else:
            hit_w = float("nan")
            auc_w = float("nan")
            denom_local = float("nan")
            n_eff = float("nan")
            binom_p = binomial_one_sided_p_ge(n_hits_local, n_local, VALIDATION_TOP_AREA_PCT / 100.0)
            post_ci_l, post_ci_u, post_prob = beta_posterior_summary(
                n_hits_local, n_local, VALIDATION_TOP_AREA_PCT / 100.0, rng
            )

        st = np.array([str(holdouts[i].get("site_type", "")).lower() for i in np.where(mask)[0]], dtype=object)
        n_burial = int(np.sum(st == "burial"))
        n_settlement = int(np.sum(st == "settlement"))
        n_hoard = int(np.sum(st == "hoard"))
        temporal_rows.append(
            {
                "scenario": label,
                "target_year_from": TARGET_YEAR_FROM,
                "target_year_to": TARGET_YEAR_TO,
                "n_validation": n_local,
                "n_hits_top15": n_hits_local,
                "hit_rate": hit_local,
                "hit_ci_lower": ci_l_local,
                "hit_ci_upper": ci_u_local,
                "auc": auc_local,
                "hit_rate_temporal_weighted": hit_w,
                "auc_temporal_weighted": auc_w,
                "n_validation_temporal_weighted": denom_local,
                "n_validation_effective_temporal": n_eff,
                "binomial_p_one_sided_vs_csr": binom_p,
                "posterior_hit_rate_ci_lower": post_ci_l,
                "posterior_hit_rate_ci_upper": post_ci_u,
                "posterior_prob_hit_rate_gt_csr": post_prob,
                "n_burial": n_burial,
                "n_settlement": n_settlement,
                "n_hoard": n_hoard,
            }
        )

    inside_mask = np.array(
        [
            temporal_pass(int(hold_from[i]), int(hold_to[i]), TARGET_YEAR_FROM, TARGET_YEAR_TO, "inside")
            for i in range(hold_from.size)
        ],
        dtype=bool,
    )
    overlap_mask = np.array(
        [
            temporal_pass(int(hold_from[i]), int(hold_to[i]), TARGET_YEAR_FROM, TARGET_YEAR_TO, "overlap")
            for i in range(hold_from.size)
        ],
        dtype=bool,
    )
    add_temporal_row("strict_inside_unweighted", inside_mask, weighted=False)
    add_temporal_row("overlap_unweighted", overlap_mask, weighted=False)
    add_temporal_row("overlap_temporal_weighted", overlap_mask, weighted=True)
    write_tsv(
        OUT_TABLES / "temporal_holdout_scenarios.tsv",
        temporal_rows,
        [
            "scenario",
            "target_year_from",
            "target_year_to",
            "n_validation",
            "n_hits_top15",
            "hit_rate",
            "hit_ci_lower",
            "hit_ci_upper",
            "auc",
            "hit_rate_temporal_weighted",
            "auc_temporal_weighted",
            "n_validation_temporal_weighted",
            "n_validation_effective_temporal",
            "binomial_p_one_sided_vs_csr",
            "posterior_hit_rate_ci_lower",
            "posterior_hit_rate_ci_upper",
            "posterior_prob_hit_rate_gt_csr",
            "n_burial",
            "n_settlement",
            "n_hoard",
        ],
    )
    with (OUT_TABLES / "leakage_check_results.txt").open("w", encoding="utf-8") as f:
        f.write("Spatial leakage check (buffer radius = 10,000 m)\n")
        f.write("\n".join(leakage_lines))
        f.write("\n")
    write_tsv(
        OUT_TABLES / "crossvalidation_summary.tsv",
        cv_summary_rows,
        ["metric", "mean", "sd", "se", "ci_lower", "ci_upper"],
    )
    write_tsv(
        OUT_TABLES / "region_holdout_transportability.tsv",
        region_transport_rows,
        [
            "scenario",
            "train_region",
            "test_region",
            "n_train_pas",
            "n_test_holdouts",
            "n_hits_top15",
            "hit_rate",
            "auc",
            "gain_vs_csr",
            "gain_vs_opportunity",
            "status",
        ],
    )
    write_tsv(
        OUT_TABLES / "spatial_autocorrelation_moransI.tsv",
        [{"moran_i": moran_i, "p_perm": moran_p, "permutations": N_PERM_MORAN}],
        ["moran_i", "p_perm", "permutations"],
    )
    write_tsv(
        OUT_TABLES / "bandwidth_sensitivity_results.tsv",
        bw_rows,
        ["bandwidth_km", "surface_correlation_h10", "pattern_preserved"],
    )
    write_tsv(
        OUT_TABLES / "kernel_comparison_results.tsv",
        kernel_rows,
        ["kernel_function", "surface_correlation_with_gaussian", "regional_pattern_fidelity"],
    )
    write_tsv(
        OUT_TABLES / "coordinate_perturbation_results.tsv",
        coord_rows,
        [
            "perturbation_magnitude_m",
            "replicates",
            "surface_correlation_mean",
            "surface_correlation_sd",
            "high_conf_zone_overlap_pct_mean",
            "high_conf_zone_overlap_pct_sd",
        ],
    )
    write_tsv(
        OUT_TABLES / "grade_perturbation_results.tsv",
        grade_rows,
        [
            "scenario",
            "replicates",
            "surface_correlation_mean",
            "surface_correlation_sd",
            "high_conf_zone_overlap_pct_mean",
            "high_conf_zone_overlap_pct_sd",
        ],
    )
    write_tsv(
        OUT_TABLES / "jackknife_results.tsv",
        jack_rows,
        ["removed_proxy", "surface_correlation_with_full", "regional_pattern_preservation"],
    )
    write_tsv(
        OUT_TABLES / "soil_chisquare_stats.tsv",
        chi_rows,
        ["test_id", "chi_square", "df", "p_value", "zone_proportion", "baseline_proportion", "enrichment_ratio"],
    )
    write_tsv(
        OUT_TABLES / "elevation_chisquare_stats.tsv",
        elev_rows,
        ["test_id", "chi_square", "df", "p_value", "zone_proportion", "baseline_proportion", "enrichment_ratio"],
    )
    write_tsv(
        OUT_TABLES / "roman_roads_proximity.tsv",
        roman_rows,
        ["mean_zone_km", "mean_baseline_km", "mean_diff", "p_value"],
    )
    write_tsv(
        OUT_TABLES / "waterways_proximity.tsv",
        water_rows,
        ["mean_zone_km", "mean_baseline_km", "mean_diff", "p_value"],
    )
    write_tsv(
        OUT_TABLES / "infrastructure_proximity_results.tsv",
        infra_rows,
        [
            "infrastructure_type",
            "distance_threshold",
            "zone_percent_within",
            "landscape_percent_within",
            "enrichment_ratio",
            "chi2_statistic",
            "p_value",
            "df",
        ],
    )
    write_tsv(
        OUT_TABLES / "permutation_tests_extended.tsv",
        perm_rows,
        ["variable", "observed_statistic", "null_mean", "null_sd", "p_perm", "replicates", "interpretation"],
    )
    # One-line permutation output expected by Supplement S2 placeholder
    elev_perm = [r for r in perm_rows if r["variable"] == "Elevation <244 m enrichment"]
    write_tsv(
        OUT_TABLES / "permutation_test_results.tsv",
        [{"p_perm": elev_perm[0]["p_perm"] if elev_perm else float("nan"), "replicates": N_PERM_ENV}],
        ["p_perm", "replicates"],
    )
    write_tsv(
        OUT_TABLES / "environmental_correlation_matrix.tsv",
        env_rows,
        [
            "variable",
            "zone_proportion_pct",
            "baseline_proportion_pct",
            "enrichment_ratio",
            "pearson_r",
            "chi_square",
            "p_unadjusted",
            "p_fdr",
            "fdr_significant",
        ],
    )
    write_tsv(
        OUT_TABLES / "fdr_adjusted_pvalues.tsv",
        fdr_rows,
        ["test_id", "p_unadjusted", "p_fdr", "significant_after_fdr"],
    )
    write_tsv(
        OUT_TABLES / "regional_stratification.tsv",
        reg_rows,
        [
            "region",
            "soil_grade1_enrichment",
            "elevation_lt244_pct",
            "roman_road_lt5km_pct",
            "waterway_lt2km_pct",
            "zone_cells_in_region",
            "total_cells_in_region",
        ],
    )
    write_tsv(
        OUT_TABLES / "interrater_kappa_ac1.tsv",
        inter_rows,
        ["metric", "value", "ci_lower", "ci_upper", "n_validation_sample", "source"],
    )
    if conf_rows:
        write_tsv(
            OUT_TABLES / "interrater_confusion_matrix.tsv",
            conf_rows,
            ["rater1_grade", "rater2_grade", "count"],
        )

    # ROC curve data
    all_scores = np.concatenate([hold_scores[np.isfinite(hold_scores)], neg_scores[np.isfinite(neg_scores)]])
    labels = np.concatenate(
        [
            np.ones(np.sum(np.isfinite(hold_scores)), dtype=np.int8),
            np.zeros(np.sum(np.isfinite(neg_scores)), dtype=np.int8),
        ]
    )
    thr = np.unique(np.quantile(all_scores, np.linspace(0, 1, 201)))
    roc_rows: List[Dict[str, object]] = []
    for t in thr:
        pred = all_scores >= t
        tp = int(np.sum(pred & (labels == 1)))
        fp = int(np.sum(pred & (labels == 0)))
        fn = int(np.sum((~pred) & (labels == 1)))
        tn = int(np.sum((~pred) & (labels == 0)))
        tpr = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
        fpr = fp / (fp + tn) if (fp + tn) > 0 else float("nan")
        roc_rows.append({"threshold": float(t), "tpr": tpr, "fpr": fpr, "tp": tp, "fp": fp, "fn": fn, "tn": tn})
    write_tsv(OUT_TABLES / "roc_curve_data.tsv", roc_rows, ["threshold", "tpr", "fpr", "tp", "fp", "fn", "tn"])

    # Dataset summary (for Table 1)
    pas_note = (
        f"Current analytical package; temporal_mode={PAS_TEMPORAL_MODE}; "
        f"temporal_weighting={PAS_TEMPORAL_WEIGHTING}; "
        f"target_window={TARGET_YEAR_FROM}-{TARGET_YEAR_TO}; "
        f"pre_temporal_n={pas_pre_temporal_n}; post_temporal_n={pas_post_temporal_n}; "
        f"temporal_weight_mean={pas_temporal_weight_mean:.4f}; "
        f"temporal_weight_mean_g1g2={pas_temporal_weight_mean_g12:.4f}; "
        f"temporal_weight_effective_n={pas_temporal_weight_effective_n:.2f}"
    )
    holdout_note = (
        f"Validation holdouts loaded from {holdout_mode}; "
        f"target_window={TARGET_YEAR_FROM}-{TARGET_YEAR_TO}; temporal_mode={HOLDOUT_TEMPORAL_MODE}; "
        f"pre_temporal_n={holdouts_pre_temporal}; post_temporal_n={len(holdouts)}; "
        f"type_counts={holdout_type_counts}; "
        f"include_sensitivity={INCLUDE_SENSITIVITY_HOLDOUTS}; "
        f"include_context={INCLUDE_CONTEXT_HOLDOUTS}"
    )
    table1_rows = [
        {"dataset": "PAS_full_northumbria_metalwork", "n": len(pas), "notes": pas_note},
        {"dataset": "KEPN_full_northumbria", "n": len(kepn), "notes": "Full KEPN export"},
        {"dataset": "CASSS_full_northumbria", "n": len(casss), "notes": "Full CASSS export"},
        {"dataset": "Validation_holdouts", "n": len(holdouts), "notes": holdout_note},
        {"dataset": "Grid_mask_cells_1km", "n": mask_n, "notes": "Analysis mask raster cells"},
        {"dataset": "High_confidence_zone_top15pct_cells", "n": zone_count, "notes": "Top 15% KDE area"},
    ]
    write_tsv(OUT_TABLES / "dataset_summary_computed.tsv", table1_rows, ["dataset", "n", "notes"])
    write_tsv(
        OUT_TABLES / "pas_temporal_weight_diagnostics.tsv",
        [
            {
                "pas_temporal_weighting": PAS_TEMPORAL_WEIGHTING,
                "pas_temporal_mode": PAS_TEMPORAL_MODE,
                "target_year_from": TARGET_YEAR_FROM,
                "target_year_to": TARGET_YEAR_TO,
                "n_pas_points": len(pas),
                "temporal_weight_mean_all": pas_temporal_weight_mean,
                "temporal_weight_mean_g1g2": pas_temporal_weight_mean_g12,
                "temporal_weight_effective_n_all": pas_temporal_weight_effective_n,
                "temporal_weight_min": float(np.min(pas_temporal_w)) if pas_temporal_w.size > 0 else float("nan"),
                "temporal_weight_p25": float(np.percentile(pas_temporal_w, 25.0))
                if pas_temporal_w.size > 0
                else float("nan"),
                "temporal_weight_median": float(np.median(pas_temporal_w)) if pas_temporal_w.size > 0 else float("nan"),
                "temporal_weight_p75": float(np.percentile(pas_temporal_w, 75.0))
                if pas_temporal_w.size > 0
                else float("nan"),
                "temporal_weight_max": float(np.max(pas_temporal_w)) if pas_temporal_w.size > 0 else float("nan"),
            }
        ],
        [
            "pas_temporal_weighting",
            "pas_temporal_mode",
            "target_year_from",
            "target_year_to",
            "n_pas_points",
            "temporal_weight_mean_all",
            "temporal_weight_mean_g1g2",
            "temporal_weight_effective_n_all",
            "temporal_weight_min",
            "temporal_weight_p25",
            "temporal_weight_median",
            "temporal_weight_p75",
            "temporal_weight_max",
        ],
    )

    # Provenance manifest
    manifest = {
        "run_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED,
        "grid_cell_size_m": CELL_SIZE_M,
        "kde_bandwidth_m": KDE_BANDWIDTH_M,
        "top_area_pct": VALIDATION_TOP_AREA_PCT,
        "bootstrap_replicates": N_BOOTSTRAP,
        "moran_permutations": N_PERM_MORAN,
        "environmental_permutations": N_PERM_ENV,
        "distance_permutations": N_PERM_DIST,
        "sensitivity_replicates": N_SENS_REPS,
        "inputs": {
            "pas_path": str(pas_path),
            "pas_sha256": sha256(pas_path),
            "kepn_path": str(kepn_path),
            "kepn_sha256": sha256(kepn_path),
            "casss_path": str(casss_path),
            "casss_sha256": sha256(casss_path),
            "boundary_path": str(boundary_path),
            "boundary_sha256": sha256(boundary_path),
            "roads_path": str(roads_path),
            "roads_sha256": sha256(roads_path),
            "waterways_path": str(waterways_path),
            "waterways_sha256": sha256(waterways_path),
            "soil_path": str(soil_path),
            "soil_sha256": sha256(soil_path),
            "dem_path": str(dem_path),
            "dem_sha256": sha256(dem_path),
            "arable_path": str(arable_path),
            "arable_sha256": sha256(arable_path),
            "pas_all_path": str(pas_all_path) if pas_all_path.exists() else "",
            "pas_all_sha256": sha256(pas_all_path) if pas_all_path.exists() else "",
            "pas_all_layer": pas_all_layer if pas_all_layer else "",
            "opportunity_model": opp_surface_mode,
            "opportunity_constraints_path": str(opp_constraints_path) if opp_constraints_path.exists() else "",
            "opportunity_constraints_sha256": sha256(opp_constraints_path) if opp_constraints_path.exists() else "",
            "opportunity_sampling_mode": opp_sampling_note,
            "holdout_path_used": str(holdout_path),
            "holdout_sha256": sha256(holdout_path),
            "holdout_mode": holdout_mode,
            "holdout_type_counts": holdout_type_counts,
            "include_sensitivity_holdouts": INCLUDE_SENSITIVITY_HOLDOUTS,
            "include_context_holdouts": INCLUDE_CONTEXT_HOLDOUTS,
            "target_year_from": TARGET_YEAR_FROM,
            "target_year_to": TARGET_YEAR_TO,
            "pas_temporal_mode": PAS_TEMPORAL_MODE,
            "pas_temporal_weighting": PAS_TEMPORAL_WEIGHTING,
            "pas_temporal_weight_mean_all": pas_temporal_weight_mean,
            "pas_temporal_weight_mean_g1g2": pas_temporal_weight_mean_g12,
            "pas_temporal_weight_effective_n_all": pas_temporal_weight_effective_n,
            "holdout_temporal_mode": HOLDOUT_TEMPORAL_MODE,
            "pas_pre_temporal_n": pas_pre_temporal_n,
            "pas_post_temporal_n": pas_post_temporal_n,
            "holdout_pre_temporal_n": holdouts_pre_temporal,
            "holdout_post_temporal_n": len(holdouts),
        },
        "limitations": [],
    }
    if holdout_mode == "vetted_gpkg":
        manifest["limitations"].append(
            "Independent validation uses provided holdout layer with role filtering "
            "(primary only by default; sensitivity/context excluded unless explicitly enabled)."
        )
    else:
        manifest["limitations"].append(
            "Independent validation uses burial TSV fallback because validation_holdout_sites.gpkg was unavailable."
        )
    if opp_surface_mode == "observed_detector_effort_baseline_pas_kde":
        manifest["limitations"].append(
            "Opportunity baseline uses observed detector-effort inferred from independent baseline PAS (grade 0/mid-Saxon) KDE intensity; "
            "this is empirical reporting intensity, not direct detector-count telemetry."
        )
    elif opp_surface_mode == "observed_detector_effort_pas_all_kde":
        manifest["limitations"].append(
            "Opportunity baseline uses observed detector-effort inferred from full PAS all-finds/all-periods KDE intensity; "
            "this is empirical reporting intensity, not direct detector-count telemetry."
        )
    elif opp_surface_mode == "observed_detector_effort_constraints_layer":
        manifest["limitations"].append(
            "Opportunity baseline uses observed detector-effort constraints from the delivered opportunity layer; "
            "weights are attenuation-based (woodland/urban penalties), not direct detector counts."
        )
    else:
        manifest["limitations"].append(
            "Opportunity baseline uses a modeled continuous covariate surface (arable+slope+road); "
            "observed constraints surface was unavailable/unusable."
        )
    if inter_source_note == "thesis_locked_value_due_empty_local_rater_columns":
        manifest["limitations"].append(
            "Inter-rater metrics are thesis-locked due empty local rater columns in SICI_validation_sample.tsv."
        )
    manifest["limitations"].append(
        "Estuarine metric uses waterway proximity proxy (no tidal-limit layer in current package)."
    )
    if HOLDOUT_TEMPORAL_MODE == "overlap":
        manifest["limitations"].append(
            "Holdout temporal inclusion uses overlap with target window; strict-inside and overlap-weighted scenarios are reported separately."
        )
    if PAS_TEMPORAL_WEIGHTING == "overlap_fraction":
        manifest["limitations"].append(
            "PAS evidence weights are temporally attenuated by overlap fraction with the target window "
            "(effective_weight = SICI_weight x overlap_fraction)."
        )
    with (OUT_LOGS / "analysis_manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    with (OUT_LOGS / "RUN_COMPLETE.txt").open("w", encoding="utf-8") as f:
        f.write("JAS analysis pipeline complete.\n")
        f.write(f"UTC: {manifest['run_timestamp_utc']}\n")
        f.write(f"Outputs directory: {OUT}\n")

    print("Complete. Outputs written to:")
    print(f"  {OUT_TABLES}")
    print(f"  {OUT_RASTERS}")
    print(f"  {OUT_LOGS}")


if __name__ == "__main__":
    main()
