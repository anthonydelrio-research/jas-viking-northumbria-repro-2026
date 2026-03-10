#!/usr/bin/env python3
"""
Derive reproducible SICI high-confidence placename site registers for Northumbria.

Definition used here:
  - PAS evidence: Scandinavian_Metalwork only, SICI grade in {1, 2}
  - Temporal filter: configurable year window with overlap logic
  - Placename filter: KEPN entries with Scandinavian language evidence
  - Site criterion: placename point within 500 m of >=1 qualifying PAS object
  - Tertiary support: CASSS point proximity within 1000 m (descriptive only)

Outputs:
  - TSV + GPKG site register for each window
  - Northumberland/Durham focus subset TSV for each window
  - PAS hotspot trace table for provided IDs (nearest placename diagnostics)
  - JSON manifest + Markdown summary
"""

from __future__ import annotations

import csv
import json
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from osgeo import ogr, osr


ROOT = Path(os.environ.get("JAS_ROOT", str(Path(__file__).resolve().parents[1]))).resolve()
INBOX = Path(os.environ.get("JAS_INBOX", str(ROOT / "incoming_arcgis" / "01_inbox"))).resolve()

PAS_PATH = INBOX / "PAS_full_northumbria_metalwork.gpkg"
PAS_LAYER = "PAS_full_northumbria_metalwork"
KEPN_PATH = INBOX / "KEPN_full_northumbria.gpkg"
KEPN_LAYER = "KEPN_full_northumbria"
CASSS_PATH = INBOX / "CASSS_full_northumbria.gpkg"
CASSS_LAYER = "CASSS_full_northumbria"

TARGET_WINDOWS: Sequence[Tuple[int, int]] = [(750, 954), (850, 1050)]
TEMPORAL_MODE = "overlap"
PAS_RADIUS_M = float(os.environ.get("SICI_PAS_RADIUS_M", "500"))
CASSS_RADIUS_M = float(os.environ.get("SICI_CASSS_RADIUS_M", "1000"))
OUTPUT_DIR = Path(
    os.environ.get("SICI_OUTPUT_DIR", str(ROOT / "outputs_submission_2026-03-05_sites"))
)

SCANDI_LANG_CODES = {"N", "ND", "NW", "NN"}
SCANDI_LANG_LABELS = {"Old Norse", "Old Danish", "Old West Scandinavian", "Old Norwegian"}

HOTSPOT_PAS_IDS = [
    "288858",
    "418825",
    "431432",
    "433983",
    "443337",
    "496711",
    "614618",
    "614626",
    "738829",
]


@dataclass
class PasPoint:
    site_id: str
    x: float
    y: float
    sici_grade: int
    source_subset: str
    date_from_ad: Optional[int]
    date_to_ad: Optional[int]
    county_hist: str
    object_class: str


@dataclass
class PlacenamePoint:
    key: str
    placeno: str
    placename: str
    county: str
    x: float
    y: float
    langcodes: Tuple[str, ...]
    languages: Tuple[str, ...]
    headwords: Tuple[str, ...]
    is_scandinavian: bool


@dataclass
class CasssPoint:
    x: float
    y: float
    name: str
    category: str
    object_type: str


def parse_year(value: object) -> Optional[int]:
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    try:
        return int(float(s))
    except ValueError:
        return None


def temporal_overlap_years(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def temporal_pass(a0: Optional[int], a1: Optional[int], win0: int, win1: int, mode: str) -> bool:
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


def open_layer(path: Path, layer_name: str):
    ds = ogr.Open(str(path))
    if ds is None:
        raise RuntimeError(f"Failed to open {path}")
    lyr = ds.GetLayerByName(layer_name)
    if lyr is None:
        raise RuntimeError(f"Layer {layer_name} not found in {path}")
    return ds, lyr


def read_pas_points(path: Path, layer_name: str, year_from: int, year_to: int, temporal_mode: str) -> List[PasPoint]:
    ds, lyr = open_layer(path, layer_name)
    rows: List[PasPoint] = []
    for feat in lyr:
        source_subset = str(feat.GetField("source_subset") or "").strip()
        if source_subset != "Scandinavian_Metalwork":
            continue
        grade_raw = feat.GetField("sici_grade")
        try:
            sici_grade = int(grade_raw)
        except Exception:
            continue
        if sici_grade not in (1, 2):
            continue
        d0 = parse_year(feat.GetField("date_from_ad"))
        d1 = parse_year(feat.GetField("date_to_ad"))
        if not temporal_pass(d0, d1, year_from, year_to, temporal_mode):
            continue
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        rows.append(
            PasPoint(
                site_id=str(feat.GetField("site_id") or "").strip(),
                x=float(geom.GetX()),
                y=float(geom.GetY()),
                sici_grade=sici_grade,
                source_subset=source_subset,
                date_from_ad=d0,
                date_to_ad=d1,
                county_hist=str(feat.GetField("county_hist") or "").strip(),
                object_class=str(feat.GetField("object_class") or "").strip(),
            )
        )
    ds = None
    return rows


def read_kepn_grouped(path: Path, layer_name: str) -> List[PlacenamePoint]:
    ds, lyr = open_layer(path, layer_name)
    grouped: Dict[str, Dict[str, object]] = {}
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        placeno = str(feat.GetField("placeno") or "").strip()
        placename = str(feat.GetField("placename") or "").strip()
        county = str(feat.GetField("county") or "").strip()
        x = float(geom.GetX())
        y = float(geom.GetY())
        key = placeno if placeno else f"{placename}|{x:.3f}|{y:.3f}"
        if key not in grouped:
            grouped[key] = {
                "key": key,
                "placeno": placeno,
                "placename": placename,
                "county": county,
                "x": x,
                "y": y,
                "langcodes": set(),
                "languages": set(),
                "headwords": set(),
            }
        langcode = str(feat.GetField("elements__langcode") or "").strip()
        language = str(feat.GetField("elements__language") or "").strip()
        headword = str(feat.GetField("elements__headword") or "").strip()
        if langcode:
            grouped[key]["langcodes"].add(langcode)
        if language:
            grouped[key]["languages"].add(language)
        if headword:
            grouped[key]["headwords"].add(headword)
    ds = None

    rows: List[PlacenamePoint] = []
    for _, rec in grouped.items():
        langcodes = tuple(sorted(rec["langcodes"]))
        languages = tuple(sorted(rec["languages"]))
        headwords = tuple(sorted(rec["headwords"]))
        is_scandi = bool(set(langcodes) & SCANDI_LANG_CODES) or bool(set(languages) & SCANDI_LANG_LABELS)
        rows.append(
            PlacenamePoint(
                key=str(rec["key"]),
                placeno=str(rec["placeno"]),
                placename=str(rec["placename"]),
                county=str(rec["county"]),
                x=float(rec["x"]),
                y=float(rec["y"]),
                langcodes=langcodes,
                languages=languages,
                headwords=headwords,
                is_scandinavian=is_scandi,
            )
        )
    return rows


def read_casss_points(path: Path, layer_name: str) -> List[CasssPoint]:
    ds, lyr = open_layer(path, layer_name)
    rows: List[CasssPoint] = []
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        rows.append(
            CasssPoint(
                x=float(geom.GetX()),
                y=float(geom.GetY()),
                name=str(feat.GetField("Name") or "").strip(),
                category=str(feat.GetField("Category") or "").strip(),
                object_type=str(feat.GetField("OBJECT_TYPE") or "").strip(),
            )
        )
    ds = None
    return rows


def county_bucket(county_hist: str) -> str:
    c = county_hist.strip().lower()
    if "northumberland" in c:
        return "Northumberland"
    if "durham" in c:
        return "Durham"
    return "Other"


def safe_ratio(num: float, den: float) -> float:
    if den == 0:
        return float("nan")
    return float(num / den)


def derive_site_rows(
    placenames: Sequence[PlacenamePoint],
    pas_points: Sequence[PasPoint],
    casss_points: Sequence[CasssPoint],
    pas_radius_m: float,
    casss_radius_m: float,
) -> List[Dict[str, object]]:
    if not pas_points:
        return []
    pas_x = np.array([p.x for p in pas_points], dtype=np.float64)
    pas_y = np.array([p.y for p in pas_points], dtype=np.float64)
    pas_grade = np.array([p.sici_grade for p in pas_points], dtype=np.int64)

    casss_x = np.array([c.x for c in casss_points], dtype=np.float64) if casss_points else np.array([], dtype=np.float64)
    casss_y = np.array([c.y for c in casss_points], dtype=np.float64) if casss_points else np.array([], dtype=np.float64)

    rows: List[Dict[str, object]] = []
    for pl in placenames:
        dx = pas_x - pl.x
        dy = pas_y - pl.y
        dist = np.sqrt(dx * dx + dy * dy)
        within = dist <= pas_radius_m
        if not np.any(within):
            continue
        idx = np.where(within)[0]
        local_dist = dist[idx]
        nearest_local = int(idx[int(np.argmin(local_dist))])
        nearest_pas = pas_points[nearest_local]

        linked_counties = [pas_points[i].county_hist for i in idx if pas_points[i].county_hist.strip()]
        county_counts = Counter(linked_counties)
        dominant_county = county_counts.most_common(1)[0][0] if county_counts else ""
        county_bucket_counts = Counter(county_bucket(c) for c in linked_counties)
        n_nd_linked = int(county_bucket_counts.get("Northumberland", 0) + county_bucket_counts.get("Durham", 0))

        casss_min = float("nan")
        casss_n = 0
        if casss_x.size > 0:
            cdx = casss_x - pl.x
            cdy = casss_y - pl.y
            cdist = np.sqrt(cdx * cdx + cdy * cdy)
            casss_n = int(np.sum(cdist <= casss_radius_m))
            casss_min = float(np.min(cdist)) if cdist.size else float("nan")

        linked_ids = sorted({pas_points[i].site_id for i in idx if pas_points[i].site_id})
        linked_classes = sorted({pas_points[i].object_class for i in idx if pas_points[i].object_class})

        rows.append(
            {
                "site_key": pl.key,
                "placeno": pl.placeno,
                "placename": pl.placename,
                "placename_county": pl.county,
                "x": pl.x,
                "y": pl.y,
                "radius_m_used": pas_radius_m,
                "langcodes": "; ".join(pl.langcodes),
                "languages": "; ".join(pl.languages),
                "headwords": "; ".join(pl.headwords),
                "is_scandinavian_toponym": 1 if pl.is_scandinavian else 0,
                "n_pas_within_500m": int(idx.size),
                "n_sici_grade2_within_500m": int(np.sum(pas_grade[idx] == 2)),
                "n_sici_grade1_within_500m": int(np.sum(pas_grade[idx] == 1)),
                "prop_grade2_within_500m": safe_ratio(float(np.sum(pas_grade[idx] == 2)), float(idx.size)),
                "min_pas_dist_m": float(np.min(local_dist)),
                "nearest_pas_site_id": nearest_pas.site_id,
                "nearest_pas_county_hist": nearest_pas.county_hist,
                "nearest_pas_object_class": nearest_pas.object_class,
                "linked_pas_county_hist_mode": dominant_county,
                "n_linked_pas_northumberland_durham": n_nd_linked,
                "is_nd_focused_site": 1 if n_nd_linked > 0 else 0,
                "linked_pas_ids": "; ".join(linked_ids),
                "linked_pas_object_classes": "; ".join(linked_classes),
                "n_casss_within_1000m": casss_n,
                "min_casss_dist_m": casss_min,
            }
        )
    rows.sort(key=lambda r: (-(int(r["n_sici_grade2_within_500m"])), float(r["min_pas_dist_m"]), str(r["placename"])))
    return rows


def derive_hotspot_trace(
    pas_points: Sequence[PasPoint],
    all_placenames: Sequence[PlacenamePoint],
    hotspot_ids: Sequence[str],
) -> List[Dict[str, object]]:
    id_map = {p.site_id: p for p in pas_points}
    all_x = np.array([p.x for p in all_placenames], dtype=np.float64)
    all_y = np.array([p.y for p in all_placenames], dtype=np.float64)
    scandi_indices = np.array([i for i, p in enumerate(all_placenames) if p.is_scandinavian], dtype=np.int64)
    rows: List[Dict[str, object]] = []

    for sid in hotspot_ids:
        if sid not in id_map:
            rows.append(
                {
                    "pas_site_id": sid,
                    "found_in_filtered_pas": 0,
                    "pas_x": "",
                    "pas_y": "",
                    "nearest_any_placename": "",
                    "nearest_any_dist_m": "",
                    "nearest_any_languages": "",
                    "nearest_scandi_placename": "",
                    "nearest_scandi_dist_m": "",
                    "nearest_scandi_languages": "",
                    "nearest_scandi_within_500m": "",
                }
            )
            continue

        p = id_map[sid]
        dx = all_x - p.x
        dy = all_y - p.y
        dist = np.sqrt(dx * dx + dy * dy)

        i_any = int(np.argmin(dist))
        any_pl = all_placenames[i_any]
        any_dist = float(dist[i_any])

        if scandi_indices.size > 0:
            sdist = dist[scandi_indices]
            j = int(np.argmin(sdist))
            i_sc = int(scandi_indices[j])
            sc_pl = all_placenames[i_sc]
            sc_dist = float(sdist[j])
            sc_within = 1 if sc_dist <= PAS_RADIUS_M else 0
            sc_name = sc_pl.placename
            sc_lang = "; ".join(sc_pl.languages)
        else:
            sc_dist = float("nan")
            sc_within = 0
            sc_name = ""
            sc_lang = ""

        rows.append(
            {
                "pas_site_id": sid,
                "found_in_filtered_pas": 1,
                "pas_x": p.x,
                "pas_y": p.y,
                "nearest_any_placename": any_pl.placename,
                "nearest_any_dist_m": any_dist,
                "nearest_any_languages": "; ".join(any_pl.languages),
                "nearest_scandi_placename": sc_name,
                "nearest_scandi_dist_m": sc_dist,
                "nearest_scandi_languages": sc_lang,
                "nearest_scandi_within_500m": sc_within,
            }
        )
    return rows


def derive_nd_distance_diagnostics(
    pas_points: Sequence[PasPoint],
    scandi_placenames: Sequence[PlacenamePoint],
) -> List[Dict[str, object]]:
    sx = np.array([p.x for p in scandi_placenames], dtype=np.float64)
    sy = np.array([p.y for p in scandi_placenames], dtype=np.float64)
    rows: List[Dict[str, object]] = []
    if sx.size == 0:
        return rows

    for p in pas_points:
        c = p.county_hist.lower()
        if ("northumberland" not in c) and ("durham" not in c):
            continue
        d = np.sqrt((sx - p.x) * (sx - p.x) + (sy - p.y) * (sy - p.y))
        i = int(np.argmin(d))
        nearest = scandi_placenames[i]
        rows.append(
            {
                "pas_site_id": p.site_id,
                "county_hist": p.county_hist,
                "sici_grade": p.sici_grade,
                "object_class": p.object_class,
                "pas_x": p.x,
                "pas_y": p.y,
                "nearest_scandi_placename": nearest.placename,
                "nearest_scandi_language": "; ".join(nearest.languages),
                "nearest_scandi_dist_m": float(d[i]),
                "radius_m_used": PAS_RADIUS_M,
                "within_radius_scandi_toponym": 1 if float(d[i]) <= PAS_RADIUS_M else 0,
            }
        )
    rows.sort(key=lambda r: float(r["nearest_scandi_dist_m"]))
    return rows


def write_tsv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with path.open("w", encoding="utf-8", newline="") as f:
            f.write("")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_sites_gpkg(path: Path, layer_name: str, rows: Sequence[Dict[str, object]]) -> None:
    if path.exists():
        path.unlink()
    drv = ogr.GetDriverByName("GPKG")
    ds = drv.CreateDataSource(str(path))
    if ds is None:
        raise RuntimeError(f"Failed to create {path}")
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(27700)
    lyr = ds.CreateLayer(layer_name, srs, ogr.wkbPoint)
    field_defs = [
        ("site_key", ogr.OFTString),
        ("placeno", ogr.OFTString),
        ("placename", ogr.OFTString),
        ("placename_county", ogr.OFTString),
        ("radius_m_used", ogr.OFTReal),
        ("langcodes", ogr.OFTString),
        ("languages", ogr.OFTString),
        ("headwords", ogr.OFTString),
        ("is_scandinavian_toponym", ogr.OFTInteger),
        ("n_pas_within_500m", ogr.OFTInteger),
        ("n_sici_grade2_within_500m", ogr.OFTInteger),
        ("n_sici_grade1_within_500m", ogr.OFTInteger),
        ("prop_grade2_within_500m", ogr.OFTReal),
        ("min_pas_dist_m", ogr.OFTReal),
        ("nearest_pas_site_id", ogr.OFTString),
        ("nearest_pas_county_hist", ogr.OFTString),
        ("nearest_pas_object_class", ogr.OFTString),
        ("linked_pas_county_hist_mode", ogr.OFTString),
        ("n_linked_pas_northumberland_durham", ogr.OFTInteger),
        ("is_nd_focused_site", ogr.OFTInteger),
        ("linked_pas_ids", ogr.OFTString),
        ("linked_pas_object_classes", ogr.OFTString),
        ("n_casss_within_1000m", ogr.OFTInteger),
        ("min_casss_dist_m", ogr.OFTReal),
    ]
    for name, ftype in field_defs:
        lyr.CreateField(ogr.FieldDefn(name, ftype))
    defn = lyr.GetLayerDefn()
    for r in rows:
        feat = ogr.Feature(defn)
        for name, _ in field_defs:
            val = r.get(name)
            if val is None or val == "":
                continue
            feat.SetField(name, val)
        g = ogr.Geometry(ogr.wkbPoint)
        g.AddPoint(float(r["x"]), float(r["y"]))
        feat.SetGeometry(g)
        lyr.CreateFeature(feat)
    ds = None


def write_summary_md(path: Path, summaries: Sequence[Dict[str, object]]) -> None:
    lines: List[str] = []
    lines.append("# SICI Site Register Summary")
    lines.append("")
    radius_label = f"{int(PAS_RADIUS_M)}" if float(PAS_RADIUS_M).is_integer() else f"{PAS_RADIUS_M}"
    for s in summaries:
        lines.append(f"## Window {s['window']}")
        lines.append("")
        lines.append(f"- PAS Scandinavian G1/G2 points used: {s['n_pas_points']}")
        lines.append(f"- KEPN placename points (grouped): {s['n_kepn_grouped']}")
        lines.append(f"- KEPN Scandinavian toponym points: {s['n_kepn_scandi']}")
        lines.append(f"- High-confidence site candidates (<= {radius_label} m): {s['n_sites']}")
        lines.append(f"- ND-focused candidates: {s['n_sites_nd']}")
        lines.append("")
        if s.get("top_sites"):
            lines.append(f"Top sites by SICI grade-2 density then distance (<= {radius_label} m):")
            for row in s["top_sites"]:
                lines.append(
                    "- "
                    f"{row['placename']} | n<=radius={row['n_pas_within_500m']} "
                    f"(g2={row['n_sici_grade2_within_500m']}, g1={row['n_sici_grade1_within_500m']}), "
                    f"min_dist={row['min_pas_dist_m']:.1f} m, "
                    f"ND_focus={row['is_nd_focused_site']}"
                )
            lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    out_dir = OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    all_kepn = read_kepn_grouped(KEPN_PATH, KEPN_LAYER)
    scandi_kepn = [k for k in all_kepn if k.is_scandinavian]
    casss_points = read_casss_points(CASSS_PATH, CASSS_LAYER)

    summaries: List[Dict[str, object]] = []
    manifest: Dict[str, object] = {
        "created_utc": "",
        "input_paths": {
            "pas": str(PAS_PATH),
            "kepn": str(KEPN_PATH),
            "casss": str(CASSS_PATH),
        },
        "parameters": {
            "temporal_mode": TEMPORAL_MODE,
            "pas_radius_m": PAS_RADIUS_M,
            "casss_radius_m": CASSS_RADIUS_M,
            "scandi_langcodes": sorted(SCANDI_LANG_CODES),
            "scandi_languages": sorted(SCANDI_LANG_LABELS),
        },
        "windows": {},
    }

    for y0, y1 in TARGET_WINDOWS:
        pas_points = read_pas_points(PAS_PATH, PAS_LAYER, y0, y1, TEMPORAL_MODE)
        site_rows = derive_site_rows(scandi_kepn, pas_points, casss_points, PAS_RADIUS_M, CASSS_RADIUS_M)
        nd_rows = [r for r in site_rows if int(r["is_nd_focused_site"]) == 1]

        stem = f"sici_site_register_northumbria_{y0}_{y1}"
        write_tsv(out_dir / f"{stem}.tsv", site_rows)
        write_sites_gpkg(out_dir / f"{stem}.gpkg", stem, site_rows)
        write_tsv(out_dir / f"{stem}_nd_focus.tsv", nd_rows)

        trace_rows = derive_hotspot_trace(pas_points, all_kepn, HOTSPOT_PAS_IDS)
        write_tsv(out_dir / f"northumberland_hotspot_trace_{y0}_{y1}.tsv", trace_rows)
        nd_diag_rows = derive_nd_distance_diagnostics(pas_points, scandi_kepn)
        write_tsv(out_dir / f"northumberland_durham_scandi_distance_diagnostics_{y0}_{y1}.tsv", nd_diag_rows)

        summaries.append(
            {
                "window": f"{y0}-{y1}",
                "n_pas_points": len(pas_points),
                "n_kepn_grouped": len(all_kepn),
                "n_kepn_scandi": len(scandi_kepn),
                "n_sites": len(site_rows),
                "n_sites_nd": len(nd_rows),
                "top_sites": site_rows[:10],
            }
        )
        manifest["windows"][f"{y0}_{y1}"] = {
            "n_pas_points": len(pas_points),
            "n_sites": len(site_rows),
            "n_sites_nd": len(nd_rows),
                "outputs": {
                    "register_tsv": str(out_dir / f"{stem}.tsv"),
                    "register_gpkg": str(out_dir / f"{stem}.gpkg"),
                    "nd_focus_tsv": str(out_dir / f"{stem}_nd_focus.tsv"),
                    "hotspot_trace_tsv": str(out_dir / f"northumberland_hotspot_trace_{y0}_{y1}.tsv"),
                    "nd_distance_diagnostics_tsv": str(
                        out_dir / f"northumberland_durham_scandi_distance_diagnostics_{y0}_{y1}.tsv"
                    ),
                },
            }

    write_summary_md(out_dir / "SICI_SITE_REGISTER_SUMMARY.md", summaries)
    manifest["created_utc"] = np.datetime_as_string(np.datetime64("now"), timezone="UTC")
    (out_dir / "SICI_SITE_REGISTER_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"Wrote site register outputs to: {out_dir}")
    for s in summaries:
        print(
            f"Window {s['window']}: PAS={s['n_pas_points']}, "
            f"Sites={s['n_sites']}, ND-focus={s['n_sites_nd']}"
        )


if __name__ == "__main__":
    main()
