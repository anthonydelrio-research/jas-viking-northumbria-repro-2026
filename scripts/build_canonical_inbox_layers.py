#!/usr/bin/env python3
"""
Build canonical KEPN-derived SICI site layers and archive legacy inbox tables.

Canonical policy used here:
  - PAS evidence pool: Scandinavian_Metalwork, SICI grade 1/2, temporal overlap 750-954
  - Scandinavian toponym layer: radius 1000 m
  - Anglo-Saxon/Anglian contextual layer: radius 500 m (plus optional 1000 m contextual export)
"""

from __future__ import annotations

import csv
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from osgeo import ogr, osr


ROOT = Path(os.environ.get("JAS_ROOT", str(Path(__file__).resolve().parents[1]))).resolve()
INBOX = Path(os.environ.get("JAS_INBOX", str(ROOT / "incoming_arcgis" / "01_inbox"))).resolve()
ARCHIVE_ROOT = ROOT / "incoming_arcgis" / "03_validated" / "2026-03-05_legacy_layer_archive"

PAS_PATH = INBOX / "PAS_full_northumbria_metalwork.gpkg"
PAS_LAYER = "PAS_full_northumbria_metalwork"
KEPN_PATH = INBOX / "KEPN_full_northumbria.gpkg"
KEPN_LAYER = "KEPN_full_northumbria"

WINDOW_FROM = 750
WINDOW_TO = 954
TEMPORAL_MODE = "overlap"

SCANDI_RADIUS_M = 1000.0
ANGLO_RADIUS_M = 500.0
ANGLO_RADIUS_CONTEXT_M = 1000.0

SCANDI_CODES = {"N", "ND", "NW", "NN"}
SCANDI_LANGS = {"Old Norse", "Old Danish", "Old West Scandinavian", "Old Norwegian"}
ANGLO_CODES = {"O", "OA"}
ANGLO_LANGS = {"Old English", "Anglian"}

SCANDI_HEADER = [
    "placeno",
    "placename",
    "easting",
    "northing",
    "county",
    "pre74countycode",
    "hundred",
    "parish",
    "ptypecode",
    "inuse",
    "etymology",
    "comment",
    "meetknowncode",
    "elements__headword",
    "elements__hword",
    "elements__hversion",
    "elements__language",
    "elements__langcode",
    "elements__note",
    "NAME_Element",
    "Join_Count",
    "ALC_GRADE",
    "ALC_Grades__Provisional____ADAS___Defra_ALC_GRADE",
    "RomanRoadNEAR_FID",
    "RomanRoadNEAR_DIST",
    "WatercourseNEAR_FID",
    "WatercourseNEAR_DIST",
    "Near_SICI_DIST",
    "NEAR_FID",
    "PAS_RADIUS_M",
    "TEMPORAL_WINDOW",
    "TEMPORAL_MODE",
    "CANONICAL_VERSION",
    "CANONICAL_TIMESTAMP_UTC",
]

ANGLO_HEADER = [
    "Mound",
    "F1",
    "placeno",
    "placename",
    "easting",
    "northing",
    "lat",
    "lng",
    "east_res",
    "north_res",
    "county",
    "pre74countycode",
    "hundred",
    "parish",
    "ptypecode",
    "inuse",
    "etymology",
    "comment",
    "meetknowncode",
    "elements__headword",
    "elements__hword",
    "elements__hversion",
    "elements__language",
    "elements__langcode",
    "elements__note",
    "NEAR_DIST",
    "NEAR_FID",
    "PAS_RADIUS_M",
    "TEMPORAL_WINDOW",
    "TEMPORAL_MODE",
    "CANONICAL_VERSION",
    "CANONICAL_TIMESTAMP_UTC",
]


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


def parse_year(v: object) -> Optional[int]:
    if v is None:
        return None
    s = str(v).strip()
    if not s:
        return None
    try:
        return int(float(s))
    except ValueError:
        return None


def open_layer(path: Path, layer: str):
    ds = ogr.Open(str(path))
    if ds is None:
        raise RuntimeError(f"Failed to open {path}")
    lyr = ds.GetLayerByName(layer)
    if lyr is None:
        raise RuntimeError(f"Layer {layer} not found in {path}")
    return ds, lyr


def load_pas_points() -> Tuple[np.ndarray, np.ndarray, List[int]]:
    ds, lyr = open_layer(PAS_PATH, PAS_LAYER)
    xs: List[float] = []
    ys: List[float] = []
    fids: List[int] = []
    for feat in lyr:
        if str(feat.GetField("source_subset") or "").strip() != "Scandinavian_Metalwork":
            continue
        try:
            grade = int(feat.GetField("sici_grade"))
        except Exception:
            continue
        if grade not in (1, 2):
            continue
        d0 = parse_year(feat.GetField("date_from_ad"))
        d1 = parse_year(feat.GetField("date_to_ad"))
        if not temporal_pass(d0, d1, WINDOW_FROM, WINDOW_TO, TEMPORAL_MODE):
            continue
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        xs.append(float(geom.GetX()))
        ys.append(float(geom.GetY()))
        fids.append(int(feat.GetFID()))
    ds = None
    if not xs:
        raise RuntimeError("No PAS points loaded under canonical criteria.")
    return np.array(xs, dtype=np.float64), np.array(ys, dtype=np.float64), fids


def load_kepn_features() -> List[Dict[str, object]]:
    ds, lyr = open_layer(KEPN_PATH, KEPN_LAYER)
    out: List[Dict[str, object]] = []
    for feat in lyr:
        def get_field_text(name: str) -> str:
            # Some KEPN exports omit fields (e.g., "inuse"); treat as optional.
            try:
                idx = feat.GetFieldIndex(name)
            except Exception:
                idx = -1
            if idx is not None and idx >= 0:
                try:
                    v = feat.GetField(idx)
                    return str(v or "").strip()
                except Exception:
                    return ""
            try:
                v = feat.GetField(name)
                return str(v or "").strip()
            except Exception:
                return ""

        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        out.append(
            {
                "placeno": get_field_text("placeno"),
                "placename": get_field_text("placename"),
                "county": get_field_text("county"),
                "pre74countycode": get_field_text("pre74countycode"),
                "hundred": get_field_text("hundred"),
                "parish": get_field_text("parish"),
                "ptypecode": get_field_text("ptypecode"),
                "inuse": get_field_text("inuse"),
                "etymology": get_field_text("etymology"),
                "comment": get_field_text("comment"),
                "meetknowncode": get_field_text("meetknowncode"),
                "elements__headword": get_field_text("elements__headword"),
                "elements__hword": get_field_text("elements__hword"),
                "elements__hversion": get_field_text("elements__hversion"),
                "elements__language": get_field_text("elements__language"),
                "elements__langcode": get_field_text("elements__langcode"),
                "elements__note": get_field_text("elements__note"),
                "x": float(geom.GetX()),
                "y": float(geom.GetY()),
            }
        )
    ds = None
    return out


def keep_scandi(rec: Dict[str, object]) -> bool:
    code = str(rec.get("elements__langcode", ""))
    lang = str(rec.get("elements__language", ""))
    return (code in SCANDI_CODES) or (lang in SCANDI_LANGS)


def keep_anglo(rec: Dict[str, object]) -> bool:
    code = str(rec.get("elements__langcode", ""))
    lang = str(rec.get("elements__language", ""))
    return (code in ANGLO_CODES) or (lang in ANGLO_LANGS)


def build_rows(
    kepn: Sequence[Dict[str, object]],
    px: np.ndarray,
    py: np.ndarray,
    pfids: Sequence[int],
    radius_m: float,
    keep_fn,
    header: Sequence[str],
    canonical_version: str,
) -> List[Dict[str, str]]:
    ts = datetime.now(timezone.utc).isoformat()
    rows: List[Dict[str, str]] = []
    src = osr.SpatialReference()
    src.ImportFromEPSG(27700)
    dst = osr.SpatialReference()
    dst.ImportFromEPSG(4326)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    dst.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    ct = osr.CoordinateTransformation(src, dst)

    for idx, rec in enumerate(kepn, start=1):
        if not keep_fn(rec):
            continue
        x = float(rec["x"])
        y = float(rec["y"])
        d = np.sqrt((px - x) * (px - x) + (py - y) * (py - y))
        if d.size == 0:
            continue
        i = int(np.argmin(d))
        min_d = float(d[i])
        if min_d > radius_m:
            continue
        near_fid = str(pfids[i])

        row: Dict[str, str] = {k: "" for k in header}
        row["placeno"] = str(rec.get("placeno", ""))
        row["placename"] = str(rec.get("placename", ""))
        row["easting"] = f"{x:.15f}"
        row["northing"] = f"{y:.15f}"
        row["county"] = str(rec.get("county", ""))
        row["pre74countycode"] = str(rec.get("pre74countycode", ""))
        row["hundred"] = str(rec.get("hundred", ""))
        row["parish"] = str(rec.get("parish", ""))
        row["ptypecode"] = str(rec.get("ptypecode", ""))
        row["inuse"] = str(rec.get("inuse", ""))
        row["etymology"] = str(rec.get("etymology", ""))
        row["comment"] = str(rec.get("comment", ""))
        row["meetknowncode"] = str(rec.get("meetknowncode", ""))
        row["elements__headword"] = str(rec.get("elements__headword", ""))
        row["elements__hword"] = str(rec.get("elements__hword", ""))
        row["elements__hversion"] = str(rec.get("elements__hversion", ""))
        row["elements__language"] = str(rec.get("elements__language", ""))
        row["elements__langcode"] = str(rec.get("elements__langcode", ""))
        row["elements__note"] = str(rec.get("elements__note", ""))

        if "Near_SICI_DIST" in row:
            row["Near_SICI_DIST"] = f"{min_d:.15f}"
        if "NEAR_DIST" in row:
            row["NEAR_DIST"] = f"{min_d:.15f}"
        row["NEAR_FID"] = near_fid
        row["PAS_RADIUS_M"] = f"{radius_m:.1f}"
        row["TEMPORAL_WINDOW"] = f"{WINDOW_FROM}-{WINDOW_TO}"
        row["TEMPORAL_MODE"] = TEMPORAL_MODE
        row["CANONICAL_VERSION"] = canonical_version
        row["CANONICAL_TIMESTAMP_UTC"] = ts

        if "Mound" in row:
            row["Mound"] = ""
        if "F1" in row:
            row["F1"] = f"{float(idx):.15f}"
        if "east_res" in row:
            row["east_res"] = "100.000000000000000"
        if "north_res" in row:
            row["north_res"] = "100.000000000000000"
        if "lat" in row or "lng" in row:
            pt = ogr.Geometry(ogr.wkbPoint)
            pt.AddPoint(x, y)
            pt.Transform(ct)
            lng = float(pt.GetX())
            lat = float(pt.GetY())
            if "lat" in row:
                row["lat"] = f"{lat:.15f}"
            if "lng" in row:
                row["lng"] = f"{lng:.15f}"

        rows.append(row)

    def sort_key(r: Dict[str, str]) -> Tuple[float, str]:
        dist_field = "Near_SICI_DIST" if "Near_SICI_DIST" in r else "NEAR_DIST"
        d = float(r.get(dist_field, "inf") or "inf")
        return (d, str(r.get("placename", "")))

    rows.sort(key=sort_key)
    return rows


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(header), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_gpkg(path: Path, layer_name: str, header: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
    if path.exists():
        path.unlink()
    drv = ogr.GetDriverByName("GPKG")
    ds = drv.CreateDataSource(str(path))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(27700)
    lyr = ds.CreateLayer(layer_name, srs, ogr.wkbPoint)
    for h in header:
        lyr.CreateField(ogr.FieldDefn(h, ogr.OFTString))
    defn = lyr.GetLayerDefn()
    for r in rows:
        try:
            x = float(r.get("easting", ""))
            y = float(r.get("northing", ""))
        except Exception:
            continue
        feat = ogr.Feature(defn)
        for h in header:
            feat.SetField(h, str(r.get(h, "")))
        g = ogr.Geometry(ogr.wkbPoint)
        g.AddPoint(x, y)
        feat.SetGeometry(g)
        lyr.CreateFeature(feat)
    ds = None


def archive_if_exists(path: Path, archive_dir: Path) -> None:
    if not path.exists():
        return
    archive_dir.mkdir(parents=True, exist_ok=True)
    dst = archive_dir / path.name
    shutil.copy2(path, dst)


def main() -> None:
    px, py, pfids = load_pas_points()
    kepn = load_kepn_features()

    scandi_rows_1000 = build_rows(
        kepn,
        px,
        py,
        pfids,
        SCANDI_RADIUS_M,
        keep_scandi,
        SCANDI_HEADER,
        canonical_version="2026-03-05-scandi1000-v1",
    )
    scandi_rows_500 = build_rows(
        kepn,
        px,
        py,
        pfids,
        500.0,
        keep_scandi,
        SCANDI_HEADER,
        canonical_version="2026-03-05-scandi500-v1",
    )

    anglo_rows_500 = build_rows(
        kepn,
        px,
        py,
        pfids,
        ANGLO_RADIUS_M,
        keep_anglo,
        ANGLO_HEADER,
        canonical_version="2026-03-05-anglo500-v1",
    )
    anglo_rows_1000 = build_rows(
        kepn,
        px,
        py,
        pfids,
        ANGLO_RADIUS_CONTEXT_M,
        keep_anglo,
        ANGLO_HEADER,
        canonical_version="2026-03-05-anglo1000-v1",
    )

    # Archive legacy files before replacement.
    archive_if_exists(INBOX / "SICI_High_Confidence_Scandinavian_Sites.tsv", ARCHIVE_ROOT)
    archive_if_exists(INBOX / "SICI_High_Confidence_Scandinavian_Sites.tsv.xml", ARCHIVE_ROOT)
    archive_if_exists(INBOX / "SICI_High_Confidence_Anglo_Saxon_Sites_500m.tsv", ARCHIVE_ROOT)
    archive_if_exists(INBOX / "SICI_High_Confidence_Anglo_Saxon_Sites_500m.tsv.xml", ARCHIVE_ROOT)
    archive_if_exists(INBOX / "SICI_High_Confidence_Anglo_Saxon_Sites_500m_updated.gpkg", ARCHIVE_ROOT)

    # Replace working inbox layer names with canonical tables.
    write_tsv(INBOX / "SICI_High_Confidence_Scandinavian_Sites.tsv", SCANDI_HEADER, scandi_rows_1000)
    write_tsv(INBOX / "SICI_High_Confidence_Anglo_Saxon_Sites_500m.tsv", ANGLO_HEADER, anglo_rows_500)

    # Canonical explicit variants + gpkg companions.
    write_tsv(INBOX / "SICI_High_Confidence_Scandinavian_Sites_1000m_CANONICAL.tsv", SCANDI_HEADER, scandi_rows_1000)
    write_tsv(INBOX / "SICI_High_Confidence_Scandinavian_Sites_500m_CANONICAL.tsv", SCANDI_HEADER, scandi_rows_500)
    write_tsv(INBOX / "SICI_Context_Anglo_Saxon_Sites_500m_CANONICAL.tsv", ANGLO_HEADER, anglo_rows_500)
    write_tsv(INBOX / "SICI_Context_Anglo_Saxon_Sites_1000m_CANONICAL.tsv", ANGLO_HEADER, anglo_rows_1000)

    write_gpkg(
        INBOX / "SICI_High_Confidence_Scandinavian_Sites_1000m_CANONICAL.gpkg",
        "SICI_High_Confidence_Scandinavian_Sites_1000m_CANONICAL",
        SCANDI_HEADER,
        scandi_rows_1000,
    )
    write_gpkg(
        INBOX / "SICI_High_Confidence_Scandinavian_Sites_500m_CANONICAL.gpkg",
        "SICI_High_Confidence_Scandinavian_Sites_500m_CANONICAL",
        SCANDI_HEADER,
        scandi_rows_500,
    )
    write_gpkg(
        INBOX / "SICI_Context_Anglo_Saxon_Sites_500m_CANONICAL.gpkg",
        "SICI_Context_Anglo_Saxon_Sites_500m_CANONICAL",
        ANGLO_HEADER,
        anglo_rows_500,
    )
    write_gpkg(
        INBOX / "SICI_Context_Anglo_Saxon_Sites_1000m_CANONICAL.gpkg",
        "SICI_Context_Anglo_Saxon_Sites_1000m_CANONICAL",
        ANGLO_HEADER,
        anglo_rows_1000,
    )

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "window": f"{WINDOW_FROM}-{WINDOW_TO}",
        "temporal_mode": TEMPORAL_MODE,
        "pas_pool": "Scandinavian_Metalwork SICI 1/2",
        "counts": {
            "scandinavian_1000m": len(scandi_rows_1000),
            "scandinavian_500m": len(scandi_rows_500),
            "anglo_500m": len(anglo_rows_500),
            "anglo_1000m": len(anglo_rows_1000),
        },
        "archive_dir": str(ARCHIVE_ROOT),
        "replaced_files": [
            str(INBOX / "SICI_High_Confidence_Scandinavian_Sites.tsv"),
            str(INBOX / "SICI_High_Confidence_Anglo_Saxon_Sites_500m.tsv"),
        ],
    }
    (INBOX / "SICI_CANONICAL_LAYER_MANIFEST_2026-03-05.json").write_text(
        __import__("json").dumps(manifest, indent=2), encoding="utf-8"
    )

    print("Canonical layers built.")
    print(f"  Scandinavian 1000m: {len(scandi_rows_1000)}")
    print(f"  Scandinavian 500m : {len(scandi_rows_500)}")
    print(f"  Anglo-Saxon 500m  : {len(anglo_rows_500)}")
    print(f"  Anglo-Saxon 1000m : {len(anglo_rows_1000)}")
    print(f"  Archive: {ARCHIVE_ROOT}")


if __name__ == "__main__":
    main()
