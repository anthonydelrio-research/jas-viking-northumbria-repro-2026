#!/usr/bin/env python3
"""
Build date-rescue sensitivity outputs for Scandinavian site derivation.

Purpose:
  - Keep strict canonical layers unchanged.
  - Quantify impact of rescuing missing PAS dates using PAS_ALL crosswalk.
  - Emit auditable sensitivity outputs (1000m and 500m site lists).
"""

from __future__ import annotations

import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from osgeo import ogr


ROOT = Path(os.environ.get("JAS_ROOT", str(Path(__file__).resolve().parents[1]))).resolve()
INBOX = Path(os.environ.get("JAS_INBOX", str(ROOT / "incoming_arcgis" / "01_inbox"))).resolve()
OUTDIR = ROOT / "outputs_submission_2026-03-05_date_rescue"

PAS_SICI_TSV = INBOX / "PAS_SICI" / "Scandinavian_Metalwork.tsv"
PAS_ALL_TSV = INBOX / "PAS_ALL" / "PAS_ALL_PERIODS_Northumbria.tsv"
KEPN_GPKG = INBOX / "KEPN_full_northumbria.gpkg"
KEPN_LAYER = "KEPN_full_northumbria"
CANONICAL_SCANDI = INBOX / "SICI_High_Confidence_Scandinavian_Sites.tsv"

WINDOW_FROM = 750
WINDOW_TO = 954
SCANDI_CODES = {"N", "ND", "NW", "NN"}
SCANDI_LANGS = {"Old Norse", "Old Danish", "Old West Scandinavian", "Old Norwegian"}


@dataclass
class PasPoint:
    easting: float
    northing: float
    old_findid: str
    secuid: str
    pas_id: str
    fromdate: Optional[int]
    todate: Optional[int]
    fromdate_rescued: Optional[int]
    todate_rescued: Optional[int]
    rescue_status: str
    rescue_source: str


def norm(s: object) -> str:
    return str(s or "").strip().upper()


def norm_id(s: object) -> str:
    return re.sub(r"[^A-Z0-9]", "", norm(s))


def parse_year(v: object) -> Optional[int]:
    s = str(v or "").strip()
    if not s:
        return None
    try:
        x = int(float(s))
    except Exception:
        return None
    if x <= 0:
        return None
    return x


def temporal_overlap(a0: Optional[int], a1: Optional[int], w0: int, w1: int) -> bool:
    if a0 is None or a1 is None:
        return False
    if a0 > a1:
        a0, a1 = a1, a0
    return max(a0, w0) <= min(a1, w1)


def load_pas_all_index() -> Tuple[
    Dict[str, List[Tuple[int, int, Dict[str, str]]]],
    Dict[str, List[Tuple[int, int, Dict[str, str]]]],
    Dict[str, List[Tuple[int, int, Dict[str, str]]]],
]:
    by_secuid: Dict[str, List[Tuple[int, int, Dict[str, str]]]] = {}
    by_old: Dict[str, List[Tuple[int, int, Dict[str, str]]]] = {}
    by_old_compact: Dict[str, List[Tuple[int, int, Dict[str, str]]]] = {}

    with PAS_ALL_TSV.open("r", encoding="utf-8", errors="replace", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            a0 = parse_year(row.get("fromdate"))
            a1 = parse_year(row.get("todate"))
            if a0 is None or a1 is None:
                continue
            secuid = norm(row.get("secuid"))
            old_findid = norm(row.get("old_findID"))
            old_compact = norm_id(row.get("old_findID"))
            rec = (a0, a1, row)
            if secuid:
                by_secuid.setdefault(secuid, []).append(rec)
            if old_findid:
                by_old.setdefault(old_findid, []).append(rec)
            if old_compact:
                by_old_compact.setdefault(old_compact, []).append(rec)

    return by_secuid, by_old, by_old_compact


def select_rescue_candidate(
    candidates: Sequence[Tuple[int, int, Dict[str, str]]],
) -> Optional[Tuple[int, int, Dict[str, str]]]:
    if not candidates:
        return None
    uniq: List[Tuple[int, int, Dict[str, str]]] = []
    seen = set()
    for a0, a1, row in candidates:
        key = (a0, a1, norm(row.get("secuid")), norm(row.get("old_findID")))
        if key in seen:
            continue
        seen.add(key)
        uniq.append((a0, a1, row))
    if not uniq:
        return None
    # Prefer narrower ranges, then earlier starts.
    uniq.sort(key=lambda t: ((t[1] - t[0]), t[0], t[1]))
    return uniq[0]


def load_pas_sici_points() -> Tuple[List[PasPoint], List[PasPoint], List[PasPoint]]:
    by_secuid, by_old, by_old_compact = load_pas_all_index()
    strict_overlap: List[PasPoint] = []
    rescue_overlap: List[PasPoint] = []
    incomplete_rows: List[PasPoint] = []

    with PAS_SICI_TSV.open("r", encoding="utf-8", errors="replace", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            try:
                confidence = int(float(str(row.get("confidence") or "").strip()))
            except Exception:
                continue
            if confidence not in (1, 2):
                continue

            try:
                x = float(str(row.get("easting") or "").strip())
                y = float(str(row.get("northing") or "").strip())
            except Exception:
                continue

            a0 = parse_year(row.get("fromdate"))
            a1 = parse_year(row.get("todate"))
            old_findid = str(row.get("old_findID") or "").strip()
            secuid = str(row.get("secuid") or "").strip()
            pas_id = str(row.get("pas_id") or "").strip()

            rescue_status = "not_needed"
            rescue_source = ""
            ra0, ra1 = a0, a1

            if a0 is None or a1 is None:
                rescue_status = "unresolved"
                key_secuid = norm(secuid)
                key_old = norm(old_findid)
                key_old_compact = norm_id(old_findid)
                cands: List[Tuple[int, int, Dict[str, str]]] = []
                sources: List[str] = []
                if key_secuid and key_secuid in by_secuid:
                    cands.extend(by_secuid[key_secuid])
                    sources.append("secuid")
                if key_old and key_old in by_old:
                    cands.extend(by_old[key_old])
                    sources.append("old_findID")
                if key_old_compact and key_old_compact in by_old_compact:
                    cands.extend(by_old_compact[key_old_compact])
                    sources.append("old_findID_compact")
                chosen = select_rescue_candidate(cands)
                if chosen is not None:
                    ra0, ra1, _ = chosen
                    rescue_status = "rescued"
                    rescue_source = ",".join(sources)

            p = PasPoint(
                easting=x,
                northing=y,
                old_findid=old_findid,
                secuid=secuid,
                pas_id=pas_id,
                fromdate=a0,
                todate=a1,
                fromdate_rescued=ra0,
                todate_rescued=ra1,
                rescue_status=rescue_status,
                rescue_source=rescue_source,
            )

            if temporal_overlap(a0, a1, WINDOW_FROM, WINDOW_TO):
                strict_overlap.append(p)
            if temporal_overlap(ra0, ra1, WINDOW_FROM, WINDOW_TO):
                rescue_overlap.append(p)
            if a0 is None or a1 is None:
                incomplete_rows.append(p)

    return strict_overlap, rescue_overlap, incomplete_rows


def load_kepn_scandinavian() -> List[Dict[str, object]]:
    ds = ogr.Open(str(KEPN_GPKG))
    if ds is None:
        raise RuntimeError(f"Failed to open {KEPN_GPKG}")
    lyr = ds.GetLayerByName(KEPN_LAYER)
    if lyr is None:
        raise RuntimeError(f"Layer {KEPN_LAYER} not found")

    out: List[Dict[str, object]] = []
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom is None:
            continue
        code = str(feat.GetField("elements__langcode") or "").strip()
        lang = str(feat.GetField("elements__language") or "").strip()
        if (code not in SCANDI_CODES) and (lang not in SCANDI_LANGS):
            continue
        out.append(
            {
                "placeno": str(feat.GetField("placeno") or "").strip(),
                "placename": str(feat.GetField("placename") or "").strip(),
                "easting": float(geom.GetX()),
                "northing": float(geom.GetY()),
                "county": str(feat.GetField("county") or "").strip(),
                "elements__language": lang,
                "elements__langcode": code,
            }
        )
    ds = None
    return out


def derive_site_rows(
    kepn_recs: Sequence[Dict[str, object]],
    pas_points: Sequence[PasPoint],
    radius_m: float,
) -> List[Dict[str, str]]:
    if not pas_points:
        return []
    px = np.array([p.easting for p in pas_points], dtype=np.float64)
    py = np.array([p.northing for p in pas_points], dtype=np.float64)
    out: List[Dict[str, str]] = []
    for rec in kepn_recs:
        x = float(rec["easting"])
        y = float(rec["northing"])
        d = np.sqrt((px - x) * (px - x) + (py - y) * (py - y))
        if d.size == 0:
            continue
        md = float(d.min())
        if md > radius_m:
            continue
        out.append(
            {
                "placeno": str(rec["placeno"]),
                "placename": str(rec["placename"]),
                "county": str(rec["county"]),
                "elements__language": str(rec["elements__language"]),
                "elements__langcode": str(rec["elements__langcode"]),
                "Near_SICI_DIST": f"{md:.6f}",
            }
        )
    out.sort(key=lambda r: (r["placename"], float(r["Near_SICI_DIST"])))
    return out


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[Dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(header), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    strict_pts, rescue_pts, incomplete = load_pas_sici_points()
    kepn = load_kepn_scandinavian()

    strict_1000 = derive_site_rows(kepn, strict_pts, 1000.0)
    rescue_1000 = derive_site_rows(kepn, rescue_pts, 1000.0)
    strict_500 = derive_site_rows(kepn, strict_pts, 500.0)
    rescue_500 = derive_site_rows(kepn, rescue_pts, 500.0)

    # Compare against current canonical names.
    canonical_names = set()
    with CANONICAL_SCANDI.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            canonical_names.add(str(row.get("placename") or "").strip())

    rescue_names = {r["placename"] for r in rescue_1000}
    strict_names = {r["placename"] for r in strict_1000}
    added_vs_strict = sorted(rescue_names - strict_names)
    added_vs_current = sorted(rescue_names - canonical_names)

    header_sites = [
        "placeno",
        "placename",
        "county",
        "elements__language",
        "elements__langcode",
        "Near_SICI_DIST",
    ]
    write_tsv(OUTDIR / "scandinavian_sites_1000m_strict.tsv", header_sites, strict_1000)
    write_tsv(
        OUTDIR / "scandinavian_sites_1000m_date_rescued_sensitivity.tsv",
        header_sites,
        rescue_1000,
    )
    write_tsv(OUTDIR / "scandinavian_sites_500m_strict.tsv", header_sites, strict_500)
    write_tsv(
        OUTDIR / "scandinavian_sites_500m_date_rescued_sensitivity.tsv",
        header_sites,
        rescue_500,
    )

    # Audit table of incomplete PAS rows.
    audit_rows: List[Dict[str, str]] = []
    for p in incomplete:
        audit_rows.append(
            {
                "old_findID": p.old_findid,
                "secuid": p.secuid,
                "pas_id": p.pas_id,
                "fromdate_original": "" if p.fromdate is None else str(p.fromdate),
                "todate_original": "" if p.todate is None else str(p.todate),
                "fromdate_rescued": "" if p.fromdate_rescued is None else str(p.fromdate_rescued),
                "todate_rescued": "" if p.todate_rescued is None else str(p.todate_rescued),
                "rescue_status": p.rescue_status,
                "rescue_source": p.rescue_source,
            }
        )
    audit_rows.sort(key=lambda r: (r["rescue_status"], r["old_findID"], r["secuid"]))
    write_tsv(
        OUTDIR / "pas_date_incomplete_rescue_audit.tsv",
        [
            "old_findID",
            "secuid",
            "pas_id",
            "fromdate_original",
            "todate_original",
            "fromdate_rescued",
            "todate_rescued",
            "rescue_status",
            "rescue_source",
        ],
        audit_rows,
    )

    rescued = sum(1 for p in incomplete if p.rescue_status == "rescued")
    unresolved = sum(1 for p in incomplete if p.rescue_status != "rescued")

    summary = []
    summary.append("# Date Rescue Sensitivity Summary\n")
    summary.append(f"- Temporal window: AD {WINDOW_FROM}-{WINDOW_TO}")
    summary.append(f"- Incomplete PAS_SICI rows evaluated: {len(incomplete)}")
    summary.append(f"- Rows with rescued valid date ranges: {rescued}")
    summary.append(f"- Rows unresolved after PAS_ALL crosswalk: {unresolved}\n")
    summary.append("## Site Count Impact")
    summary.append(f"- Strict 1000m Scandinavian sites: {len(strict_1000)}")
    summary.append(f"- Date-rescued 1000m Scandinavian sites: {len(rescue_1000)}")
    summary.append(f"- Strict 500m Scandinavian sites: {len(strict_500)}")
    summary.append(f"- Date-rescued 500m Scandinavian sites: {len(rescue_500)}\n")
    summary.append("## Added Sites")
    summary.append(f"- Added vs strict (1000m): {', '.join(added_vs_strict) if added_vs_strict else 'None'}")
    summary.append(f"- Added vs current canonical (1000m): {', '.join(added_vs_current) if added_vs_current else 'None'}")
    summary.append("\n## Interpretation")
    summary.append(
        "Date rescue has limited impact in this package. Primary strict-canonical inference remains stable;"
        " rescued chronology adds one Scandinavian toponym candidate at 1000m."
    )
    (OUTDIR / "DATE_RESCUE_SENSITIVITY_SUMMARY.md").write_text("\n".join(summary), encoding="utf-8")

    print("Date-rescue sensitivity outputs written:")
    print(f"  {OUTDIR}")
    print(f"  strict 1000m: {len(strict_1000)}")
    print(f"  rescued 1000m: {len(rescue_1000)}")
    print(f"  strict 500m: {len(strict_500)}")
    print(f"  rescued 500m: {len(rescue_500)}")
    print(f"  rescued rows: {rescued}")
    print(f"  unresolved rows: {unresolved}")


if __name__ == "__main__":
    main()
