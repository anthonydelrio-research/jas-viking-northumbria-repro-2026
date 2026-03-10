#!/usr/bin/env python3
"""
Harvest open-source archaeological records to build provisional validation holdouts.

Sources currently supported:
1) ADS Data Catalogue API (JSON)
2) York HER monument search (public HTML)

Outputs are written under:
JAS/incoming_arcgis/01_inbox/open_harvest_YYYY-MM-DD/
"""

from __future__ import annotations

import csv
import datetime as dt
import html
import json
import os
import re
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ROOT = Path(os.environ.get("JAS_ROOT", str(Path(__file__).resolve().parents[1]))).resolve()
INBOX = Path(os.environ.get("JAS_INBOX", str(ROOT / "incoming_arcgis" / "01_inbox"))).resolve()
TODAY = dt.date.today().isoformat()
OUTDIR = INBOX / f"open_harvest_{TODAY}"

ADS_API = "https://archaeologydataservice.ac.uk/data-catalogue-api/api/search"
YORK_SEARCH = "https://her.york.gov.uk/monuments/search"
UA = (
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0 Safari/537.36"
)

ADS_QUERY_TERMS = [
    "viking yorkshire",
    "viking cumbria",
    "viking lancashire",
    "scandinavian yorkshire",
    "scandinavian cumbria",
    "norse yorkshire",
    "anglo-scandinavian yorkshire",
    "early medieval settlement yorkshire",
    "early medieval settlement cumbria",
]
MAX_ADS_PAGES_PER_TERM = 8

COUNTY_TOKENS = [
    "YORKSHIRE",
    "NORTH YORKSHIRE",
    "WEST YORKSHIRE",
    "EAST RIDING",
    "LANCASHIRE",
    "CUMBRIA",
    "WESTMORLAND",
    "CUMBERLAND",
    "NORTHUMBERLAND",
    "DURHAM",
]

SETTLEMENT_KEYWORDS = [
    "settlement",
    "farm",
    "farmstead",
    "village",
    "hamlet",
    "toft",
    "croft",
    "grubenhaus",
    "sunken-featured",
    "sunken featured",
    "hall",
    "occupation",
    "domestic",
    "enclosure",
]

EXCAVATION_KEYWORDS = [
    "excavat",
    "trench",
    "phase",
    "post-hole",
    "post hole",
    "pit",
    "hearth",
    "structure",
]

TARGET_CULTURE_TERMS = [
    "anglian",
    "anglo-saxon",
    "anglo saxon",
    "anglo-danish",
    "anglo danish",
    "anglo-scandinavian",
    "anglo scandinavian",
    "scandinavian",
    "viking",
    "early medieval",
    "pre-conquest",
]

LIKELY_NON_SETTLEMENT = [
    "church",
    "chapel",
    "bridge",
    "hospital",
    "abbey",
    "priory",
    "cross",
]

ARTIFACT_TERMS = [
    "hoard",
    "clasp",
    "hooked tag",
    "brooch",
    "coin",
    "strap",
    "mount",
    "pin",
    "sword",
]


def ensure_outdir() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)


def fetch_text(url: str, retries: int = 4, sleep_s: float = 1.0) -> str:
    last_err: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(
            url,
            headers={
                "User-Agent": UA,
                "Accept": "application/json,text/html,*/*",
            },
        )
        try:
            with urllib.request.urlopen(req, timeout=45) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as err:  # noqa: BLE001
            last_err = err
            if attempt == retries:
                break
            time.sleep(sleep_s * attempt)
    raise RuntimeError(f"Failed URL after retries: {url} :: {last_err}")


def parse_int(value: object) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    m = re.search(r"-?\d+", text)
    if not m:
        return None
    try:
        return int(m.group(0))
    except ValueError:
        return None


def normalize_text(value: object) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def any_keyword(text: str, keywords: Iterable[str]) -> List[str]:
    t = text.lower()
    return [k for k in keywords if k in t]


def county_hit(place_name: str) -> Optional[str]:
    upper = place_name.upper()
    for token in COUNTY_TOKENS:
        if token in upper:
            return token
    return None


def temporal_overlap(temporal_items: object, y0: int = 750, y1: int = 1100) -> bool:
    if not isinstance(temporal_items, list):
        return False
    for item in temporal_items:
        if not isinstance(item, dict):
            continue
        f = parse_int(item.get("from"))
        u = parse_int(item.get("until"))
        name = normalize_text(item.get("periodName")).lower()
        if "early medieval" in name or "saxon" in name or "viking" in name:
            return True
        if f is None or u is None:
            continue
        lo, hi = min(f, u), max(f, u)
        if hi >= y0 and lo <= y1:
            return True
    return False


def transform_wgs84_to_bng(lon: float, lat: float) -> Tuple[Optional[float], Optional[float]]:
    try:
        proc = subprocess.run(
            [
                "gdaltransform",
                "-s_srs",
                "EPSG:4326",
                "-t_srs",
                "EPSG:27700",
            ],
            input=f"{lon} {lat}\n".encode("utf-8"),
            capture_output=True,
            check=True,
        )
        out = proc.stdout.decode("utf-8", errors="replace").strip()
        if not out:
            return None, None
        parts = out.split()
        if len(parts) < 2:
            return None, None
        return float(parts[0]), float(parts[1])
    except Exception:  # noqa: BLE001
        return None, None


def harvest_ads() -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    raw_rows: List[Dict[str, object]] = []
    candidates: List[Dict[str, object]] = []
    seen_keys = set()

    for term in ADS_QUERY_TERMS:
        print(f"[ADS] term='{term}'", flush=True)
        page = 0
        while True:
            if page >= MAX_ADS_PAGES_PER_TERM:
                print(f"[ADS] capped at {MAX_ADS_PAGES_PER_TERM} pages for term '{term}'", flush=True)
                break
            params = {
                "q": term,
                "country": "England",
                "size": 50,
                "page": page,
            }
            url = ADS_API + "?" + urllib.parse.urlencode(params)
            text = fetch_text(url)
            payload = json.loads(text)
            hits = payload.get("hits", [])
            if not hits:
                break

            for hit in hits:
                data = hit.get("data", {})
                spatial = data.get("spatial", [])
                first_spatial = spatial[0] if spatial else {}
                place_name = normalize_text(first_spatial.get("placeName"))
                county_token = county_hit(place_name)

                centroid = first_spatial.get("centroid", {}) or {}
                geopoint = first_spatial.get("geopoint", {}) or {}
                lat = centroid.get("lat", geopoint.get("lat"))
                lon = centroid.get("lon", geopoint.get("lon"))

                title = normalize_text((data.get("title") or {}).get("text"))
                desc = normalize_text((data.get("description") or {}).get("text"))
                native_subject = " | ".join(
                    normalize_text(x.get("prefLabel"))
                    for x in data.get("nativeSubject", [])
                    if isinstance(x, dict)
                )
                derived_subject = " | ".join(
                    normalize_text(x.get("prefLabel"))
                    for x in data.get("derivedSubject", [])
                    if isinstance(x, dict)
                )
                period_names = " | ".join(
                    normalize_text(x.get("periodName"))
                    for x in data.get("temporal", [])
                    if isinstance(x, dict)
                )

                text_blob = " ".join([title, desc, native_subject, derived_subject, period_names]).lower()
                settlement_hits = any_keyword(text_blob, SETTLEMENT_KEYWORDS)
                excavation_hits = any_keyword(text_blob, EXCAVATION_KEYWORDS)
                non_settlement_hits = any_keyword(text_blob, LIKELY_NON_SETTLEMENT)
                target_terms = any_keyword(text_blob, TARGET_CULTURE_TERMS)
                artifact_terms = any_keyword(title.lower(), ARTIFACT_TERMS)
                temporal_ok = temporal_overlap(data.get("temporal", []))
                likely_settlement = bool(
                    settlement_hits
                    and (not non_settlement_hits or excavation_hits)
                    and target_terms
                    and (not artifact_terms)
                )
                false_positive = "viking way" in title.lower()

                row = {
                    "source": "ADS",
                    "query_term": term,
                    "ads_hit_id": normalize_text(hit.get("id")),
                    "original_id": normalize_text(data.get("originalId")),
                    "title": title,
                    "place_name": place_name,
                    "county_token": county_token or "",
                    "lat": lat,
                    "lon": lon,
                    "period_names": period_names,
                    "native_subject": native_subject,
                    "derived_subject": derived_subject,
                    "description": desc,
                    "settlement_hits": ",".join(settlement_hits),
                    "excavation_hits": ",".join(excavation_hits),
                    "non_settlement_hits": ",".join(non_settlement_hits),
                    "target_terms": ",".join(target_terms),
                    "artifact_terms": ",".join(artifact_terms),
                    "temporal_ok": temporal_ok,
                    "likely_settlement": likely_settlement,
                    "false_positive": false_positive,
                    "record_url": f"https://archaeologydataservice.ac.uk/archsearch/record?titleId={normalize_text(data.get('originalId'))}",
                }
                raw_rows.append(row)

                key = (
                    normalize_text(data.get("originalId")) or normalize_text(hit.get("id")),
                    title.lower(),
                )
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                if county_token and temporal_ok and likely_settlement and (lat is not None) and (lon is not None) and not false_positive:
                    y0, y1 = infer_temporal_bounds(data.get("temporal", []))
                    e, n = transform_wgs84_to_bng(float(lon), float(lat))
                    candidates.append(
                        {
                            "site_id": f"ADS_{normalize_text(data.get('originalId') or hit.get('id'))}",
                            "site_type": "settlement_open_candidate",
                            "date_from_ad": y0 if y0 is not None else 850,
                            "date_to_ad": y1 if y1 is not None else 954,
                            "precision_m": 1000,
                            "easting": e,
                            "northing": n,
                            "lat": lat,
                            "lon": lon,
                            "source": "ADS",
                            "source_record": row["record_url"],
                            "title": title,
                            "county_token": county_token,
                            "notes": "Open-source candidate; manual archaeological validation required.",
                        }
                    )

            total = int((payload.get("total") or {}).get("value", 0))
            if (page + 1) * 50 >= total:
                break
            page += 1
            time.sleep(0.25)

    return raw_rows, candidates


def infer_temporal_bounds(temporal_items: object) -> Tuple[Optional[int], Optional[int]]:
    years_from: List[int] = []
    years_to: List[int] = []
    if isinstance(temporal_items, list):
        for item in temporal_items:
            if not isinstance(item, dict):
                continue
            f = parse_int(item.get("from"))
            u = parse_int(item.get("until"))
            if f is not None:
                years_from.append(f)
            if u is not None:
                years_to.append(u)
    y0 = min(years_from) if years_from else None
    y1 = max(years_to) if years_to else None
    return y0, y1


def strip_tags(text: str) -> str:
    text = re.sub(r"<[^>]+>", " ", text)
    return normalize_text(html.unescape(text))


def harvest_york_her() -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    all_rows: List[Dict[str, object]] = []
    candidates: List[Dict[str, object]] = []
    seen = set()

    for pn in range(1, 30):
        params = {
            "RecordType": "mon",
            "periodfrom": "EM",
            "periodto": "EM",
            "PeriodIncludeOverlaps": "true",
            "SortOrder": "1",
            "ps": "50",
            "pn": str(pn),
        }
        url = YORK_SEARCH + "?" + urllib.parse.urlencode(params)
        text = fetch_text(url)

        # Parse each record block.
        for block in re.findall(r"<li itemscope=.*?</li>", text, flags=re.S):
            href_m = re.search(r'href="(/Monument/[^"]+)"', block)
            if not href_m:
                continue
            href = href_m.group(1)
            rec_id = href.rsplit("/", 1)[-1]
            if rec_id in seen:
                continue
            seen.add(rec_id)

            wkt_m = re.search(r'data-wkt="([^"]+)"', block)
            wkt = wkt_m.group(1) if wkt_m else ""
            xy_m = re.search(r"POINT\(([-0-9.]+)\s+([-0-9.]+)\)", wkt)
            easting = float(xy_m.group(1)) if xy_m else None
            northing = float(xy_m.group(2)) if xy_m else None

            name_m = re.search(r'<span itemprop="name">\s*(.*?)\s*</span>', block, flags=re.S)
            name = strip_tags(name_m.group(1)) if name_m else ""
            class_m = re.search(r"</span>\s*\(([^)]+)\)\s*</a>", block, flags=re.S)
            rec_class = strip_tags(class_m.group(1)) if class_m else ""
            desc_m = re.search(r'<div itemprop="description">\s*(.*?)\s*</div>', block, flags=re.S)
            desc = strip_tags(desc_m.group(1)) if desc_m else ""

            blob = f"{name} {desc}".lower()
            settlement_hits = any_keyword(blob, SETTLEMENT_KEYWORDS)
            excavation_hits = any_keyword(blob, EXCAVATION_KEYWORDS)
            non_settlement_hits = any_keyword(blob, LIKELY_NON_SETTLEMENT)
            target_terms = any_keyword(blob, TARGET_CULTURE_TERMS)

            likely_settlement = bool(
                settlement_hits
                and (not non_settlement_hits or excavation_hits)
                and target_terms
            )

            row = {
                "source": "York_HER",
                "record_id": rec_id,
                "record_url": f"https://her.york.gov.uk{href}",
                "name": name,
                "record_class": rec_class,
                "description": desc,
                "easting": easting,
                "northing": northing,
                "settlement_hits": ",".join(settlement_hits),
                "excavation_hits": ",".join(excavation_hits),
                "non_settlement_hits": ",".join(non_settlement_hits),
                "target_terms": ",".join(target_terms),
                "likely_settlement": likely_settlement,
            }
            all_rows.append(row)

            if likely_settlement and easting is not None and northing is not None:
                candidates.append(
                    {
                        "site_id": f"YORKHER_{rec_id}",
                        "site_type": "settlement_open_candidate",
                        "date_from_ad": 850,
                        "date_to_ad": 954,
                        "precision_m": 100,
                        "easting": easting,
                        "northing": northing,
                        "lat": "",
                        "lon": "",
                        "source": "York_HER",
                        "source_record": f"https://her.york.gov.uk{href}",
                        "title": name,
                        "county_token": "YORKSHIRE",
                        "notes": "Open-source candidate from York HER EM records; manual archaeological validation required.",
                    }
                )

        # If page has no records we can stop.
        if "search-result-listitem-details" not in text:
            break
        time.sleep(0.2)

    return all_rows, candidates


def parse_float(value: str) -> Optional[float]:
    v = normalize_text(value)
    if not v:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def standardize_local_holdouts() -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []

    # Burials (NW)
    p = INBOX / "Scandinavian_Burials.tsv"
    if p.exists():
        with p.open("r", encoding="utf-8-sig", newline="") as f:
            rdr = csv.DictReader(f, delimiter="\t")
            for i, r in enumerate(rdr, start=1):
                e = parse_float(r.get("Easting", ""))
                n = parse_float(r.get("Northing", ""))
                if e is None or n is None:
                    continue
                rows.append(
                    {
                        "site_id": f"BURIAL_NW_{i}",
                        "site_type": "burial",
                        "date_from_ad": 850,
                        "date_to_ad": 954,
                        "precision_m": 100,
                        "easting": e,
                        "northing": n,
                        "lat": "",
                        "lon": "",
                        "source": "User_Burials_NW",
                        "source_record": normalize_text(r.get("Location")),
                        "title": normalize_text(r.get("Location")),
                        "county_token": normalize_text(r.get("County")),
                        "notes": normalize_text(r.get("Type")),
                    }
                )

    # Burials (Yorkshire)
    p = INBOX / "Scandinavian_Burials_Yorkshire.tsv"
    if p.exists():
        with p.open("r", encoding="utf-8-sig", newline="") as f:
            rdr = csv.DictReader(f, delimiter="\t")
            for i, r in enumerate(rdr, start=1):
                e = parse_float(r.get("Eastings", ""))
                n = parse_float(r.get("Northings", ""))
                if e is None or n is None:
                    continue
                rows.append(
                    {
                        "site_id": f"BURIAL_YKS_{i}",
                        "site_type": "burial",
                        "date_from_ad": 850,
                        "date_to_ad": 954,
                        "precision_m": 100,
                        "easting": e,
                        "northing": n,
                        "lat": "",
                        "lon": "",
                        "source": "User_Burials_Yorkshire",
                        "source_record": normalize_text(r.get("Site_Name")),
                        "title": normalize_text(r.get("Site_Name")),
                        "county_token": normalize_text(r.get("County")),
                        "notes": normalize_text(r.get("Type")),
                    }
                )

    # Hoards
    p = INBOX / "Scandinavian_Hoards.tsv"
    if p.exists():
        with p.open("r", encoding="utf-8-sig", newline="") as f:
            rdr = csv.DictReader(f, delimiter="\t")
            for i, r in enumerate(rdr, start=1):
                e = parse_float(r.get("easting", ""))
                n = parse_float(r.get("northing", ""))
                if e is None or n is None:
                    continue
                y0 = parse_int(r.get("fromdate"))
                y1 = parse_int(r.get("todate"))
                rows.append(
                    {
                        "site_id": f"HOARD_{i}",
                        "site_type": "hoard",
                        "date_from_ad": y0 if y0 is not None else 850,
                        "date_to_ad": y1 if y1 is not None else 954,
                        "precision_m": 100,
                        "easting": e,
                        "northing": n,
                        "lat": "",
                        "lon": "",
                        "source": "User_Hoards",
                        "source_record": normalize_text(r.get("hoard_name")),
                        "title": normalize_text(r.get("hoard_name")),
                        "county_token": normalize_text(r.get("county")),
                        "notes": normalize_text(r.get("category")),
                    }
                )

    return rows


def write_tsv(path: Path, rows: List[Dict[str, object]], field_order: List[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=field_order, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in field_order})


def make_gpkg_from_holdout_tsv(tsv_path: Path, gpkg_path: Path, layer_name: str = "validation_holdout_sites") -> None:
    if gpkg_path.exists():
        gpkg_path.unlink()
    cmd = [
        "ogr2ogr",
        "-f",
        "GPKG",
        str(gpkg_path),
        str(tsv_path),
        "-oo",
        "SEPARATOR=TAB",
        "-oo",
        "X_POSSIBLE_NAMES=easting",
        "-oo",
        "Y_POSSIBLE_NAMES=northing",
        "-a_srs",
        "EPSG:27700",
        "-nln",
        layer_name,
    ]
    subprocess.run(cmd, check=True)


def write_summary(
    ads_raw: List[Dict[str, object]],
    ads_candidates: List[Dict[str, object]],
    york_raw: List[Dict[str, object]],
    york_candidates: List[Dict[str, object]],
    local_holdouts: List[Dict[str, object]],
    final_holdouts: List[Dict[str, object]],
) -> None:
    report = OUTDIR / "OPEN_HARVEST_SUMMARY.md"
    lines = []
    lines.append(f"# Open Harvest Summary ({TODAY})")
    lines.append("")
    lines.append("## Sources harvested")
    lines.append("- ADS Data Catalogue API")
    lines.append("- York HER public monument search (Early Medieval)")
    lines.append("- User-provided local holdouts (burials + hoards)")
    lines.append("")
    lines.append("## Row counts")
    lines.append(f"- ADS raw screened rows: {len(ads_raw)}")
    lines.append(f"- ADS settlement candidates (auto-filtered): {len(ads_candidates)}")
    lines.append(f"- York HER raw EM records: {len(york_raw)}")
    lines.append(f"- York HER settlement candidates (auto-filtered): {len(york_candidates)}")
    lines.append(f"- Local user holdouts (burial + hoard): {len(local_holdouts)}")
    lines.append(f"- Final provisional holdout rows: {len(final_holdouts)}")
    lines.append("")
    lines.append("## Important caveats")
    lines.append("- Open-source settlement candidates require manual archaeological vetting before use as final independent validation.")
    lines.append("- County HER systems outside York do not provide straightforward bulk API exports on public pages; direct HER data request is still required for defensible regional completeness.")
    lines.append("- This provisional holdout is suitable for sensitivity tests, not as final publication-grade independent settlement control.")
    lines.append("")
    report.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_outdir()

    local_holdouts = standardize_local_holdouts()
    ads_raw, ads_candidates = harvest_ads()
    york_raw, york_candidates = harvest_york_her()

    # Merge and deduplicate by site_id.
    merged = {}
    for row in local_holdouts + ads_candidates + york_candidates:
        merged[row["site_id"]] = row
    final_holdouts = list(merged.values())

    # Write source outputs.
    write_tsv(
        OUTDIR / "ADS_raw_screened.tsv",
        ads_raw,
        [
            "source",
            "query_term",
            "ads_hit_id",
            "original_id",
            "title",
            "place_name",
            "county_token",
            "lat",
            "lon",
            "period_names",
            "native_subject",
            "derived_subject",
            "description",
            "settlement_hits",
            "excavation_hits",
            "non_settlement_hits",
            "target_terms",
            "artifact_terms",
            "temporal_ok",
            "likely_settlement",
            "false_positive",
            "record_url",
        ],
    )

    write_tsv(
        OUTDIR / "ADS_open_settlement_candidates.tsv",
        ads_candidates,
        [
            "site_id",
            "site_type",
            "date_from_ad",
            "date_to_ad",
            "precision_m",
            "easting",
            "northing",
            "lat",
            "lon",
            "source",
            "source_record",
            "title",
            "county_token",
            "notes",
        ],
    )

    write_tsv(
        OUTDIR / "YorkHER_early_medieval_records.tsv",
        york_raw,
        [
            "source",
            "record_id",
            "record_url",
            "name",
            "record_class",
            "description",
            "easting",
            "northing",
            "settlement_hits",
            "excavation_hits",
            "non_settlement_hits",
            "target_terms",
            "likely_settlement",
        ],
    )

    write_tsv(
        OUTDIR / "YorkHER_settlement_candidates.tsv",
        york_candidates,
        [
            "site_id",
            "site_type",
            "date_from_ad",
            "date_to_ad",
            "precision_m",
            "easting",
            "northing",
            "lat",
            "lon",
            "source",
            "source_record",
            "title",
            "county_token",
            "notes",
        ],
    )

    final_tsv = OUTDIR / "validation_holdout_sites_open_provisional.tsv"
    write_tsv(
        final_tsv,
        final_holdouts,
        [
            "site_id",
            "site_type",
            "date_from_ad",
            "date_to_ad",
            "precision_m",
            "easting",
            "northing",
            "lat",
            "lon",
            "source",
            "source_record",
            "title",
            "county_token",
            "notes",
        ],
    )

    final_gpkg = OUTDIR / "validation_holdout_sites_open_provisional.gpkg"
    make_gpkg_from_holdout_tsv(final_tsv, final_gpkg, layer_name="validation_holdout_sites")

    write_summary(
        ads_raw=ads_raw,
        ads_candidates=ads_candidates,
        york_raw=york_raw,
        york_candidates=york_candidates,
        local_holdouts=local_holdouts,
        final_holdouts=final_holdouts,
    )

    print(f"Open harvest complete. Outputs in: {OUTDIR}")
    print(f"Final provisional holdout rows: {len(final_holdouts)}")


if __name__ == "__main__":
    main()
