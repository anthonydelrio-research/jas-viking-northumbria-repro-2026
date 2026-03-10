# JAS Reproducibility Package (JASC26-162)

Code and reproducibility materials for:

**Del Rio, Anthony Brian.**
*Systematic Confidence-Weighted Integration of Heterogeneous Archaeological Databases: Application to Viking-Age Northumbria*

This repository is prepared for Journal of Archaeological Science reproducibility screening.

## What is included

- `scripts/`:
  - `run_jas_analysis.py` (main analysis pipeline)
  - `derive_sici_sites_register.py`
  - `harvest_open_holdouts.py`
  - `build_canonical_inbox_layers.py`
  - `build_date_rescue_sensitivity.py`
  - `blocker_mitigation_stress_tests.py`
  - `build_figure5_validation_panels.R`
  - `repro_check_report.py` (quick verifier for manuscript headline metrics)
- `data/derived/primary/tables/*.tsv` (primary run derived tables)
- `data/derived/sensitivity/tables/*.tsv` (sensitivity run derived tables)
- `REPRODUCIBILITY_NOTE.md` (original run context)
- `docs/TSV_DATA_DICTIONARY.md` (field-level guide to key derived TSV outputs)
- `CITATION.cff` (preferred citation metadata)
- `requirements-gis.txt` (Python dependency baseline)

## Quick reproducibility check (no restricted source data required)

This validates that packaged derived outputs reproduce manuscript headline metrics.

```bash
python3 scripts/repro_check_report.py
cat docs/REPRO_CHECK_REPORT.md
```

Expected headline values include:
- Primary Hit@15% = `0.24`
- Primary AUC = `0.617026`
- Primary opportunity-null AUC = `0.531501`
- Sensitivity Hit@15% = `0.222222...`
- Sensitivity AUC = `0.5947518519`

## Full analytical rerun (requires licensed source datasets)

### 1) Environment

Python 3.11+ recommended.

```bash
pip install -r requirements-gis.txt
```

R is required only for `build_figure5_validation_panels.R`.
That script uses base R only (no additional R package installation required).
By default it reads redistributed tables under `data/derived/*` and writes figure outputs to `outputs/figures/`.
If holdout point coordinates are unavailable, panel C is rendered as a documented placeholder while panels A-B are still fully reproducible.

### 2) Required input layers

The main pipeline expects an inbox directory containing the following layers:

- `PAS_full_northumbria_metalwork.gpkg` (layer `PAS_full_northumbria_metalwork`)
- `KEPN_full_northumbria.gpkg` (layer `KEPN_full_northumbria`)
- `CASSS_full_northumbria.gpkg` (layer `CASSS_full_northumbria`)
- `northumbria_boundary.gpkg`
- `roman_roads.gpkg`
- `major_waterways.gpkg`
- `ALC_soil_grades.tif`
- `elevation_dem.tif`
- `landuse_arable.tif`
- `Yorkshire_Polygon.shp`
- `Northwest_Polygon.shp`
- holdout layers (e.g., `validation_holdout_sites_vetted_primary_2026-03-03.gpkg`)

### 3) Run

```bash
export JAS_ROOT=$(pwd)
export JAS_INBOX=/absolute/path/to/inbox_with_layers
export JAS_OUTPUT_DIR=$(pwd)/outputs_run

python3 scripts/run_jas_analysis.py
```

The script is deterministic with fixed seeds.

## Data licensing and access

Raw input datasets are not redistributed here when controlled by third-party licensing/terms:

- PAS (Portable Antiquities Scheme)
- KEPN (EPNS/Institute for Name-Studies)
- CASSS corpus data
- OS-derived and other licensed geospatial layers

This repository therefore provides:
- full analysis code,
- exact derived tables used for manuscript metrics,
- explicit rerun requirements and file contracts.

## Contact

Anthony Brian Del Rio  
ORCID: `0009-0002-3799-8225`  
Email: `contact@anthonydelrio.com`
