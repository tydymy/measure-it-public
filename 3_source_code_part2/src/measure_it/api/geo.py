"""Simplified county / state boundaries as GeoJSON for display (API `/geojson/{level}` and the dashboard maps).

    boundaries_geojson("county")            -> FeatureCollection (dict), cached per (level, tolerance, state)
    boundaries_geojson_bytes("state")       -> the same, serialised once (what the API returns)

Source: `data/processed/{county,state}_boundaries_2024.geoparquet` (Census 2024 cartographic boundaries, 1:5M, built by
`measure_it.geography.crosswalk`; SOURCE_REGISTRY id `census_geography`). Geometries are simplified per feature with
shapely (topology preserved within a feature, not across neighbours) and coordinates snapped to a 0.001-degree grid, so
neighbouring polygons can leave hairline gaps. They are for drawing maps only: never use them for spatial joins or
point-in-polygon assignment (the pipeline uses the unsimplified boundaries for that).

Every feature has ``id`` = FIPS (state 2 chars, county 5 chars) and properties ``geo_id``, ``name``, ``state_abbr``,
``object_id`` (``geo:<fips>``) and ``in_ranking_universe`` (50 states + DC, the scoring universe). The collection carries
a ``metadata`` foreign member with the source, simplification and caveat.
"""
from __future__ import annotations

import hashlib
import json
from functools import lru_cache

from ..config import PROCESSED

LEVELS = ("county", "state")
DEFAULT_TOLERANCE_DEG = {"county": 0.01, "state": 0.02}
GRID_DEG = 0.001
MAX_TOLERANCE_DEG = 0.2
NON_STATE_FIPS = {"60", "66", "69", "72", "78"}   # territories: outside the 50 states + DC ranking universe
SOURCE = {"source_id": "census_geography",
          "name": "Census cartographic boundaries (2024, 1:5,000,000)",
          "tables": ["county_boundaries_2024.geoparquet", "state_boundaries_2024.geoparquet"],
          "producer": "measure_it.api.geo"}
CAVEAT = ("Display geometry only: simplified per feature and snapped to a 0.001-degree grid; neighbouring polygons can "
          "leave hairline gaps. Not for spatial joins.")


def _path(level: str):
    return PROCESSED / f"{level}_boundaries_2024.geoparquet"


def _mtime(level: str) -> float:
    p = _path(level)
    return p.stat().st_mtime if p.exists() else 0.0


def _round(coords, nd: int = 3):
    if isinstance(coords, (list, tuple)):
        if coords and isinstance(coords[0], (int, float)):
            return [round(float(c), nd) for c in coords]
        return [_round(c, nd) for c in coords]
    return coords


@lru_cache(maxsize=16)
def _build(level: str, tolerance: float, state: str | None, mtime: float) -> tuple[dict, bytes, str]:
    import geopandas as gpd
    import shapely
    from shapely.geometry import mapping

    g = gpd.read_parquet(_path(level))
    if g.crs is not None and g.crs.to_epsg() != 4326:
        g = g.to_crs(4326)   # NAD83 -> WGS84 is a sub-metre shift; irrelevant at display precision
    g = g.sort_values("geo_id").reset_index(drop=True)
    if state:
        g = g[g["geo_id"].str[:2] == state].reset_index(drop=True)
    geoms = shapely.simplify(g.geometry.values, tolerance, preserve_topology=True) if tolerance > 0 else \
        g.geometry.values
    geoms = shapely.set_precision(geoms, GRID_DEG)
    feats = []
    for gid, name, st, geom in zip(g["geo_id"], g["name"], g["state_abbr"], geoms):
        if geom is None or geom.is_empty:
            continue
        m = mapping(geom)
        feats.append({"type": "Feature", "id": gid,
                      "properties": {"geo_id": gid, "name": name, "state_abbr": st, "object_id": f"geo:{gid}",
                                     "in_ranking_universe": gid[:2] not in NON_STATE_FIPS},
                      "geometry": {"type": m["type"], "coordinates": _round(m["coordinates"])}})
    fc = {"type": "FeatureCollection", "features": feats,
          "metadata": {"level": level, "n_features": len(feats), "state_filter": state,
                       "simplify_tolerance_deg": tolerance, "coordinate_grid_deg": GRID_DEG, "crs": "EPSG:4326",
                       "source": SOURCE, "caveat": CAVEAT}}
    raw = json.dumps(fc, separators=(",", ":")).encode()
    etag = hashlib.sha1(raw).hexdigest()[:20]
    return fc, raw, etag


def _args(level: str, tolerance: float | None, state: str | None) -> tuple[str, float, str | None]:
    if level not in LEVELS:
        raise ValueError(f"level must be one of {LEVELS}")
    tol = DEFAULT_TOLERANCE_DEG[level] if tolerance is None else float(tolerance)
    if not 0 <= tol <= MAX_TOLERANCE_DEG:
        raise ValueError(f"tolerance must be between 0 and {MAX_TOLERANCE_DEG} degrees")
    st = None
    if state is not None and str(state).strip():
        st = str(state).strip().zfill(2)
        if not (len(st) == 2 and st.isdigit()):
            raise ValueError("state must be a 2-digit state FIPS code, e.g. '06'")
    return level, round(tol, 6), st


def boundaries_geojson(level: str = "county", tolerance: float | None = None, state: str | None = None) -> dict:
    """Simplified boundaries as a GeoJSON FeatureCollection (dict; do not mutate: it is cached)."""
    lv, tol, st = _args(level, tolerance, state)
    return _build(lv, tol, st, _mtime(lv))[0]


def boundaries_geojson_bytes(level: str = "county", tolerance: float | None = None,
                             state: str | None = None) -> tuple[bytes, str]:
    """(serialised FeatureCollection, etag) - serialised once per argument set."""
    lv, tol, st = _args(level, tolerance, state)
    _, raw, etag = _build(lv, tol, st, _mtime(lv))
    return raw, etag
