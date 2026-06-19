#!/usr/bin/env python3
"""Descarga precipitación diaria de GCMs CMIP6 *crudos* (sin downscaling) en
el punto de grilla más cercano a Quinta Normal, desde el archivo público
Pangeo CMIP6 en Google Cloud (zarr, sin autenticación).

A diferencia de los productos HighResMIP de Open-Meteo (ya downscaled a
~10 km, ver docs/gcm-validation.md), estos son los GCMs originales a
~200–250 km: la celda que contiene Santiago promedia océano Pacífico, valle
central y cordillera, de modo que su sesgo respecto a una estación puntual
es estructural y grande — el caso de uso más exigente del motor.

Salida: data/gcm_raw_quinta_normal.csv (una columna por modelo) + la
resolución de grilla de cada uno.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
LAT, LON = -33.45, -70.68

MODELS = {
    "IPSL-CM6A-LR": "gs://cmip6/CMIP6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/day/pr/gr/v20180803/",
    "MPI-ESM1-2-LR": "gs://cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/day/pr/gn/v20190710/",
    "CanESM5": "gs://cmip6/CMIP6/CMIP/CCCma/CanESM5/historical/r1i1p1f1/day/pr/gn/v20190429/",
}


def extract(store):
    ds = xr.open_zarr(store, consolidated=True, storage_options={"token": "anon"})
    # Longitud del GCM suele ir 0–360.
    lon_q = LON % 360 if float(ds.lon.max()) > 180 else LON
    cell = ds.pr.sel(lat=LAT, lon=lon_q, method="nearest")
    cell = cell.sel(time=slice("1950-01-01", "2014-12-31")).load()
    # Resolución de grilla aproximada (grados → km en el ecuador).
    dlat = float(np.abs(np.diff(ds.lat.values)).mean())
    dlon = float(np.abs(np.diff(ds.lon.values)).mean())
    res_km = 0.5 * (dlat + dlon) * 111.0
    years = cell["time"].dt.year.values
    pr_mm = np.asarray(cell.values) * 86400.0  # kg m⁻² s⁻¹ → mm/día
    return years, pr_mm, res_km, float(cell.lat), float(cell.lon)


def main():
    out = {}
    meta = []
    for name, store in MODELS.items():
        print(f"  {name} ...", flush=True)
        years, pr, res, glat, glon = extract(store)
        out[name] = pr
        out.setdefault("year", years)
        meta.append((name, res, glat, glon, len(pr)))
    # Las grillas difieren entre modelos → guardamos por modelo con su año.
    n = min(len(v) for v in out.values())
    df = pd.DataFrame({k: v[:n] for k, v in out.items()})
    dest = DATA / "gcm_raw_quinta_normal.csv"
    df.to_csv(dest, index=False, float_format="%.4f")
    print(f"\nescrito {dest} ({n} días)")
    print(f"{'modelo':<16}{'resolución':>12}{'celda (lat,lon)':>22}")
    for name, res, glat, glon in [(m[0], m[1], m[2], m[3]) for m in meta]:
        glon_s = glon - 360 if glon > 180 else glon
        print(f"{name:<16}{res:>9.0f} km   ({glat:+.1f}, {glon_s:+.1f})")


if __name__ == "__main__":
    main()
