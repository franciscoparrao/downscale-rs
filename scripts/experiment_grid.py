#!/usr/bin/env python3
"""Corrección de sesgo de un campo grillado real: precipitación de un GCM
CMIP6 crudo (MPI-ESM1-2-LR, ~200 km) sobre Chile, corregida celda por celda
hacia ERA5 — el caso de uso raster del motor.

A diferencia de las series puntuales de los demás experimentos, aquí entra
y sale un campo `[time, lat, lon]`: el núcleo Rust corrige cada celda y el
I/O grillado lo maneja xarray. Genera un mapa de la precipitación media
(GCM crudo / ERA5 / corregido) en `site/grid_correction.png`.

Requiere: wheel downscale_rs + gcsfs, zarr, xarray, matplotlib; acceso a
Pangeo CMIP6 (Google Cloud) y Open-Meteo.
"""

import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import xarray as xr

import downscale_rs as ds

ROOT = Path(__file__).resolve().parent.parent
STORE = "gs://cmip6/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/day/pr/gn/v20190710/"
# ventana Chile central-sur (pocas celdas a ~200 km → bulk ERA5 manejable).
LAT0, LAT1, LON0, LON1 = -42.0, -30.0, -73.5, -70.0


def load_gcm():
    g = xr.open_zarr(STORE, consolidated=True, storage_options={"token": "anon"})
    lon0, lon1 = LON0 % 360, LON1 % 360
    sub = g.pr.sel(lat=slice(LAT0, LAT1), lon=slice(lon0, lon1),
                   time=slice("1980-01-01", "2014-12-31")).load()
    lats = sub.lat.values
    lons = np.where(sub.lon.values > 180, sub.lon.values - 360, sub.lon.values)
    field = np.asarray(sub.values) * 86400.0  # kg/m²/s → mm/d
    return field, lats, lons


def fetch_era5_grid(lats, lons, n_time_hint):
    # Bulk: todos los centroides en una sola request (evita rate limit).
    glat, glon = np.meshgrid(lats, lons, indexing="ij")
    flat_lat = ",".join(f"{v:.4f}" for v in glat.ravel())
    flat_lon = ",".join(f"{v:.4f}" for v in glon.ravel())
    params = urllib.parse.urlencode({
        "latitude": flat_lat, "longitude": flat_lon,
        "start_date": "1980-01-01", "end_date": "2014-12-31",
        "daily": "precipitation_sum", "timezone": "GMT",
    })
    url = f"https://archive-api.open-meteo.com/v1/archive?{params}"
    print("  ERA5 grilla (bulk) ...", flush=True)
    with urllib.request.urlopen(url, timeout=300) as r:
        d = json.load(r)
    locs = d if isinstance(d, list) else [d]
    nt = len(locs[0]["daily"]["precipitation_sum"])
    grid = np.full((nt, len(lats), len(lons)), np.nan)
    for k, loc in enumerate(locs):
        i, j = divmod(k, len(lons))
        grid[:, i, j] = loc["daily"]["precipitation_sum"]
    return grid


def main():
    print("Corrección de un campo grillado: GCM crudo MPI → ERA5, sobre Chile\n")
    gcm, lats, lons = load_gcm()
    print(f"  GCM recorte: {gcm.shape[0]} días × {len(lats)}×{len(lons)} celdas "
          f"(~200 km), lat {lats.min():.1f}..{lats.max():.1f}")
    era = fetch_era5_grid(lats, lons, gcm.shape[0])

    # Corrección celda por celda (distribucional: el GCM no se parea en el tiempo).
    era = np.ascontiguousarray(era, dtype=np.float64)
    gcm = np.ascontiguousarray(gcm, dtype=np.float64)
    corr = ds.correct_grid(era, gcm, n_quantiles=100, kind="mult")

    mean_mm = lambda f: np.nanmean(f, axis=0) * 365.25
    g_map, e_map, c_map = mean_mm(gcm), mean_mm(era), mean_mm(corr)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vmax = np.nanpercentile(np.concatenate([g_map.ravel(), e_map.ravel()]), 98)
    fig, axes = plt.subplots(1, 3, figsize=(11, 5.2), constrained_layout=True)
    ext = [lons.min(), lons.max(), lats.min(), lats.max()]
    for ax, m, t in zip(axes, [g_map, e_map, c_map],
                        ["GCM crudo (MPI, ~200 km)", "ERA5 (referencia)",
                         "GCM corregido (celda a celda)"]):
        im = ax.imshow(m, origin="lower", extent=ext, aspect="auto",
                       cmap="YlGnBu", vmin=0, vmax=vmax)
        ax.set_title(t, fontsize=11)
        ax.set_xlabel("lon")
    axes[0].set_ylabel("lat")
    fig.colorbar(im, ax=axes, shrink=0.8, label="precipitación media (mm/año)")
    fig.suptitle("downscale-rs · corrección de sesgo de un campo grillado",
                 fontsize=13, weight="bold")
    dest = ROOT / "site" / "grid_correction.png"
    fig.savefig(dest, dpi=130)
    print(f"\n  mapa → {dest}")

    # Sesgo areal medio antes/después (vs ERA5, sobre celdas válidas).
    valid = np.isfinite(g_map) & np.isfinite(e_map)
    b_raw = np.nanmean((g_map - e_map)[valid])
    b_cor = np.nanmean((c_map - e_map)[valid])
    print(f"  sesgo areal medio vs ERA5: crudo {b_raw:+.0f} mm/año → "
          f"corregido {b_cor:+.1f} mm/año")


if __name__ == "__main__":
    main()
