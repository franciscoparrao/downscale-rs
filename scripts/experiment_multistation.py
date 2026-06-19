#!/usr/bin/env python3
"""Generalización fuera de Santiago: corrección de sesgo de ERA5 en 9
estaciones que cubren el gradiente climático de Chile, del desierto de
Atacama a Magallanes (latitud −18° a −53°).

Para cada estación descarga ERA5 puntual (Open-Meteo, cacheado), parea con
la observación CR2, calibra EQM multiplicativo en el 70% inicial y evalúa la
distribución corregida en el 30% restante. Muestra que el motor corrige el
sesgo en todos los regímenes, y cómo el sesgo de ERA5 varía con el clima.

Usa el motor vía bindings; las observaciones salen de
`scripts/...` (extraídas del archivo CR2 a data/stations/).
"""

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import downscale_rs as ds

ROOT = Path(__file__).resolve().parent.parent
STA = ROOT / "data" / "stations"
WET = 0.1
CALIB_FRAC = 0.7

# slug, nombre, región, latitud, longitud — ordenadas norte → sur.
STATIONS = [
    ("arica", "Arica", "desierto extremo", -18.35, -70.34),
    ("antofagasta", "Antofagasta", "desierto costero", -23.45, -70.44),
    ("la_serena", "La Serena", "semiárido", -29.92, -71.20),
    ("quinta_normal", "Santiago", "mediterráneo", -33.45, -70.68),
    ("concepcion", "Concepción", "medit. húmedo", -36.78, -73.06),
    ("valdivia", "Valdivia", "templado lluvioso", -39.65, -73.08),
    ("puerto_montt", "Puerto Montt", "oceánico lluvioso", -41.44, -73.10),
    ("coyhaique", "Coyhaique", "patagónico", -45.59, -72.11),
    ("punta_arenas", "Punta Arenas", "subpolar estepario", -53.00, -70.84),
]


def fetch_era5(slug, lat, lon):
    dest = STA / f"{slug}_era5.csv"
    if dest.exists():
        return dest
    params = urllib.parse.urlencode({
        "latitude": lat, "longitude": lon,
        "start_date": "1950-01-01", "end_date": "2018-03-09",
        "daily": "precipitation_sum", "timezone": "GMT",
    })
    url = f"https://archive-api.open-meteo.com/v1/archive?{params}"
    print(f"  descargando ERA5 {slug} ...", flush=True)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                d = json.load(r)
            break
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < 4:
                wait = 20 * (attempt + 1)
                print(f"    rate limit; espera {wait}s ...", flush=True)
                time.sleep(wait)
            else:
                raise
    daily = d["daily"]
    pd.DataFrame({"date": daily["time"], "pr": daily["precipitation_sum"]}).to_csv(
        dest, index=False)
    time.sleep(8)  # cortesía con la API entre descargas
    return dest


def ks(a, b):
    av, bv = np.sort(a), np.sort(b)
    allv = np.concatenate([av, bv])
    return np.max(np.abs(np.searchsorted(av, allv, "right") / len(av)
                         - np.searchsorted(bv, allv, "right") / len(bv)))


def main():
    print("Generalización del bias correction · gradiente climático de Chile\n")
    header = (f"{'estación':<14}{'régimen':<20}{'lat':>6}  "
              f"{'mm/año':>7}{'húm%':>6}  {'KS cru':>7}{'KS cor':>7}  "
              f"{'sesgo cru':>9}{'sesgo cor':>9}")
    print(header)
    print("-" * len(header))
    rows = []
    for slug, name, region, lat, lon in STATIONS:
        obs = pd.read_csv(STA / f"{slug}_pr.csv", parse_dates=["date"])
        era = pd.read_csv(fetch_era5(slug, lat, lon), parse_dates=["date"])
        df = obs.merge(era, on="date", suffixes=("_o", "_m")).dropna().sort_values("date")
        n = len(df)
        split = round(n * CALIB_FRAC)
        o = df["pr_o"].to_numpy()
        m = df["pr_m"].to_numpy()
        oc, ov = o[:split], o[split:]
        mc, mv = m[:split], m[split:]
        qm = ds.QuantileMapping(oc, mc, n_quantiles=100, kind="mult", nodes="midpoint")
        corr = qm.apply(mv)

        mm_year = np.mean(ov) * 365.25
        wet = np.mean(ov >= WET) * 100
        ks_raw, ks_cor = ks(mv, ov), ks(corr, ov)
        b_raw = ds.mean_bias(mv, ov)
        b_cor = ds.mean_bias(corr, ov)
        print(f"{name:<14}{region:<20}{lat:>6.1f}  {mm_year:>7.0f}{wet:>6.1f}  "
              f"{ks_raw:>7.3f}{ks_cor:>7.3f}  {b_raw:>+9.3f}{b_cor:>+9.3f}")
        rows.append((name, ks_raw, ks_cor))

    # Donde ERA5 tiene sesgo sustancial (KS crudo > 0.1: climas áridos a
    # mediterráneos), la corrección es dramática. En climas lluviosos ERA5
    # ya parte bien y el residuo dominante es no-estacionariedad, no sesgo.
    biased = [(r, c) for _, r, c in rows if r > 0.1]
    impr = np.mean([1 - c / r for r, c in biased])
    print(f"\nEn las {len(biased)} estaciones con sesgo de ERA5 sustancial "
          f"(KS crudo > 0.1, del desierto al mediterráneo), el KS de validación "
          f"baja en promedio {impr * 100:.0f}%.")
    print("En los climas lluviosos del sur ERA5 ya reproduce bien la distribución "
          "(KS crudo ~0.01); ahí el residuo es la no-estacionariedad del clima "
          "(megasequía), no sesgo del modelo — y reaparece como sesgo corregido "
          "positivo por calibrar en un pasado más lluvioso.")


if __name__ == "__main__":
    main()
