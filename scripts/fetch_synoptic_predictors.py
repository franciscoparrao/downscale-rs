#!/usr/bin/env python3
"""Descarga predictores sinópticos de superficie (ERA5 vía Open-Meteo) para
Quinta Normal y los agrega a media diaria, alineados con la estación.

Predictores de gran escala (perfect-prognosis) que modulan la precipitación
en Santiago: presión a nivel del mar (circulación sinóptica), humedad y
punto de rocío, déficit de presión de vapor, nubosidad, viento (u, v) y
temperatura. El viento se descompone a componentes u/v ANTES de promediar
(la dirección es circular y no se puede promediar directamente).

Salida: data/predictors_quinta_normal.csv con columnas
    date, pmsl, rh, dewpoint, cloud, vpd, u, v, t2m
agregadas a media diaria. La precipitación corregida y la observada se
toman de los CSV existentes en el ejemplo Rust.

Requisitos: pandas, numpy. Periodo 1950-2018 (mismo que el caso de paridad).
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
import urllib.request
import urllib.parse
import json

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
LAT, LON = -33.445, -70.6828
START, END = 1950, 2018
CHUNK_YEARS = 5

HOURLY_VARS = [
    "pressure_msl",
    "relative_humidity_2m",
    "dew_point_2m",
    "cloud_cover",
    "vapour_pressure_deficit",
    "wind_speed_10m",
    "wind_direction_10m",
    "temperature_2m",
]


def fetch_chunk(y0: int, y1: int) -> pd.DataFrame:
    params = urllib.parse.urlencode(
        {
            "latitude": LAT,
            "longitude": LON,
            "start_date": f"{y0}-01-01",
            "end_date": f"{y1}-12-31",
            "hourly": ",".join(HOURLY_VARS),
            "timezone": "GMT",
        }
    )
    url = f"https://archive-api.open-meteo.com/v1/archive?{params}"
    with urllib.request.urlopen(url, timeout=300) as r:
        d = json.load(r)
    if "error" in d:
        raise RuntimeError(d.get("reason", "API error"))
    h = d["hourly"]
    df = pd.DataFrame(h)
    df["time"] = pd.to_datetime(df["time"])
    return df


def main():
    chunks = []
    y = START
    while y <= END:
        y1 = min(y + CHUNK_YEARS - 1, END)
        print(f"  descargando {y}-{y1} ...", flush=True)
        for attempt in range(3):
            try:
                chunks.append(fetch_chunk(y, y1))
                break
            except Exception as e:
                print(f"    intento {attempt + 1} falló: {e}", flush=True)
                time.sleep(10)
        else:
            raise SystemExit(f"no se pudo descargar {y}-{y1}")
        time.sleep(3)  # cortesía con la API
        y = y1 + 1

    df = pd.concat(chunks, ignore_index=True)
    df = df.set_index("time").sort_index()

    # Descomponer viento a u/v (convención meteorológica: dir = de donde
    # viene). u: oeste->este positivo; v: sur->norte positivo.
    spd = df["wind_speed_10m"].to_numpy()
    rad = np.deg2rad(df["wind_direction_10m"].to_numpy())
    df["u"] = -spd * np.sin(rad)
    df["v"] = -spd * np.cos(rad)

    # Media diaria de cada predictor.
    daily = df.resample("1D").mean(numeric_only=True)
    out = pd.DataFrame(
        {
            "date": daily.index.strftime("%Y-%m-%d"),
            "pmsl": daily["pressure_msl"],
            "rh": daily["relative_humidity_2m"],
            "dewpoint": daily["dew_point_2m"],
            "cloud": daily["cloud_cover"],
            "vpd": daily["vapour_pressure_deficit"],
            "u": daily["u"],
            "v": daily["v"],
            "t2m": daily["temperature_2m"],
        }
    ).reset_index(drop=True)
    out = out.dropna()
    dest = DATA / "predictors_quinta_normal.csv"
    out.to_csv(dest, index=False, float_format="%.4f")
    print(f"escrito {dest} ({len(out)} días, {out.shape[1] - 1} predictores)")


if __name__ == "__main__":
    main()
