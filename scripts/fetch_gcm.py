#!/usr/bin/env python3
"""Descarga precipitación diaria de modelos CMIP6 (Open-Meteo Climate API)
en el punto de Quinta Normal, para validar la corrección de sesgo de un
GCM — el caso de uso declarado del motor (GCM → escala local).

A diferencia de ERA5 (reanálisis, sincronizado día a día con la realidad),
un GCM produce su propio clima cuyas fechas NO corresponden a los días
meteorológicos reales. Por eso el bias correction de GCM es DISTRIBUCIONAL
y se evalúa con métricas de distribución (KS, sesgo medio, cuantiles), no
con RMSE día a día pareado.

Salida: data/gcm_quinta_normal.csv con una columna por modelo.
"""

import json
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
LAT, LON = -33.445, -70.6828
MODELS = ["MRI_AGCM3_2_S", "EC_Earth3P_HR", "MPI_ESM1_2_XR", "CMCC_CM2_VHR4"]


def main():
    params = urllib.parse.urlencode(
        {
            "latitude": LAT,
            "longitude": LON,
            "start_date": "1950-01-01",
            "end_date": "2018-12-31",
            "models": ",".join(MODELS),
            "daily": "precipitation_sum",
        }
    )
    url = f"https://climate-api.open-meteo.com/v1/climate?{params}"
    print("descargando CMIP6 (4 modelos, 1950-2018) ...", flush=True)
    with urllib.request.urlopen(url, timeout=300) as r:
        d = json.load(r)
    if "error" in d:
        raise SystemExit(d.get("reason", "API error"))

    daily = d["daily"]
    out = pd.DataFrame({"date": daily["time"]})
    for m in MODELS:
        col = f"precipitation_sum_{m}"
        out[m] = daily[col]
    out = out.dropna()
    dest = DATA / "gcm_quinta_normal.csv"
    out.to_csv(dest, index=False, float_format="%.3f")
    print(f"escrito {dest} ({len(out)} días, modelos: {', '.join(MODELS)})")


if __name__ == "__main__":
    main()
