#!/usr/bin/env python3
"""Impacto de la corrección de sesgo en el modelado hidrológico aguas abajo
— la cadena downscale-rs → rainflow sobre cuencas CAMELS-CL reales.

Pregunta: si en vez de la precipitación de referencia CR2MET (grillada,
regional, costosa) se fuerza un modelo lluvia-escorrentía con ERA5 cruda
(global, gratis, sesgada), ¿cuánto se degrada? ¿Y cuánto recupera
downscale-rs corrigiendo ERA5 contra CR2MET?

Para cada cuenca pluvial se arman tres forzantes —CR2MET (referencia), ERA5
cruda, ERA5 corregida (QM multiplicativo)— y se calibra GR4J con el
split-sample de Klemeš (rainflow), comparando el KGE de validación.

Requiere: wheel downscale_rs, binario rainflow, CAMELS-CL en rainflow,
acceso a Open-Meteo. ERA5 se cachea en data/basins/.
"""

import json
import re
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import downscale_rs as ds

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "basins"
CAMELS = Path.home() / "proyectos/rainflow/data/camels-cl"
RAINFLOW = Path.home() / "proyectos/rainflow/target/release/rainflow"

# cuencas pluviales near-natural (GR4J puro), gauge lat/lon aproximado.
BASINS = [
    ("8123001", "Río Itata en Cholguán", -37.15, -71.93),
    ("7330001", "Río Perquilauquén en San Manuel", -36.37, -71.65),
]


def fetch_era5(gid, lat, lon):
    dest = OUT / f"{gid}_era5_dl.csv"  # distinto del escenario de forzante "era5"
    if dest.exists():
        return pd.read_csv(dest, parse_dates=["date"])
    params = urllib.parse.urlencode({
        "latitude": lat, "longitude": lon,
        "start_date": "1979-01-01", "end_date": "2016-12-31",
        "daily": "precipitation_sum", "timezone": "GMT",
    })
    url = f"https://archive-api.open-meteo.com/v1/archive?{params}"
    print(f"  ERA5 {gid} ...", flush=True)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                d = json.load(r)
            break
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < 4:
                time.sleep(20 * (attempt + 1))
            else:
                raise
    df = pd.DataFrame({"date": pd.to_datetime(d["daily"]["time"]),
                       "era5": d["daily"]["precipitation_sum"]})
    df.to_csv(dest, index=False)
    time.sleep(8)
    return df


def write_forcing(path, df, pcol):
    out = df[["date", pcol, "pet", "qobs"]].copy()
    out.columns = ["date", "p", "pet", "qobs"]
    out["date"] = out["date"].dt.strftime("%Y-%m-%d")
    out.to_csv(path, index=False)


def split_sample_kge(forcing):
    r = subprocess.run([str(RAINFLOW), "split-sample", "--forcing", str(forcing),
                        "--model", "gr4j", "--objective", "kge"],
                       capture_output=True, text=True, check=True)
    vals = [float(m) for m in re.findall(r"val [AB] .*?:\s*([0-9.]+)", r.stdout)]
    return np.mean(vals) if vals else float("nan")


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    def ks(a, b):
        av, bv = np.sort(a), np.sort(b)
        allv = np.concatenate([av, bv])
        return np.max(np.abs(np.searchsorted(av, allv, "right") / len(av)
                             - np.searchsorted(bv, allv, "right") / len(bv)))

    print("Cadena downscale-rs → rainflow sobre cuencas CAMELS-CL pluviales\n")
    print("(A) Sesgo de la precipitación de cuenca vs CR2MET (KS):")
    print(f"    {'cuenca':<32}{'ERA5 cru':>10}{'ERA5 cor':>10}")
    rows = []
    for gid, name, lat, lon in BASINS:
        cam = pd.read_csv(CAMELS / f"{gid}.csv", parse_dates=["date"])
        era = fetch_era5(gid, lat, lon)
        df = cam.merge(era, on="date").sort_values("date").reset_index(drop=True)
        ref = df["p"].to_numpy()
        wet = df["era5"].to_numpy()
        qm = ds.QuantileMapping(ref, wet, n_quantiles=100, kind="mult", nodes="midpoint")
        df["era5_corr"] = qm.apply(wet)
        ks_raw = ks(wet, ref)
        ks_cor = ks(df["era5_corr"].to_numpy(), ref)
        print(f"    {name:<32}{ks_raw:>10.3f}{ks_cor:>10.3f}")

        kge = {}
        for tag, pcol in [("cr2met", "p"), ("era5", "era5"), ("era5_corr", "era5_corr")]:
            f = OUT / f"{gid}_{tag}.csv"
            write_forcing(f, df, pcol)
            kge[tag] = split_sample_kge(f)
        rows.append((name, kge, ks_raw, ks_cor))

    print("\n(B) Impacto en GR4J (KGE de validación, split-sample de Klemeš):")
    print(f"    {'cuenca':<32}{'CR2MET':>9}{'ERA5 cru':>10}{'ERA5 cor':>10}")
    for name, kge, _, _ in rows:
        print(f"    {name:<32}{kge['cr2met']:>9.3f}{kge['era5']:>10.3f}{kge['era5_corr']:>10.3f}")


if __name__ == "__main__":
    main()
