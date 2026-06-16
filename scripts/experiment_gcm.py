#!/usr/bin/env python3
"""Validación de la corrección de sesgo de GCMs reales (CMIP6) contra la
estación Quinta Normal — el caso de uso declarado del motor.

Metodología DISTRIBUCIONAL: un GCM no está sincronizado día a día con la
realidad, así que las series no se parean en el tiempo. Se separa por
PERÍODO (calibración 1950–1999 / validación 2000–2018), se calibra el
quantile mapping comparando las DISTRIBUCIONES de cada período, y se evalúa
la distribución corregida con KS, sesgo medio y frecuencia de días húmedos
(no RMSE, que requeriría pareo temporal inexistente).

Compara ERA5 (reanálisis, sesgo suave de llovizna) con 4 GCMs CMIP6 cuyo
sesgo es típicamente más fuerte y estructural. Usa el motor vía bindings.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import downscale_rs as ds

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
WET = 0.1  # umbral día húmedo (mm)
SPLIT_YEAR = 2000


def load():
    obs = pd.read_csv(DATA / "quinta_normal_330020_pr.csv", parse_dates=["date"])
    obs = obs.rename(columns={"pr_mm": "obs"})
    era = pd.read_csv(DATA / "era5_quinta_normal_pr.csv", skiprows=3)
    era.columns = ["date", "ERA5"]
    era["date"] = pd.to_datetime(era["date"])
    gcm = pd.read_csv(DATA / "gcm_quinta_normal.csv", parse_dates=["date"])
    df = obs.merge(era, on="date").merge(gcm, on="date")
    return df


def wet_freq(x):
    return float(np.mean(x >= WET))


def describe(model, obs):
    return {
        "ks": ds.ks_statistic(model, obs),
        "bias": ds.mean_bias(model, obs),
        "ratio": float(np.mean(model) / np.mean(obs)),
        "wet": wet_freq(model),
    }


def main():
    df = load()
    sources = ["ERA5", "MRI_AGCM3_2_S", "EC_Earth3P_HR", "MPI_ESM1_2_XR", "CMCC_CM2_VHR4"]
    year = df["date"].dt.year
    cal, val = df[year < SPLIT_YEAR], df[year >= SPLIT_YEAR]
    obs_cal = cal["obs"].to_numpy()
    obs_val = val["obs"].to_numpy()

    print(f"Quinta Normal — calibración {cal['date'].dt.year.min()}–{SPLIT_YEAR - 1} "
          f"({len(cal)} d) / validación {SPLIT_YEAR}–{val['date'].dt.year.max()} ({len(val)} d)")
    print(f"obs: precipitación media {np.mean(obs_val):.2f} mm/d, "
          f"días húmedos {wet_freq(obs_val) * 100:.1f}%\n")

    header = (f"{'fuente':<15} {'KS crudo':>8} {'KS corr':>8} "
              f"{'sesgo cru':>9} {'sesgo cor':>9} {'húmed cru':>9} {'húmed cor':>9}")
    print(header)
    print("-" * len(header))
    lines = [header, "-" * len(header)]
    obs_wet = wet_freq(obs_val) * 100
    for s in sources:
        mc, mv = cal[s].to_numpy(), val[s].to_numpy()
        raw = describe(mv, obs_val)
        # QM multiplicativo con corrección de frecuencia de días húmedos.
        qm = ds.QuantileMapping(obs_cal, mc, n_quantiles=100, kind="mult", nodes="midpoint")
        corr = qm.apply(mv)
        cor = describe(corr, obs_val)
        line = (f"{s:<15} {raw['ks']:>8.3f} {cor['ks']:>8.3f} "
                f"{raw['bias']:>+9.3f} {cor['bias']:>+9.3f} "
                f"{raw['wet'] * 100:>8.1f}% {cor['wet'] * 100:>8.1f}%")
        print(line)
        lines.append(line)
    print(f"\n(referencia: días húmedos observados = {obs_wet:.1f}%)")
    return lines


if __name__ == "__main__":
    main()
