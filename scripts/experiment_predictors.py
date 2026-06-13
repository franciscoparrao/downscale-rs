#!/usr/bin/env python3
"""Experimento: ¿mejoran análogos/regresión con predictores sinópticos?

Compara métodos de downscaling de precipitación diaria en Quinta Normal
(estación DMC 330020) usando ERA5 como modelo de gran escala, sobre un
split temporal 70/30. Usa el motor Rust vía los bindings `downscale_rs`.

Conjuntos de predictores para análogos/regresión (perfect-prognosis):
  - pr        : precipitación ERA5 (el único predictor del baseline actual)
  - pr+season : pr + armónicos del día del año
  - synoptic  : estado atmosférico de gran escala SIN precipitación
                (presión msl, humedad, punto de rocío, nubosidad, déficit
                 de presión de vapor, viento u/v, temperatura) + armónicos
  - all       : synoptic + pr

EQM y QDM (que solo usan pr) van como referencia de bias correction.

Requisitos: el wheel downscale_rs instalado, data/ con la serie de la
estación, ERA5 pr y predictors_quinta_normal.csv (ver
fetch_synoptic_predictors.py). Salida: tabla por stdout + docs/predictors.md
"""

from pathlib import Path

import numpy as np
import pandas as pd

import downscale_rs as ds

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CALIB_FRAC = 0.7


def load_aligned():
    obs = pd.read_csv(DATA / "quinta_normal_330020_pr.csv", parse_dates=["date"])
    obs = obs.rename(columns={"pr_mm": "obs"})

    era = pd.read_csv(DATA / "era5_quinta_normal_pr.csv", skiprows=3)
    era.columns = ["date", "pr"]
    era["date"] = pd.to_datetime(era["date"])

    syn = pd.read_csv(DATA / "predictors_quinta_normal.csv", parse_dates=["date"])

    df = obs.merge(era, on="date").merge(syn, on="date").sort_values("date")
    df = df.dropna().reset_index(drop=True)
    # Armónicos del día del año.
    doy = df["date"].dt.dayofyear.to_numpy()
    df["sin"] = np.sin(2 * np.pi * doy / 365.25)
    df["cos"] = np.cos(2 * np.pi * doy / 365.25)
    return df


PREDICTOR_SETS = {
    "pr": ["pr"],
    "pr+season": ["pr", "sin", "cos"],
    "synoptic": ["pmsl", "rh", "dewpoint", "cloud", "vpd", "u", "v", "t2m", "sin", "cos"],
    "all": ["pmsl", "rh", "dewpoint", "cloud", "vpd", "u", "v", "t2m", "sin", "cos", "pr"],
}


def metrics(sim, obs):
    return ds.rmse(sim, obs), ds.ks_statistic(sim, obs), ds.mean_bias(sim, obs)


def main():
    df = load_aligned()
    n = len(df)
    split = round(n * CALIB_FRAC)
    cal, val = df.iloc[:split], df.iloc[split:]
    obs_cal = cal["obs"].to_numpy()
    obs_val = val["obs"].to_numpy()
    pr_cal = cal["pr"].to_numpy()
    pr_val = val["pr"].to_numpy()

    print(f"Quinta Normal 330020 vs ERA5 — {n} días pareados "
          f"({df['date'].iloc[0].date()}..{df['date'].iloc[-1].date()})")
    print(f"calibración {split} días / validación {n - split} días\n")

    rows = []
    rows.append(("raw ERA5", "-", *metrics(pr_val, obs_val)))

    eqm = ds.QuantileMapping(obs_cal, pr_cal, n_quantiles=100, kind="mult")
    rows.append(("EQM", "pr", *metrics(eqm.apply(pr_val), obs_val)))

    qdm = ds.QuantileDeltaMapping(obs_cal, pr_cal, n_quantiles=100, kind="mult")
    rows.append(("QDM", "pr", *metrics(qdm.apply(pr_val), obs_val)))

    for set_name, cols in PREDICTOR_SETS.items():
        xcal = cal[cols].to_numpy()
        xval = val[cols].to_numpy()
        for k in (10, 1):
            ad = ds.AnalogDownscaling(xcal, obs_cal, k=k)
            pred = ad.predict(xval)
            rows.append((f"análogos k={k}", set_name, *metrics(pred, obs_val)))
        try:
            lm = ds.LinearDownscaling(xcal, obs_cal)
            pred = np.maximum(lm.predict(xval), 0.0)
            rows.append((f"regresión (R²={lm.r2:.2f})", set_name, *metrics(pred, obs_val)))
        except ValueError as e:
            rows.append(("regresión", set_name, float("nan"), float("nan"), float("nan")))

    header = f"{'método':<18} {'predictores':<12} {'RMSE':>7} {'KS':>8} {'sesgo':>8}"
    print(header)
    print("-" * len(header))
    lines = [header, "-" * len(header)]
    for name, pset, rmse, ks, bias in rows:
        line = f"{name:<18} {pset:<12} {rmse:>7.3f} {ks:>8.4f} {bias:>+8.3f}"
        print(line)
        lines.append(line)

    return lines


if __name__ == "__main__":
    main()
