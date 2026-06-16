#!/usr/bin/env python3
"""Corrección multivariada (MBCn) de precipitación + temperatura en Quinta
Normal, contra ERA5. Muestra lo que la corrección univariada NO hace:
preservar la estructura de dependencia entre variables.

En clima mediterráneo (Santiago) la precipitación y la temperatura
correlacionan negativamente — los días lluviosos son más fríos. El quantile
mapping aplicado a cada variable por separado corrige las marginales pero
deja intacta la correlación del modelo. MBCn (Cannon 2018) corrige ambas.

ERA5 es reanálisis (sincronizado día a día con la realidad), por lo que la
correlación pr–temp día a día es significativa tanto en obs como en modelo
(a diferencia de un GCM). Split temporal 70/30; ambos métodos corrigen hacia
el clima de calibración y se evalúan en validación.

Usa el motor vía bindings: QuantileMapping (univariado) y mbcn.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import downscale_rs as ds

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CALIB_FRAC = 0.7


def load():
    pr_o = pd.read_csv(DATA / "quinta_normal_330020_pr.csv", parse_dates=["date"]).rename(
        columns={"pr_mm": "pr_obs"}
    )
    t_o = pd.read_csv(DATA / "quinta_normal_330020_temp.csv", parse_dates=["date"]).rename(
        columns={"temp": "t_obs"}
    )
    pr_e = pd.read_csv(DATA / "era5_quinta_normal_pr.csv", skiprows=3)
    pr_e.columns = ["date", "pr_mod"]
    pr_e["date"] = pd.to_datetime(pr_e["date"])
    t_e = pd.read_csv(DATA / "era5_quinta_normal_tmean.csv", skiprows=3)
    t_e.columns = ["date", "t_mod"]
    t_e["date"] = pd.to_datetime(t_e["date"])
    df = pr_o.merge(t_o, on="date").merge(pr_e, on="date").merge(t_e, on="date")
    return df.dropna().sort_values("date").reset_index(drop=True)


def corr(m):
    return np.corrcoef(m[:, 0], m[:, 1])[0, 1]


def main():
    df = load()
    n = len(df)
    split = round(n * CALIB_FRAC)
    cal, val = df.iloc[:split], df.iloc[split:]

    obs_cal = cal[["pr_obs", "t_obs"]].to_numpy()
    obs_val = val[["pr_obs", "t_obs"]].to_numpy()
    mod_val = val[["pr_mod", "t_mod"]].to_numpy()

    print(f"Quinta Normal pr+temp vs ERA5 — {n} días pareados "
          f"({df['date'].iloc[0].date()}..{df['date'].iloc[-1].date()})")
    print(f"calibración {split} / validación {n - split}\n")

    # Univariado: QM de cada variable hacia el clima de calibración.
    qm_pr = ds.QuantileMapping(obs_cal[:, 0], cal["pr_mod"].to_numpy(), kind="mult")
    qm_t = ds.QuantileMapping(obs_cal[:, 1], cal["t_mod"].to_numpy(), kind="add")
    uni = np.column_stack([qm_pr.apply(mod_val[:, 0]), qm_t.apply(mod_val[:, 1])])

    # Multivariado: MBCn corrige model_val hacia el clima de calibración.
    mv = ds.mbcn(obs_cal, mod_val, n_iterations=30, seed=42)

    def ks2(a, b):
        av, bv = np.sort(a), np.sort(b)
        allv = np.concatenate([av, bv])
        return np.max(np.abs(np.searchsorted(av, allv, "right") / len(av)
                             - np.searchsorted(bv, allv, "right") / len(bv)))

    print("correlación pr–temp (validación):")
    print(f"  observado      {corr(obs_val):+.3f}")
    print(f"  ERA5 crudo     {corr(mod_val):+.3f}")
    print(f"  QM univariado  {corr(uni):+.3f}   (≈ la del modelo: no corrige la dependencia)")
    print(f"  MBCn           {corr(mv):+.3f}   (≈ la observada: corrige la dependencia)")

    print("\nKS marginal vs observado (ambos corrigen las marginales):")
    for j, name in enumerate(["pr", "temp"]):
        print(f"  {name:5s}  QM univariado {ks2(uni[:, j], obs_val[:, j]):.3f}   "
              f"MBCn {ks2(mv[:, j], obs_val[:, j]):.3f}   "
              f"crudo {ks2(mod_val[:, j], obs_val[:, j]):.3f}")


if __name__ == "__main__":
    main()
