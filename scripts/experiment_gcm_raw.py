#!/usr/bin/env python3
"""Corrección de sesgo de GCMs CMIP6 *crudos* (~200–300 km, sin downscaling)
contra la estación Quinta Normal — el caso de uso más exigente del motor.

Metodología distribucional (un GCM no se sincroniza con la realidad: split
por período, no pareo día a día; ver docs/gcm-validation.md). Muestra el
gradiente resolución → sesgo: reanálisis ERA5 (~25 km) y HighResMIP (~10 km)
ya parten relativamente bien; un GCM crudo a ~200 km tiene sesgo estructural
grande porque su celda mezcla océano, valle y cordillera. El quantile
mapping lo corrige en todos los casos.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import downscale_rs as ds

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
WET = 0.1
SPLIT_YEAR = 2000

RES = {"IPSL-CM6A-LR": "209 km", "MPI-ESM1-2-LR": "208 km", "CanESM5": "311 km"}


def ks(a, b):
    av, bv = np.sort(a), np.sort(b)
    allv = np.concatenate([av, bv])
    return np.max(np.abs(np.searchsorted(av, allv, "right") / len(av)
                         - np.searchsorted(bv, allv, "right") / len(bv)))


def split_stats(name, res, mc, mv, oc, ov):
    qm = ds.QuantileMapping(oc, mc, n_quantiles=100, kind="mult", nodes="midpoint")
    corr = qm.apply(mv)
    return (name, res,
            np.mean(mv) * 365.25,            # mm/año crudo
            np.mean(mv >= WET) * 100,         # días húmedos crudo
            ks(mv, ov), ks(corr, ov),
            ds.mean_bias(mv, ov), ds.mean_bias(corr, ov))


def main():
    obs = pd.read_csv(DATA / "quinta_normal_330020_pr.csv", parse_dates=["date"])
    obs["year"] = obs["date"].dt.year
    oc = obs[obs.year < SPLIT_YEAR]["pr_mm"].to_numpy()
    ov = obs[obs.year >= SPLIT_YEAR]["pr_mm"].to_numpy()
    obs_mm = np.mean(ov) * 365.25
    obs_wet = np.mean(ov >= WET) * 100

    raw = pd.read_csv(DATA / "gcm_raw_quinta_normal.csv")
    rows = []
    for m in ["IPSL-CM6A-LR", "MPI-ESM1-2-LR", "CanESM5"]:
        mc = raw[raw.year < SPLIT_YEAR][m].to_numpy()
        mv = raw[raw.year >= SPLIT_YEAR][m].to_numpy()
        rows.append(split_stats(m, RES[m], mc, mv, oc, ov))

    print("GCMs CMIP6 crudos (~200–300 km) → estación Quinta Normal "
          f"(obs: {obs_mm:.0f} mm/año, {obs_wet:.0f}% días húmedos)\n")
    h = (f"{'modelo':<15}{'resol.':>8}{'mm/año':>8}{'húm%':>6}  "
         f"{'KS cru':>7}{'KS cor':>7}  {'sesgo cru':>10}{'sesgo cor':>10}")
    print(h); print("-" * len(h))
    for name, res, mm, wet, kr, kc, br, bc in rows:
        print(f"{name:<15}{res:>8}{mm:>8.0f}{wet:>6.0f}  "
              f"{kr:>7.3f}{kc:>7.3f}  {br:>+10.2f}{bc:>+10.3f}")

    # Gradiente de resolución (familia MPI: crudo LR vs HighResMIP XR vs ERA5).
    print("\nGradiente resolución → sesgo (KS crudo vs estación):")
    era = pd.read_csv(DATA / "era5_quinta_normal_pr.csv", skiprows=3)
    era.columns = ["date", "pr"]; era["date"] = pd.to_datetime(era["date"])
    era["year"] = era["date"].dt.year
    era_v = era[era.year >= SPLIT_YEAR]["pr"].to_numpy()
    hr = pd.read_csv(DATA / "gcm_quinta_normal.csv", parse_dates=["date"])
    hr["year"] = hr["date"].dt.year
    hr_v = hr[hr.year >= SPLIT_YEAR]["MPI_ESM1_2_XR"].to_numpy()
    mpi_raw_v = raw[raw.year >= SPLIT_YEAR]["MPI-ESM1-2-LR"].to_numpy()
    print(f"  ERA5 reanálisis    ~25 km   KS {ks(era_v, ov):.3f}")
    print(f"  MPI HighResMIP-XR  ~10 km   KS {ks(hr_v, ov):.3f}")
    print(f"  MPI crudo (LR)    ~200 km   KS {ks(mpi_raw_v, ov):.3f}")


if __name__ == "__main__":
    main()
