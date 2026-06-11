#!/usr/bin/env python3
"""Cross-check de paridad: downscale (Rust) vs xsdba vs cmethods.

EQM multiplicativo, 100 cuantiles, split temporal 70/30, sobre
Quinta Normal 330020 (CR2) vs ERA5 puntual (Open-Meteo).

Requisitos: .venv con xsdba, python-cmethods, pandas; datos según
data/README.md; binario `cargo build --release`.
Resultados de referencia: docs/parity.md.
"""

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
PARITY = DATA / "parity"
N_QUANTILES = 100
CALIB_FRAC = 0.7


def build_split():
    PARITY.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(DATA / "quinta_normal_330020_pr.csv", parse_dates=["date"])
    era = pd.read_csv(DATA / "era5_quinta_normal_pr.csv", skiprows=3)
    era.columns = ["date", "pr_mm"]
    era["date"] = pd.to_datetime(era["date"])
    m = obs.merge(era, on="date", suffixes=("_obs", "_mod")).sort_values("date")
    split = round(len(m) * CALIB_FRAC)  # mismo criterio que split_temporal (round)
    cal, val = m.iloc[:split], m.iloc[split:]
    for name, df, col in [
        ("obs_cal", cal, "pr_mm_obs"),
        ("model_cal", cal, "pr_mm_mod"),
        ("model_val", val, "pr_mm_mod"),
        ("obs_val", val, "pr_mm_obs"),
    ]:
        out = df[["date", col]].copy()
        out["date"] = out["date"].dt.strftime("%Y-%m-%d")
        out.columns = ["date", "value"]
        out.to_csv(PARITY / f"{name}.csv", index=False)
    print(f"split: {split} calibración / {len(val)} validación")


def run_rust():
    subprocess.run(
        [
            ROOT / "target/release/downscale", "correct",
            "--obs", PARITY / "obs_cal.csv",
            "--model", PARITY / "model_cal.csv",
            "--target", PARITY / "model_val.csv",
            "--kind", "mult",
            "--quantiles", str(N_QUANTILES),
            "-o", PARITY / "rust_corrected.csv",
        ],
        check=True,
    )


def load_da(name):
    df = pd.read_csv(PARITY / f"{name}.csv", parse_dates=["date"])
    da = xr.DataArray(
        df["value"].values, dims="time", coords={"time": df["date"].values}, name="pr"
    )
    da.attrs["units"] = "mm/d"
    return da


def run_references():
    import xsdba
    from cmethods import adjust

    ref, hist, sim = load_da("obs_cal"), load_da("model_cal"), load_da("model_val")

    eqm = xsdba.EmpiricalQuantileMapping.train(
        ref, hist, nquantiles=N_QUANTILES, kind="*", group="time"
    )
    out = eqm.adjust(sim, interp="linear", extrapolation="constant")
    save(sim, np.asarray(out), "xclim_corrected.csv")

    cm = adjust(
        method="quantile_mapping", obs=ref, simh=hist, simp=sim,
        n_quantiles=N_QUANTILES, kind="*",
    )
    vals = np.asarray(cm["pr"]) if hasattr(cm, "data_vars") else np.asarray(cm)
    save(sim, vals, "cmethods_corrected.csv")


def save(sim, values, name):
    pd.DataFrame({"date": sim.time.values, "value": values}).to_csv(
        PARITY / name, index=False
    )


def ks(a, b):
    av, bv = np.sort(a), np.sort(b)
    allv = np.concatenate([av, bv])
    fa = np.searchsorted(av, allv, side="right") / len(av)
    fb = np.searchsorted(bv, allv, side="right") / len(bv)
    return np.max(np.abs(fa - fb))


def compare():
    col = lambda n: pd.read_csv(PARITY / f"{n}.csv")["value"].values
    rust, xc, cm = col("rust_corrected"), col("xclim_corrected"), col("cmethods_corrected")
    obs, raw = col("obs_val"), col("model_val")

    print("\n== diferencias entre implementaciones ==")
    for name, a, b in [("rust-xsdba", rust, xc), ("rust-cmethods", rust, cm)]:
        d = np.abs(a - b)
        print(
            f"{name:15s} med={np.median(d):.4f} P90={np.percentile(d, 90):.4f} "
            f"P99={np.percentile(d, 99):.4f} rmsd={np.sqrt(np.mean(d**2)):.4f} "
            f"max={d.max():.4f}"
        )

    print("\n== calidad vs obs (holdout) ==")
    for name, a in [("raw", raw), ("rust", rust), ("xsdba", xc), ("cmethods", cm)]:
        print(f"{name:10s} KS={ks(a, obs):.4f} sesgo={np.mean(a) - np.mean(obs):+.4f}")


if __name__ == "__main__":
    build_split()
    run_rust()
    run_references()
    compare()
