#!/usr/bin/env python3
"""Benchmark cross-language: downscale-rs (bindings) vs xsdba vs cmethods.

Mide el wall-time de una corrección completa (calibrar + aplicar) de EQM y
QDM sobre el dataset Quinta Normal (data/parity/, ~17 430 cal / 7 470 val),
con 100 cuantiles, multiplicativo. Para xsdba/xclim se fuerza el cómputo
(`.values`) porque son perezosos (dask). Reporta la mediana de N repeticiones
y el speedup de Rust.

Requisitos: wheel downscale_rs + xsdba + python-cmethods en el venv, y
data/parity/ generado por scripts/parity_quinta_normal.py.
"""

import statistics as stats
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import downscale_rs as ds

ROOT = Path(__file__).resolve().parent.parent
PARITY = ROOT / "data" / "parity"
NQ = 100
REPS = 50


def load(name):
    return pd.read_csv(PARITY / f"{name}.csv")["value"].to_numpy()


def as_da(values, name="pr"):
    n = len(values)
    time_axis = pd.date_range("1950-01-01", periods=n, freq="D")
    da = xr.DataArray(values, dims="time", coords={"time": time_axis}, name=name)
    da.attrs["units"] = "mm/d"
    return da


def timeit(fn, reps=REPS, warmup=3):
    for _ in range(warmup):
        fn()
    samples = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - t0)
    return stats.median(samples)


def main():
    obs_cal, mod_cal = load("obs_cal"), load("model_cal")
    mod_val = load("model_val")
    ref, hist, sim = as_da(obs_cal), as_da(mod_cal), as_da(mod_val)

    import xsdba
    from cmethods import adjust

    def rust_eqm():
        ds.QuantileMapping(obs_cal, mod_cal, n_quantiles=NQ, kind="mult").apply(mod_val)

    def xsdba_eqm():
        eqm = xsdba.EmpiricalQuantileMapping.train(ref, hist, nquantiles=NQ, kind="*", group="time")
        eqm.adjust(sim, interp="linear", extrapolation="constant").values

    def cmethods_eqm():
        adjust(method="quantile_mapping", obs=ref, simh=hist, simp=sim, n_quantiles=NQ, kind="*")

    def rust_qdm():
        ds.QuantileDeltaMapping(obs_cal, mod_cal, n_quantiles=NQ, kind="mult").apply(mod_val)

    def xsdba_qdm():
        qdm = xsdba.QuantileDeltaMapping.train(ref, hist, nquantiles=NQ, kind="*", group="time")
        qdm.adjust(sim, interp="linear").values

    bench = {
        ("EQM", "downscale-rs"): rust_eqm,
        ("EQM", "xsdba"): xsdba_eqm,
        ("EQM", "cmethods"): cmethods_eqm,
        ("QDM", "downscale-rs"): rust_qdm,
        ("QDM", "xsdba"): xsdba_qdm,
    }
    results = {key: timeit(fn) for key, fn in bench.items()}

    print(f"Corrección completa (calibrar+aplicar), {len(obs_cal)}+{len(mod_val)} días, "
          f"{NQ} cuantiles, mediana de {REPS} reps\n")
    header = f"{'método':<6} {'implementación':<14} {'tiempo':>11} {'speedup':>9}"
    print(header)
    print("-" * len(header))
    lines = [header, "-" * len(header)]
    for method in ("EQM", "QDM"):
        rust_t = results[(method, "downscale-rs")]
        for impl in ("downscale-rs", "xsdba", "cmethods"):
            if (method, impl) not in results:
                continue
            t = results[(method, impl)]
            ms = t * 1e3
            speed = "1.0× (ref)" if impl == "downscale-rs" else f"{t / rust_t:.1f}×"
            line = f"{method:<6} {impl:<14} {ms:>8.3f} ms {speed:>9}"
            print(line)
            lines.append(line)
    return lines


if __name__ == "__main__":
    main()
