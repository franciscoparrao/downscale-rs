# downscale-rs

Statistical downscaling and bias correction of climate variables
(GCM/RCM/reanalysis → local scale) in Rust: deterministic, fast,
single-binary, no GDAL/Python chain.

## Methods (v0.1)

**Bias correction**
- Empirical quantile mapping (EQM), additive and multiplicative, with
  endpoint or midpoint quantile nodes (the latter matches xclim/xsdba).
- Parametric quantile mapping: normal (temperature) and zero-inflated
  gamma mixture (precipitation — corrects wet-day frequency and intensity
  simultaneously; gamma fitted by maximum likelihood).
- Delta change (additive/multiplicative perturbation of observations).
- Wet-day frequency correction by threshold adaptation (Themeßl et al. 2012).

**Downscaling**
- Analog method: standardized k-NN over large-scale predictors,
  inverse-distance-weighted mean (k=1 → pure resampling).
- Multiple linear regression (OLS with calibration R²).

**Validation**
- Temporal split (calibration/holdout) with raw-model baselines: RMSE,
  mean bias, two-sample Kolmogorov–Smirnov, per-quantile bias.

**Hydrological forcing interface**
- Validated contiguous daily axis, multi-variable alignment, and the
  canonical wide CSV (`date,pr,pet[,tmean]`) consumed by
  [rainflow](../rainflow) (GR4J/HBV). See `docs/forcing-interface.md`.

## Numerical parity

Cross-checked against **xclim/xsdba** and a pure-NumPy replication on a
real case (DMC station Quinta Normal 330020 from the CR2 daily
precipitation dataset vs point-sampled ERA5, 24 900 paired days
1950–2018). With midpoint nodes, max |Δ| vs xsdba is 0.087 mm over the
7 470-day holdout (RMSD 0.012 mm), with identical KS. The comparison also
documents that cmethods' histogram-binned EQM fails to correct wet-day
frequency on zero-inflated precipitation. Details and declared
tolerances: `docs/parity.md`; reproduce with
`scripts/parity_quinta_normal.py`.

## Performance

A full EQM correction (calibrate + apply) over the 24 900-day Quinta Normal
case runs in **1.5 ms** — about **11× faster than cmethods** (pure NumPy,
algorithm-to-algorithm) and **~150× faster than the xclim/xsdba** xarray
flow. QDM is similar (1.9 ms). Details, microbenchmarks and methodology:
`docs/performance.md` (`cargo bench` + `scripts/benchmark_vs_python.py`).

## Layout

- `crates/downscale-core` — the engine. Pure `f64` slices, no I/O, one
  dependency (`thiserror`). Special functions (incomplete gamma, normal
  PPF, digamma) implemented in-crate and tested against SciPy values.
- `crates/downscale-cli` — the `downscale` binary: `validate`, `correct`,
  `forcing` over `date,value` CSV files, plus `analog`/`regress`
  downscaling over `date,col1,col2,...` predictor matrices (tolerant
  reader: headers, NA, DGA/CR2 `-9999` sentinels).
- `crates/downscale-python` — PyO3 bindings (`downscale_rs` module):
  numpy-in/numpy-out for every method, built with maturin.
- `crates/downscale-wasm` — WebAssembly bindings (~74 KB) with a
  browser demo: bias correction with no server. See its README.

## Quick start

```bash
cargo build --release

# Validate quantile mapping with a 70/30 temporal split
downscale validate --obs station_pr.csv --model era5_pr.csv \
  --kind mult --calib-frac 0.7

# Calibrate and correct a target series
downscale correct --obs station_pr.csv --model era5_pr.csv \
  --kind mult --nodes midpoint -o pr_corrected.csv

# Assemble forcings for a rainfall-runoff model
downscale forcing --pr pr_corrected.csv --pet era5_pet.csv \
  --temp era5_tmean.csv -o forcing.csv
```

Input CSVs are `date,value` with ISO dates; series are paired by date
(distribution-based methods do not require pairing for calibration).

## Python

```bash
cd crates/downscale-python && maturin build --release
pip install ../../target/wheels/downscale_rs-*.whl
```

```python
import numpy as np
import downscale_rs as ds

qm = ds.QuantileMapping(obs, era5, n_quantiles=100, kind="mult", nodes="midpoint")
corrected = qm.apply(era5_future)              # np.ndarray

pqm = ds.ParametricQuantileMapping(obs, era5, dist="gamma", wet_threshold=0.1)
report = ds.validate_split(obs, era5, calib_frac=0.7, kind="mult")
print(report["ks"], report["ks_raw"])          # corrected vs raw baseline
```

Classes mirror the Rust API: `QuantileMapping`, `ParametricQuantileMapping`,
`DeltaChange`, `WetDayCorrection`, `AnalogDownscaling` (2-D predictor
arrays), `LinearDownscaling`; functions `rmse`, `mean_bias`,
`ks_statistic`, `quantile_bias`, `validate_split`. Tests:
`crates/downscale-python/tests/test_downscale.py`.

## Data for the test case

Not versioned; provenance and regeneration commands in `data/README.md`
(CR2 `cr2_prDaily_2018` + Open-Meteo ERA5 archive API).

## License

MIT OR Apache-2.0, at your option.
