# downscale-rs — Downscaling estadístico y bias correction de clima en Rust

> **Estado:** MVP en desarrollo. Creado 2026-06-10; primer código 2026-06-10.
> Familia de motores Rust del autor: SurtGIS, Hydroflux, Smelt, Anvil, Cantus, Criterium.
> Doc madre: `~/proyectos/ideas-motores-rust.md` (idea C2).

## Qué es
Motor para corrección de sesgo y downscaling estadístico de variables
climáticas (GCM/RCM → escala local), reproducible y rápido.

## El gap que llena
Tienes carpetas `downscaling` y `super-resolution-dem` pero no un motor
operacional. El campo es scripts Python dispersos (xclim, cmethods). Rust lo
hace determinista y batch-friendly para muchos puntos/grillas.

## Forzantes multi-sitio para rainflow (2026-06-15)
`forcing.rs` gana `areal_average(sets, weights)`: promedio areal ponderado
(Thiessen/área/uniforme) de N forzantes → forzante de cuenca, recortado al
período común, mismas variables. CLI `downscale areal --forcing a.csv
--forcing b.csv --weight 0.6 --weight 0.4 -o cuenca.csv`. Verificado e2e:
2 sitios → areal → `rainflow run` GR4J 24.905 pasos. Dos caminos multi-sitio
documentados en forcing-interface.md: (1) cuencas independientes →
`rainflow batch`; (2) areal cuando una cuenca tiene varias estaciones. La
banda de elevación NO usa esto (rainflow hbv-bands deriva bandas con
gradientes TCALT/PCALT desde una forzante). 74 tests lib.

## k-d tree para análogos (2026-06-15)
`analog.rs` ahora usa un k-d tree (búsqueda k-NN exacta, sin deps) en vez de
fuerza bruta. Microbenchmark análogos k=10 predict: 1.14 s → 283 ms (~4× en
el peor caso: predictores uniformes 4-D; mayor con datos reales correlacionados
o k menor). Test `kdtree_knn_matches_bruteforce` + `predict_matches_bruteforce_idw`
garantizan resultados idénticos a la fuerza bruta. API pública sin cambios.
95 tests Rust. docs/performance.md actualizado.

## Rendimiento (2026-06-13) — `docs/performance.md`
Microbenchmarks criterion (`cargo bench`) + comparación cross-language
(`scripts/benchmark_vs_python.py`). Corrección EQM completa del caso Quinta
Normal: **1.5 ms** en downscale-rs vs **16.8 ms cmethods (11×)** y
**233 ms xsdba/xclim (156×)**. QDM 1.9 ms vs 240 ms xsdba (125×). El 11× vs
cmethods es algoritmo-a-algoritmo (NumPy puro); el 125–156× vs xsdba incluye
el overhead real del stack xarray/dask. Demuestra el sello del portafolio
("determinista y rápido") — argumento directo para el paper EMS.

## Métodos v0.2 (2026-06-11)
- [x] **QDM** (`qdm.rs`, Cannon et al. 2015): preserva la señal de cambio
  cuantil a cuantil; CLI `correct --method qdm`. Paridad vs
  xsdba.QuantileDeltaMapping en Quinta Normal: misma fórmula, difiere el
  estimador de `p` (rankdata vs interpolación de nodos); KS holdout
  rust 0.0147 vs xsdba 0.0301 (docs/parity.md).
- [x] **Schaake shuffle** (`multivariate.rs`, Clark et al. 2004):
  restaura dependencia entre variables tras corrección univariada;
  marginales exactas + rangos de la plantilla.
- [x] **PET Hargreaves** (`pet.rs`, FAO-56): Ra de latitud+doy (test vs
  ejemplo 8 FAO), `hargreaves_from_epoch_days` para el eje de forcing.
- [x] **Predictores ricos para análogos** (`docs/predictors.md`,
  2026-06-13): 8 predictores sinópticos de superficie ERA5 (pmsl, rh,
  dewpoint, cloud, vpd, u, v, t2m) + armónicos, experimento vía bindings
  Python. Hallazgos: enriquecer baja RMSE de análogos k=10 de 3.61 a 3.14
  (−13%); análogos k=1 + synoptic da el MEJOR KS de la tabla (0.0058 <
  EQM 0.0179); regresión mejor RMSE (3.03) pero colapsa varianza (KS 0.54).
  Confirma trade-off RMSE↔distribución, accionable. Tabla lista para paper.
- [ ] Pendiente v0.2: MBCn si se necesita multivariado iterativo;
  aplicación por ventanas de QDM.

## Alcance MVP (v0.1) — COMPLETO
- [x] Bias correction: quantile mapping empírico (aditivo/multiplicativo, nodos endpoints/midpoint), paramétrico (normal y gamma mixta con masa en cero), delta change, adaptación de umbral seco/húmedo.
- [x] Downscaling por análogos (k-NN estandarizado, media ponderada por distancia inversa) y regresión lineal múltiple (OLS).
- [x] Validación: split temporal, métricas (RMSE, KS, sesgo de cuantiles) con líneas base raw vs corregido.
- [ ] (v0.2) QDM (quantile delta mapping), multivariado.

## Arquitectura tentativa
- `downscale-core`: mapeos de cuantiles, interpolación, validación.
- Targets: native (Rayon) + Python (PyO3) + CLI.
- I/O: series CSV / NetCDF mínimo / GeoTIFF (vía SurtGIS).

## Validación / paridad numérica
Cross-check contra **cmethods**/xclim sobre estaciones DGA/CR2 (Chile).

## Venue objetivo
**Environmental Modelling & Software** o revista de clima aplicado.

## Conexiones con tu ecosistema
- **rainflow / Hydroflux**: provee forzantes corregidas a los modelos hidro.
- **Smelt**: downscaling vía ML como variante.

## Estado del código (2026-06-10)
Workspace Cargo (edition 2024, resolver 3) con `crates/downscale-core` (sin I/O,
solo `thiserror`):
- `qm.rs` — `QuantileMapping::fit/fit_with_nodes/apply` (EQM, cuantiles tipo 7,
  `NodePlacement::{Endpoints,Midpoint}`, extrapolación aditiva/multiplicativa).
- `parametric.rs` — `ParametricQuantileMapping` normal y gamma mixta con masa
  en cero (MLE gamma vía Thom+Newton; corrige frecuencia de días húmedos).
- `special.rs` — ln-gamma (Lanczos), digamma/trigamma, gamma incompleta P y
  P⁻¹ (NR), normal CDF/PPF (Acklam). Verificadas contra SciPy en tests.
- `delta.rs` — `DeltaChange` aditivo/multiplicativo (perturbación de obs).
- `wetday.rs` — `WetDayCorrection` (adaptación de umbral, Themeßl 2012).
- `analog.rs` — `AnalogDownscaling` (k-NN z-score, euclidiana, media
  ponderada 1/d; matriz aplanada fila-por-día).
- `regression.rs` — `LinearDownscaling` (OLS, ecuaciones normales con
  pivoteo parcial, R² de calibración).
- `metrics.rs` — `rmse`, `mean_bias`, `ks_statistic` (2 muestras),
  `quantile_bias` (deciles + P5/P95/P99).
- `validation.rs` — `split_temporal`, `validate_split[_with]` + `QmOptions`
  → `ValidationReport` con métricas corregido vs raw.
- `forcing.rs` — interfaz de forzantes hacia rainflow: `Variable`
  (pr/pet/tmean canónicos), `ForcingSeries` (eje diario contiguo validado,
  fechas civiles Hinnant sin chrono), `ForcingSet` (alineación a período
  común + CSV canónico). Contrato: `docs/forcing-interface.md`.
- 63 tests (unitarios + integración + doctests), clippy default limpio.
- Ejemplo `quinta_normal_downscaling` (cargo example): compara raw/EQM/
  análogos/regresión en el holdout real. Resultado: EQM domina con un solo
  predictor (KS 0.018); regresión gana RMSE (3.05) pero colapsa varianza
  (KS 0.63); análogos k=1 preserva distribución (KS 0.054) y k=10 minimiza
  error puntual — trade-off clásico, documentado para el paper. Análogos/
  regresión rinden con predictores ricos (campos de presión/humedad), v0.2.

`crates/downscale-cli` (binario `downscale`; deps clap + anyhow):
- `validate`: parea dos CSV `fecha,valor` por fecha y reporta
  corregido-vs-raw (RMSE, sesgo medio, KS, sesgo por cuantil).
- `correct`: calibra obs+model y escribe CSV corregido (`--target` opcional).
- `forcing`: ensambla pr+pet[+tmean] corregidos en el CSV ancho que
  consume rainflow, validando contigüidad y alineando al período común.
- Flags `--nodes endpoints|midpoint` y `--wet-day-threshold <mm>`.
- Lector CSV tolerante (headers, NA, centinela -9999 de DGA/CR2).

Integración verificada (2026-06-11): ERA5 corregido (Quinta Normal) →
`downscale forcing` → `rainflow run` GR4J, 24.905 pasos sin errores.

`crates/downscale-python` (módulo `downscale_rs`, PyO3 0.23 + numpy 0.23,
maturin; convención SurtGIS): clases espejo de la API Rust con
numpy-in/numpy-out (QuantileMapping, ParametricQuantileMapping,
DeltaChange, WetDayCorrection, AnalogDownscaling con predictores 2D,
LinearDownscaling) + funciones (rmse, mean_bias, ks_statistic,
quantile_bias, validate_split → dict). `[lib] test = false` (los tests
viven en Python: `tests/test_downscale.py`, 9 tests incl. paridad
bindings ≡ CLI en Quinta Normal). Build: `maturin build --release` +
pip install del wheel.

Hallazgo en Quinta Normal: el sesgo residual P95 (+2.0 mm) NO es llovizna
(la adaptación de umbral no lo cambia) sino **no-estacionariedad**: la
calibración 1950–1997 es más lluviosa que la validación 1997–2018
(megasequía). Relevante para la discusión del paper.

## Prueba sobre datos reales (2026-06-10)
Quinta Normal 330020 (CR2 prDaily) vs ERA5 puntual (Open-Meteo), 24.900 días
pareados 1950–2018, QM multiplicativo, split 70/30. Resultado: KS 0.149 → 0.018,
sesgo medio 0.572 → 0.052 mm/día, RMSE −17.8%. ERA5 muestra el sesgo húmedo
clásico (llovizna excesiva); QM lo remueve. Sesgo residual en P95 (+2.0 mm)
→ motiva corrección de frecuencia de días húmedos en v0.2.
Datos y procedencia: `data/README.md` (no versionados).

## Paridad numérica (2026-06-10) — ver `docs/parity.md`
Cross-check EQM multiplicativo (100 cuantiles, split 70/30) sobre Quinta
Normal vs ERA5, reproducible con `scripts/parity_quinta_normal.py`:
- **Rust ≡ réplica NumPy del algoritmo** (Δmax 5e-4 = redondeo CSV).
- **Rust ≈ xsdba** (KS idéntico 0.0179; P99 |Δ| ≤ 0.28 mm; diferencia de
  cola por nodos `i/(n-1)` vs `(i+0.5)/n` — candidato `NodePlacement` v0.2).
- **cmethods diverge**: su CDF binned por histograma no corrige frecuencia
  de días húmedos (KS queda en 0.1486 = raw). Argumento para el paper.

## Repo (2026-06-11)
Git local con commit raíz c825282 (29 archivos, 4.402 líneas). CI GitHub
Actions (`.github/workflows/ci.yml`): fmt --check, clippy -D warnings,
test, build release — secuencia verificada localmente. LICENSE-MIT +
LICENSE-APACHE, README.md en inglés (cara pública / paper EMS).
**Sin remote aún** — crear con `gh repo create franciscoparrao/downscale-rs`.

## Pulido técnico (2026-06-14)
- [x] CLI `analog`/`regress`: downscaling por análogos/regresión sobre
  matriz de predictores (`date,col1,...`); lector `Matrix` + `pair_matrix_series`
  en series.rs; flag `--non-negative` para regresión de precipitación.
- [x] Property tests (`tests/properties.rs`, proptest): QM monótono, KS en
  [0,1] y auto-cero, sesgo medio antisimétrico, Schaake preserva marginales,
  parse/format de fecha roundtrip. 93 tests Rust en total.
- [x] Job `python` en CI: maturin build + pip install wheel + pytest
  (12 tests Python ahora en CI; el de datos reales hace pytest.skip).

## Validación con GCMs reales CMIP6 (2026-06-15) — `docs/gcm-validation.md`
4 GCMs CMIP6 (MRI-AGCM3-2-S, EC-Earth3P-HR, MPI-ESM1-2-XR, CMCC-CM2-VHR4)
vía Open-Meteo Climate API. Cierra la brecha: el motor probado con su caso
de uso declarado (GCM→local, no solo reanálisis). **Metodología
distribucional** (un GCM no asimila datos → no se parea día a día → solo
KS/sesgo/frecuencia, no RMSE; split por período 1950-99/2000-18). Hallazgos:
(1) drizzle bias universal — todos sobreestiman días húmedos (14-25% vs 10%
obs), el QM corrige a ~10%; (2) el motor corrige los 4 GCMs tan bien como
ERA5 (KS val 0.008-0.018); (3) sesgo corregido levemente negativo = misma
no-estacionariedad (megasequía post-2010). Limitación honesta: los CMIP6 de
Open-Meteo (HighResMIP) ya vienen downscaled ~10km → sesgo comparable a ERA5,
no el de un GCM crudo (ESGF, pendiente). scripts/fetch_gcm.py + experiment_gcm.py.

## WASM (2026-06-14) — patrón multi-target COMPLETO
`crates/downscale-wasm` (wasm-bindgen, convención surtgis-wasm): expone
QuantileMapping, QuantileDeltaMapping y métricas (rmse/ksStatistic/meanBias)
a JS; series como Float64Array. Binario ~74 KB. Demo en navegador
`www/index.html` (genera clima sintético + sesgo, corrige con EQM, grafica
las 3 CDFs empíricas + métricas raw vs corregido en canvas puro, sin libs).
Verificada por screenshot: KS 0.287→0.001, sesgo 2.98→0.00. Build:
`wasm-pack build --target web --out-dir pkg crates/downscale-wasm` (pkg/
gitignoreado). CI: job `wasm` con wasm-pack. **Superficies: core+CLI+Python+WASM.**

## Próximos pasos al retomar
1. **Software paper EMS** (próximo gran hito): outline calibrado con
   `/paper-review-ems`; figuras (paridad vs xsdba, caso Quinta Normal,
   tabla de predictores, performance, cadena hazard, demo WASM); draft.
   Todo el material ya existe en docs/.
2. Pulido restante menor: paridad QM paramétrico vs xsdba/scipy; actualizar
   GH actions a Node24.
3. v0.2 extensiones: MBCn multivariado iterativo; QDM por ventanas;
   forzantes multi-sitio para rainflow semi-distribuido; validación con GCM
   real CMIP6 (Open-Meteo Climate API).
