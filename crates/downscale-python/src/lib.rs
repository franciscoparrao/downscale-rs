//! Bindings Python de downscale-rs vía PyO3.
//!
//! Las series son arrays 1D de numpy (f64); los predictores de
//! análogos/regresión, arrays 2D (filas = días, columnas = features).
//! `kind` ∈ {"add", "mult"}, `nodes` ∈ {"endpoints", "midpoint"},
//! `dist` ∈ {"normal", "gamma"}.

use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use downscale_core as core;

fn err(e: core::DownscaleError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn parse_kind(kind: &str) -> PyResult<core::Kind> {
    match kind {
        "add" | "additive" => Ok(core::Kind::Additive),
        "mult" | "multiplicative" => Ok(core::Kind::Multiplicative),
        other => Err(PyValueError::new_err(format!(
            "kind inválido: '{other}' (esperado 'add' o 'mult')"
        ))),
    }
}

fn parse_nodes(nodes: &str) -> PyResult<core::NodePlacement> {
    match nodes {
        "endpoints" => Ok(core::NodePlacement::Endpoints),
        "midpoint" => Ok(core::NodePlacement::Midpoint),
        other => Err(PyValueError::new_err(format!(
            "nodes inválido: '{other}' (esperado 'endpoints' o 'midpoint')"
        ))),
    }
}

/// Aplana un array 2D (filas = días) y devuelve (flat, n_features).
fn flatten_2d(arr: &PyReadonlyArray2<'_, f64>) -> (Vec<f64>, usize) {
    let view = arr.as_array();
    (view.iter().copied().collect(), view.ncols())
}

/// Quantile mapping empírico (EQM).
#[pyclass]
struct QuantileMapping {
    inner: core::QuantileMapping,
}

#[pymethods]
impl QuantileMapping {
    /// Calibra con observaciones y modelo del mismo período.
    #[new]
    #[pyo3(signature = (obs, model, n_quantiles=100, kind="add", nodes="endpoints"))]
    fn new(
        obs: PyReadonlyArray1<'_, f64>,
        model: PyReadonlyArray1<'_, f64>,
        n_quantiles: usize,
        kind: &str,
        nodes: &str,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: core::QuantileMapping::fit_with_nodes(
                obs.as_slice()?,
                model.as_slice()?,
                n_quantiles,
                parse_kind(kind)?,
                parse_nodes(nodes)?,
            )
            .map_err(err)?,
        })
    }

    /// Corrige una serie del modelo.
    fn apply<'py>(
        &self,
        py: Python<'py>,
        series: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(self
            .inner
            .apply(series.as_slice()?)
            .map_err(err)?
            .into_pyarray(py))
    }

    /// Corrige un único valor.
    fn correct_one(&self, x: f64) -> f64 {
        self.inner.correct_one(x)
    }
}

/// Quantile mapping paramétrico: normal o gamma mixta con masa en cero.
#[pyclass]
struct ParametricQuantileMapping {
    inner: core::ParametricQuantileMapping,
}

#[pymethods]
impl ParametricQuantileMapping {
    #[new]
    #[pyo3(signature = (obs, model, dist="normal", wet_threshold=0.1))]
    fn new(
        obs: PyReadonlyArray1<'_, f64>,
        model: PyReadonlyArray1<'_, f64>,
        dist: &str,
        wet_threshold: f64,
    ) -> PyResult<Self> {
        let dist = match dist {
            "normal" => core::Distribution::Normal,
            "gamma" => core::Distribution::Gamma { wet_threshold },
            other => {
                return Err(PyValueError::new_err(format!(
                    "dist inválida: '{other}' (esperado 'normal' o 'gamma')"
                )));
            }
        };
        Ok(Self {
            inner: core::ParametricQuantileMapping::fit(obs.as_slice()?, model.as_slice()?, dist)
                .map_err(err)?,
        })
    }

    fn apply<'py>(
        &self,
        py: Python<'py>,
        series: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(self
            .inner
            .apply(series.as_slice()?)
            .map_err(err)?
            .into_pyarray(py))
    }

    fn correct_one(&self, x: f64) -> f64 {
        self.inner.correct_one(x)
    }
}

/// Quantile delta mapping (Cannon et al. 2015): corrige preservando la
/// señal de cambio de la serie objetivo cuantil a cuantil.
#[pyclass]
struct QuantileDeltaMapping {
    inner: core::QuantileDeltaMapping,
}

#[pymethods]
impl QuantileDeltaMapping {
    #[new]
    #[pyo3(signature = (obs, model_hist, n_quantiles=100, kind="add", nodes="endpoints"))]
    fn new(
        obs: PyReadonlyArray1<'_, f64>,
        model_hist: PyReadonlyArray1<'_, f64>,
        n_quantiles: usize,
        kind: &str,
        nodes: &str,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: core::QuantileDeltaMapping::fit_with_nodes(
                obs.as_slice()?,
                model_hist.as_slice()?,
                n_quantiles,
                parse_kind(kind)?,
                parse_nodes(nodes)?,
            )
            .map_err(err)?,
        })
    }

    /// Corrige una serie de proyección (su CDF se estima de ella misma).
    fn apply<'py>(
        &self,
        py: Python<'py>,
        proj: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(self
            .inner
            .apply(proj.as_slice()?)
            .map_err(err)?
            .into_pyarray(py))
    }
}

/// Delta change: perturba observaciones con la señal de cambio del modelo.
#[pyclass]
struct DeltaChange {
    inner: core::DeltaChange,
}

#[pymethods]
impl DeltaChange {
    #[new]
    #[pyo3(signature = (model_hist, model_fut, kind="add"))]
    fn new(
        model_hist: PyReadonlyArray1<'_, f64>,
        model_fut: PyReadonlyArray1<'_, f64>,
        kind: &str,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: core::DeltaChange::fit(
                model_hist.as_slice()?,
                model_fut.as_slice()?,
                parse_kind(kind)?,
            )
            .map_err(err)?,
        })
    }

    /// Delta calibrado (diferencia o razón de medias).
    #[getter]
    fn delta(&self) -> f64 {
        self.inner.delta()
    }

    fn apply<'py>(
        &self,
        py: Python<'py>,
        obs: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(self
            .inner
            .apply(obs.as_slice()?)
            .map_err(err)?
            .into_pyarray(py))
    }
}

/// Adaptación de umbral seco/húmedo (corrección de frecuencia de días húmedos).
#[pyclass]
struct WetDayCorrection {
    inner: core::WetDayCorrection,
}

#[pymethods]
impl WetDayCorrection {
    #[new]
    #[pyo3(signature = (obs, model, obs_wet_threshold=0.1))]
    fn new(
        obs: PyReadonlyArray1<'_, f64>,
        model: PyReadonlyArray1<'_, f64>,
        obs_wet_threshold: f64,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: core::WetDayCorrection::fit(
                obs.as_slice()?,
                model.as_slice()?,
                obs_wet_threshold,
            )
            .map_err(err)?,
        })
    }

    #[getter]
    fn model_threshold(&self) -> f64 {
        self.inner.model_threshold()
    }

    fn transform<'py>(
        &self,
        py: Python<'py>,
        series: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        Ok(self.inner.transform(series.as_slice()?).into_pyarray(py))
    }
}

/// Downscaling por análogos (k-NN estandarizado).
#[pyclass]
struct AnalogDownscaling {
    inner: core::AnalogDownscaling,
}

#[pymethods]
impl AnalogDownscaling {
    /// `predictors`: array 2D (filas = días, columnas = features).
    #[new]
    #[pyo3(signature = (predictors, obs, k=5))]
    fn new(
        predictors: PyReadonlyArray2<'_, f64>,
        obs: PyReadonlyArray1<'_, f64>,
        k: usize,
    ) -> PyResult<Self> {
        let (flat, n_features) = flatten_2d(&predictors);
        Ok(Self {
            inner: core::AnalogDownscaling::fit(&flat, n_features, obs.as_slice()?, k)
                .map_err(err)?,
        })
    }

    fn predict<'py>(
        &self,
        py: Python<'py>,
        queries: PyReadonlyArray2<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let (flat, _) = flatten_2d(&queries);
        Ok(self.inner.predict(&flat).map_err(err)?.into_pyarray(py))
    }
}

/// Downscaling por regresión lineal múltiple (OLS).
#[pyclass]
struct LinearDownscaling {
    inner: core::LinearDownscaling,
}

#[pymethods]
impl LinearDownscaling {
    #[new]
    fn new(
        predictors: PyReadonlyArray2<'_, f64>,
        obs: PyReadonlyArray1<'_, f64>,
    ) -> PyResult<Self> {
        let (flat, n_features) = flatten_2d(&predictors);
        Ok(Self {
            inner: core::LinearDownscaling::fit(&flat, n_features, obs.as_slice()?).map_err(err)?,
        })
    }

    #[getter]
    fn intercept(&self) -> f64 {
        self.inner.intercept()
    }

    #[getter]
    fn coefs<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.inner.coefs().to_vec().into_pyarray(py)
    }

    #[getter]
    fn r2(&self) -> f64 {
        self.inner.r2()
    }

    fn predict<'py>(
        &self,
        py: Python<'py>,
        queries: PyReadonlyArray2<'_, f64>,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let (flat, _) = flatten_2d(&queries);
        Ok(self.inner.predict(&flat).map_err(err)?.into_pyarray(py))
    }
}

/// Schaake shuffle (Clark et al. 2004): reordena `corrected` (2D,
/// filas = días) para reproducir la estructura de rangos de `template`,
/// preservando exactamente las marginales corregidas.
#[pyfunction]
fn schaake_shuffle<'py>(
    py: Python<'py>,
    template: PyReadonlyArray2<'_, f64>,
    corrected: PyReadonlyArray2<'_, f64>,
) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
    use numpy::ndarray::Array2;
    let (flat_t, n_vars) = flatten_2d(&template);
    let (flat_c, n_vars_c) = flatten_2d(&corrected);
    if n_vars != n_vars_c {
        return Err(PyValueError::new_err(format!(
            "template tiene {n_vars} columnas y corrected {n_vars_c}"
        )));
    }
    let out = core::multivariate::schaake_shuffle(&flat_t, &flat_c, n_vars).map_err(err)?;
    let rows = out.len() / n_vars;
    let arr = Array2::from_shape_vec((rows, n_vars), out)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(arr.into_pyarray(py))
}

/// PET de Hargreaves (mm/día) con radiación extraterrestre FAO-56.
/// `days_of_year`: día del año (1..=366) por paso; `latitude_deg`: sur negativo.
#[pyfunction]
fn hargreaves<'py>(
    py: Python<'py>,
    tmin: PyReadonlyArray1<'_, f64>,
    tmax: PyReadonlyArray1<'_, f64>,
    days_of_year: Vec<u32>,
    latitude_deg: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    Ok(core::pet::hargreaves(
        tmin.as_slice()?,
        tmax.as_slice()?,
        &days_of_year,
        latitude_deg,
    )
    .map_err(err)?
    .into_pyarray(py))
}

/// Raíz del error cuadrático medio entre series pareadas.
#[pyfunction]
fn rmse(sim: PyReadonlyArray1<'_, f64>, obs: PyReadonlyArray1<'_, f64>) -> PyResult<f64> {
    core::metrics::rmse(sim.as_slice()?, obs.as_slice()?).map_err(err)
}

/// Sesgo medio `mean(sim) - mean(obs)`.
#[pyfunction]
fn mean_bias(sim: PyReadonlyArray1<'_, f64>, obs: PyReadonlyArray1<'_, f64>) -> PyResult<f64> {
    core::metrics::mean_bias(sim.as_slice()?, obs.as_slice()?).map_err(err)
}

/// Estadístico KS de dos muestras.
#[pyfunction]
fn ks_statistic(sim: PyReadonlyArray1<'_, f64>, obs: PyReadonlyArray1<'_, f64>) -> PyResult<f64> {
    core::metrics::ks_statistic(sim.as_slice()?, obs.as_slice()?).map_err(err)
}

/// Sesgo por cuantil. Devuelve lista de tuplas `(prob, sim, obs, bias)`.
#[pyfunction]
#[pyo3(signature = (sim, obs, probs=None))]
fn quantile_bias(
    sim: PyReadonlyArray1<'_, f64>,
    obs: PyReadonlyArray1<'_, f64>,
    probs: Option<PyReadonlyArray1<'_, f64>>,
) -> PyResult<Vec<(f64, f64, f64, f64)>> {
    let default = core::metrics::DEFAULT_QUANTILE_PROBS;
    let probs_vec: Vec<f64> = match &probs {
        Some(p) => p.as_slice()?.to_vec(),
        None => default.to_vec(),
    };
    Ok(
        core::metrics::quantile_bias(sim.as_slice()?, obs.as_slice()?, &probs_vec)
            .map_err(err)?
            .into_iter()
            .map(|q| (q.prob, q.sim, q.obs, q.bias))
            .collect(),
    )
}

/// Validación con split temporal. Devuelve un dict con las métricas del
/// período de validación (corregido y línea base `*_raw`).
#[pyfunction]
#[pyo3(signature = (obs, model, calib_frac=0.7, n_quantiles=100, kind="add", nodes="endpoints", wet_day_threshold=None))]
#[allow(clippy::too_many_arguments)]
fn validate_split<'py>(
    py: Python<'py>,
    obs: PyReadonlyArray1<'_, f64>,
    model: PyReadonlyArray1<'_, f64>,
    calib_frac: f64,
    n_quantiles: usize,
    kind: &str,
    nodes: &str,
    wet_day_threshold: Option<f64>,
) -> PyResult<Bound<'py, PyDict>> {
    let opts = core::QmOptions {
        n_quantiles,
        kind: parse_kind(kind)?,
        placement: parse_nodes(nodes)?,
        wet_day_threshold,
    };
    let report = core::validate_split_with(obs.as_slice()?, model.as_slice()?, calib_frac, &opts)
        .map_err(err)?;

    let dict = PyDict::new(py);
    dict.set_item("split_index", report.split_index)?;
    dict.set_item("rmse", report.rmse)?;
    dict.set_item("rmse_raw", report.rmse_raw)?;
    dict.set_item("mean_bias", report.mean_bias)?;
    dict.set_item("mean_bias_raw", report.mean_bias_raw)?;
    dict.set_item("ks", report.ks)?;
    dict.set_item("ks_raw", report.ks_raw)?;
    let qb: Vec<(f64, f64, f64, f64)> = report
        .quantile_bias
        .into_iter()
        .map(|q| (q.prob, q.sim, q.obs, q.bias))
        .collect();
    dict.set_item("quantile_bias", qb)?;
    Ok(dict)
}

/// Módulo Python `downscale_rs`.
#[pymodule]
fn downscale_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<QuantileMapping>()?;
    m.add_class::<QuantileDeltaMapping>()?;
    m.add_class::<ParametricQuantileMapping>()?;
    m.add_class::<DeltaChange>()?;
    m.add_class::<WetDayCorrection>()?;
    m.add_class::<AnalogDownscaling>()?;
    m.add_class::<LinearDownscaling>()?;
    m.add_function(wrap_pyfunction!(schaake_shuffle, m)?)?;
    m.add_function(wrap_pyfunction!(hargreaves, m)?)?;
    m.add_function(wrap_pyfunction!(rmse, m)?)?;
    m.add_function(wrap_pyfunction!(mean_bias, m)?)?;
    m.add_function(wrap_pyfunction!(ks_statistic, m)?)?;
    m.add_function(wrap_pyfunction!(quantile_bias, m)?)?;
    m.add_function(wrap_pyfunction!(validate_split, m)?)?;
    Ok(())
}
