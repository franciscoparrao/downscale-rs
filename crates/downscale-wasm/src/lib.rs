//! Bindings WebAssembly de downscale-rs vía wasm-bindgen.
//!
//! Expone la corrección de sesgo (EQM, QDM) y las métricas a JavaScript;
//! las series son `Float64Array` (entran como `&[f64]`, salen como
//! `Vec<f64>`). Permite correr el motor en el navegador sin servidor —
//! la cara de demostración del patrón multi-target del portafolio.
//!
//! `kind` ∈ {"add", "mult"}, `nodes` ∈ {"endpoints", "midpoint"}.

use wasm_bindgen::prelude::*;

use downscale_core as core;

fn err(e: core::DownscaleError) -> JsValue {
    JsValue::from_str(&e.to_string())
}

fn parse_kind(kind: &str) -> Result<core::Kind, JsValue> {
    match kind {
        "add" | "additive" => Ok(core::Kind::Additive),
        "mult" | "multiplicative" => Ok(core::Kind::Multiplicative),
        other => Err(JsValue::from_str(&format!(
            "kind inválido: '{other}' (esperado 'add' o 'mult')"
        ))),
    }
}

fn parse_nodes(nodes: &str) -> Result<core::NodePlacement, JsValue> {
    match nodes {
        "endpoints" => Ok(core::NodePlacement::Endpoints),
        "midpoint" => Ok(core::NodePlacement::Midpoint),
        other => Err(JsValue::from_str(&format!(
            "nodes inválido: '{other}' (esperado 'endpoints' o 'midpoint')"
        ))),
    }
}

/// Quantile mapping empírico (EQM).
#[wasm_bindgen]
pub struct QuantileMapping {
    inner: core::QuantileMapping,
}

#[wasm_bindgen]
impl QuantileMapping {
    /// Calibra con observaciones y modelo del mismo período.
    #[wasm_bindgen(constructor)]
    pub fn new(
        obs: &[f64],
        model: &[f64],
        n_quantiles: usize,
        kind: &str,
        nodes: &str,
    ) -> Result<QuantileMapping, JsValue> {
        Ok(Self {
            inner: core::QuantileMapping::fit_with_nodes(
                obs,
                model,
                n_quantiles,
                parse_kind(kind)?,
                parse_nodes(nodes)?,
            )
            .map_err(err)?,
        })
    }

    /// Corrige una serie del modelo.
    pub fn apply(&self, series: &[f64]) -> Result<Vec<f64>, JsValue> {
        self.inner.apply(series).map_err(err)
    }

    /// Corrige un único valor.
    #[wasm_bindgen(js_name = correctOne)]
    pub fn correct_one(&self, x: f64) -> f64 {
        self.inner.correct_one(x)
    }
}

/// Quantile delta mapping (Cannon et al. 2015).
#[wasm_bindgen]
pub struct QuantileDeltaMapping {
    inner: core::QuantileDeltaMapping,
}

#[wasm_bindgen]
impl QuantileDeltaMapping {
    /// Calibra con observaciones y modelo del período histórico.
    #[wasm_bindgen(constructor)]
    pub fn new(
        obs: &[f64],
        model_hist: &[f64],
        n_quantiles: usize,
        kind: &str,
        nodes: &str,
    ) -> Result<QuantileDeltaMapping, JsValue> {
        Ok(Self {
            inner: core::QuantileDeltaMapping::fit_with_nodes(
                obs,
                model_hist,
                n_quantiles,
                parse_kind(kind)?,
                parse_nodes(nodes)?,
            )
            .map_err(err)?,
        })
    }

    /// Corrige una serie de proyección (su CDF se estima de ella misma).
    pub fn apply(&self, proj: &[f64]) -> Result<Vec<f64>, JsValue> {
        self.inner.apply(proj).map_err(err)
    }
}

/// Raíz del error cuadrático medio entre series pareadas.
#[wasm_bindgen]
pub fn rmse(sim: &[f64], obs: &[f64]) -> Result<f64, JsValue> {
    core::metrics::rmse(sim, obs).map_err(err)
}

/// Sesgo medio `mean(sim) - mean(obs)`.
#[wasm_bindgen(js_name = meanBias)]
pub fn mean_bias(sim: &[f64], obs: &[f64]) -> Result<f64, JsValue> {
    core::metrics::mean_bias(sim, obs).map_err(err)
}

/// Estadístico de Kolmogorov–Smirnov de dos muestras.
#[wasm_bindgen(js_name = ksStatistic)]
pub fn ks_statistic(sim: &[f64], obs: &[f64]) -> Result<f64, JsValue> {
    core::metrics::ks_statistic(sim, obs).map_err(err)
}
