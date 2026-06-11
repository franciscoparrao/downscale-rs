//! Quantile mapping empírico (EQM) para corrección de sesgo.
//!
//! Implementa la transformación clásica `x_corr = F_obs⁻¹(F_mod(x))` sobre
//! cuantiles empíricos estimados en un período de calibración común entre
//! observaciones y modelo (GCM/RCM). Dentro del rango calibrado se interpola
//! linealmente; fuera del rango se extrapola con corrección constante
//! (aditiva o multiplicativa según [`Kind`]), siguiendo la práctica de
//! `cmethods`/`xclim`.

use crate::error::{DownscaleError, Result, check_series};

/// Número mínimo de puntos para calibrar un mapeo.
const MIN_FIT_LEN: usize = 2;

/// Tipo de corrección en las colas (fuera del rango calibrado).
///
/// Dentro del rango calibrado ambas variantes aplican la misma transformación
/// cuantil-a-cuantil; difieren solo en la extrapolación.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Kind {
    /// Corrección aditiva (`x + Δ`). Usual para temperatura.
    #[default]
    Additive,
    /// Corrección multiplicativa (`x · δ`). Usual para precipitación
    /// (preserva el cero y evita valores negativos).
    Multiplicative,
}

/// Colocación de los nodos de probabilidad de la CDF empírica.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NodePlacement {
    /// `i/(n−1)`, incluye 0 y 1 (min y max de la muestra). Default.
    #[default]
    Endpoints,
    /// `(i+0.5)/n`, sin extremos — convención de xclim/xsdba. Útil para
    /// paridad bit-cercana con esa referencia.
    Midpoint,
}

/// Mapeo de cuantiles empírico calibrado.
///
/// Se construye con [`QuantileMapping::fit`] sobre un período de calibración
/// y luego se aplica a cualquier serie del modelo con
/// [`QuantileMapping::apply`].
///
/// # Ejemplo
///
/// ```
/// use downscale_core::qm::{Kind, QuantileMapping};
///
/// // Modelo con sesgo aditivo constante de +2.0 respecto a lo observado.
/// let obs: Vec<f64> = (0..100).map(|i| f64::from(i) * 0.1).collect();
/// let model: Vec<f64> = obs.iter().map(|v| v + 2.0).collect();
///
/// let qm = QuantileMapping::fit(&obs, &model, 50, Kind::Additive).unwrap();
/// let corrected = qm.apply(&model).unwrap();
///
/// // La corrección remueve el sesgo de la media.
/// let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
/// assert!((mean(&corrected) - mean(&obs)).abs() < 1e-9);
/// ```
#[derive(Debug, Clone)]
pub struct QuantileMapping {
    kind: Kind,
    /// Probabilidades comunes en \[0, 1\] (crecientes).
    probs: Vec<f64>,
    /// Cuantiles observados en `probs`.
    obs_q: Vec<f64>,
    /// Cuantiles del modelo (período histórico) en `probs`.
    mod_q: Vec<f64>,
}

impl QuantileMapping {
    /// Calibra el mapeo con observaciones y modelo del mismo período.
    ///
    /// `n_quantiles` controla la resolución de la CDF empírica (típico:
    /// 100). Las series no necesitan estar pareadas en el tiempo ni tener
    /// el mismo largo: solo se comparan sus distribuciones.
    ///
    /// # Errors
    ///
    /// - [`DownscaleError::SeriesTooShort`] si alguna serie tiene menos de
    ///   2 puntos.
    /// - [`DownscaleError::NonFinite`] si alguna serie contiene NaN/inf.
    /// - [`DownscaleError::InvalidParameter`] si `n_quantiles < 2`.
    pub fn fit(obs: &[f64], model: &[f64], n_quantiles: usize, kind: Kind) -> Result<Self> {
        Self::fit_with_nodes(obs, model, n_quantiles, kind, NodePlacement::Endpoints)
    }

    /// Como [`QuantileMapping::fit`], con colocación de nodos explícita.
    ///
    /// # Errors
    ///
    /// Igual que [`QuantileMapping::fit`].
    pub fn fit_with_nodes(
        obs: &[f64],
        model: &[f64],
        n_quantiles: usize,
        kind: Kind,
        placement: NodePlacement,
    ) -> Result<Self> {
        check_series("obs", obs, MIN_FIT_LEN)?;
        check_series("model", model, MIN_FIT_LEN)?;
        if n_quantiles < 2 {
            return Err(DownscaleError::InvalidParameter {
                name: "n_quantiles",
                value: n_quantiles as f64,
                expected: ">= 2",
            });
        }

        let probs = node_probs(n_quantiles, placement);
        let obs_q = empirical_quantiles(obs, &probs);
        let mod_q = empirical_quantiles(model, &probs);

        Ok(Self {
            kind,
            probs,
            obs_q,
            mod_q,
        })
    }

    /// Corrige una serie del modelo (histórica o futura).
    ///
    /// # Errors
    ///
    /// [`DownscaleError::NonFinite`] si la serie contiene NaN/inf.
    pub fn apply(&self, series: &[f64]) -> Result<Vec<f64>> {
        check_series("series", series, 0)?;
        Ok(series.iter().map(|&x| self.correct_one(x)).collect())
    }

    /// Corrige un único valor.
    #[must_use]
    pub fn correct_one(&self, x: f64) -> f64 {
        let lo = self.mod_q[0];
        let hi = self.mod_q[self.mod_q.len() - 1];

        if x < lo {
            return self.extrapolate(x, 0);
        }
        if x > hi {
            return self.extrapolate(x, self.mod_q.len() - 1);
        }
        // F_mod(x): probabilidad empírica de x en la CDF del modelo.
        let p = interp(&self.mod_q, &self.probs, x);
        // F_obs⁻¹(p): cuantil observado en esa probabilidad.
        interp(&self.probs, &self.obs_q, p)
    }

    /// Corrección constante en la cola (`idx` = 0 o último).
    fn extrapolate(&self, x: f64, idx: usize) -> f64 {
        match self.kind {
            Kind::Additive => x + (self.obs_q[idx] - self.mod_q[idx]),
            Kind::Multiplicative => {
                if self.mod_q[idx] == 0.0 {
                    // Cola degenerada (p. ej. precipitación con muchos ceros):
                    // sin razón definida, se conserva el valor.
                    x
                } else {
                    x * (self.obs_q[idx] / self.mod_q[idx])
                }
            }
        }
    }

    /// Tipo de corrección usado en las colas.
    #[must_use]
    pub fn kind(&self) -> Kind {
        self.kind
    }

    /// Probabilidades de los cuantiles calibrados.
    #[must_use]
    pub fn probs(&self) -> &[f64] {
        &self.probs
    }
}

/// Nodos de probabilidad según la colocación elegida.
pub(crate) fn node_probs(n_quantiles: usize, placement: NodePlacement) -> Vec<f64> {
    match placement {
        NodePlacement::Endpoints => (0..n_quantiles)
            .map(|i| i as f64 / (n_quantiles - 1) as f64)
            .collect(),
        NodePlacement::Midpoint => (0..n_quantiles)
            .map(|i| (i as f64 + 0.5) / n_quantiles as f64)
            .collect(),
    }
}

/// Cuantiles empíricos con interpolación lineal (tipo 7 de Hyndman & Fan,
/// el default de NumPy/R), sobre una copia ordenada de la serie.
pub(crate) fn empirical_quantiles(series: &[f64], probs: &[f64]) -> Vec<f64> {
    let mut sorted = series.to_vec();
    sorted.sort_unstable_by(|a, b| a.partial_cmp(b).expect("series validada sin NaN"));
    probs.iter().map(|&p| quantile_sorted(&sorted, p)).collect()
}

/// Cuantil tipo 7 sobre una serie ya ordenada. `p` se satura a \[0, 1\].
pub(crate) fn quantile_sorted(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    debug_assert!(n >= 1);
    let h = p.clamp(0.0, 1.0) * (n - 1) as f64;
    let lo = h.floor() as usize;
    let hi = h.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        sorted[lo] + (h - lo as f64) * (sorted[hi] - sorted[lo])
    }
}

/// Interpolación lineal de `y(x)` sobre nodos crecientes `xs` → `ys`.
///
/// Para `x` fuera del rango devuelve el extremo (clamp). Tramos con `xs`
/// repetidos (CDF plana) devuelven el primer `y` del tramo.
pub(crate) fn interp(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    debug_assert_eq!(xs.len(), ys.len());
    let n = xs.len();
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[n - 1] {
        return ys[n - 1];
    }
    // partition_point: primer índice con xs[i] > x; como x < xs[n-1], i < n.
    let i = xs.partition_point(|&v| v <= x);
    let (x0, x1) = (xs[i - 1], xs[i]);
    let (y0, y1) = (ys[i - 1], ys[i]);
    if x1 == x0 {
        y0
    } else {
        y0 + (x - x0) / (x1 - x0) * (y1 - y0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Generador determinista (LCG) para series sintéticas reproducibles.
    fn lcg(seed: u64, n: usize) -> Vec<f64> {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (state >> 11) as f64 / (1u64 << 53) as f64
            })
            .collect()
    }

    fn mean(s: &[f64]) -> f64 {
        s.iter().sum::<f64>() / s.len() as f64
    }

    #[test]
    fn identity_when_distributions_match() {
        let obs = lcg(42, 500);
        let qm = QuantileMapping::fit(&obs, &obs, 100, Kind::Additive).unwrap();
        let corrected = qm.apply(&obs).unwrap();
        for (c, o) in corrected.iter().zip(&obs) {
            assert!((c - o).abs() < 1e-9, "corrected={c}, original={o}");
        }
    }

    #[test]
    fn removes_additive_bias() {
        let obs: Vec<f64> = lcg(1, 1000).iter().map(|v| 10.0 + 5.0 * v).collect();
        let model: Vec<f64> = obs.iter().map(|v| v + 3.0).collect();
        let qm = QuantileMapping::fit(&obs, &model, 100, Kind::Additive).unwrap();
        let corrected = qm.apply(&model).unwrap();
        assert!((mean(&corrected) - mean(&obs)).abs() < 1e-6);
    }

    #[test]
    fn removes_multiplicative_bias() {
        let obs: Vec<f64> = lcg(7, 1000).iter().map(|v| 4.0 * v).collect();
        let model: Vec<f64> = obs.iter().map(|v| v * 1.5).collect();
        let qm = QuantileMapping::fit(&obs, &model, 100, Kind::Multiplicative).unwrap();
        let corrected = qm.apply(&model).unwrap();
        assert!((mean(&corrected) - mean(&obs)).abs() < 1e-6);
    }

    #[test]
    fn additive_extrapolation_beyond_calibration_range() {
        let obs: Vec<f64> = (0..100).map(f64::from).collect();
        let model: Vec<f64> = obs.iter().map(|v| v + 10.0).collect();
        let qm = QuantileMapping::fit(&obs, &model, 50, Kind::Additive).unwrap();
        // Valor futuro sobre el máximo calibrado del modelo (109).
        assert!((qm.correct_one(200.0) - 190.0).abs() < 1e-9);
        // Valor bajo el mínimo calibrado del modelo (10).
        assert!((qm.correct_one(0.0) - (-10.0)).abs() < 1e-9);
    }

    #[test]
    fn multiplicative_extrapolation_scales_tail() {
        let obs: Vec<f64> = (1..=100).map(f64::from).collect();
        let model: Vec<f64> = obs.iter().map(|v| v * 2.0).collect();
        let qm = QuantileMapping::fit(&obs, &model, 50, Kind::Multiplicative).unwrap();
        // Sobre el máximo del modelo (200): razón obs/mod en la cola = 0.5.
        assert!((qm.correct_one(400.0) - 200.0).abs() < 1e-9);
    }

    #[test]
    fn multiplicative_preserves_zeros() {
        // Precipitación: muchos ceros en ambas series.
        let obs = [0.0, 0.0, 0.0, 1.0, 2.0, 5.0];
        let model = [0.0, 0.0, 0.0, 2.0, 4.0, 10.0];
        let qm = QuantileMapping::fit(&obs, &model, 20, Kind::Multiplicative).unwrap();
        assert_eq!(qm.correct_one(0.0), 0.0);
        assert!(qm.correct_one(10.0) <= 5.0 + 1e-9);
    }

    #[test]
    fn midpoint_nodes_match_xsdba_convention() {
        let obs: Vec<f64> = (0..200).map(f64::from).collect();
        let qm = QuantileMapping::fit_with_nodes(
            &obs,
            &obs,
            10,
            Kind::Additive,
            NodePlacement::Midpoint,
        )
        .unwrap();
        let expected: Vec<f64> = (0..10).map(|i| (f64::from(i) + 0.5) / 10.0).collect();
        for (p, e) in qm.probs().iter().zip(&expected) {
            assert!((p - e).abs() < 1e-15);
        }
        // Sigue siendo ~identidad cuando las distribuciones coinciden.
        let c = qm.apply(&[50.0, 100.0]).unwrap();
        assert!((c[0] - 50.0).abs() < 1e-9 && (c[1] - 100.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_short_series() {
        let err = QuantileMapping::fit(&[1.0], &[1.0, 2.0], 10, Kind::Additive).unwrap_err();
        assert_eq!(
            err,
            DownscaleError::SeriesTooShort {
                name: "obs",
                len: 1,
                min: 2
            }
        );
    }

    #[test]
    fn rejects_nan() {
        let err =
            QuantileMapping::fit(&[1.0, f64::NAN], &[1.0, 2.0], 10, Kind::Additive).unwrap_err();
        assert_eq!(
            err,
            DownscaleError::NonFinite {
                name: "obs",
                index: 1
            }
        );
    }

    #[test]
    fn rejects_invalid_n_quantiles() {
        let err = QuantileMapping::fit(&[1.0, 2.0], &[1.0, 2.0], 1, Kind::Additive).unwrap_err();
        assert!(matches!(
            err,
            DownscaleError::InvalidParameter {
                name: "n_quantiles",
                ..
            }
        ));
    }

    #[test]
    fn quantile_sorted_matches_numpy_type7() {
        let sorted = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(quantile_sorted(&sorted, 0.0), 1.0);
        assert_eq!(quantile_sorted(&sorted, 1.0), 4.0);
        assert_eq!(quantile_sorted(&sorted, 0.5), 2.5);
        // numpy.quantile([1,2,3,4], 0.25) == 1.75
        assert!((quantile_sorted(&sorted, 0.25) - 1.75).abs() < 1e-12);
    }

    #[test]
    fn interp_handles_flat_segments() {
        let xs = [0.0, 1.0, 1.0, 2.0];
        let ys = [0.0, 10.0, 20.0, 30.0];
        // En el tramo plano devuelve el primer y del tramo siguiente al salto.
        let v = interp(&xs, &ys, 1.0);
        assert!((10.0..=20.0).contains(&v));
        assert!((interp(&xs, &ys, 1.5) - 25.0).abs() < 1e-12);
    }
}
