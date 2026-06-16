//! Quantile Delta Mapping (QDM, Cannon et al. 2015): corrige el sesgo
//! preservando la señal de cambio del modelo en cada cuantil.
//!
//! Para cada valor `x` de la serie a corregir (típicamente la proyección
//! futura) se calcula su probabilidad de no-excedencia `p` en la CDF de la
//! **propia serie**; el delta del modelo en ese cuantil
//! (`x / F_hist⁻¹(p)` multiplicativo o `x − F_hist⁻¹(p)` aditivo) se
//! aplica sobre el cuantil observado `F_obs⁻¹(p)`. A diferencia del EQM,
//! el cambio proyectado por el modelo (p. ej. intensificación de extremos)
//! no se "descorrige" hacia el clima histórico.

use crate::error::{DownscaleError, Result, check_series};
use crate::qm::{Kind, NodePlacement, empirical_quantiles, interp, node_probs};

/// Mínimo de puntos por serie.
const MIN_FIT_LEN: usize = 2;

/// QDM calibrado con observaciones y modelo del período histórico.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::qdm::QuantileDeltaMapping;
/// use downscale_core::qm::Kind;
///
/// // Modelo histórico con sesgo +2; proyección con +3 °C de señal de cambio.
/// let obs: Vec<f64> = (0..500).map(|i| 10.0 + (f64::from(i) * 0.13).sin() * 5.0).collect();
/// let hist: Vec<f64> = obs.iter().map(|v| v + 2.0).collect();
/// let proj: Vec<f64> = hist.iter().map(|v| v + 3.0).collect();
///
/// let qdm = QuantileDeltaMapping::fit(&obs, &hist, 100, Kind::Additive).unwrap();
/// let corrected = qdm.apply(&proj).unwrap();
///
/// // Remueve el sesgo (+2) pero conserva la señal de cambio (+3).
/// let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
/// assert!((mean(&corrected) - (mean(&obs) + 3.0)).abs() < 1e-9);
/// ```
#[derive(Debug, Clone)]
pub struct QuantileDeltaMapping {
    kind: Kind,
    n_quantiles: usize,
    placement: NodePlacement,
    probs: Vec<f64>,
    obs_q: Vec<f64>,
    hist_q: Vec<f64>,
}

impl QuantileDeltaMapping {
    /// Calibra con observaciones y modelo del mismo período histórico.
    ///
    /// # Errors
    ///
    /// Igual que [`crate::qm::QuantileMapping::fit`].
    pub fn fit(obs: &[f64], model_hist: &[f64], n_quantiles: usize, kind: Kind) -> Result<Self> {
        Self::fit_with_nodes(obs, model_hist, n_quantiles, kind, NodePlacement::Endpoints)
    }

    /// Como [`QuantileDeltaMapping::fit`], con colocación de nodos explícita.
    ///
    /// # Errors
    ///
    /// Igual que [`QuantileDeltaMapping::fit`].
    pub fn fit_with_nodes(
        obs: &[f64],
        model_hist: &[f64],
        n_quantiles: usize,
        kind: Kind,
        placement: NodePlacement,
    ) -> Result<Self> {
        check_series("obs", obs, MIN_FIT_LEN)?;
        check_series("model_hist", model_hist, MIN_FIT_LEN)?;
        if n_quantiles < 2 {
            return Err(DownscaleError::InvalidParameter {
                name: "n_quantiles",
                value: n_quantiles as f64,
                expected: ">= 2",
            });
        }
        let probs = node_probs(n_quantiles, placement);
        let obs_q = empirical_quantiles(obs, &probs);
        let hist_q = empirical_quantiles(model_hist, &probs);
        Ok(Self {
            kind,
            n_quantiles,
            placement,
            probs,
            obs_q,
            hist_q,
        })
    }

    /// Corrige una serie de proyección completa.
    ///
    /// La CDF de la proyección se estima de la serie entrante, por lo que
    /// QDM opera sobre series (no valor a valor). Para escenarios largos,
    /// la práctica usual es aplicar por ventanas (p. ej. 30 años) para que
    /// `F_proj` sea representativa del período.
    ///
    /// # Errors
    ///
    /// [`DownscaleError::SeriesTooShort`] / [`DownscaleError::NonFinite`].
    pub fn apply(&self, proj: &[f64]) -> Result<Vec<f64>> {
        check_series("proj", proj, MIN_FIT_LEN)?;
        let proj_q = empirical_quantiles(proj, &self.probs);
        Ok(proj
            .iter()
            .map(|&x| self.correct_value(x, &proj_q))
            .collect())
    }

    /// Corrige por **ventanas no solapadas** de `window` días: la CDF de la
    /// proyección se estima dentro de cada bloque, capturando el cambio
    /// temporal de la distribución (recomendado para escenarios largos, p.
    /// ej. ventanas de 30 años, donde la CDF global mezclaría climas muy
    /// distintos). Un bloque final con menos de 2 días reusa la CDF global.
    ///
    /// # Errors
    ///
    /// [`DownscaleError::SeriesTooShort`] / [`DownscaleError::NonFinite`], o
    /// `window < 2`.
    pub fn apply_windowed(&self, proj: &[f64], window: usize) -> Result<Vec<f64>> {
        check_series("proj", proj, MIN_FIT_LEN)?;
        if window < MIN_FIT_LEN {
            return Err(DownscaleError::InvalidParameter {
                name: "window",
                value: window as f64,
                expected: ">= 2",
            });
        }
        let global_q = empirical_quantiles(proj, &self.probs);
        let mut out = Vec::with_capacity(proj.len());
        for block in proj.chunks(window) {
            let block_q = if block.len() >= MIN_FIT_LEN {
                empirical_quantiles(block, &self.probs)
            } else {
                global_q.clone()
            };
            out.extend(block.iter().map(|&x| self.correct_value(x, &block_q)));
        }
        Ok(out)
    }

    /// Aplica la transformación QDM a un valor dada la CDF `proj_q` de la
    /// proyección (estimada global o por ventana).
    fn correct_value(&self, x: f64, proj_q: &[f64]) -> f64 {
        // Probabilidad de no-excedencia en la CDF de la proyección.
        let p = interp(proj_q, &self.probs, x);
        let obs_at_p = interp(&self.probs, &self.obs_q, p);
        let hist_at_p = interp(&self.probs, &self.hist_q, p);
        match self.kind {
            Kind::Additive => obs_at_p + (x - hist_at_p),
            Kind::Multiplicative => {
                if hist_at_p == 0.0 {
                    // Cola seca degenerada: sin razón definida, se conserva
                    // el cuantil observado (delta = 1).
                    obs_at_p
                } else {
                    obs_at_p * (x / hist_at_p)
                }
            }
        }
    }

    /// Número de cuantiles de la CDF empírica.
    #[must_use]
    pub fn n_quantiles(&self) -> usize {
        self.n_quantiles
    }

    /// Colocación de nodos usada.
    #[must_use]
    pub fn placement(&self) -> NodePlacement {
        self.placement
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::ks_statistic;
    use crate::qm::QuantileMapping;

    /// LCG determinista en (0, 1).
    fn uniform(seed: u64, n: usize) -> Vec<f64> {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((state >> 11) as f64 + 0.5) / (1u64 << 53) as f64
            })
            .collect()
    }

    #[test]
    fn identity_when_no_bias_and_no_change() {
        let obs = uniform(1, 1000);
        let qdm = QuantileDeltaMapping::fit(&obs, &obs, 100, Kind::Additive).unwrap();
        let corrected = qdm.apply(&obs).unwrap();
        for (c, o) in corrected.iter().zip(&obs) {
            assert!((c - o).abs() < 1e-9);
        }
    }

    #[test]
    fn preserves_additive_change_signal() {
        let obs: Vec<f64> = uniform(3, 2000).iter().map(|u| 10.0 + 8.0 * u).collect();
        let hist: Vec<f64> = obs.iter().map(|v| v + 2.5).collect(); // sesgo
        let proj: Vec<f64> = hist.iter().map(|v| v + 4.0).collect(); // señal

        let qdm = QuantileDeltaMapping::fit(&obs, &hist, 100, Kind::Additive).unwrap();
        let corrected = qdm.apply(&proj).unwrap();

        let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
        // Corregido = obs + señal de cambio, sin el sesgo.
        assert!((mean(&corrected) - (mean(&obs) + 4.0)).abs() < 1e-6);
    }

    #[test]
    fn preserves_multiplicative_change_in_extremes() {
        // Proyección: el modelo intensifica 50% todos los cuantiles.
        let obs: Vec<f64> = uniform(7, 4000).iter().map(|u| 10.0 * u * u).collect();
        let hist: Vec<f64> = obs.iter().map(|v| v * 1.3).collect();
        let proj: Vec<f64> = hist.iter().map(|v| v * 1.5).collect();

        let qdm = QuantileDeltaMapping::fit(&obs, &hist, 100, Kind::Multiplicative).unwrap();
        let corrected = qdm.apply(&proj).unwrap();

        // P99 corregido ≈ P99 obs × 1.5 (señal preservada, sesgo 1.3 removido).
        let q99 = |s: &[f64]| {
            let mut v = s.to_vec();
            v.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            crate::qm::quantile_sorted(&v, 0.99)
        };
        let expected = q99(&obs) * 1.5;
        assert!(
            (q99(&corrected) - expected).abs() / expected < 0.02,
            "P99 corregido = {}, esperado = {expected}",
            q99(&corrected)
        );
    }

    #[test]
    fn eqm_shrinks_change_signal_but_qdm_does_not() {
        // La motivación de QDM: con EQM la señal de cambio fuera del rango
        // calibrado se trata con corrección constante; QDM la preserva
        // cuantil a cuantil.
        let obs: Vec<f64> = uniform(11, 3000).iter().map(|u| 5.0 + 10.0 * u).collect();
        let hist: Vec<f64> = obs.iter().map(|v| v * 1.2 + 1.0).collect();
        let proj: Vec<f64> = hist.iter().map(|v| v + 6.0).collect();

        let qdm = QuantileDeltaMapping::fit(&obs, &hist, 100, Kind::Additive).unwrap();
        let qdm_out = qdm.apply(&proj).unwrap();
        let eqm = QuantileMapping::fit(&obs, &hist, 100, Kind::Additive).unwrap();
        let eqm_out = eqm.apply(&proj).unwrap();

        let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
        let qdm_signal = mean(&qdm_out) - mean(&obs);
        let eqm_signal = mean(&eqm_out) - mean(&obs);
        // Señal real del modelo: +6 (aditiva pura sobre hist).
        assert!((qdm_signal - 6.0).abs() < 0.05, "QDM señal = {qdm_signal}");
        // Ambas son correcciones válidas; la distribución corregida por QDM
        // de la proyección histórica equivale a obs (sanity).
        let hist_corr = qdm.apply(&hist).unwrap();
        assert!(ks_statistic(&hist_corr, &obs).unwrap() < 0.02);
        // El test clave: QDM clava la señal; EQM aquí la distorsiona más.
        assert!((qdm_signal - 6.0).abs() <= (eqm_signal - 6.0).abs() + 1e-12);
    }

    #[test]
    fn multiplicative_dry_tail_is_safe() {
        // Precipitación con muchos ceros: hist_q(p) = 0 en gran parte de
        // la CDF; el delta degenera a 1 y se devuelve el cuantil observado.
        let u = uniform(13, 3000);
        let obs: Vec<f64> = u
            .iter()
            .map(|&x| if x < 0.7 { 0.0 } else { x * 10.0 })
            .collect();
        let hist: Vec<f64> = obs.iter().map(|v| v * 1.4).collect();
        let qdm = QuantileDeltaMapping::fit(&obs, &hist, 100, Kind::Multiplicative).unwrap();
        let corrected = qdm.apply(&hist).unwrap();
        assert!(corrected.iter().all(|v| v.is_finite() && *v >= 0.0));
        assert!(ks_statistic(&corrected, &obs).unwrap() < 0.02);
    }

    #[test]
    fn rejects_invalid_input() {
        let s = [1.0, 2.0, 3.0];
        assert!(QuantileDeltaMapping::fit(&s, &s, 1, Kind::Additive).is_err());
        assert!(QuantileDeltaMapping::fit(&[1.0], &s, 10, Kind::Additive).is_err());
        let qdm = QuantileDeltaMapping::fit(&s, &s, 3, Kind::Additive).unwrap();
        assert!(qdm.apply(&[f64::NAN, 1.0]).is_err());
    }

    #[test]
    fn windowed_with_full_window_equals_global() {
        let obs: Vec<f64> = uniform(2, 1000).iter().map(|u| 10.0 + 8.0 * u).collect();
        let hist: Vec<f64> = obs.iter().map(|v| v + 2.0).collect();
        let proj: Vec<f64> = hist.iter().map(|v| v + 3.0).collect();
        let qdm = QuantileDeltaMapping::fit(&obs, &hist, 100, Kind::Additive).unwrap();
        let global = qdm.apply(&proj).unwrap();
        let windowed = qdm.apply_windowed(&proj, proj.len()).unwrap();
        for (a, b) in global.iter().zip(&windowed) {
            assert!((a - b).abs() < 1e-12);
        }
        assert!(qdm.apply_windowed(&proj, 1).is_err());
    }

    #[test]
    fn windowed_tracks_drifting_distribution() {
        // Proyección no estacionaria: primera mitad fría, segunda caliente.
        // La ventana captura cada régimen; la global los mezcla.
        let obs: Vec<f64> = uniform(5, 2000).iter().map(|u| 10.0 + 5.0 * u).collect();
        let hist: Vec<f64> = obs.iter().map(|v| v + 2.0).collect(); // sesgo +2
        let qdm = QuantileDeltaMapping::fit(&obs, &hist, 100, Kind::Additive).unwrap();

        let mut proj: Vec<f64> = (0..1000).map(|i| hist[i] + 1.0).collect(); // +1 de cambio
        proj.extend((1000..2000).map(|i| hist[i] + 9.0)); // +9 de cambio

        let w = qdm.apply_windowed(&proj, 1000).unwrap();
        let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
        // Cada ventana preserva su señal de cambio (1 y 9) sin el sesgo (2).
        assert!((mean(&w[..1000]) - (mean(&obs) + 1.0)).abs() < 0.2);
        assert!((mean(&w[1000..]) - (mean(&obs) + 9.0)).abs() < 0.2);
    }
}
