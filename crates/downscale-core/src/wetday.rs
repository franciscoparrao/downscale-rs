//! Corrección de frecuencia de días húmedos por adaptación de umbral
//! (Themeßl et al. 2012). Preprocesamiento para EQM en precipitación:
//! encuentra el umbral del modelo cuyo percentil iguala la fracción de
//! días secos observada, y trunca a cero todo lo que quede bajo él.
//! Remueve la llovizna espuria típica de reanálisis/GCM.

use crate::error::{DownscaleError, Result, check_series};
use crate::qm::quantile_sorted;

/// Mínimo de puntos para calibrar el umbral.
const MIN_FIT_LEN: usize = 10;

/// Umbral de adaptación seco/húmedo calibrado.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::wetday::WetDayCorrection;
///
/// // Obs: 60% de días secos. Modelo: llovizna en casi todos los días.
/// let obs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 9.0, 20.0];
/// let model = [0.2, 0.3, 0.1, 0.4, 0.2, 0.5, 3.0, 6.0, 10.0, 22.0];
/// let wd = WetDayCorrection::fit(&obs, &model, 0.1).unwrap();
/// let t = wd.transform(&model);
/// let dry = t.iter().filter(|&&v| v == 0.0).count();
/// assert_eq!(dry, 6); // misma fracción seca que obs
/// ```
#[derive(Debug, Clone, Copy)]
pub struct WetDayCorrection {
    model_threshold: f64,
}

impl WetDayCorrection {
    /// Calibra el umbral del modelo en el período común.
    ///
    /// `obs_wet_threshold` define qué cuenta como día seco en las
    /// observaciones (típico 0.1 mm). El umbral del modelo es el cuantil
    /// de la serie del modelo en la fracción seca observada.
    ///
    /// # Errors
    ///
    /// Series cortas/no finitas o `obs_wet_threshold < 0`.
    pub fn fit(obs: &[f64], model: &[f64], obs_wet_threshold: f64) -> Result<Self> {
        check_series("obs", obs, MIN_FIT_LEN)?;
        check_series("model", model, MIN_FIT_LEN)?;
        if obs_wet_threshold < 0.0 {
            return Err(DownscaleError::InvalidParameter {
                name: "obs_wet_threshold",
                value: obs_wet_threshold,
                expected: ">= 0",
            });
        }
        let dry_frac =
            obs.iter().filter(|&&v| v < obs_wet_threshold).count() as f64 / obs.len() as f64;
        let mut sorted = model.to_vec();
        sorted.sort_unstable_by(|a, b| a.partial_cmp(b).expect("serie validada sin NaN"));
        Ok(Self {
            model_threshold: quantile_sorted(&sorted, dry_frac),
        })
    }

    /// Trunca a cero los valores bajo el umbral calibrado.
    #[must_use]
    pub fn transform(&self, series: &[f64]) -> Vec<f64> {
        series
            .iter()
            .map(|&v| if v <= self.model_threshold { 0.0 } else { v })
            .collect()
    }

    /// Umbral calibrado en las unidades del modelo.
    #[must_use]
    pub fn model_threshold(&self) -> f64 {
        self.model_threshold
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_obs_dry_fraction_on_drizzle_model() {
        // Obs: 80% seco. Modelo: nunca seco (llovizna 0.1–0.5 mm).
        let mut obs = vec![0.0; 80];
        obs.extend((1..=20).map(f64::from));
        let model: Vec<f64> = (0..100)
            .map(|i| 0.1 + 0.004 * f64::from(i) * 100.0)
            .collect();

        let wd = WetDayCorrection::fit(&obs, &model, 0.1).unwrap();
        let t = wd.transform(&model);
        let dry = t.iter().filter(|&&v| v == 0.0).count();
        assert!((78..=82).contains(&dry), "secos = {dry}");
    }

    #[test]
    fn no_op_when_model_already_dry_enough() {
        let obs = [0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 4.0, 5.0, 6.0];
        let model = [0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let wd = WetDayCorrection::fit(&obs, &model, 0.1).unwrap();
        // El umbral cae en la zona de ceros del modelo → solo trunca ceros.
        let t = wd.transform(&model);
        assert_eq!(t.iter().filter(|&&v| v > 0.0).count(), 6);
    }

    #[test]
    fn rejects_negative_threshold() {
        let s = [0.0; 20];
        let m = [1.0; 20];
        assert!(matches!(
            WetDayCorrection::fit(&s, &m, -0.1).unwrap_err(),
            DownscaleError::InvalidParameter { .. }
        ));
    }
}
