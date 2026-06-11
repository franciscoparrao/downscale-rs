//! Validación por split temporal: calibrar en el primer tramo de la serie
//! y evaluar la corrección en el tramo restante (holdout).

use crate::error::{DownscaleError, Result, check_same_len, check_series};
use crate::metrics::{
    DEFAULT_QUANTILE_PROBS, QuantileBias, ks_statistic, mean_bias, quantile_bias, rmse,
};
use crate::qm::{Kind, NodePlacement, QuantileMapping};
use crate::wetday::WetDayCorrection;

/// Configuración del quantile mapping a validar.
#[derive(Debug, Clone, Copy)]
pub struct QmOptions {
    /// Número de cuantiles de la CDF empírica.
    pub n_quantiles: usize,
    /// Tipo de corrección en colas.
    pub kind: Kind,
    /// Colocación de nodos de probabilidad.
    pub placement: NodePlacement,
    /// Si está presente, aplica adaptación de umbral seco/húmedo
    /// ([`WetDayCorrection`]) con este umbral observado antes del EQM.
    pub wet_day_threshold: Option<f64>,
}

impl Default for QmOptions {
    fn default() -> Self {
        Self {
            n_quantiles: 100,
            kind: Kind::Additive,
            placement: NodePlacement::Endpoints,
            wet_day_threshold: None,
        }
    }
}

/// Divide una serie temporal en calibración y validación.
///
/// `calib_frac` es la fracción inicial (cronológica) destinada a
/// calibración; el resto queda para validación. Ambos tramos deben quedar
/// no vacíos.
///
/// # Errors
///
/// [`DownscaleError::InvalidParameter`] si `calib_frac` no deja al menos un
/// punto en cada tramo.
///
/// # Ejemplo
///
/// ```
/// let serie = [1.0, 2.0, 3.0, 4.0];
/// let (cal, val) = downscale_core::validation::split_temporal(&serie, 0.75).unwrap();
/// assert_eq!(cal, &[1.0, 2.0, 3.0]);
/// assert_eq!(val, &[4.0]);
/// ```
pub fn split_temporal(series: &[f64], calib_frac: f64) -> Result<(&[f64], &[f64])> {
    check_series("series", series, 2)?;
    if !(0.0..=1.0).contains(&calib_frac) {
        return Err(DownscaleError::InvalidParameter {
            name: "calib_frac",
            value: calib_frac,
            expected: "en (0, 1)",
        });
    }
    let split = (series.len() as f64 * calib_frac).round() as usize;
    if split == 0 || split == series.len() {
        return Err(DownscaleError::InvalidParameter {
            name: "calib_frac",
            value: calib_frac,
            expected: "que deje >= 1 punto en calibración y validación",
        });
    }
    Ok(series.split_at(split))
}

/// Resultado de la validación holdout de un quantile mapping.
#[derive(Debug, Clone)]
pub struct ValidationReport {
    /// Índice donde se cortó la serie (largo del tramo de calibración).
    pub split_index: usize,
    /// RMSE de la serie corregida vs observada (período de validación).
    pub rmse: f64,
    /// RMSE del modelo crudo vs observado, como línea base.
    pub rmse_raw: f64,
    /// Sesgo medio de la serie corregida.
    pub mean_bias: f64,
    /// Sesgo medio del modelo crudo, como línea base.
    pub mean_bias_raw: f64,
    /// Estadístico KS de la serie corregida vs observada.
    pub ks: f64,
    /// Estadístico KS del modelo crudo vs observado, como línea base.
    pub ks_raw: f64,
    /// Sesgo por cuantil de la serie corregida vs observada.
    pub quantile_bias: Vec<QuantileBias>,
}

/// Valida un quantile mapping con split temporal sobre series pareadas.
///
/// Calibra el mapeo con el primer tramo (`calib_frac`) de `obs` y `model`,
/// corrige el tramo de validación del modelo y lo compara contra las
/// observaciones del mismo tramo. Reporta cada métrica junto a su línea
/// base sin corrección, para cuantificar la ganancia del método.
///
/// # Errors
///
/// Series de largos distintos, demasiado cortas, con NaN/inf, o
/// parámetros fuera de rango.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::qm::Kind;
/// use downscale_core::validation::validate_split;
///
/// let obs: Vec<f64> = (0..200).map(|i| 10.0 + f64::from(i % 50)).collect();
/// let model: Vec<f64> = obs.iter().map(|v| v + 3.0).collect();
/// let report = validate_split(&obs, &model, 0.7, 100, Kind::Additive).unwrap();
/// assert!(report.rmse < report.rmse_raw);
/// ```
pub fn validate_split(
    obs: &[f64],
    model: &[f64],
    calib_frac: f64,
    n_quantiles: usize,
    kind: Kind,
) -> Result<ValidationReport> {
    validate_split_with(
        obs,
        model,
        calib_frac,
        &QmOptions {
            n_quantiles,
            kind,
            ..QmOptions::default()
        },
    )
}

/// Como [`validate_split`], con configuración completa ([`QmOptions`]):
/// colocación de nodos y adaptación opcional de umbral seco/húmedo
/// (calibrada solo con el tramo de calibración, como corresponde).
///
/// # Errors
///
/// Igual que [`validate_split`].
pub fn validate_split_with(
    obs: &[f64],
    model: &[f64],
    calib_frac: f64,
    opts: &QmOptions,
) -> Result<ValidationReport> {
    check_same_len("obs", obs, "model", model)?;
    let (obs_cal, obs_val) = split_temporal(obs, calib_frac)?;
    let (mod_cal_raw, mod_val_raw) = split_temporal(model, calib_frac)?;

    // Adaptación de umbral seco/húmedo (opcional), calibrada en el
    // primer tramo y aplicada a ambos tramos del modelo. Las métricas
    // "raw" siempre se calculan contra el modelo sin transformar.
    let (mod_cal, mod_val): (Vec<f64>, Vec<f64>) = match opts.wet_day_threshold {
        Some(thr) => {
            let wd = WetDayCorrection::fit(obs_cal, mod_cal_raw, thr)?;
            (wd.transform(mod_cal_raw), wd.transform(mod_val_raw))
        }
        None => (mod_cal_raw.to_vec(), mod_val_raw.to_vec()),
    };

    let qm = QuantileMapping::fit_with_nodes(
        obs_cal,
        &mod_cal,
        opts.n_quantiles,
        opts.kind,
        opts.placement,
    )?;
    let corrected = qm.apply(&mod_val)?;

    Ok(ValidationReport {
        split_index: obs_cal.len(),
        rmse: rmse(&corrected, obs_val)?,
        rmse_raw: rmse(mod_val_raw, obs_val)?,
        mean_bias: mean_bias(&corrected, obs_val)?,
        mean_bias_raw: mean_bias(mod_val_raw, obs_val)?,
        ks: ks_statistic(&corrected, obs_val)?,
        ks_raw: ks_statistic(mod_val_raw, obs_val)?,
        quantile_bias: quantile_bias(&corrected, obs_val, &DEFAULT_QUANTILE_PROBS)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_respects_fraction() {
        let s: Vec<f64> = (0..10).map(f64::from).collect();
        let (cal, val) = split_temporal(&s, 0.7).unwrap();
        assert_eq!(cal.len(), 7);
        assert_eq!(val.len(), 3);
        assert_eq!(cal[6], 6.0);
        assert_eq!(val[0], 7.0);
    }

    #[test]
    fn split_rejects_degenerate_fractions() {
        let s = [1.0, 2.0, 3.0];
        assert!(split_temporal(&s, 0.0).is_err());
        assert!(split_temporal(&s, 1.0).is_err());
        assert!(split_temporal(&s, 1.5).is_err());
    }

    #[test]
    fn validate_rejects_mismatched_series() {
        let err =
            validate_split(&[1.0, 2.0], &[1.0, 2.0, 3.0], 0.5, 10, Kind::Additive).unwrap_err();
        assert!(matches!(err, DownscaleError::LengthMismatch { .. }));
    }

    #[test]
    fn validation_improves_biased_model() {
        // Observaciones con estacionalidad sintética; modelo con sesgo
        // aditivo de +3 y leve amplificación.
        let obs: Vec<f64> = (0..730)
            .map(|i| 15.0 + 8.0 * (f64::from(i) * std::f64::consts::TAU / 365.0).sin())
            .collect();
        let model: Vec<f64> = obs.iter().map(|v| v * 1.1 + 3.0).collect();

        let report = validate_split(&obs, &model, 0.5, 100, Kind::Additive).unwrap();

        assert_eq!(report.split_index, 365);
        assert!(report.rmse < report.rmse_raw, "QM debe reducir el RMSE");
        assert!(
            report.mean_bias.abs() < report.mean_bias_raw.abs(),
            "QM debe reducir el sesgo medio"
        );
        assert!(report.ks <= report.ks_raw, "QM no debe empeorar el KS");
        // El sesgo de la mediana corregida debe ser pequeño.
        let median = report
            .quantile_bias
            .iter()
            .find(|q| (q.prob - 0.5).abs() < 1e-12)
            .unwrap();
        assert!(median.bias.abs() < 0.5, "sesgo mediana = {}", median.bias);
    }
}
