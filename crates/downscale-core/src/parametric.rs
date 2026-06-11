//! Quantile mapping paramétrico: normal (temperatura) y gamma mixta con
//! masa en cero (precipitación).
//!
//! La variante gamma modela cada serie como mezcla `p_dry · δ₀ +
//! (1−p_dry) · Gamma(k, θ)` ajustada por máxima verosimilitud sobre los
//! días húmedos, y mapea `x_corr = F_obs⁻¹(F_mod(x))` entre mezclas. Esto
//! corrige simultáneamente intensidad **y frecuencia de días húmedos**
//! (la llovizna espuria del modelo cae bajo la masa seca observada).

use crate::error::{DownscaleError, Result, check_series};
use crate::special::{digamma, gamma_p, gamma_p_inv, norm_cdf, norm_ppf, trigamma};

/// Mínimo de puntos (totales o húmedos) para ajustar una distribución.
const MIN_FIT_LEN: usize = 10;

/// Familia de distribución para el QM paramétrico.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Distribution {
    /// Normal. Para variables ~simétricas (temperatura).
    Normal,
    /// Gamma con masa en cero. Para precipitación; `wet_threshold` define
    /// el límite seco/húmedo en las unidades de la serie (típico 0.1 mm).
    Gamma {
        /// Valores < `wet_threshold` se tratan como secos (masa en cero).
        wet_threshold: f64,
    },
}

/// Parámetros ajustados de una serie.
#[derive(Debug, Clone, Copy)]
enum Fitted {
    Normal { mean: f64, sd: f64 },
    MixedGamma { p_dry: f64, shape: f64, scale: f64 },
}

impl Fitted {
    /// CDF de la serie ajustada.
    fn cdf(&self, x: f64) -> f64 {
        match *self {
            Fitted::Normal { mean, sd } => norm_cdf((x - mean) / sd),
            Fitted::MixedGamma {
                p_dry,
                shape,
                scale,
            } => {
                if x <= 0.0 {
                    p_dry
                } else {
                    p_dry + (1.0 - p_dry) * gamma_p(shape, x / scale)
                }
            }
        }
    }

    /// Inversa de la CDF.
    fn ppf(&self, p: f64) -> f64 {
        match *self {
            Fitted::Normal { mean, sd } => mean + sd * norm_ppf(p.clamp(1e-12, 1.0 - 1e-12)),
            Fitted::MixedGamma {
                p_dry,
                shape,
                scale,
            } => {
                if p <= p_dry {
                    0.0
                } else {
                    let q = (p - p_dry) / (1.0 - p_dry);
                    scale * gamma_p_inv(shape, q)
                }
            }
        }
    }
}

/// Mapeo de cuantiles paramétrico calibrado.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::parametric::{Distribution, ParametricQuantileMapping};
///
/// // Temperatura: modelo con +3 °C de sesgo y varianza inflada.
/// let obs: Vec<f64> = (0..400).map(|i| 12.0 + 5.0 * (f64::from(i) * 0.7).sin()).collect();
/// let model: Vec<f64> = obs.iter().map(|v| (v - 12.0) * 1.4 + 15.0).collect();
///
/// let pqm = ParametricQuantileMapping::fit(&obs, &model, Distribution::Normal).unwrap();
/// let corrected = pqm.apply(&model).unwrap();
/// let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
/// assert!((mean(&corrected) - mean(&obs)).abs() < 0.05);
/// ```
#[derive(Debug, Clone)]
pub struct ParametricQuantileMapping {
    obs: Fitted,
    model: Fitted,
}

impl ParametricQuantileMapping {
    /// Ajusta la misma familia a observaciones y modelo (período común).
    ///
    /// # Errors
    ///
    /// - [`DownscaleError::SeriesTooShort`] si una serie (o su subconjunto
    ///   húmedo, en gamma) tiene menos de 10 puntos.
    /// - [`DownscaleError::NonFinite`] si hay NaN/inf.
    /// - [`DownscaleError::InvalidParameter`] si la serie es degenerada
    ///   (desviación cero, sin días húmedos) o `wet_threshold < 0`.
    pub fn fit(obs: &[f64], model: &[f64], dist: Distribution) -> Result<Self> {
        check_series("obs", obs, MIN_FIT_LEN)?;
        check_series("model", model, MIN_FIT_LEN)?;
        Ok(Self {
            obs: fit_one("obs", obs, dist)?,
            model: fit_one("model", model, dist)?,
        })
    }

    /// Corrige una serie del modelo.
    ///
    /// # Errors
    ///
    /// [`DownscaleError::NonFinite`] si la serie contiene NaN/inf.
    pub fn apply(&self, series: &[f64]) -> Result<Vec<f64>> {
        check_series("series", series, 0)?;
        Ok(series.iter().map(|&x| self.correct_one(x)).collect())
    }

    /// Corrige un único valor: `F_obs⁻¹(F_mod(x))`.
    #[must_use]
    pub fn correct_one(&self, x: f64) -> f64 {
        self.obs.ppf(self.model.cdf(x))
    }
}

fn fit_one(name: &'static str, series: &[f64], dist: Distribution) -> Result<Fitted> {
    match dist {
        Distribution::Normal => {
            let n = series.len() as f64;
            let mean = series.iter().sum::<f64>() / n;
            let var = series.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
            let sd = var.sqrt();
            if sd <= 0.0 {
                return Err(DownscaleError::InvalidParameter {
                    name,
                    value: sd,
                    expected: "desviación estándar > 0",
                });
            }
            Ok(Fitted::Normal { mean, sd })
        }
        Distribution::Gamma { wet_threshold } => {
            if wet_threshold < 0.0 {
                return Err(DownscaleError::InvalidParameter {
                    name: "wet_threshold",
                    value: wet_threshold,
                    expected: ">= 0",
                });
            }
            let wet: Vec<f64> = series
                .iter()
                .copied()
                .filter(|&v| v >= wet_threshold && v > 0.0)
                .collect();
            if wet.len() < MIN_FIT_LEN {
                return Err(DownscaleError::SeriesTooShort {
                    name,
                    len: wet.len(),
                    min: MIN_FIT_LEN,
                });
            }
            let p_dry = 1.0 - wet.len() as f64 / series.len() as f64;
            let (shape, scale) = gamma_mle(&wet);
            Ok(Fitted::MixedGamma {
                p_dry,
                shape,
                scale,
            })
        }
    }
}

/// MLE de Gamma(k, θ) sobre valores estrictamente positivos.
///
/// Inicio con la aproximación de Thom y refinamiento Newton sobre
/// `f(k) = ln k − ψ(k) − s`, con `s = ln(media) − media(ln x)`.
fn gamma_mle(wet: &[f64]) -> (f64, f64) {
    let n = wet.len() as f64;
    let mean = wet.iter().sum::<f64>() / n;
    let mean_ln = wet.iter().map(|v| v.ln()).sum::<f64>() / n;
    let s = (mean.ln() - mean_ln).max(1e-9);

    let mut k = (3.0 - s + ((s - 3.0).powi(2) + 24.0 * s).sqrt()) / (12.0 * s);
    for _ in 0..20 {
        let f = k.ln() - digamma(k) - s;
        let fp = 1.0 / k - trigamma(k);
        let step = f / fp;
        let next = (k - step).max(k * 0.1);
        if (next - k).abs() < 1e-12 * k {
            k = next;
            break;
        }
        k = next;
    }
    (k, mean / k)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// LCG determinista uniforme en (0, 1).
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

    /// Muestra Gamma(k, θ) determinista vía inversa de la CDF.
    fn gamma_sample(seed: u64, n: usize, shape: f64, scale: f64) -> Vec<f64> {
        uniform(seed, n)
            .iter()
            .map(|&u| scale * gamma_p_inv(shape, u))
            .collect()
    }

    fn mean(s: &[f64]) -> f64 {
        s.iter().sum::<f64>() / s.len() as f64
    }

    #[test]
    fn gamma_mle_recovers_parameters() {
        let x = gamma_sample(7, 4000, 2.5, 3.0);
        let (k, theta) = gamma_mle(&x);
        assert!((k - 2.5).abs() < 0.15, "shape = {k}");
        assert!((theta - 3.0).abs() < 0.2, "scale = {theta}");
    }

    #[test]
    fn normal_qm_corrects_mean_and_variance() {
        let u = uniform(3, 4000);
        // Normal aproximada por suma de uniformes (CLT, 12 términos).
        let obs: Vec<f64> = u
            .chunks(4)
            .map(|c| 10.0 + 2.0 * (c.iter().sum::<f64>() - 2.0) * (3.0f64).sqrt())
            .collect();
        let model: Vec<f64> = obs.iter().map(|v| (v - 10.0) * 1.5 + 13.0).collect();

        let pqm = ParametricQuantileMapping::fit(&obs, &model, Distribution::Normal).unwrap();
        let corrected = pqm.apply(&model).unwrap();

        assert!((mean(&corrected) - mean(&obs)).abs() < 1e-6);
        let sd = |s: &[f64]| {
            let m = mean(s);
            (s.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (s.len() as f64 - 1.0)).sqrt()
        };
        assert!((sd(&corrected) - sd(&obs)).abs() < 1e-6);
    }

    #[test]
    fn mixed_gamma_corrects_wet_day_frequency() {
        // Obs: 70% seco, lluvia Gamma(0.8, 8). Modelo: solo 40% seco
        // (llovizna espuria) y lluvia más débil Gamma(0.8, 5).
        let u_obs = uniform(11, 6000);
        let obs: Vec<f64> = u_obs
            .iter()
            .map(|&u| {
                if u < 0.7 {
                    0.0
                } else {
                    8.0 * gamma_p_inv(0.8, (u - 0.7) / 0.3)
                }
            })
            .collect();
        let u_mod = uniform(13, 6000);
        let model: Vec<f64> = u_mod
            .iter()
            .map(|&u| {
                if u < 0.4 {
                    0.0
                } else {
                    5.0 * gamma_p_inv(0.8, (u - 0.4) / 0.6)
                }
            })
            .collect();

        let pqm = ParametricQuantileMapping::fit(
            &obs,
            &model,
            Distribution::Gamma { wet_threshold: 0.1 },
        )
        .unwrap();
        let corrected = pqm.apply(&model).unwrap();

        let dry_frac = |s: &[f64]| s.iter().filter(|&&v| v < 0.1).count() as f64 / s.len() as f64;
        // La frecuencia de días secos corregida debe acercarse a la observada.
        assert!(
            (dry_frac(&corrected) - dry_frac(&obs)).abs() < 0.02,
            "dry corregido = {}, dry obs = {}",
            dry_frac(&corrected),
            dry_frac(&obs)
        );
        // Y la media también.
        assert!((mean(&corrected) - mean(&obs)).abs() < 0.15 * mean(&obs));
        // Sin valores negativos.
        assert!(corrected.iter().all(|&v| v >= 0.0));
    }

    #[test]
    fn rejects_degenerate_series() {
        let flat = vec![5.0; 50];
        let varied: Vec<f64> = (0..50).map(f64::from).collect();
        assert!(matches!(
            ParametricQuantileMapping::fit(&flat, &varied, Distribution::Normal).unwrap_err(),
            DownscaleError::InvalidParameter { .. }
        ));
        // Gamma sin días húmedos suficientes.
        let dry = vec![0.0; 50];
        assert!(matches!(
            ParametricQuantileMapping::fit(
                &dry,
                &varied,
                Distribution::Gamma { wet_threshold: 0.1 }
            )
            .unwrap_err(),
            DownscaleError::SeriesTooShort { .. }
        ));
    }
}
