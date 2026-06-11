//! Métricas de validación: RMSE, sesgo medio, estadístico KS de dos
//! muestras y sesgo por cuantil.

use crate::error::{Result, check_same_len, check_series};
use crate::qm::quantile_sorted;

/// Raíz del error cuadrático medio entre series pareadas en el tiempo.
///
/// # Errors
///
/// Series vacías, con NaN/inf o de largos distintos.
///
/// # Ejemplo
///
/// ```
/// let r = downscale_core::metrics::rmse(&[1.0, 2.0], &[1.0, 4.0]).unwrap();
/// assert!((r - 2.0_f64.sqrt()).abs() < 1e-12);
/// ```
pub fn rmse(sim: &[f64], obs: &[f64]) -> Result<f64> {
    check_series("sim", sim, 1)?;
    check_series("obs", obs, 1)?;
    check_same_len("sim", sim, "obs", obs)?;
    let sse: f64 = sim.iter().zip(obs).map(|(s, o)| (s - o).powi(2)).sum();
    Ok((sse / sim.len() as f64).sqrt())
}

/// Sesgo medio: `mean(sim) - mean(obs)`. No requiere series pareadas.
///
/// # Errors
///
/// Series vacías o con NaN/inf.
pub fn mean_bias(sim: &[f64], obs: &[f64]) -> Result<f64> {
    check_series("sim", sim, 1)?;
    check_series("obs", obs, 1)?;
    let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
    Ok(mean(sim) - mean(obs))
}

/// Estadístico de Kolmogorov–Smirnov de dos muestras:
/// `D = sup |F_sim(x) - F_obs(x)|` sobre las CDFs empíricas.
///
/// `D = 0` indica distribuciones idénticas; `D = 1`, soportes disjuntos.
///
/// # Errors
///
/// Series vacías o con NaN/inf.
pub fn ks_statistic(sim: &[f64], obs: &[f64]) -> Result<f64> {
    check_series("sim", sim, 1)?;
    check_series("obs", obs, 1)?;

    let mut a = sim.to_vec();
    let mut b = obs.to_vec();
    a.sort_unstable_by(|x, y| x.partial_cmp(y).expect("serie validada sin NaN"));
    b.sort_unstable_by(|x, y| x.partial_cmp(y).expect("serie validada sin NaN"));

    let (na, nb) = (a.len() as f64, b.len() as f64);
    let (mut i, mut j) = (0usize, 0usize);
    let mut d_max = 0.0f64;
    while i < a.len() && j < b.len() {
        let x = a[i].min(b[j]);
        while i < a.len() && a[i] <= x {
            i += 1;
        }
        while j < b.len() && b[j] <= x {
            j += 1;
        }
        d_max = d_max.max((i as f64 / na - j as f64 / nb).abs());
    }
    Ok(d_max)
}

/// Sesgo en un cuantil dado.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuantileBias {
    /// Probabilidad del cuantil, en \[0, 1\].
    pub prob: f64,
    /// Cuantil de la serie simulada/corregida.
    pub sim: f64,
    /// Cuantil de la serie observada.
    pub obs: f64,
    /// Diferencia `sim - obs`.
    pub bias: f64,
}

/// Probabilidades por defecto para [`quantile_bias`]: deciles más colas
/// (P5, P95, P99), relevantes para extremos hidroclimáticos.
pub const DEFAULT_QUANTILE_PROBS: [f64; 12] = [
    0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99,
];

/// Sesgo por cuantil entre dos series (no necesitan estar pareadas).
///
/// # Errors
///
/// Series vacías, con NaN/inf, o alguna probabilidad fuera de \[0, 1\].
pub fn quantile_bias(sim: &[f64], obs: &[f64], probs: &[f64]) -> Result<Vec<QuantileBias>> {
    check_series("sim", sim, 1)?;
    check_series("obs", obs, 1)?;
    check_series("probs", probs, 1)?;
    if let Some(&p) = probs.iter().find(|p| !(0.0..=1.0).contains(*p)) {
        return Err(crate::error::DownscaleError::InvalidParameter {
            name: "probs",
            value: p,
            expected: "en [0, 1]",
        });
    }

    let mut sim_sorted = sim.to_vec();
    let mut obs_sorted = obs.to_vec();
    sim_sorted.sort_unstable_by(|x, y| x.partial_cmp(y).expect("serie validada sin NaN"));
    obs_sorted.sort_unstable_by(|x, y| x.partial_cmp(y).expect("serie validada sin NaN"));

    Ok(probs
        .iter()
        .map(|&prob| {
            let s = quantile_sorted(&sim_sorted, prob);
            let o = quantile_sorted(&obs_sorted, prob);
            QuantileBias {
                prob,
                sim: s,
                obs: o,
                bias: s - o,
            }
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::DownscaleError;

    #[test]
    fn rmse_known_value() {
        let r = rmse(&[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0]).unwrap();
        assert_eq!(r, 0.0);
        let r = rmse(&[0.0, 0.0], &[3.0, 4.0]).unwrap();
        // sqrt((9+16)/2) = sqrt(12.5)
        assert!((r - 12.5_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn rmse_rejects_length_mismatch() {
        let err = rmse(&[1.0], &[1.0, 2.0]).unwrap_err();
        assert!(matches!(err, DownscaleError::LengthMismatch { .. }));
    }

    #[test]
    fn mean_bias_known_value() {
        let b = mean_bias(&[3.0, 5.0], &[1.0, 1.0]).unwrap();
        assert_eq!(b, 3.0);
    }

    #[test]
    fn ks_identical_is_zero() {
        let s = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(ks_statistic(&s, &s).unwrap(), 0.0);
    }

    #[test]
    fn ks_disjoint_is_one() {
        let a = [1.0, 2.0, 3.0];
        let b = [10.0, 11.0, 12.0];
        assert_eq!(ks_statistic(&a, &b).unwrap(), 1.0);
    }

    #[test]
    fn ks_known_half() {
        // F_a salta a 1 en 1.0; F_b en 1.0 vale 0.5 → D = 0.5.
        let a = [1.0, 1.0];
        let b = [1.0, 2.0];
        assert!((ks_statistic(&a, &b).unwrap() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn ks_different_lengths_allowed() {
        let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let b = [1.0, 3.0, 5.0];
        let d = ks_statistic(&a, &b).unwrap();
        assert!((0.0..=1.0).contains(&d));
    }

    #[test]
    fn quantile_bias_detects_shift() {
        let obs: Vec<f64> = (0..101).map(f64::from).collect();
        let sim: Vec<f64> = obs.iter().map(|v| v + 5.0).collect();
        let qb = quantile_bias(&sim, &obs, &[0.5]).unwrap();
        assert_eq!(qb.len(), 1);
        assert!((qb[0].bias - 5.0).abs() < 1e-9);
    }

    #[test]
    fn quantile_bias_rejects_bad_prob() {
        let err = quantile_bias(&[1.0], &[1.0], &[1.5]).unwrap_err();
        assert!(matches!(
            err,
            DownscaleError::InvalidParameter { name: "probs", .. }
        ));
    }
}
