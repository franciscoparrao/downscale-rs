//! Corrección de sesgo multivariada N-dimensional (MBCn, Cannon 2018).
//!
//! El Schaake shuffle ([`crate::multivariate`]) restaura la dependencia de
//! rangos respecto a una plantilla observada, pero no corrige la estructura
//! de dependencia *completa* (correlaciones entre variables). MBCn sí:
//! alterna rotaciones ortogonales aleatorias con quantile mapping marginal
//! en el espacio rotado, de modo que tras varias iteraciones la muestra
//! corregida adquiere la copula multivariada de las observaciones. El paso
//! final reordena las series con marginales QM-corregidas (unidades
//! correctas) según el orden de rangos de la muestra iterada (dependencia
//! correcta) — un Schaake shuffle.
//!
//! Es determinista: las rotaciones provienen de un PRNG sembrado.

use crate::analog::check_matrix;
use crate::error::{DownscaleError, Result};
use crate::multivariate::schaake_shuffle;
use crate::qm::{Kind, QuantileMapping};

/// Parámetros de [`mbcn`].
#[derive(Debug, Clone, Copy)]
pub struct MbcnOptions {
    /// Iteraciones rotación + corrección marginal (default 30).
    pub n_iterations: usize,
    /// Cuantiles del quantile mapping interno (default 100).
    pub n_quantiles: usize,
    /// Semilla del PRNG de rotaciones (reproducibilidad).
    pub seed: u64,
}

impl Default for MbcnOptions {
    fn default() -> Self {
        Self {
            n_iterations: 30,
            n_quantiles: 100,
            seed: 42,
        }
    }
}

/// Corrige `model` hacia `obs` en marginales y estructura de dependencia
/// multivariada (MBCn). Matrices aplanadas fila-por-día (`n × n_vars`);
/// `obs` y `model` pueden tener distinto número de filas. Devuelve `model`
/// corregido con el mismo tamaño que la entrada.
///
/// # Errors
///
/// Dimensiones inconsistentes, NaN/inf, `n_vars == 0`, o series demasiado
/// cortas para el quantile mapping.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::mbcn::{mbcn, MbcnOptions};
///
/// // 2 variables correlacionadas; modelo con sesgo y dependencia distinta.
/// let obs = [0.0, 0.0, 1.0, 0.9, 2.0, 2.1, 3.0, 2.9, 4.0, 4.1];
/// let model = [5.0, 0.0, 6.0, 1.0, 7.0, 2.0, 8.0, 3.0, 9.0, 4.0];
/// let out = mbcn(&obs, &model, 2, &MbcnOptions::default()).unwrap();
/// assert_eq!(out.len(), model.len());
/// ```
pub fn mbcn(obs: &[f64], model: &[f64], n_vars: usize, opts: &MbcnOptions) -> Result<Vec<f64>> {
    // Valida ambas matrices (reusa la validación dimensional de análogos
    // con una pseudo-observación del largo de cada una).
    let n_obs = validate_matrix(obs, n_vars)?;
    let n_mod = validate_matrix(model, n_vars)?;
    let d = n_vars;

    // 1. Quantile mapping marginal por variable (unidades originales).
    let mut model_qdm = model.to_vec();
    for j in 0..d {
        let obs_col = column(obs, n_obs, d, j);
        let mod_col = column(model, n_mod, d, j);
        let qm = QuantileMapping::fit(&obs_col, &mod_col, opts.n_quantiles, Kind::Additive)?;
        let corr = qm.apply(&mod_col)?;
        for (i, &v) in corr.iter().enumerate() {
            model_qdm[i * d + j] = v;
        }
    }

    // 2. Estandariza obs y la muestra actual para la iteración (las
    //    rotaciones requieren escalas comparables).
    let (om, osd) = col_stats(obs, n_obs, d);
    let obs_std = standardize(obs, n_obs, d, &om, &osd);
    let (mm, msd) = col_stats(&model_qdm, n_mod, d);
    let mut cur = standardize(&model_qdm, n_mod, d, &mm, &msd);

    // 3. Itera: rota, corrige cada eje marginalmente, rota de vuelta.
    let mut rng = Pcg::new(opts.seed);
    for _ in 0..opts.n_iterations {
        let r = random_orthogonal(d, &mut rng);
        let obs_r = matmul(&obs_std, &r, n_obs, d, false);
        let mut cur_r = matmul(&cur, &r, n_mod, d, false);
        for j in 0..d {
            let obs_col = column(&obs_r, n_obs, d, j);
            let mod_col = column(&cur_r, n_mod, d, j);
            let qm = QuantileMapping::fit(&obs_col, &mod_col, opts.n_quantiles, Kind::Additive)?;
            let corr = qm.apply(&mod_col)?;
            for (i, &v) in corr.iter().enumerate() {
                cur_r[i * d + j] = v;
            }
        }
        // Rota de vuelta con Rᵀ (matriz ortogonal: R⁻¹ = Rᵀ).
        cur = matmul(&cur_r, &r, n_mod, d, true);
    }

    // 4. Schaake: marginales de model_qdm con el orden de rangos de la
    //    muestra iterada (que porta la dependencia de obs).
    schaake_shuffle(&cur, &model_qdm, d)
}

/// Valida una matriz aplanada y devuelve su número de filas.
fn validate_matrix(matrix: &[f64], n_vars: usize) -> Result<usize> {
    if n_vars == 0 || !matrix.len().is_multiple_of(n_vars.max(1)) {
        return Err(DownscaleError::InvalidParameter {
            name: "n_vars",
            value: n_vars as f64,
            expected: "n_vars >= 1 y largo múltiplo de n_vars",
        });
    }
    let rows = matrix.len() / n_vars;
    let pseudo = vec![0.0; rows];
    check_matrix(matrix, n_vars, &pseudo)
}

/// Columna `j` de una matriz aplanada `n × d`.
fn column(m: &[f64], n: usize, d: usize, j: usize) -> Vec<f64> {
    (0..n).map(|i| m[i * d + j]).collect()
}

/// Media y desviación estándar por columna.
fn col_stats(a: &[f64], n: usize, d: usize) -> (Vec<f64>, Vec<f64>) {
    let mut mean = vec![0.0; d];
    let mut sd = vec![0.0; d];
    for j in 0..d {
        let m = (0..n).map(|i| a[i * d + j]).sum::<f64>() / n as f64;
        let var =
            (0..n).map(|i| (a[i * d + j] - m).powi(2)).sum::<f64>() / (n as f64 - 1.0).max(1.0);
        mean[j] = m;
        sd[j] = if var > 0.0 { var.sqrt() } else { 1.0 };
    }
    (mean, sd)
}

/// Estandariza (z-score) una matriz con estadísticas dadas por columna.
fn standardize(a: &[f64], n: usize, d: usize, mean: &[f64], sd: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; n * d];
    for i in 0..n {
        for j in 0..d {
            out[i * d + j] = (a[i * d + j] - mean[j]) / sd[j];
        }
    }
    out
}

/// Producto `A · R` (o `A · Rᵀ` si `transpose`), con `A` de `n × d` y `R`
/// de `d × d`, todo fila-mayor.
fn matmul(a: &[f64], r: &[f64], n: usize, d: usize, transpose: bool) -> Vec<f64> {
    let mut out = vec![0.0; n * d];
    for i in 0..n {
        for j in 0..d {
            let mut s = 0.0;
            for k in 0..d {
                let rkj = if transpose {
                    r[j * d + k]
                } else {
                    r[k * d + j]
                };
                s += a[i * d + k] * rkj;
            }
            out[i * d + j] = s;
        }
    }
    out
}

/// Matriz ortogonal `d × d` uniforme-ish: Gram–Schmidt sobre columnas
/// gaussianas. Fila-mayor.
fn random_orthogonal(d: usize, rng: &mut Pcg) -> Vec<f64> {
    let mut m = vec![0.0; d * d];
    for col in 0..d {
        let mut v: Vec<f64> = (0..d).map(|_| rng.gaussian()).collect();
        // Ortogonaliza contra las columnas ya fijadas.
        for prev in 0..col {
            let dot: f64 = (0..d).map(|r| v[r] * m[r * d + prev]).sum();
            for (r, vr) in v.iter_mut().enumerate() {
                *vr -= dot * m[r * d + prev];
            }
        }
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        let norm = if norm < 1e-12 { 1.0 } else { norm };
        for (r, &vr) in v.iter().enumerate() {
            m[r * d + col] = vr / norm;
        }
    }
    m
}

/// PRNG PCG-XSH-RR de 64 bits — determinista, sin dependencias.
struct Pcg(u64);

impl Pcg {
    fn new(seed: u64) -> Self {
        let mut p = Pcg(0);
        p.0 = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        p.next_u32(); // mezcla inicial
        p
    }

    fn next_u32(&mut self) -> u32 {
        let old = self.0;
        self.0 = old
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let xorshifted = (((old >> 18) ^ old) >> 27) as u32;
        let rot = (old >> 59) as u32;
        xorshifted.rotate_right(rot)
    }

    /// Uniforme en (0, 1).
    fn next_f64(&mut self) -> f64 {
        (f64::from(self.next_u32()) + 0.5) / f64::from(u32::MAX)
    }

    /// Normal estándar por Box–Muller.
    fn gaussian(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::ks_statistic;

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

    fn gaussians(seed: u64, n: usize) -> Vec<f64> {
        let u = uniform(seed, 2 * n);
        (0..n)
            .map(|i| {
                let u1 = u[2 * i].max(1e-12);
                (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u[2 * i + 1]).cos()
            })
            .collect()
    }

    fn correlation(a: &[f64], b: &[f64]) -> f64 {
        let n = a.len() as f64;
        let ma = a.iter().sum::<f64>() / n;
        let mb = b.iter().sum::<f64>() / n;
        let cov: f64 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum();
        let va: f64 = a.iter().map(|x| (x - ma).powi(2)).sum();
        let vb: f64 = b.iter().map(|y| (y - mb).powi(2)).sum();
        cov / (va * vb).sqrt()
    }

    fn cols(m: &[f64], d: usize) -> Vec<Vec<f64>> {
        (0..d).map(|j| column(m, m.len() / d, d, j)).collect()
    }

    /// Construye obs con correlación ~+0.8 y un modelo con correlación ~−0.3
    /// y sesgo; MBCn debe recuperar la correlación de obs y sus marginales.
    fn make_case() -> (Vec<f64>, Vec<f64>) {
        let n = 3000;
        let (z1, z2) = (gaussians(1, n), gaussians(2, n));
        let (w1, w2) = (gaussians(3, n), gaussians(4, n));
        let mut obs = Vec::with_capacity(n * 2);
        let mut model = Vec::with_capacity(n * 2);
        for i in 0..n {
            // obs: corr ~0.8.
            obs.push(10.0 + 2.0 * z1[i]);
            obs.push(20.0 + 0.8 * 3.0 * z1[i] + 0.6 * 3.0 * z2[i]);
            // model: corr ~−0.3, con sesgo en media y escala.
            model.push(13.0 + 3.0 * w1[i]);
            model.push(17.0 - 0.3 * 4.0 * w1[i] + 0.95 * 4.0 * w2[i]);
        }
        (obs, model)
    }

    #[test]
    fn recovers_dependence_and_marginals() {
        let (obs, model) = make_case();
        let out = mbcn(&obs, &model, 2, &MbcnOptions::default()).unwrap();

        let oc = cols(&obs, 2);
        let mc = cols(&model, 2);
        let rc = cols(&out, 2);

        let corr_obs = correlation(&oc[0], &oc[1]);
        let corr_model = correlation(&mc[0], &mc[1]);
        let corr_out = correlation(&rc[0], &rc[1]);
        // La correlación corregida se acerca a la observada, lejos de la del modelo.
        assert!(
            (corr_out - corr_obs).abs() < 0.1,
            "corr obs={corr_obs:.3} model={corr_model:.3} out={corr_out:.3}"
        );
        assert!((corr_out - corr_model).abs() > 0.3);

        // Marginales: cada variable corregida ≈ la observada (KS bajo).
        for j in 0..2 {
            let ks = ks_statistic(&rc[j], &oc[j]).unwrap();
            assert!(ks < 0.06, "KS marginal var {j} = {ks}");
        }
    }

    #[test]
    fn is_deterministic_with_seed() {
        let (obs, model) = make_case();
        let a = mbcn(&obs, &model, 2, &MbcnOptions::default()).unwrap();
        let b = mbcn(&obs, &model, 2, &MbcnOptions::default()).unwrap();
        assert_eq!(a, b);
        // Semilla distinta → resultado distinto (pero válido).
        let opts2 = MbcnOptions {
            seed: 7,
            ..MbcnOptions::default()
        };
        let c = mbcn(&obs, &model, 2, &opts2).unwrap();
        assert_ne!(a, c);
    }

    #[test]
    fn random_orthogonal_is_orthonormal() {
        let mut rng = Pcg::new(123);
        let d = 4;
        let r = random_orthogonal(d, &mut rng);
        // Rᵀ·R = I.
        for i in 0..d {
            for j in 0..d {
                let dot: f64 = (0..d).map(|k| r[k * d + i] * r[k * d + j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dot - expected).abs() < 1e-10, "({i},{j}) = {dot}");
            }
        }
    }

    #[test]
    fn rejects_bad_input() {
        assert!(mbcn(&[1.0, 2.0, 3.0], &[1.0, 2.0], 2, &MbcnOptions::default()).is_err());
        assert!(mbcn(&[1.0, 2.0], &[1.0, 2.0], 0, &MbcnOptions::default()).is_err());
    }
}
