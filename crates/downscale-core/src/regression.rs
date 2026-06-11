//! Downscaling por regresión lineal múltiple (OLS): predictores de gran
//! escala → predictando local. Resuelve las ecuaciones normales
//! `(XᵀX)β = Xᵀy` con eliminación gaussiana con pivoteo parcial
//! (p pequeño: pocos predictores, matriz (p+1)×(p+1)).

use crate::analog::check_matrix;
use crate::error::{DownscaleError, Result};

/// Modelo de regresión lineal calibrado.
///
/// Predictores como matriz aplanada por filas (igual que
/// [`crate::analog::AnalogDownscaling`]).
///
/// # Ejemplo
///
/// ```
/// use downscale_core::regression::LinearDownscaling;
///
/// // y = 1 + 2·x
/// let predictors = [0.0, 1.0, 2.0, 3.0];
/// let obs = [1.0, 3.0, 5.0, 7.0];
/// let lm = LinearDownscaling::fit(&predictors, 1, &obs).unwrap();
/// assert!((lm.intercept() - 1.0).abs() < 1e-9);
/// assert!((lm.coefs()[0] - 2.0).abs() < 1e-9);
/// assert!((lm.predict(&[10.0]).unwrap()[0] - 21.0).abs() < 1e-9);
/// ```
#[derive(Debug, Clone)]
pub struct LinearDownscaling {
    n_features: usize,
    intercept: f64,
    coefs: Vec<f64>,
    /// R² sobre el período de calibración.
    r2: f64,
}

impl LinearDownscaling {
    /// Ajusta OLS sobre el período de calibración.
    ///
    /// # Errors
    ///
    /// - Dimensiones inconsistentes o NaN/inf (igual que análogos).
    /// - [`DownscaleError::SeriesTooShort`] si hay menos filas que
    ///   parámetros (`n_features + 1`).
    /// - [`DownscaleError::InvalidParameter`] si `XᵀX` es singular
    ///   (predictores colineales o constantes).
    pub fn fit(predictors: &[f64], n_features: usize, obs: &[f64]) -> Result<Self> {
        let rows = check_matrix(predictors, n_features, obs)?;
        let p = n_features + 1; // + intercepto
        if rows < p {
            return Err(DownscaleError::SeriesTooShort {
                name: "obs",
                len: rows,
                min: p,
            });
        }

        // Ecuaciones normales con columna de unos al frente.
        // xtx[(a, b)] = Σ_i X[i,a]·X[i,b], X[i,0] = 1.
        let xval = |i: usize, a: usize| {
            if a == 0 {
                1.0
            } else {
                predictors[i * n_features + (a - 1)]
            }
        };
        let mut xtx = vec![0.0; p * p];
        let mut xty = vec![0.0; p];
        for (i, &y) in obs.iter().enumerate() {
            for a in 0..p {
                let xa = xval(i, a);
                xty[a] += xa * y;
                for b in a..p {
                    xtx[a * p + b] += xa * xval(i, b);
                }
            }
        }
        for a in 0..p {
            for b in 0..a {
                xtx[a * p + b] = xtx[b * p + a];
            }
        }

        let beta = solve(&mut xtx, &mut xty, p)?;

        // R² de calibración.
        let mean_y = obs.iter().sum::<f64>() / rows as f64;
        let mut ss_res = 0.0;
        let mut ss_tot = 0.0;
        for (i, &y) in obs.iter().enumerate() {
            let pred: f64 = (0..p).map(|a| beta[a] * xval(i, a)).sum();
            ss_res += (y - pred).powi(2);
            ss_tot += (y - mean_y).powi(2);
        }
        let r2 = if ss_tot > 0.0 {
            1.0 - ss_res / ss_tot
        } else {
            0.0
        };

        Ok(Self {
            n_features,
            intercept: beta[0],
            coefs: beta[1..].to_vec(),
            r2,
        })
    }

    /// Predice para una secuencia de días (matriz aplanada por filas).
    ///
    /// # Errors
    ///
    /// Largo no múltiplo de `n_features` o valores no finitos.
    pub fn predict(&self, queries: &[f64]) -> Result<Vec<f64>> {
        if !queries.len().is_multiple_of(self.n_features) {
            return Err(DownscaleError::InvalidParameter {
                name: "queries",
                value: queries.len() as f64,
                expected: "largo múltiplo de n_features",
            });
        }
        crate::error::check_series("queries", queries, 0)?;
        Ok(queries
            .chunks_exact(self.n_features)
            .map(|row| {
                self.intercept + row.iter().zip(&self.coefs).map(|(x, c)| x * c).sum::<f64>()
            })
            .collect())
    }

    /// Intercepto ajustado.
    #[must_use]
    pub fn intercept(&self) -> f64 {
        self.intercept
    }

    /// Coeficientes por feature.
    #[must_use]
    pub fn coefs(&self) -> &[f64] {
        &self.coefs
    }

    /// R² sobre el período de calibración.
    #[must_use]
    pub fn r2(&self) -> f64 {
        self.r2
    }
}

/// Resuelve `A·x = b` (A de n×n aplanada, se modifica in place) por
/// eliminación gaussiana con pivoteo parcial.
fn solve(a: &mut [f64], b: &mut [f64], n: usize) -> Result<Vec<f64>> {
    for col in 0..n {
        // Pivoteo parcial.
        let pivot_row = (col..n)
            .max_by(|&i, &j| {
                a[i * n + col]
                    .abs()
                    .partial_cmp(&a[j * n + col].abs())
                    .expect("matriz finita")
            })
            .expect("rango no vacío");
        if a[pivot_row * n + col].abs() < 1e-12 {
            return Err(DownscaleError::InvalidParameter {
                name: "predictors",
                value: col as f64,
                expected: "matriz XᵀX no singular (sin predictores colineales/constantes)",
            });
        }
        if pivot_row != col {
            for k in 0..n {
                a.swap(col * n + k, pivot_row * n + k);
            }
            b.swap(col, pivot_row);
        }
        // Eliminación.
        for row in (col + 1)..n {
            let f = a[row * n + col] / a[col * n + col];
            for k in col..n {
                a[row * n + k] -= f * a[col * n + k];
            }
            b[row] -= f * b[col];
        }
    }
    // Sustitución hacia atrás.
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = ((row + 1)..n).map(|k| a[row * n + k] * x[k]).sum();
        x[row] = (b[row] - s) / a[row * n + row];
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn recovers_known_coefficients_with_noise() {
        // y = 2 + 3·x1 − 1·x2 + ruido pequeño.
        let u = uniform(5, 900);
        let mut predictors = Vec::new();
        let mut obs = Vec::new();
        for c in u.chunks(3) {
            let (x1, x2, noise) = (c[0] * 10.0, c[1] * 4.0, (c[2] - 0.5) * 0.01);
            predictors.extend([x1, x2]);
            obs.push(2.0 + 3.0 * x1 - x2 + noise);
        }
        let lm = LinearDownscaling::fit(&predictors, 2, &obs).unwrap();
        assert!((lm.intercept() - 2.0).abs() < 0.01);
        assert!((lm.coefs()[0] - 3.0).abs() < 0.01);
        assert!((lm.coefs()[1] + 1.0).abs() < 0.01);
        assert!(lm.r2() > 0.999);
    }

    #[test]
    fn rejects_collinear_predictors() {
        // x2 = 2·x1 exacto → XᵀX singular.
        let mut predictors = Vec::new();
        let obs: Vec<f64> = (0..50).map(f64::from).collect();
        for i in 0..50 {
            predictors.extend([f64::from(i), 2.0 * f64::from(i)]);
        }
        assert!(matches!(
            LinearDownscaling::fit(&predictors, 2, &obs).unwrap_err(),
            DownscaleError::InvalidParameter { .. }
        ));
    }

    #[test]
    fn rejects_more_params_than_rows() {
        let predictors = [1.0, 2.0, 3.0, 4.0]; // 2 filas × 2 features
        let obs = [1.0, 2.0];
        assert!(matches!(
            LinearDownscaling::fit(&predictors, 2, &obs).unwrap_err(),
            DownscaleError::SeriesTooShort { .. }
        ));
    }
}
