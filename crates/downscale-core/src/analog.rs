//! Downscaling por análogos: para cada día objetivo se buscan los `k`
//! días más parecidos en el archivo de predictores de gran escala
//! (período de calibración) y se predice el valor local como media de las
//! observaciones de esos análogos, ponderada por distancia inversa.
//!
//! Los predictores se estandarizan (z-score) con las estadísticas del
//! archivo y se comparan con distancia euclidiana. Búsqueda por fuerza
//! bruta: O(n) por consulta, suficiente para series diarias (~10⁴ días).

use crate::error::{DownscaleError, Result, check_series};

/// Mínimo de filas del archivo.
const MIN_FIT_ROWS: usize = 2;

/// Modelo de análogos calibrado.
///
/// Los predictores se pasan como matriz aplanada por filas: el día `i`
/// ocupa `predictors[i*n_features .. (i+1)*n_features]`.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::analog::AnalogDownscaling;
///
/// // Archivo: 1 predictor; el predictando es ~2× el predictor.
/// let predictors = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let obs = [2.1, 3.9, 6.2, 8.0, 9.8];
/// let ad = AnalogDownscaling::fit(&predictors, 1, &obs, 2).unwrap();
///
/// // Consulta entre 2.0 y 3.0 → mezcla de sus análogos.
/// let y = ad.predict_one(&[2.5]).unwrap();
/// assert!((3.9..=6.2).contains(&y));
/// ```
#[derive(Debug, Clone)]
pub struct AnalogDownscaling {
    n_features: usize,
    k: usize,
    /// Media por feature (para estandarizar consultas).
    means: Vec<f64>,
    /// Desviación por feature; features constantes quedan con 1.0
    /// (no discriminan, distancia 0).
    sds: Vec<f64>,
    /// Archivo estandarizado, aplanado por filas.
    archive: Vec<f64>,
    /// Observaciones locales del archivo.
    obs: Vec<f64>,
}

impl AnalogDownscaling {
    /// Calibra con el archivo de predictores y observaciones pareadas.
    ///
    /// # Errors
    ///
    /// - [`DownscaleError::InvalidParameter`] si `n_features == 0`,
    ///   `predictors.len()` no es múltiplo de `n_features`, el número de
    ///   filas no coincide con `obs.len()`, o `k` está fuera de
    ///   `1..=filas`.
    /// - [`DownscaleError::SeriesTooShort`] / [`DownscaleError::NonFinite`]
    ///   sobre las series de entrada.
    pub fn fit(predictors: &[f64], n_features: usize, obs: &[f64], k: usize) -> Result<Self> {
        let rows = check_matrix(predictors, n_features, obs)?;
        if k == 0 || k > rows {
            return Err(DownscaleError::InvalidParameter {
                name: "k",
                value: k as f64,
                expected: "1 <= k <= filas del archivo",
            });
        }

        // Estadísticas por feature para el z-score.
        let mut means = vec![0.0; n_features];
        let mut sds = vec![0.0; n_features];
        for (j, mean) in means.iter_mut().enumerate() {
            *mean = (0..rows)
                .map(|i| predictors[i * n_features + j])
                .sum::<f64>()
                / rows as f64;
        }
        for (j, sd) in sds.iter_mut().enumerate() {
            let var = (0..rows)
                .map(|i| (predictors[i * n_features + j] - means[j]).powi(2))
                .sum::<f64>()
                / (rows as f64 - 1.0);
            *sd = if var > 0.0 { var.sqrt() } else { 1.0 };
        }

        let archive: Vec<f64> = predictors
            .iter()
            .enumerate()
            .map(|(idx, &v)| {
                let j = idx % n_features;
                (v - means[j]) / sds[j]
            })
            .collect();

        Ok(Self {
            n_features,
            k,
            means,
            sds,
            archive,
            obs: obs.to_vec(),
        })
    }

    /// Predice el valor local para un vector de predictores.
    ///
    /// # Errors
    ///
    /// [`DownscaleError::InvalidParameter`] si `query.len() != n_features`;
    /// [`DownscaleError::NonFinite`] si contiene NaN/inf.
    pub fn predict_one(&self, query: &[f64]) -> Result<f64> {
        check_series("query", query, 0)?;
        if query.len() != self.n_features {
            return Err(DownscaleError::InvalidParameter {
                name: "query",
                value: query.len() as f64,
                expected: "largo == n_features",
            });
        }
        let q: Vec<f64> = query
            .iter()
            .zip(self.means.iter().zip(&self.sds))
            .map(|(&v, (&m, &s))| (v - m) / s)
            .collect();

        // Distancias a todo el archivo y selección de los k menores.
        let rows = self.obs.len();
        let mut dist: Vec<(f64, f64)> = (0..rows)
            .map(|i| {
                let row = &self.archive[i * self.n_features..(i + 1) * self.n_features];
                let d2: f64 = row.iter().zip(&q).map(|(a, b)| (a - b).powi(2)).sum();
                (d2, self.obs[i])
            })
            .collect();
        dist.select_nth_unstable_by(self.k - 1, |a, b| {
            a.0.partial_cmp(&b.0).expect("distancias finitas")
        });

        // Media ponderada por distancia inversa (eps evita división por 0).
        let mut num = 0.0;
        let mut den = 0.0;
        for &(d2, y) in &dist[..self.k] {
            let w = 1.0 / (d2.sqrt() + 1e-12);
            num += w * y;
            den += w;
        }
        Ok(num / den)
    }

    /// Predice una secuencia de días (matriz aplanada por filas).
    ///
    /// # Errors
    ///
    /// Igual que [`AnalogDownscaling::predict_one`], más largo no múltiplo
    /// de `n_features`.
    pub fn predict(&self, queries: &[f64]) -> Result<Vec<f64>> {
        if !queries.len().is_multiple_of(self.n_features) {
            return Err(DownscaleError::InvalidParameter {
                name: "queries",
                value: queries.len() as f64,
                expected: "largo múltiplo de n_features",
            });
        }
        queries
            .chunks_exact(self.n_features)
            .map(|q| self.predict_one(q))
            .collect()
    }

    /// Número de análogos usados.
    #[must_use]
    pub fn k(&self) -> usize {
        self.k
    }
}

/// Valida la matriz aplanada y devuelve el número de filas.
pub(crate) fn check_matrix(predictors: &[f64], n_features: usize, obs: &[f64]) -> Result<usize> {
    if n_features == 0 {
        return Err(DownscaleError::InvalidParameter {
            name: "n_features",
            value: 0.0,
            expected: ">= 1",
        });
    }
    if !predictors.len().is_multiple_of(n_features) {
        return Err(DownscaleError::InvalidParameter {
            name: "predictors",
            value: predictors.len() as f64,
            expected: "largo múltiplo de n_features",
        });
    }
    check_series("predictors", predictors, MIN_FIT_ROWS * n_features)?;
    check_series("obs", obs, MIN_FIT_ROWS)?;
    let rows = predictors.len() / n_features;
    if rows != obs.len() {
        return Err(DownscaleError::LengthMismatch {
            left_name: "predictors (filas)",
            left: rows,
            right_name: "obs",
            right: obs.len(),
        });
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_match_with_k1_returns_archive_obs() {
        let predictors = [0.0, 10.0, 20.0, 30.0];
        let obs = [1.0, 2.0, 3.0, 4.0];
        let ad = AnalogDownscaling::fit(&predictors, 1, &obs, 1).unwrap();
        for (p, o) in predictors.iter().zip(&obs) {
            assert!((ad.predict_one(&[*p]).unwrap() - o).abs() < 1e-9);
        }
    }

    #[test]
    fn recovers_smooth_relationship() {
        // y = sin(x) con archivo denso; consulta interpolada.
        let predictors: Vec<f64> = (0..500).map(|i| f64::from(i) * 0.01).collect();
        let obs: Vec<f64> = predictors.iter().map(|x| x.sin()).collect();
        let ad = AnalogDownscaling::fit(&predictors, 1, &obs, 3).unwrap();
        for &x in &[0.5, 1.7, 3.3, 4.9] {
            let y = ad.predict_one(&[x]).unwrap();
            assert!((y - x.sin()).abs() < 0.02, "x={x}: {y} vs {}", x.sin());
        }
    }

    #[test]
    fn multivariate_standardization_balances_scales() {
        // Feature 2 con escala 1000×: sin estandarizar dominaría.
        // y depende solo de la feature 1.
        let mut predictors = Vec::new();
        let mut obs = Vec::new();
        for i in 0..200 {
            let x1 = f64::from(i % 20);
            let x2 = f64::from(i % 7) * 1000.0;
            predictors.extend([x1, x2]);
            obs.push(x1 * 2.0);
        }
        let ad = AnalogDownscaling::fit(&predictors, 2, &obs, 5).unwrap();
        let y = ad.predict_one(&[10.0, 3000.0]).unwrap();
        assert!((y - 20.0).abs() < 2.0, "y = {y}");
    }

    #[test]
    fn constant_feature_is_ignored() {
        let predictors = [1.0, 5.0, 2.0, 5.0, 3.0, 5.0, 4.0, 5.0];
        let obs = [10.0, 20.0, 30.0, 40.0];
        let ad = AnalogDownscaling::fit(&predictors, 2, &obs, 1).unwrap();
        assert!((ad.predict_one(&[2.0, 5.0]).unwrap() - 20.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_bad_dimensions() {
        let obs = [1.0, 2.0];
        assert!(AnalogDownscaling::fit(&[1.0, 2.0, 3.0], 2, &obs, 1).is_err());
        assert!(AnalogDownscaling::fit(&[1.0, 2.0], 0, &obs, 1).is_err());
        assert!(AnalogDownscaling::fit(&[1.0, 2.0], 1, &obs, 3).is_err());
        let ad = AnalogDownscaling::fit(&[1.0, 2.0], 1, &obs, 1).unwrap();
        assert!(ad.predict_one(&[1.0, 2.0]).is_err());
    }
}
