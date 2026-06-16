//! Downscaling por análogos: para cada día objetivo se buscan los `k`
//! días más parecidos en el archivo de predictores de gran escala
//! (período de calibración) y se predice el valor local como media de las
//! observaciones de esos análogos, ponderada por distancia inversa.
//!
//! Los predictores se estandarizan (z-score) con las estadísticas del
//! archivo y se comparan con distancia euclidiana. La búsqueda de los `k`
//! vecinos usa un k-d tree (O(log n) promedio por consulta), exacto: da
//! el mismo conjunto de análogos que la fuerza bruta.

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
    /// k-d tree sobre el archivo estandarizado.
    tree: KdTree,
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

        let standardized: Vec<f64> = predictors
            .iter()
            .enumerate()
            .map(|(idx, &v)| {
                let j = idx % n_features;
                (v - means[j]) / sds[j]
            })
            .collect();

        let tree = KdTree::build(&standardized, obs, n_features);

        Ok(Self {
            n_features,
            k,
            means,
            sds,
            tree,
        })
    }

    /// Estandariza una consulta con las estadísticas del archivo.
    fn standardize(&self, query: &[f64]) -> Vec<f64> {
        query
            .iter()
            .zip(self.means.iter().zip(&self.sds))
            .map(|(&v, (&m, &s))| (v - m) / s)
            .collect()
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
        let q = self.standardize(query);
        let neighbors = self.tree.knn(&q, self.k);
        Ok(idw_mean(&neighbors))
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

/// Media ponderada por distancia inversa sobre vecinos `(dist², obs)`.
/// El `eps` evita la división por cero cuando un análogo es exacto.
fn idw_mean(neighbors: &[(f64, f64)]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for &(d2, y) in neighbors {
        let w = 1.0 / (d2.sqrt() + 1e-12);
        num += w * y;
        den += w;
    }
    num / den
}

/// k-d tree de búsqueda de vecinos más cercanos sobre puntos
/// `n_features`-dimensionales, con el valor observado asociado a cada punto.
#[derive(Debug, Clone)]
struct KdTree {
    nodes: Vec<KdNode>,
    root: Option<usize>,
}

#[derive(Debug, Clone)]
struct KdNode {
    /// Coordenadas estandarizadas del punto.
    point: Vec<f64>,
    /// Observación local asociada.
    obs: f64,
    /// Eje de partición de este nodo.
    axis: usize,
    left: Option<usize>,
    right: Option<usize>,
}

impl KdTree {
    /// Construye el árbol particionando recursivamente por la mediana del
    /// eje que rota con la profundidad (quickselect, O(n log n)).
    fn build(points: &[f64], obs: &[f64], n_features: usize) -> Self {
        let n = obs.len();
        let mut indices: Vec<usize> = (0..n).collect();
        let mut nodes = Vec::with_capacity(n);
        let root = build_node(&mut indices, points, obs, n_features, 0, &mut nodes);
        Self { nodes, root }
    }

    /// Devuelve los `k` vecinos más cercanos a `query` como pares
    /// `(dist², obs)` ordenados por distancia ascendente.
    fn knn(&self, query: &[f64], k: usize) -> Vec<(f64, f64)> {
        let mut best: Vec<(f64, f64)> = Vec::with_capacity(k + 1);
        self.search(self.root, query, k, &mut best);
        best
    }

    fn search(&self, node: Option<usize>, query: &[f64], k: usize, best: &mut Vec<(f64, f64)>) {
        let Some(ni) = node else {
            return;
        };
        let nd = &self.nodes[ni];
        let d2: f64 = nd
            .point
            .iter()
            .zip(query)
            .map(|(a, b)| (a - b).powi(2))
            .sum();
        // Inserta manteniendo `best` ordenado ascendente, tamaño <= k.
        if best.len() < k {
            let pos = best.partition_point(|&(dd, _)| dd <= d2);
            best.insert(pos, (d2, nd.obs));
        } else if d2 < best[k - 1].0 {
            let pos = best.partition_point(|&(dd, _)| dd <= d2);
            best.insert(pos, (d2, nd.obs));
            best.truncate(k);
        }

        let diff = query[nd.axis] - nd.point[nd.axis];
        let (near, far) = if diff <= 0.0 {
            (nd.left, nd.right)
        } else {
            (nd.right, nd.left)
        };
        self.search(near, query, k, best);
        // Visita el subárbol lejano solo si puede contener un vecino mejor.
        let worst = best.last().map_or(f64::INFINITY, |&(d, _)| d);
        if best.len() < k || diff * diff < worst {
            self.search(far, query, k, best);
        }
    }

    /// k-NN por fuerza bruta (referencia para el test de equivalencia).
    #[cfg(test)]
    fn knn_bruteforce(&self, query: &[f64], k: usize) -> Vec<(f64, f64)> {
        let mut all: Vec<(f64, f64)> = self
            .nodes
            .iter()
            .map(|nd| {
                let d2 = nd
                    .point
                    .iter()
                    .zip(query)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum();
                (d2, nd.obs)
            })
            .collect();
        all.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("sin NaN"));
        all.truncate(k);
        all
    }
}

/// Construcción recursiva: elige la mediana del `axis` actual como nodo y
/// reparte el resto en los subárboles. Devuelve el índice del nodo creado.
fn build_node(
    indices: &mut [usize],
    points: &[f64],
    obs: &[f64],
    n_features: usize,
    depth: usize,
    nodes: &mut Vec<KdNode>,
) -> Option<usize> {
    if indices.is_empty() {
        return None;
    }
    let axis = depth % n_features;
    let mid = indices.len() / 2;
    indices.select_nth_unstable_by(mid, |&a, &b| {
        points[a * n_features + axis]
            .partial_cmp(&points[b * n_features + axis])
            .expect("archivo validado sin NaN")
    });
    let median = indices[mid];
    let (left_idx, rest) = indices.split_at_mut(mid);
    let right_idx = &mut rest[1..];
    let left = build_node(left_idx, points, obs, n_features, depth + 1, nodes);
    let right = build_node(right_idx, points, obs, n_features, depth + 1, nodes);
    let point = points[median * n_features..(median + 1) * n_features].to_vec();
    nodes.push(KdNode {
        point,
        obs: obs[median],
        axis,
        left,
        right,
    });
    Some(nodes.len() - 1)
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
    fn kdtree_knn_matches_bruteforce() {
        // Para datos continuos (sin empates) el k-d tree devuelve exactamente
        // el mismo conjunto de vecinos que la fuerza bruta.
        let n_features = 4;
        let n = 600;
        let points = uniform(1, n * n_features);
        let obs = uniform(2, n);
        let tree = KdTree::build(&points, &obs, n_features);

        for q in 0..50 {
            let query = uniform(100 + q as u64, n_features);
            for &k in &[1usize, 5, 10, 25] {
                let fast = tree.knn(&query, k);
                let slow = tree.knn_bruteforce(&query, k);
                assert_eq!(fast.len(), k);
                for (a, b) in fast.iter().zip(&slow) {
                    assert!((a.0 - b.0).abs() < 1e-12, "dist² difiere: {a:?} vs {b:?}");
                    assert!((a.1 - b.1).abs() < 1e-12, "obs difiere: {a:?} vs {b:?}");
                }
            }
        }
    }

    #[test]
    fn predict_matches_bruteforce_idw() {
        // El downscaling completo coincide con la referencia por fuerza bruta.
        let n_features = 3;
        let n = 400;
        let predictors = uniform(7, n * n_features);
        let obs = uniform(8, n);
        let k = 8;
        let ad = AnalogDownscaling::fit(&predictors, n_features, &obs, k).unwrap();

        for q in 0..40 {
            let query = uniform(500 + q as u64, n_features);
            let got = ad.predict_one(&query).unwrap();
            // Referencia: fuerza bruta sobre el archivo estandarizado.
            let qs = ad.standardize(&query);
            let want = idw_mean(&ad.tree.knn_bruteforce(&qs, k));
            assert!((got - want).abs() < 1e-9, "got {got}, want {want}");
        }
    }
}
