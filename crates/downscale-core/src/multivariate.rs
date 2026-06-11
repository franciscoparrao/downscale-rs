//! Corrección multivariada por reordenamiento de rangos (Schaake shuffle,
//! Clark et al. 2004).
//!
//! La corrección univariada (QM/QDM por variable) destruye la estructura
//! de dependencia entre variables (p. ej. correlación precipitación–
//! temperatura). El Schaake shuffle la restaura: reordena los valores
//! corregidos de cada variable para que sus **rangos** reproduzcan los de
//! una plantilla observada, preservando exactamente las distribuciones
//! marginales corregidas.

use crate::analog::check_matrix;
use crate::error::{DownscaleError, Result};

/// Reordena `corrected` para imitar la estructura de rangos de `template`.
///
/// Ambas son matrices aplanadas fila-por-día (`n_days × n_vars`, mismo
/// layout que análogos/regresión). Para cada variable `j`:
/// `salida[i][j] = sorted(corrected[:,j])[rango de template[i][j]]`.
///
/// Garantías: las marginales de la salida son exactamente las de
/// `corrected` (mismo multiconjunto por columna) y la correlación de
/// rangos entre columnas es la de `template`.
///
/// # Errors
///
/// Dimensiones inconsistentes, NaN/inf, o número de filas distinto entre
/// plantilla y corregida.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::multivariate::schaake_shuffle;
///
/// // Plantilla obs: las dos variables suben juntas (rango idéntico).
/// let template = [1.0, 10.0, 2.0, 20.0, 3.0, 30.0];
/// // Corregida: marginales nuevas, orden arbitrario.
/// let corrected = [0.7, 200.0, 0.1, 100.0, 0.4, 300.0];
/// let out = schaake_shuffle(&template, &corrected, 2).unwrap();
/// // Fila a fila, ambas variables quedan co-rankeadas como la plantilla.
/// assert_eq!(out, vec![0.1, 100.0, 0.4, 200.0, 0.7, 300.0]);
/// ```
pub fn schaake_shuffle(template: &[f64], corrected: &[f64], n_vars: usize) -> Result<Vec<f64>> {
    // Reusa la validación de matriz con una pseudo-serie de obs del largo
    // correcto en filas.
    let rows_t = validate(template, n_vars, "template")?;
    let rows_c = validate(corrected, n_vars, "corrected")?;
    if rows_t != rows_c {
        return Err(DownscaleError::LengthMismatch {
            left_name: "template (filas)",
            left: rows_t,
            right_name: "corrected (filas)",
            right: rows_c,
        });
    }

    let n = rows_t;
    let mut out = vec![0.0; corrected.len()];
    let mut col_t: Vec<(f64, usize)> = Vec::with_capacity(n);
    let mut col_c: Vec<f64> = Vec::with_capacity(n);

    for j in 0..n_vars {
        // Rango de cada día en la plantilla (argsort estable).
        col_t.clear();
        col_t.extend((0..n).map(|i| (template[i * n_vars + j], i)));
        col_t.sort_by(|a, b| a.0.partial_cmp(&b.0).expect("validado sin NaN"));

        // Valores corregidos ordenados.
        col_c.clear();
        col_c.extend((0..n).map(|i| corrected[i * n_vars + j]));
        col_c.sort_unstable_by(|a, b| a.partial_cmp(b).expect("validado sin NaN"));

        // El día con rango r en la plantilla recibe el r-ésimo valor corregido.
        for (rank, &(_, day)) in col_t.iter().enumerate() {
            out[day * n_vars + j] = col_c[rank];
        }
    }
    Ok(out)
}

/// Valida matriz aplanada (delegando en la validación de análogos) y
/// devuelve el número de filas.
fn validate(matrix: &[f64], n_vars: usize, _name: &'static str) -> Result<usize> {
    // check_matrix exige obs del mismo largo en filas; construimos una
    // pseudo-obs solo para la validación dimensional.
    if n_vars == 0 || !matrix.len().is_multiple_of(n_vars.max(1)) {
        return Err(DownscaleError::InvalidParameter {
            name: "n_vars",
            value: n_vars as f64,
            expected: "n_vars >= 1 y largo múltiplo de n_vars",
        });
    }
    let rows = matrix.len() / n_vars;
    let pseudo_obs = vec![0.0; rows];
    check_matrix(matrix, n_vars, &pseudo_obs)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ranks(col: &[f64]) -> Vec<usize> {
        let mut idx: Vec<usize> = (0..col.len()).collect();
        idx.sort_by(|&a, &b| col[a].partial_cmp(&col[b]).unwrap());
        let mut r = vec![0; col.len()];
        for (rank, &i) in idx.iter().enumerate() {
            r[i] = rank;
        }
        r
    }

    fn column(matrix: &[f64], n_vars: usize, j: usize) -> Vec<f64> {
        matrix.chunks_exact(n_vars).map(|row| row[j]).collect()
    }

    #[test]
    fn output_ranks_match_template_and_marginals_match_corrected() {
        // LCG determinista.
        let mut state = 99u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let n = 200;
        let n_vars = 3;
        let template: Vec<f64> = (0..n * n_vars).map(|_| next()).collect();
        let corrected: Vec<f64> = (0..n * n_vars).map(|_| next() * 50.0).collect();

        let out = schaake_shuffle(&template, &corrected, n_vars).unwrap();

        for j in 0..n_vars {
            // Rangos de la salida == rangos de la plantilla.
            assert_eq!(
                ranks(&column(&out, n_vars, j)),
                ranks(&column(&template, n_vars, j)),
                "columna {j}"
            );
            // Marginal preservada: mismo multiconjunto que corrected.
            let mut a = column(&out, n_vars, j);
            let mut b = column(&corrected, n_vars, j);
            a.sort_unstable_by(|x, y| x.partial_cmp(y).unwrap());
            b.sort_unstable_by(|x, y| x.partial_cmp(y).unwrap());
            assert_eq!(a, b, "columna {j}");
        }
    }

    #[test]
    fn rejects_mismatched_rows() {
        let t = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]; // 3 filas × 2 vars
        let c = [1.0, 2.0, 3.0, 4.0]; // 2 filas × 2 vars
        assert!(matches!(
            schaake_shuffle(&t, &c, 2).unwrap_err(),
            DownscaleError::LengthMismatch { .. }
        ));
    }
}
