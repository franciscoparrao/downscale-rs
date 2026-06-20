//! Corrección de sesgo de campos grillados, celda por celda.
//!
//! Una grilla temporal se pasa aplanada en orden row-major `[n_time ×
//! n_cells]`: el elemento del paso `t` y la celda `c` está en
//! `t*n_cells + c`. Cada celda se corrige de forma independiente con
//! quantile mapping (reusa [`crate::qm`]). Las celdas con algún valor no
//! finito —mar, máscara, dato faltante— se devuelven como `NaN` sin tocar,
//! de modo que la máscara espacial se preserva.
//!
//! El núcleo no hace I/O: las grillas reales (NetCDF, GeoTIFF) se leen y
//! escriben en las superficies (p. ej. xarray vía los bindings Python).
//! Las celdas son independientes, así que el cómputo es trivialmente
//! paralelizable aguas arriba.

use crate::error::{DownscaleError, Result};
use crate::qm::{Kind, NodePlacement, QuantileMapping};

/// Configuración de [`correct_grid`].
#[derive(Debug, Clone, Copy)]
pub struct GridOptions {
    /// Número de cuantiles de la CDF empírica.
    pub n_quantiles: usize,
    /// Tipo de corrección en colas.
    pub kind: Kind,
    /// Colocación de nodos de probabilidad.
    pub placement: NodePlacement,
}

impl Default for GridOptions {
    fn default() -> Self {
        Self {
            n_quantiles: 100,
            kind: Kind::Multiplicative,
            placement: NodePlacement::Midpoint,
        }
    }
}

/// Corrige el sesgo de `model` hacia `obs`, celda por celda.
///
/// `obs` es `[n_obs_time × n_cells]` y `model` `[n_model_time × n_cells]`,
/// ambos aplanados row-major y con el mismo número de celdas (la misma
/// grilla espacial); el número de pasos temporales puede diferir. Devuelve
/// el campo corregido con la forma de `model`.
///
/// Una celda se enmascara (sale toda `NaN`) si su serie de observaciones o
/// de modelo contiene algún valor no finito, o es demasiado corta para
/// calibrar.
///
/// # Errors
///
/// [`DownscaleError::InvalidParameter`] si `n_cells == 0`, los largos no son
/// múltiplos de `n_cells`, o `n_quantiles < 2`.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::grid::{correct_grid, GridOptions};
/// use downscale_core::qm::Kind;
///
/// // 2 pasos × 2 celdas. La celda 0 tiene sesgo +2 (aditivo).
/// let obs =   [10.0, 5.0,  20.0, 8.0];
/// let model = [12.0, 5.0,  22.0, 8.0];
/// let opts = GridOptions { kind: Kind::Additive, ..Default::default() };
/// let out = correct_grid(&obs, &model, 2, &opts).unwrap();
/// // La celda 0 se corrige hacia obs; la celda 1 ya coincide.
/// assert!((out[0] - 10.0).abs() < 1e-9 && (out[2] - 20.0).abs() < 1e-9);
/// ```
pub fn correct_grid(
    obs: &[f64],
    model: &[f64],
    n_cells: usize,
    opts: &GridOptions,
) -> Result<Vec<f64>> {
    if n_cells == 0 {
        return Err(DownscaleError::InvalidParameter {
            name: "n_cells",
            value: 0.0,
            expected: ">= 1",
        });
    }
    if !obs.len().is_multiple_of(n_cells) || !model.len().is_multiple_of(n_cells) {
        return Err(DownscaleError::InvalidParameter {
            name: "grid",
            value: model.len() as f64,
            expected: "largos múltiplos de n_cells",
        });
    }
    if opts.n_quantiles < 2 {
        return Err(DownscaleError::InvalidParameter {
            name: "n_quantiles",
            value: opts.n_quantiles as f64,
            expected: ">= 2",
        });
    }

    let n_obs = obs.len() / n_cells;
    let n_mod = model.len() / n_cells;
    let mut out = vec![f64::NAN; model.len()];

    let mut obs_col = vec![0.0; n_obs];
    let mut mod_col = vec![0.0; n_mod];
    for c in 0..n_cells {
        // Extrae la serie temporal de la celda (stride n_cells).
        for (t, v) in obs_col.iter_mut().enumerate() {
            *v = obs[t * n_cells + c];
        }
        for (t, v) in mod_col.iter_mut().enumerate() {
            *v = model[t * n_cells + c];
        }
        // Celda enmascarada (mar/dato faltante): se deja NaN.
        if obs_col.iter().chain(&mod_col).any(|v| !v.is_finite()) {
            continue;
        }
        let qm = match QuantileMapping::fit_with_nodes(
            &obs_col,
            &mod_col,
            opts.n_quantiles,
            opts.kind,
            opts.placement,
        ) {
            Ok(qm) => qm,
            Err(_) => continue, // celda demasiado corta → enmascarada
        };
        for (t, &x) in mod_col.iter().enumerate() {
            out[t * n_cells + c] = qm.correct_one(x);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::mean_bias;

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
    fn corrects_each_cell_independently() {
        // Grilla 3 celdas con sesgos aditivos distintos por celda.
        let n_cells = 3;
        let n_time = 500;
        let biases = [1.0, 4.0, -2.0];
        let base = uniform(1, n_time * n_cells);
        let mut obs = vec![0.0; n_time * n_cells];
        let mut model = vec![0.0; n_time * n_cells];
        for t in 0..n_time {
            for c in 0..n_cells {
                let v = 10.0 + 8.0 * base[t * n_cells + c];
                obs[t * n_cells + c] = v;
                model[t * n_cells + c] = v + biases[c];
            }
        }
        let opts = GridOptions {
            kind: Kind::Additive,
            ..Default::default()
        };
        let out = correct_grid(&obs, &model, n_cells, &opts).unwrap();

        // Cada celda corregida ≈ obs de esa celda (sesgo removido).
        for c in 0..n_cells {
            let obs_c: Vec<f64> = (0..n_time).map(|t| obs[t * n_cells + c]).collect();
            let out_c: Vec<f64> = (0..n_time).map(|t| out[t * n_cells + c]).collect();
            assert!(mean_bias(&out_c, &obs_c).unwrap().abs() < 1e-6, "celda {c}");
        }
    }

    #[test]
    fn masks_nan_cells() {
        // Celda 1 es "mar" (NaN en obs); debe salir toda NaN.
        let n_cells = 2;
        let obs = [10.0, f64::NAN, 12.0, f64::NAN, 11.0, f64::NAN];
        let model = [11.0, 3.0, 13.0, 4.0, 12.0, 5.0];
        let out = correct_grid(&obs, &model, n_cells, &GridOptions::default()).unwrap();
        // Celda 0 corregida (finita), celda 1 enmascarada (NaN).
        for t in 0..3 {
            assert!(out[t * n_cells].is_finite());
            assert!(out[t * n_cells + 1].is_nan());
        }
    }

    #[test]
    fn matches_per_cell_quantile_mapping() {
        let n_cells = 2;
        let n_time = 300;
        let obs = uniform(5, n_time * n_cells);
        let model: Vec<f64> = uniform(6, n_time * n_cells)
            .iter()
            .map(|v| v * 1.3)
            .collect();
        let out = correct_grid(&obs, &model, n_cells, &GridOptions::default()).unwrap();

        for c in 0..n_cells {
            let obs_c: Vec<f64> = (0..n_time).map(|t| obs[t * n_cells + c]).collect();
            let mod_c: Vec<f64> = (0..n_time).map(|t| model[t * n_cells + c]).collect();
            let qm = QuantileMapping::fit_with_nodes(
                &obs_c,
                &mod_c,
                100,
                Kind::Multiplicative,
                NodePlacement::Midpoint,
            )
            .unwrap();
            let ref_c = qm.apply(&mod_c).unwrap();
            for t in 0..n_time {
                assert!((out[t * n_cells + c] - ref_c[t]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn rejects_bad_dimensions() {
        assert!(correct_grid(&[1.0, 2.0, 3.0], &[1.0, 2.0], 2, &GridOptions::default()).is_err());
        assert!(correct_grid(&[1.0, 2.0], &[1.0, 2.0], 0, &GridOptions::default()).is_err());
    }
}
