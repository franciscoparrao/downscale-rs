//! Errores del crate.

/// Errores de `downscale-core`.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum DownscaleError {
    /// Una serie de entrada está vacía o es demasiado corta para la operación.
    #[error("serie '{name}' demasiado corta: largo {len}, mínimo {min}")]
    SeriesTooShort {
        /// Nombre de la serie (p. ej. "obs", "model").
        name: &'static str,
        /// Largo recibido.
        len: usize,
        /// Largo mínimo requerido.
        min: usize,
    },

    /// Dos series que deben tener el mismo largo difieren.
    #[error("largos incompatibles: {left_name} ({left}) vs {right_name} ({right})")]
    LengthMismatch {
        /// Nombre de la serie izquierda.
        left_name: &'static str,
        /// Largo de la serie izquierda.
        left: usize,
        /// Nombre de la serie derecha.
        right_name: &'static str,
        /// Largo de la serie derecha.
        right: usize,
    },

    /// Una serie contiene valores no finitos (NaN o ±inf).
    #[error("serie '{name}' contiene un valor no finito en el índice {index}")]
    NonFinite {
        /// Nombre de la serie.
        name: &'static str,
        /// Índice del primer valor no finito.
        index: usize,
    },

    /// Parámetro fuera de rango.
    #[error("parámetro '{name}' fuera de rango: {value} (esperado {expected})")]
    InvalidParameter {
        /// Nombre del parámetro.
        name: &'static str,
        /// Valor recibido (formateado).
        value: f64,
        /// Descripción del rango esperado.
        expected: &'static str,
    },

    /// Fecha no parseable como `YYYY-MM-DD`.
    #[error("fecha inválida: '{value}' (formato esperado YYYY-MM-DD)")]
    InvalidDate {
        /// Texto recibido.
        value: String,
    },

    /// Eje temporal con huecos o desordenado.
    #[error("eje temporal no contiguo en '{name}': salto de {gap_days} días después de {date}")]
    NonContiguous {
        /// Nombre de la serie.
        name: &'static str,
        /// Última fecha válida antes del salto.
        date: String,
        /// Tamaño del salto en días (1 = contiguo).
        gap_days: i64,
    },

    /// Las series de forzantes no comparten ningún período común.
    #[error("las series de forzantes no comparten período común")]
    NoOverlap,
}

/// Alias de `Result` con [`DownscaleError`].
pub type Result<T> = std::result::Result<T, DownscaleError>;

/// Valida que la serie no tenga NaN/inf y cumpla un largo mínimo.
pub(crate) fn check_series(name: &'static str, series: &[f64], min_len: usize) -> Result<()> {
    if series.len() < min_len {
        return Err(DownscaleError::SeriesTooShort {
            name,
            len: series.len(),
            min: min_len,
        });
    }
    if let Some(index) = series.iter().position(|v| !v.is_finite()) {
        return Err(DownscaleError::NonFinite { name, index });
    }
    Ok(())
}

/// Valida que dos series tengan el mismo largo.
pub(crate) fn check_same_len(
    left_name: &'static str,
    left: &[f64],
    right_name: &'static str,
    right: &[f64],
) -> Result<()> {
    if left.len() != right.len() {
        return Err(DownscaleError::LengthMismatch {
            left_name,
            left: left.len(),
            right_name,
            right: right.len(),
        });
    }
    Ok(())
}
