//! Delta change (perturbación de observaciones): aplica a la serie
//! observada el cambio que el modelo proyecta entre dos períodos, en vez
//! de corregir el modelo. Complementario al quantile mapping cuando solo
//! interesa la señal de cambio (práctica estándar en estudios de impacto).

use crate::error::{DownscaleError, Result, check_series};
use crate::qm::Kind;

/// Mínimo de puntos por período del modelo.
const MIN_FIT_LEN: usize = 2;

/// Factor de cambio calibrado entre período histórico y futuro del modelo.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::delta::DeltaChange;
/// use downscale_core::qm::Kind;
///
/// let hist = [10.0, 12.0, 14.0];
/// let fut = [12.0, 14.0, 16.0]; // el modelo proyecta +2
/// let dc = DeltaChange::fit(&hist, &fut, Kind::Additive).unwrap();
/// assert_eq!(dc.delta(), 2.0);
/// assert_eq!(dc.apply(&[20.0, 21.0]).unwrap(), vec![22.0, 23.0]);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct DeltaChange {
    kind: Kind,
    delta: f64,
}

impl DeltaChange {
    /// Calibra el delta con las medias del modelo en ambos períodos:
    /// aditivo `mean(fut) − mean(hist)`, multiplicativo `mean(fut) / mean(hist)`.
    ///
    /// # Errors
    ///
    /// Series cortas/no finitas, o `mean(hist) == 0` en multiplicativo.
    pub fn fit(model_hist: &[f64], model_fut: &[f64], kind: Kind) -> Result<Self> {
        check_series("model_hist", model_hist, MIN_FIT_LEN)?;
        check_series("model_fut", model_fut, MIN_FIT_LEN)?;
        let mean = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
        let (mh, mf) = (mean(model_hist), mean(model_fut));
        let delta = match kind {
            Kind::Additive => mf - mh,
            Kind::Multiplicative => {
                if mh == 0.0 {
                    return Err(DownscaleError::InvalidParameter {
                        name: "model_hist",
                        value: mh,
                        expected: "media != 0 para delta multiplicativo",
                    });
                }
                mf / mh
            }
        };
        Ok(Self { kind, delta })
    }

    /// Perturba la serie observada con el delta calibrado.
    ///
    /// # Errors
    ///
    /// [`DownscaleError::NonFinite`] si la serie contiene NaN/inf.
    pub fn apply(&self, obs: &[f64]) -> Result<Vec<f64>> {
        check_series("obs", obs, 0)?;
        Ok(obs
            .iter()
            .map(|&v| match self.kind {
                Kind::Additive => v + self.delta,
                Kind::Multiplicative => v * self.delta,
            })
            .collect())
    }

    /// Delta calibrado (diferencia o razón de medias según [`Kind`]).
    #[must_use]
    pub fn delta(&self) -> f64 {
        self.delta
    }

    /// Tipo de perturbación.
    #[must_use]
    pub fn kind(&self) -> Kind {
        self.kind
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiplicative_delta_scales_obs() {
        let hist = [2.0, 4.0];
        let fut = [3.0, 6.0]; // razón 1.5
        let dc = DeltaChange::fit(&hist, &fut, Kind::Multiplicative).unwrap();
        assert_eq!(dc.delta(), 1.5);
        assert_eq!(dc.apply(&[10.0, 0.0]).unwrap(), vec![15.0, 0.0]);
    }

    #[test]
    fn multiplicative_rejects_zero_mean_hist() {
        let err = DeltaChange::fit(&[0.0, 0.0], &[1.0, 2.0], Kind::Multiplicative).unwrap_err();
        assert!(matches!(err, DownscaleError::InvalidParameter { .. }));
    }
}
