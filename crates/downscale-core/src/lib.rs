//! # downscale-core
//!
//! Núcleo de corrección de sesgo y downscaling estadístico de variables
//! climáticas (GCM/RCM → escala local). Sin I/O: opera sobre slices de
//! `f64`; la lectura de CSV/NetCDF/GeoTIFF vive en los crates de superficie
//! (CLI, PyO3).
//!
//! ## Métodos
//!
//! - [`qm`]: quantile mapping empírico (aditivo y multiplicativo, nodos
//!   con extremos o puntos medios).
//! - [`qdm`]: quantile delta mapping (Cannon et al. 2015) — preserva la
//!   señal de cambio del modelo cuantil a cuantil (v0.2).
//! - [`multivariate`]: Schaake shuffle (Clark et al. 2004) — restaura la
//!   dependencia entre variables tras la corrección univariada (v0.2).
//! - [`pet`]: PET de Hargreaves con radiación extraterrestre FAO-56 (v0.2).
//! - [`parametric`]: QM paramétrico — normal (temperatura) y gamma mixta
//!   con masa en cero (precipitación; corrige frecuencia de días húmedos).
//! - [`analog`]: downscaling por análogos (k-NN estandarizado).
//! - [`regression`]: downscaling por regresión lineal múltiple (OLS).
//! - [`delta`]: delta change (perturbación de observaciones).
//! - [`wetday`]: adaptación de umbral seco/húmedo para EQM.
//! - [`validation`]: split temporal + reporte con líneas base.
//! - [`metrics`]: RMSE, sesgo medio, KS de dos muestras, sesgo por cuantil.
//! - [`forcing`]: interfaz de forzantes hacia motores hidrológicos
//!   (rainflow) — eje diario contiguo validado + CSV canónico.
//!
//! ## Ejemplo end-to-end
//!
//! ```
//! use downscale_core::qm::{Kind, QuantileMapping};
//! use downscale_core::validation::validate_split;
//!
//! // Series pareadas obs/modelo del período histórico.
//! let obs: Vec<f64> = (0..365).map(|i| 12.0 + (f64::from(i) * 0.017).sin() * 6.0).collect();
//! let model: Vec<f64> = obs.iter().map(|v| v + 2.5).collect();
//!
//! // 1. Validar el método con split temporal 70/30.
//! let report = validate_split(&obs, &model, 0.7, 100, Kind::Additive).unwrap();
//! assert!(report.rmse < report.rmse_raw);
//!
//! // 2. Calibrar con todo el período histórico y corregir un escenario.
//! let qm = QuantileMapping::fit(&obs, &model, 100, Kind::Additive).unwrap();
//! let futuro: Vec<f64> = model.iter().map(|v| v + 1.0).collect();
//! let corregido = qm.apply(&futuro).unwrap();
//! assert_eq!(corregido.len(), futuro.len());
//! ```

#![warn(missing_docs)]

pub mod analog;
pub mod delta;
pub mod error;
pub mod forcing;
pub mod metrics;
pub mod multivariate;
pub mod parametric;
pub mod pet;
pub mod qdm;
pub mod qm;
pub mod regression;
mod special;
pub mod validation;
pub mod wetday;

pub use analog::AnalogDownscaling;
pub use delta::DeltaChange;
pub use error::{DownscaleError, Result};
pub use forcing::{ForcingSeries, ForcingSet, Variable};
pub use multivariate::schaake_shuffle;
pub use parametric::{Distribution, ParametricQuantileMapping};
pub use pet::hargreaves;
pub use qdm::QuantileDeltaMapping;
pub use qm::{Kind, NodePlacement, QuantileMapping};
pub use regression::LinearDownscaling;
pub use validation::{QmOptions, ValidationReport, validate_split, validate_split_with};
pub use wetday::WetDayCorrection;
