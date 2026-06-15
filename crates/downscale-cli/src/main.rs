//! `downscale` — bias correction de series climáticas desde CSV.

mod series;

use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand, ValueEnum};
use downscale_core::analog::AnalogDownscaling;
use downscale_core::forcing::{ForcingSeries, ForcingSet, Variable};
use downscale_core::qdm::QuantileDeltaMapping;
use downscale_core::qm::{Kind, NodePlacement, QuantileMapping};
use downscale_core::regression::LinearDownscaling;
use downscale_core::validation::{QmOptions, validate_split_with};
use downscale_core::wetday::WetDayCorrection;

use series::{Matrix, Series, pair_by_date, pair_matrix_series};

#[derive(Parser)]
#[command(
    name = "downscale",
    version,
    about = "Bias correction y downscaling estadístico de clima"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Clone, Copy, ValueEnum)]
enum KindArg {
    /// Corrección aditiva (temperatura).
    Add,
    /// Corrección multiplicativa (precipitación).
    Mult,
}

impl From<KindArg> for Kind {
    fn from(k: KindArg) -> Self {
        match k {
            KindArg::Add => Kind::Additive,
            KindArg::Mult => Kind::Multiplicative,
        }
    }
}

#[derive(Clone, Copy, ValueEnum)]
enum MethodArg {
    /// Quantile mapping empírico (corrige hacia el clima observado).
    Eqm,
    /// Quantile delta mapping (Cannon 2015: corrige preservando la señal
    /// de cambio de la serie objetivo cuantil a cuantil).
    Qdm,
}

#[derive(Clone, Copy, ValueEnum)]
enum NodesArg {
    /// Nodos i/(n−1) con extremos (default).
    Endpoints,
    /// Nodos (i+0.5)/n, convención xclim/xsdba.
    Midpoint,
}

impl From<NodesArg> for NodePlacement {
    fn from(n: NodesArg) -> Self {
        match n {
            NodesArg::Endpoints => NodePlacement::Endpoints,
            NodesArg::Midpoint => NodePlacement::Midpoint,
        }
    }
}

#[derive(Subcommand)]
enum Command {
    /// Valida quantile mapping con split temporal sobre series pareadas por fecha.
    Validate {
        /// CSV observado (fecha,valor).
        #[arg(long)]
        obs: PathBuf,
        /// CSV del modelo/reanálisis (fecha,valor).
        #[arg(long)]
        model: PathBuf,
        /// Tipo de corrección en colas.
        #[arg(long, value_enum, default_value = "add")]
        kind: KindArg,
        /// Fracción inicial para calibración.
        #[arg(long, default_value_t = 0.7)]
        calib_frac: f64,
        /// Número de cuantiles de la CDF empírica.
        #[arg(long, default_value_t = 100)]
        quantiles: usize,
        /// Colocación de nodos de probabilidad.
        #[arg(long, value_enum, default_value = "endpoints")]
        nodes: NodesArg,
        /// Umbral seco/húmedo observado (mm) para adaptación de umbral
        /// (corrección de frecuencia de días húmedos) antes del EQM.
        #[arg(long)]
        wet_day_threshold: Option<f64>,
    },
    /// Calibra con obs+model y corrige una serie objetivo; escribe CSV corregido.
    Correct {
        /// CSV observado del período de calibración (fecha,valor).
        #[arg(long)]
        obs: PathBuf,
        /// CSV del modelo en el período de calibración (fecha,valor).
        #[arg(long)]
        model: PathBuf,
        /// CSV a corregir (fecha,valor). Si se omite, se corrige `--model` completo.
        #[arg(long)]
        target: Option<PathBuf>,
        /// Método de corrección.
        #[arg(long, value_enum, default_value = "eqm")]
        method: MethodArg,
        /// Tipo de corrección en colas.
        #[arg(long, value_enum, default_value = "add")]
        kind: KindArg,
        /// Número de cuantiles de la CDF empírica.
        #[arg(long, default_value_t = 100)]
        quantiles: usize,
        /// Colocación de nodos de probabilidad.
        #[arg(long, value_enum, default_value = "endpoints")]
        nodes: NodesArg,
        /// Umbral seco/húmedo observado (mm) para adaptación de umbral
        /// antes del EQM.
        #[arg(long)]
        wet_day_threshold: Option<f64>,
        /// Ruta del CSV de salida.
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Ensambla forzantes corregidas en el CSV que consume rainflow
    /// (`date,pr,pet[,tmean]`), validando eje diario contiguo y alineando
    /// al período común.
    Forcing {
        /// CSV de precipitación corregida (fecha,valor; mm/día).
        #[arg(long)]
        pr: PathBuf,
        /// CSV de evapotranspiración potencial (fecha,valor; mm/día).
        #[arg(long)]
        pet: PathBuf,
        /// CSV de temperatura media (fecha,valor; °C), opcional.
        #[arg(long)]
        temp: Option<PathBuf>,
        /// Ruta del CSV de forzantes de salida.
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Downscaling por análogos: k-NN sobre predictores de gran escala
    /// (CSV `date,col1,col2,...`) calibrado contra observaciones locales.
    Analog {
        /// CSV de predictores de calibración (`date,col1,...`).
        #[arg(long)]
        predictors: PathBuf,
        /// CSV de observaciones locales (fecha,valor).
        #[arg(long)]
        obs: PathBuf,
        /// CSV de predictores a predecir. Si se omite, se usa `--predictors`.
        #[arg(long)]
        target: Option<PathBuf>,
        /// Número de análogos.
        #[arg(long, default_value_t = 10)]
        k: usize,
        /// Ruta del CSV de salida (date,value).
        #[arg(short, long)]
        output: PathBuf,
    },
    /// Downscaling por regresión lineal múltiple (OLS) sobre predictores
    /// de gran escala calibrada contra observaciones locales.
    Regress {
        /// CSV de predictores de calibración (`date,col1,...`).
        #[arg(long)]
        predictors: PathBuf,
        /// CSV de observaciones locales (fecha,valor).
        #[arg(long)]
        obs: PathBuf,
        /// CSV de predictores a predecir. Si se omite, se usa `--predictors`.
        #[arg(long)]
        target: Option<PathBuf>,
        /// Trunca las predicciones a >= 0 (p. ej. precipitación).
        #[arg(long)]
        non_negative: bool,
        /// Ruta del CSV de salida (date,value).
        #[arg(short, long)]
        output: PathBuf,
    },
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Validate {
            obs,
            model,
            kind,
            calib_frac,
            quantiles,
            nodes,
            wet_day_threshold,
        } => {
            let obs_s = Series::read_csv(&obs)?;
            let model_s = Series::read_csv(&model)?;
            let (dates, o, m) = pair_by_date(&obs_s, &model_s)?;
            eprintln!(
                "series pareadas: {} días ({} .. {})",
                dates.len(),
                dates.first().unwrap(),
                dates.last().unwrap()
            );
            let opts = QmOptions {
                n_quantiles: quantiles,
                kind: kind.into(),
                placement: nodes.into(),
                wet_day_threshold,
            };
            let report = validate_split_with(&o, &m, calib_frac, &opts)?;
            print_report(&report, &dates);
            Ok(())
        }
        Command::Correct {
            obs,
            model,
            target,
            method,
            kind,
            quantiles,
            nodes,
            wet_day_threshold,
            output,
        } => {
            let obs_s = Series::read_csv(&obs)?;
            let model_s = Series::read_csv(&model)?;
            let (_, o, m) = pair_by_date(&obs_s, &model_s)?;

            // Adaptación opcional de umbral seco/húmedo, calibrada en el
            // período común y aplicada también a la serie objetivo.
            let wd = wet_day_threshold
                .map(|thr| WetDayCorrection::fit(&o, &m, thr))
                .transpose()?;
            let m = wd.as_ref().map_or(m.clone(), |w| w.transform(&m));

            let target_s = match &target {
                Some(p) => Series::read_csv(p)?,
                None => model_s,
            };
            let target_values = wd
                .as_ref()
                .map_or(target_s.values.clone(), |w| w.transform(&target_s.values));

            let corrected = match method {
                MethodArg::Eqm => {
                    QuantileMapping::fit_with_nodes(&o, &m, quantiles, kind.into(), nodes.into())?
                        .apply(&target_values)?
                }
                MethodArg::Qdm => QuantileDeltaMapping::fit_with_nodes(
                    &o,
                    &m,
                    quantiles,
                    kind.into(),
                    nodes.into(),
                )?
                .apply(&target_values)?,
            };

            let mut out = String::from("date,value\n");
            for (d, v) in target_s.dates.iter().zip(&corrected) {
                out.push_str(&format!("{d},{v:.3}\n"));
            }
            std::fs::write(&output, out)
                .with_context(|| format!("no se pudo escribir {}", output.display()))?;
            eprintln!("escrito {} ({} filas)", output.display(), corrected.len());
            Ok(())
        }
        Command::Forcing {
            pr,
            pet,
            temp,
            output,
        } => {
            let load = |path: &PathBuf, var: Variable| -> Result<ForcingSeries> {
                let s = Series::read_csv(path)?;
                ForcingSeries::from_dates(var, &s.dates, &s.values)
                    .with_context(|| format!("{} ({})", path.display(), var.column_name()))
            };
            let mut series = vec![
                load(&pr, Variable::Precipitation)?,
                load(&pet, Variable::Pet)?,
            ];
            if let Some(t) = &temp {
                series.push(load(t, Variable::TemperatureMean)?);
            }
            let set = ForcingSet::align(series)?;
            std::fs::write(&output, set.to_csv())
                .with_context(|| format!("no se pudo escribir {}", output.display()))?;
            eprintln!(
                "escrito {} ({} días desde {}, columnas: {})",
                output.display(),
                set.len(),
                downscale_core::forcing::format_date(set.start_day()),
                set.variables()
                    .iter()
                    .map(|v| v.column_name())
                    .collect::<Vec<_>>()
                    .join(",")
            );
            Ok(())
        }
        Command::Analog {
            predictors,
            obs,
            target,
            k,
            output,
        } => {
            let (cal_data, n_features, obs_cal, tgt) =
                downscale_inputs(&predictors, &obs, target.as_deref())?;
            let model = AnalogDownscaling::fit(&cal_data, n_features, &obs_cal, k)?;
            let pred = model.predict(&tgt.data)?;
            write_predictions(&output, &tgt.dates, &pred)?;
            eprintln!(
                "escrito {} ({} predicciones, k={k}, {n_features} predictores)",
                output.display(),
                pred.len()
            );
            Ok(())
        }
        Command::Regress {
            predictors,
            obs,
            target,
            non_negative,
            output,
        } => {
            let (cal_data, n_features, obs_cal, tgt) =
                downscale_inputs(&predictors, &obs, target.as_deref())?;
            let model = LinearDownscaling::fit(&cal_data, n_features, &obs_cal)?;
            let mut pred = model.predict(&tgt.data)?;
            if non_negative {
                for v in &mut pred {
                    *v = v.max(0.0);
                }
            }
            write_predictions(&output, &tgt.dates, &pred)?;
            eprintln!(
                "escrito {} ({} predicciones, R²cal={:.3}, {n_features} predictores)",
                output.display(),
                pred.len(),
                model.r2()
            );
            Ok(())
        }
    }
}

/// Carga predictores + observaciones, los parea por fecha, y resuelve el
/// target (otra matriz o los propios predictores), verificando que sus
/// columnas coincidan con las de calibración.
fn downscale_inputs(
    predictors: &std::path::Path,
    obs: &std::path::Path,
    target: Option<&std::path::Path>,
) -> Result<(Vec<f64>, usize, Vec<f64>, Matrix)> {
    let pred_m = Matrix::read_csv(predictors)?;
    let obs_s = Series::read_csv(obs)?;
    let (cal_data, n_features, obs_cal) = pair_matrix_series(&pred_m, &obs_s)?;

    let tgt = match target {
        Some(p) => {
            let t = Matrix::read_csv(p)?;
            if t.columns != pred_m.columns {
                bail!(
                    "las columnas del target {:?} no coinciden con las de calibración {:?}",
                    t.columns,
                    pred_m.columns
                );
            }
            t
        }
        None => pred_m,
    };
    Ok((cal_data, n_features, obs_cal, tgt))
}

/// Escribe un CSV `date,value` con 3 decimales.
fn write_predictions(output: &std::path::Path, dates: &[String], values: &[f64]) -> Result<()> {
    let mut out = String::from("date,value\n");
    for (d, v) in dates.iter().zip(values) {
        out.push_str(&format!("{d},{v:.3}\n"));
    }
    std::fs::write(output, out).with_context(|| format!("no se pudo escribir {}", output.display()))
}

fn print_report(report: &downscale_core::ValidationReport, dates: &[String]) {
    let pct = |new: f64, raw: f64| {
        if raw == 0.0 {
            0.0
        } else {
            100.0 * (raw - new) / raw
        }
    };
    println!("== Validación split temporal ==");
    println!(
        "calibración: {} días | validación: {} días (desde {})",
        report.split_index,
        dates.len() - report.split_index,
        dates[report.split_index]
    );
    println!("                 corregido      raw       mejora");
    println!(
        "RMSE         {:>10.3} {:>10.3} {:>10.1}%",
        report.rmse,
        report.rmse_raw,
        pct(report.rmse, report.rmse_raw)
    );
    println!(
        "sesgo medio  {:>10.3} {:>10.3} {:>10.1}%",
        report.mean_bias,
        report.mean_bias_raw,
        pct(report.mean_bias.abs(), report.mean_bias_raw.abs())
    );
    println!(
        "KS           {:>10.4} {:>10.4} {:>10.1}%",
        report.ks,
        report.ks_raw,
        pct(report.ks, report.ks_raw)
    );
    println!("\nsesgo por cuantil (corregido - obs):");
    for q in &report.quantile_bias {
        println!(
            "  P{:<4} sim={:>9.3}  obs={:>9.3}  sesgo={:>8.3}",
            (q.prob * 100.0).round(),
            q.sim,
            q.obs,
            q.bias
        );
    }
}
