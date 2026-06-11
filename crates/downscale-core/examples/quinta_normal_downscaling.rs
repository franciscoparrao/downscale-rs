//! Comparación de métodos sobre Quinta Normal (CR2) vs ERA5 (Open-Meteo):
//! raw, EQM, análogos y regresión, evaluados en el holdout 1997–2018.
//!
//! Requiere los CSV de `data/parity/` (ver `scripts/parity_quinta_normal.py`
//! y `data/README.md`). Ejecutar desde la raíz del workspace:
//!
//! ```text
//! cargo run --release -p downscale-core --example quinta_normal_downscaling
//! ```

use downscale_core::analog::AnalogDownscaling;
use downscale_core::metrics::{ks_statistic, mean_bias, rmse};
use downscale_core::qm::{Kind, QuantileMapping};
use downscale_core::regression::LinearDownscaling;

/// Lee un CSV `date,value` (formato de data/parity/) → (fechas, valores).
fn read_csv(path: &str) -> (Vec<String>, Vec<f64>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| {
        panic!("no se pudo leer {path}: {e} (¿corriste el script de paridad?)")
    });
    let mut dates = Vec::new();
    let mut values = Vec::new();
    for line in text.lines().skip(1) {
        let (d, v) = line.split_once(',').expect("formato date,value");
        dates.push(d.to_string());
        values.push(v.parse::<f64>().expect("valor numérico"));
    }
    (dates, values)
}

/// Día del año aproximado (1..=366) desde fecha ISO `YYYY-MM-DD`.
fn day_of_year(date: &str) -> f64 {
    const CUM: [u32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let month: usize = date[5..7].parse().expect("mes");
    let day: u32 = date[8..10].parse().expect("día");
    f64::from(CUM[month - 1] + day)
}

/// Predictores por día: \[pr ERA5, sin(doy), cos(doy)\], aplanados.
fn build_predictors(dates: &[String], era5: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(dates.len() * 3);
    for (d, &pr) in dates.iter().zip(era5) {
        let phase = day_of_year(d) * std::f64::consts::TAU / 365.25;
        out.extend([pr, phase.sin(), phase.cos()]);
    }
    out
}

fn main() {
    let (cal_dates, obs_cal) = read_csv("data/parity/obs_cal.csv");
    let (_, mod_cal) = read_csv("data/parity/model_cal.csv");
    let (val_dates, mod_val) = read_csv("data/parity/model_val.csv");
    let (_, obs_val) = read_csv("data/parity/obs_val.csv");

    let x_cal = build_predictors(&cal_dates, &mod_cal);
    let x_val = build_predictors(&val_dates, &mod_val);

    // EQM (referencia de bias correction).
    let qm = QuantileMapping::fit(&obs_cal, &mod_cal, 100, Kind::Multiplicative).expect("fit EQM");
    let eqm = qm.apply(&mod_val).expect("apply EQM");

    // Análogos sobre [pr, sin, cos]: k=10 (media ponderada, minimiza error
    // puntual) y k=1 (remuestreo puro, preserva la distribución).
    let ad10 = AnalogDownscaling::fit(&x_cal, 3, &obs_cal, 10).expect("fit análogos");
    let analog10 = ad10.predict(&x_val).expect("predict análogos");
    let ad1 = AnalogDownscaling::fit(&x_cal, 3, &obs_cal, 1).expect("fit análogos");
    let analog1 = ad1.predict(&x_val).expect("predict análogos");

    // Regresión lineal sobre los mismos predictores (clip a >= 0).
    let lm = LinearDownscaling::fit(&x_cal, 3, &obs_cal).expect("fit regresión");
    let regression: Vec<f64> = lm
        .predict(&x_val)
        .expect("predict regresión")
        .iter()
        .map(|&v| v.max(0.0))
        .collect();
    eprintln!(
        "regresión: R²(cal) = {:.3}, coefs = {:?}",
        lm.r2(),
        lm.coefs()
            .iter()
            .map(|c| (c * 1000.0).round() / 1000.0)
            .collect::<Vec<_>>()
    );

    println!(
        "{:<14} {:>8} {:>8} {:>8}  (validación {} días)",
        "método",
        "RMSE",
        "KS",
        "sesgo",
        obs_val.len()
    );
    for (name, series) in [
        ("raw ERA5", &mod_val),
        ("EQM", &eqm),
        ("análogos k=10", &analog10),
        ("análogos k=1", &analog1),
        ("regresión", &regression),
    ] {
        println!(
            "{:<14} {:>8.3} {:>8.4} {:>8.3}",
            name,
            rmse(series, &obs_val).expect("rmse"),
            ks_statistic(series, &obs_val).expect("ks"),
            mean_bias(series, &obs_val).expect("bias"),
        );
    }
}
