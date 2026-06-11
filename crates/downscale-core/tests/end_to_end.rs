//! Integración end-to-end: serie sintética tipo precipitación con sesgo
//! multiplicativo y temperatura con sesgo aditivo, flujo completo
//! calibración → corrección → validación.

use downscale_core::metrics::{ks_statistic, mean_bias};
use downscale_core::qm::{Kind, QuantileMapping};
use downscale_core::validation::validate_split;

/// LCG determinista en [0, 1) para reproducibilidad sin dependencias.
fn uniform(seed: u64, n: usize) -> Vec<f64> {
    let mut state = seed;
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

/// Precipitación sintética diaria: ~70% días secos, lluvia exponencial.
fn synthetic_precip(seed: u64, n: usize) -> Vec<f64> {
    uniform(seed, n)
        .iter()
        .map(|&u| {
            if u < 0.7 {
                0.0
            } else {
                // Reescala el 30% húmedo a (0,1] y aplica inversa exponencial.
                let v = (u - 0.7) / 0.3;
                -10.0 * (1.0 - v.clamp(0.0, 1.0 - 1e-9)).ln()
            }
        })
        .collect()
}

#[test]
fn temperature_additive_pipeline() {
    // Temperatura: ciclo anual + ruido determinista, modelo sesgado +2.8 °C.
    let noise = uniform(11, 3650);
    let obs: Vec<f64> = (0..3650)
        .map(|i| {
            14.0 + 9.0 * (f64::from(i) * std::f64::consts::TAU / 365.0).sin()
                + (noise[i as usize] - 0.5) * 2.0
        })
        .collect();
    let model: Vec<f64> = obs.iter().map(|v| v * 1.05 + 2.8).collect();

    let report = validate_split(&obs, &model, 0.7, 100, Kind::Additive).unwrap();

    assert!(report.rmse < report.rmse_raw);
    assert!(report.mean_bias.abs() < 0.3);
    assert!(report.ks < report.ks_raw);
    // Colas: sesgo del P95 corregido acotado.
    let p95 = report
        .quantile_bias
        .iter()
        .find(|q| (q.prob - 0.95).abs() < 1e-12)
        .unwrap();
    assert!(p95.bias.abs() < 1.0, "sesgo P95 = {}", p95.bias);
}

#[test]
fn precipitation_multiplicative_pipeline() {
    let obs = synthetic_precip(101, 5000);
    // Modelo: llueve 40% más en intensidad (sesgo multiplicativo típico).
    let model: Vec<f64> = obs.iter().map(|v| v * 1.4).collect();

    let qm = QuantileMapping::fit(&obs, &model, 200, Kind::Multiplicative).unwrap();
    let corrected = qm.apply(&model).unwrap();

    // Los días secos siguen secos.
    for (c, m) in corrected.iter().zip(&model) {
        if *m == 0.0 {
            assert_eq!(*c, 0.0);
        }
        assert!(*c >= 0.0, "precipitación corregida negativa: {c}");
    }
    // Distribución corregida ≈ observada.
    let ks = ks_statistic(&corrected, &obs).unwrap();
    assert!(ks < 0.02, "KS post-corrección = {ks}");
    let bias = mean_bias(&corrected, &obs).unwrap();
    assert!(bias.abs() < 0.05, "sesgo medio post-corrección = {bias}");
}
