//! Benchmarks de los métodos de corrección/downscaling sobre un tamaño
//! representativo del caso Quinta Normal (~17 430 días de calibración,
//! 7 470 de validación). Mide calibración (`fit`) y aplicación (`apply`)
//! por separado donde tiene sentido.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

use downscale_core::analog::AnalogDownscaling;
use downscale_core::parametric::{Distribution, ParametricQuantileMapping};
use downscale_core::qdm::QuantileDeltaMapping;
use downscale_core::qm::{Kind, QuantileMapping};
use downscale_core::regression::LinearDownscaling;

const N_CAL: usize = 17_430;
const N_VAL: usize = 7_470;
const N_QUANTILES: usize = 100;

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

/// Precipitación sintética: 70 % de días secos, lluvia exponencial.
fn precip(seed: u64, n: usize) -> Vec<f64> {
    uniform(seed, n)
        .iter()
        .map(|&u| {
            if u < 0.7 {
                0.0
            } else {
                -8.0 * (1.0 - (u - 0.7) / 0.3).ln()
            }
        })
        .collect()
}

fn bench_eqm(c: &mut Criterion) {
    let obs = precip(1, N_CAL);
    let model = precip(2, N_CAL);
    let proj = precip(3, N_VAL);
    let mut g = c.benchmark_group("eqm");
    g.bench_function("fit", |b| {
        b.iter(|| {
            QuantileMapping::fit(
                black_box(&obs),
                black_box(&model),
                N_QUANTILES,
                Kind::Multiplicative,
            )
            .unwrap()
        });
    });
    let qm = QuantileMapping::fit(&obs, &model, N_QUANTILES, Kind::Multiplicative).unwrap();
    g.bench_function("apply", |b| {
        b.iter(|| qm.apply(black_box(&proj)).unwrap());
    });
    g.finish();
}

fn bench_qdm(c: &mut Criterion) {
    let obs = precip(1, N_CAL);
    let model = precip(2, N_CAL);
    let proj = precip(3, N_VAL);
    let qdm = QuantileDeltaMapping::fit(&obs, &model, N_QUANTILES, Kind::Multiplicative).unwrap();
    c.bench_function("qdm_apply", |b| {
        b.iter(|| qdm.apply(black_box(&proj)).unwrap());
    });
}

fn bench_parametric(c: &mut Criterion) {
    let obs = precip(1, N_CAL);
    let model = precip(2, N_CAL);
    let proj = precip(3, N_VAL);
    let dist = Distribution::Gamma { wet_threshold: 0.1 };
    let mut g = c.benchmark_group("parametric_gamma");
    g.bench_function("fit", |b| {
        b.iter(|| {
            ParametricQuantileMapping::fit(black_box(&obs), black_box(&model), dist).unwrap()
        });
    });
    let pqm = ParametricQuantileMapping::fit(&obs, &model, dist).unwrap();
    g.bench_function("apply", |b| {
        b.iter(|| pqm.apply(black_box(&proj)).unwrap());
    });
    g.finish();
}

fn bench_analog(c: &mut Criterion) {
    // 4 predictores sinópticos sintéticos.
    let n_features = 4;
    let pred_cal: Vec<f64> = uniform(10, N_CAL * n_features);
    let obs = precip(1, N_CAL);
    let pred_val: Vec<f64> = uniform(20, N_VAL * n_features);
    let ad = AnalogDownscaling::fit(&pred_cal, n_features, &obs, 10).unwrap();
    c.bench_function("analog_k10_predict", |b| {
        b.iter(|| ad.predict(black_box(&pred_val)).unwrap());
    });
}

fn bench_regression(c: &mut Criterion) {
    let n_features = 4;
    let pred_cal: Vec<f64> = uniform(10, N_CAL * n_features);
    let obs = precip(1, N_CAL);
    c.bench_function("regression_fit", |b| {
        b.iter(|| {
            LinearDownscaling::fit(black_box(&pred_cal), n_features, black_box(&obs)).unwrap()
        });
    });
}

criterion_group!(
    benches,
    bench_eqm,
    bench_qdm,
    bench_parametric,
    bench_analog,
    bench_regression
);
criterion_main!(benches);
