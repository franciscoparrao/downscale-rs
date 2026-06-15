//! Property tests sobre invariantes de la API pública. Complementan los
//! tests por ejemplo de cada módulo: en vez de casos fijos, verifican que
//! ciertas propiedades se cumplan para entradas arbitrarias.

use downscale_core::forcing::{civil_from_epoch_day, format_date, parse_date};
use downscale_core::metrics::{ks_statistic, mean_bias};
use downscale_core::multivariate::schaake_shuffle;
use downscale_core::qm::{Kind, QuantileMapping};
use proptest::prelude::*;

/// Serie finita de largo razonable, valores en un rango amplio para que los
/// empates exactos sean improbables.
fn series(min_len: usize) -> impl Strategy<Value = Vec<f64>> {
    prop::collection::vec(-1000.0f64..1000.0, min_len..300)
}

proptest! {
    /// El quantile mapping es monótono: si x1 <= x2 entonces la corrección
    /// preserva el orden. Es la invariante fundamental de un mapeo de
    /// cuantiles (composición de CDFs no decrecientes).
    #[test]
    fn qm_is_monotone(
        obs in series(20),
        model in series(20),
        a in -1500.0f64..1500.0,
        b in -1500.0f64..1500.0,
    ) {
        let qm = QuantileMapping::fit(&obs, &model, 50, Kind::Additive).unwrap();
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let clo = qm.correct_one(lo);
        let chi = qm.correct_one(hi);
        prop_assert!(clo <= chi + 1e-6, "correct({lo})={clo} > correct({hi})={chi}");
    }

    /// El estadístico KS de dos muestras siempre cae en [0, 1].
    #[test]
    fn ks_in_unit_interval(a in series(1), b in series(1)) {
        let d = ks_statistic(&a, &b).unwrap();
        prop_assert!((0.0..=1.0).contains(&d), "KS = {d}");
    }

    /// KS de una serie consigo misma es exactamente 0.
    #[test]
    fn ks_self_is_zero(a in series(1)) {
        prop_assert_eq!(ks_statistic(&a, &a).unwrap(), 0.0);
    }

    /// El sesgo medio es antisimétrico: bias(a, b) == -bias(b, a).
    #[test]
    fn mean_bias_is_antisymmetric(a in series(1), b in series(1)) {
        let ab = mean_bias(&a, &b).unwrap();
        let ba = mean_bias(&b, &a).unwrap();
        prop_assert!((ab + ba).abs() < 1e-9, "bias(a,b)={ab}, bias(b,a)={ba}");
    }

    /// El Schaake shuffle preserva exactamente las marginales corregidas:
    /// cada columna de la salida es una permutación de la columna de entrada.
    #[test]
    fn schaake_preserves_marginals(
        n in 5usize..150,
        template_seed in prop::collection::vec(-100.0f64..100.0, 5 * 2..150 * 2),
        corrected_seed in prop::collection::vec(-100.0f64..100.0, 5 * 2..150 * 2),
    ) {
        let n_vars = 2;
        // Recorta ambos vectores a n filas × n_vars.
        let len = n * n_vars;
        prop_assume!(template_seed.len() >= len && corrected_seed.len() >= len);
        let template = &template_seed[..len];
        let corrected = &corrected_seed[..len];

        let out = schaake_shuffle(template, corrected, n_vars).unwrap();
        for j in 0..n_vars {
            let mut got: Vec<f64> = (0..n).map(|i| out[i * n_vars + j]).collect();
            let mut want: Vec<f64> = (0..n).map(|i| corrected[i * n_vars + j]).collect();
            got.sort_by(|x, y| x.partial_cmp(y).unwrap());
            want.sort_by(|x, y| x.partial_cmp(y).unwrap());
            prop_assert_eq!(got, want, "columna {}", j);
        }
    }

    /// parse_date ∘ format_date es la identidad sobre días de la época
    /// dentro de un rango amplio (1801–2199, cubre todo uso climático).
    #[test]
    fn date_roundtrip(epoch_day in -61_000i64..84_000) {
        let s = format_date(epoch_day);
        let back = parse_date(&s).unwrap();
        prop_assert_eq!(back, epoch_day);
        // Y la fecha civil formateada vuelve a la misma civil.
        let (y, m, d) = civil_from_epoch_day(epoch_day);
        prop_assert_eq!(s, format!("{y:04}-{m:02}-{d:02}"));
    }
}
