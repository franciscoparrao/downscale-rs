//! Funciones especiales mínimas para QM paramétrico: ln-gamma (Lanczos),
//! digamma/trigamma, gamma incompleta regularizada y su inversa, y CDF/PPF
//! de la normal estándar.
//!
//! Implementaciones clásicas (Lanczos g=7; serie/fracción continua de
//! Numerical Recipes; PPF de Acklam) verificadas contra SciPy en los tests.

// Coeficientes publicados transcritos textuales, aunque excedan f64.
#![allow(clippy::excessive_precision)]

/// Precisión objetivo de las iteraciones internas.
const EPS: f64 = 3e-14;

/// ln Γ(x) por aproximación de Lanczos (g = 7, 9 coeficientes).
pub(crate) fn ln_gamma(x: f64) -> f64 {
    const COEF: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_13,
        -176.615_029_162_140_59,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_571_6e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflexión: Γ(x)Γ(1-x) = π / sin(πx)
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = COEF[0];
    let t = x + 7.5;
    for (i, &c) in COEF.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Digamma ψ(x) = d/dx ln Γ(x): recurrencia + serie asintótica.
pub(crate) fn digamma(mut x: f64) -> f64 {
    let mut r = 0.0;
    while x < 6.0 {
        r -= 1.0 / x;
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    r + x.ln()
        - 0.5 * inv
        - inv2 * (1.0 / 12.0 - inv2 * (1.0 / 120.0 - inv2 * (1.0 / 252.0 - inv2 / 240.0)))
}

/// Trigamma ψ′(x): recurrencia + serie asintótica.
pub(crate) fn trigamma(mut x: f64) -> f64 {
    let mut r = 0.0;
    while x < 6.0 {
        r += 1.0 / (x * x);
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    r + inv
        * (1.0
            + inv
                * (0.5
                    + inv * (1.0 / 6.0 - inv2 * (1.0 / 30.0 - inv2 * (1.0 / 42.0 - inv2 / 30.0)))))
}

/// Gamma incompleta inferior regularizada P(a, x) = γ(a, x) / Γ(a).
pub(crate) fn gamma_p(a: f64, x: f64) -> f64 {
    debug_assert!(a > 0.0);
    if x <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        gamma_p_series(a, x)
    } else {
        1.0 - gamma_q_cf(a, x)
    }
}

/// Serie de P(a, x) (converge rápido para x < a + 1).
fn gamma_p_series(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..500 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * EPS {
            break;
        }
    }
    sum * (-x + a * x.ln() - gln).exp()
}

/// Fracción continua de Q(a, x) = 1 − P(a, x) (método de Lentz).
fn gamma_q_cf(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let fpmin = f64::MIN_POSITIVE / EPS;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / fpmin;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..500 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = b + an / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    h * (-x + a * x.ln() - gln).exp()
}

/// Inversa de P(a, ·): devuelve x tal que `gamma_p(a, x) = p`.
///
/// Guess inicial de Wilson–Hilferty refinado con Halley
/// (Numerical Recipes §6.2.1). `p` se satura a \[0, 1).
pub(crate) fn gamma_p_inv(a: f64, p: f64) -> f64 {
    debug_assert!(a > 0.0);
    if p <= 0.0 {
        return 0.0;
    }
    let p = p.min(1.0 - 1e-15);
    let gln = ln_gamma(a);
    let a1 = a - 1.0;

    let mut x = if a > 1.0 {
        let z = norm_ppf(p);
        let t = 1.0 - 1.0 / (9.0 * a) + z / (3.0 * a.sqrt());
        (a * t * t * t).max(1e-3)
    } else {
        let t = 1.0 - a * (0.253 + a * 0.12);
        if p < t {
            (p / t).powf(1.0 / a)
        } else {
            1.0 - (1.0 - (p - t) / (1.0 - t)).ln()
        }
    };

    for _ in 0..14 {
        if x <= 0.0 {
            return 0.0;
        }
        let err = gamma_p(a, x) - p;
        let pdf = (-x + a1 * x.ln() - gln).exp();
        if pdf == 0.0 {
            break;
        }
        let u = err / pdf;
        // Paso de Halley (Numerical Recipes invgammp).
        let t = u / (1.0 - 0.5 * (u * (a1 / x - 1.0)).min(1.0));
        x -= t;
        if x <= 0.0 {
            x = 0.5 * (x + t);
        }
        if t.abs() < EPS.max(1e-12) * x {
            break;
        }
    }
    x
}

/// CDF de la normal estándar, vía P(1/2, x²/2).
pub(crate) fn norm_cdf(x: f64) -> f64 {
    let p = gamma_p(0.5, 0.5 * x * x);
    if x >= 0.0 {
        0.5 + 0.5 * p
    } else {
        0.5 - 0.5 * p
    }
}

/// Inversa de la CDF normal estándar (algoritmo de Acklam, ~1e-9 relativo).
pub(crate) fn norm_ppf(p: f64) -> f64 {
    debug_assert!(p > 0.0 && p < 1.0);
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_690e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];
    const P_LOW: f64 = 0.02425;

    let tail = |q: f64| {
        let t = (-2.0 * q.ln()).sqrt();
        (((((C[0] * t + C[1]) * t + C[2]) * t + C[3]) * t + C[4]) * t + C[5])
            / ((((D[0] * t + D[1]) * t + D[2]) * t + D[3]) * t + 1.0)
    };

    if p < P_LOW {
        tail(p)
    } else if p <= 1.0 - P_LOW {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        -tail(1.0 - p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Referencias generadas con SciPy 2026-06-10 (ver docs/parity.md).
    #[test]
    fn ln_gamma_matches_scipy() {
        assert!((ln_gamma(0.5) - 0.572_364_942_924_7).abs() < 1e-12);
        assert!((ln_gamma(7.3) - 7.147_892_523_022_249).abs() < 1e-12);
        assert!((ln_gamma(5.0) - 24.0_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn digamma_matches_scipy() {
        assert!((digamma(0.7) - (-1.220_023_553_697_934_7)).abs() < 1e-10);
        assert!((digamma(15.0) - 2.674_346_661_660_793_6).abs() < 1e-12);
    }

    #[test]
    fn trigamma_matches_scipy() {
        assert!((trigamma(0.7) - 2.834_049_16).abs() < 1e-7);
        assert!((trigamma(15.0) - 0.068_938_23).abs() < 1e-7);
    }

    #[test]
    fn gamma_p_matches_scipy() {
        assert!((gamma_p(2.0, 3.0) - 0.800_851_726_528_544_2).abs() < 1e-12);
        assert!((gamma_p(0.5, 0.1) - 0.345_279_153_981_423_17).abs() < 1e-12);
        assert!((gamma_p(9.0, 12.0) - 0.844_972_218_232_537).abs() < 1e-12);
        assert_eq!(gamma_p(2.0, 0.0), 0.0);
    }

    #[test]
    fn gamma_p_inv_matches_scipy() {
        assert!((gamma_p_inv(2.0, 0.8) - 2.994_308_347_002_123_2).abs() < 1e-8);
        assert!((gamma_p_inv(0.5, 0.3) - 0.074_235_930_916_272_69).abs() < 1e-8);
    }

    #[test]
    fn gamma_p_inv_is_inverse_of_gamma_p() {
        for &a in &[0.3, 0.9, 1.7, 4.2, 25.0] {
            for &p in &[0.01, 0.1, 0.5, 0.9, 0.99] {
                let x = gamma_p_inv(a, p);
                assert!(
                    (gamma_p(a, x) - p).abs() < 1e-9,
                    "a={a}, p={p}: P(a, {x}) = {}",
                    gamma_p(a, x)
                );
            }
        }
    }

    #[test]
    fn norm_cdf_and_ppf_match_scipy() {
        assert!((norm_cdf(1.5) - 0.933_192_798_731_141_9).abs() < 1e-12);
        assert!((norm_ppf(0.975) - 1.959_963_984_540_054).abs() < 1e-8);
        assert!((norm_ppf(0.01) - (-2.326_347_874_040_840_8)).abs() < 1e-8);
        for &p in &[0.001, 0.05, 0.5, 0.95, 0.999] {
            assert!((norm_cdf(norm_ppf(p)) - p).abs() < 1e-9);
        }
    }
}
