//! Evapotranspiración potencial de Hargreaves–Samani (1985) con radiación
//! extraterrestre FAO-56 (Allen et al. 1998, ecs. 21–25).
//!
//! `PET = 0.0023 · (0.408·Ra) · (Tmean + 17.8) · √(Tmax − Tmin)` en
//! mm/día, donde `Ra` es la radiación extraterrestre (MJ m⁻² día⁻¹)
//! calculada de latitud y día del año. Solo requiere tmin/tmax — el caso
//! típico al construir forzantes para rainflow desde temperaturas
//! corregidas.

use crate::error::{DownscaleError, Result, check_same_len, check_series};
use crate::forcing::civil_from_epoch_day;

/// Constante solar (MJ m⁻² min⁻¹), FAO-56.
const GSC: f64 = 0.082;

/// Radiación extraterrestre diaria `Ra` (MJ m⁻² día⁻¹) para una latitud
/// (grados, sur negativo) y día del año (1..=366). FAO-56 ec. 21.
#[must_use]
pub fn extraterrestrial_radiation(latitude_deg: f64, day_of_year: u32) -> f64 {
    let phi = latitude_deg.to_radians();
    let j = f64::from(day_of_year);
    let dr = 1.0 + 0.033 * (std::f64::consts::TAU * j / 365.0).cos();
    let delta = 0.409 * (std::f64::consts::TAU * j / 365.0 - 1.39).sin();
    // Ángulo horario de puesta de sol; clamp cubre latitudes polares
    // (sol de medianoche / noche polar).
    let ws = (-phi.tan() * delta.tan()).clamp(-1.0, 1.0).acos();
    24.0 * 60.0 / std::f64::consts::PI
        * GSC
        * dr
        * (ws * phi.sin() * delta.sin() + phi.cos() * delta.cos() * ws.sin())
}

/// PET de Hargreaves (mm/día) para series diarias de tmin/tmax (°C).
///
/// `days_of_year[i]` es el día del año (1..=366) de cada paso; ver
/// [`hargreaves_from_epoch_days`] para derivarlo de fechas.
/// Si `tmax < tmin` en algún día, el término `√ΔT` se trata como 0
/// (PET = 0) en vez de fallar: ocurre en datos reales por redondeo.
///
/// # Errors
///
/// Series vacías, de largos distintos, con NaN/inf, latitud fuera de
/// \[-90, 90\] o día del año fuera de 1..=366.
pub fn hargreaves(
    tmin: &[f64],
    tmax: &[f64],
    days_of_year: &[u32],
    latitude_deg: f64,
) -> Result<Vec<f64>> {
    check_series("tmin", tmin, 1)?;
    check_series("tmax", tmax, 1)?;
    check_same_len("tmin", tmin, "tmax", tmax)?;
    if days_of_year.len() != tmin.len() {
        return Err(DownscaleError::LengthMismatch {
            left_name: "days_of_year",
            left: days_of_year.len(),
            right_name: "tmin",
            right: tmin.len(),
        });
    }
    if !(-90.0..=90.0).contains(&latitude_deg) {
        return Err(DownscaleError::InvalidParameter {
            name: "latitude_deg",
            value: latitude_deg,
            expected: "en [-90, 90]",
        });
    }
    if let Some(&bad) = days_of_year.iter().find(|d| !(1..=366).contains(*d)) {
        return Err(DownscaleError::InvalidParameter {
            name: "days_of_year",
            value: f64::from(bad),
            expected: "en 1..=366",
        });
    }

    Ok(tmin
        .iter()
        .zip(tmax)
        .zip(days_of_year)
        .map(|((&tn, &tx), &doy)| {
            let ra = extraterrestrial_radiation(latitude_deg, doy);
            let trange = (tx - tn).max(0.0);
            let tmean = 0.5 * (tx + tn);
            // 0.408 convierte MJ m⁻² día⁻¹ a mm/día equivalentes.
            (0.0023 * 0.408 * ra * (tmean + 17.8) * trange.sqrt()).max(0.0)
        })
        .collect())
}

/// Como [`hargreaves`], derivando el día del año desde días de la época
/// Unix (el eje de [`crate::forcing::ForcingSeries`]).
///
/// # Errors
///
/// Igual que [`hargreaves`].
pub fn hargreaves_from_epoch_days(
    tmin: &[f64],
    tmax: &[f64],
    start_epoch_day: i64,
    latitude_deg: f64,
) -> Result<Vec<f64>> {
    let doys: Vec<u32> = (0..tmin.len() as i64)
        .map(|i| day_of_year(start_epoch_day + i))
        .collect();
    hargreaves(tmin, tmax, &doys, latitude_deg)
}

/// Día del año (1..=366) de un día de la época Unix.
#[must_use]
pub fn day_of_year(epoch_day: i64) -> u32 {
    let (y, m, d) = civil_from_epoch_day(epoch_day);
    let leap = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
    const CUM: [u32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    CUM[(m - 1) as usize] + d + u32::from(leap && m > 2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcing::parse_date;

    #[test]
    fn ra_matches_fao56_example_8() {
        // FAO-56, ejemplo 8: lat 20°S, 3 de septiembre (J=246) → Ra = 32.2.
        let ra = extraterrestrial_radiation(-20.0, 246);
        assert!((ra - 32.2).abs() < 0.1, "Ra = {ra}");
    }

    #[test]
    fn ra_polar_night_is_zero() {
        // Invierno austral profundo en lat -80°: noche polar, Ra ≈ 0.
        let ra = extraterrestrial_radiation(-80.0, 172);
        assert!(ra.abs() < 0.5, "Ra = {ra}");
    }

    #[test]
    fn hargreaves_plausible_for_santiago_summer() {
        // Santiago (-33.45°), enero: tmin 13, tmax 30 → PET ~5-7 mm/día.
        let pet = hargreaves(&[13.0], &[30.0], &[15], -33.45).unwrap();
        assert!((4.0..=8.0).contains(&pet[0]), "PET = {} mm/día", pet[0]);
        // Invierno: día corto y frío → mucho menor.
        let pet_winter = hargreaves(&[3.0], &[14.0], &[180], -33.45).unwrap();
        assert!(pet_winter[0] < pet[0] * 0.5);
    }

    #[test]
    fn inverted_range_yields_zero_not_error() {
        let pet = hargreaves(&[10.0], &[9.5], &[100], -33.0).unwrap();
        assert_eq!(pet[0], 0.0);
    }

    #[test]
    fn day_of_year_handles_leap_years() {
        assert_eq!(day_of_year(parse_date("2023-01-01").unwrap()), 1);
        assert_eq!(day_of_year(parse_date("2023-12-31").unwrap()), 365);
        assert_eq!(day_of_year(parse_date("2024-02-29").unwrap()), 60);
        assert_eq!(day_of_year(parse_date("2024-03-01").unwrap()), 61);
        assert_eq!(day_of_year(parse_date("2024-12-31").unwrap()), 366);
    }

    #[test]
    fn rejects_bad_latitude_and_doy() {
        assert!(hargreaves(&[1.0], &[2.0], &[1], 95.0).is_err());
        assert!(hargreaves(&[1.0], &[2.0], &[0], -33.0).is_err());
        assert!(hargreaves(&[1.0], &[2.0], &[367], -33.0).is_err());
    }
}
