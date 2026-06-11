//! Interfaz de forzantes hacia motores hidrológicos (rainflow, snowmelt-rs).
//!
//! Contrato (ver `docs/forcing-interface.md`): forzantes **diarias,
//! contiguas y sin huecos**, con nombres de columna que el lector de
//! rainflow reconoce (`pr`, `pet`, `tmean`) y unidades fijas (mm/día, °C).
//! El eje temporal usa días desde la época Unix (aritmética de fechas
//! civiles propia, sin dependencias).

use std::fmt::Write as _;

use crate::error::{DownscaleError, Result, check_series};

/// Variable de forzante con nombre de columna y unidad canónicos.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Variable {
    /// Precipitación (mm/día). Columna `pr`.
    Precipitation,
    /// Evapotranspiración potencial (mm/día). Columna `pet`.
    Pet,
    /// Temperatura media (°C). Columna `tmean`.
    TemperatureMean,
}

impl Variable {
    /// Nombre de columna canónico (reconocido por el lector de rainflow).
    #[must_use]
    pub fn column_name(self) -> &'static str {
        match self {
            Variable::Precipitation => "pr",
            Variable::Pet => "pet",
            Variable::TemperatureMean => "tmean",
        }
    }

    /// Unidad del contrato.
    #[must_use]
    pub fn unit(self) -> &'static str {
        match self {
            Variable::Precipitation | Variable::Pet => "mm/d",
            Variable::TemperatureMean => "degC",
        }
    }
}

/// Días desde 1970-01-01 a partir de una fecha civil (algoritmo de
/// Hinnant, válido para todo el calendario gregoriano proléptico).
#[must_use]
pub fn epoch_day_from_civil(year: i64, month: u32, day: u32) -> i64 {
    let y = if month <= 2 { year - 1 } else { year };
    let era = if y >= 0 { y } else { y - 399 } / 400;
    let yoe = y - era * 400;
    let mp = i64::from(if month > 2 { month - 3 } else { month + 9 });
    let doy = (153 * mp + 2) / 5 + i64::from(day) - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// Fecha civil `(año, mes, día)` desde días de la época Unix.
#[must_use]
pub fn civil_from_epoch_day(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// Parsea `YYYY-MM-DD` a días de la época. Valida que la fecha exista.
pub fn parse_date(date: &str) -> Result<i64> {
    let invalid = || DownscaleError::InvalidDate {
        value: date.to_string(),
    };
    let bytes = date.as_bytes();
    if bytes.len() < 10 || bytes[4] != b'-' || bytes[7] != b'-' {
        return Err(invalid());
    }
    let year: i64 = date[..4].parse().map_err(|_| invalid())?;
    let month: u32 = date[5..7].parse().map_err(|_| invalid())?;
    let day: u32 = date[8..10].parse().map_err(|_| invalid())?;
    if !(1..=12).contains(&month) || !(1..=31).contains(&day) {
        return Err(invalid());
    }
    let z = epoch_day_from_civil(year, month, day);
    // Rechaza fechas inexistentes (p. ej. 2023-02-30) por ida y vuelta.
    if civil_from_epoch_day(z) != (year, month, day) {
        return Err(invalid());
    }
    Ok(z)
}

/// Formatea días de la época como `YYYY-MM-DD`.
#[must_use]
pub fn format_date(epoch_day: i64) -> String {
    let (y, m, d) = civil_from_epoch_day(epoch_day);
    format!("{y:04}-{m:02}-{d:02}")
}

/// Serie de forzante diaria contigua.
#[derive(Debug, Clone)]
pub struct ForcingSeries {
    variable: Variable,
    start_day: i64,
    values: Vec<f64>,
}

impl ForcingSeries {
    /// Construye desde fechas ISO y valores pareados, validando que el
    /// eje sea diario contiguo (sin huecos ni desorden) y los valores
    /// finitos.
    ///
    /// # Errors
    ///
    /// [`DownscaleError::InvalidDate`], [`DownscaleError::NonContiguous`],
    /// [`DownscaleError::LengthMismatch`], [`DownscaleError::NonFinite`] o
    /// serie vacía.
    pub fn from_dates(variable: Variable, dates: &[String], values: &[f64]) -> Result<Self> {
        if dates.len() != values.len() {
            return Err(DownscaleError::LengthMismatch {
                left_name: "dates",
                left: dates.len(),
                right_name: "values",
                right: values.len(),
            });
        }
        check_series(variable.column_name(), values, 1)?;
        let start_day = parse_date(&dates[0])?;
        let mut prev = start_day;
        for date in &dates[1..] {
            let day = parse_date(date)?;
            if day != prev + 1 {
                return Err(DownscaleError::NonContiguous {
                    name: variable.column_name(),
                    date: format_date(prev),
                    gap_days: day - prev,
                });
            }
            prev = day;
        }
        Ok(Self {
            variable,
            start_day,
            values: values.to_vec(),
        })
    }

    /// Variable de la serie.
    #[must_use]
    pub fn variable(&self) -> Variable {
        self.variable
    }

    /// Primer día (época Unix).
    #[must_use]
    pub fn start_day(&self) -> i64 {
        self.start_day
    }

    /// Último día (época Unix, inclusivo).
    #[must_use]
    pub fn end_day(&self) -> i64 {
        self.start_day + self.values.len() as i64 - 1
    }

    /// Valores diarios.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }
}

/// Conjunto de forzantes alineadas a un período común — la estructura que
/// consume un motor hidrológico.
///
/// # Ejemplo
///
/// ```
/// use downscale_core::forcing::{ForcingSeries, ForcingSet, Variable};
///
/// let dates: Vec<String> = ["2000-01-01", "2000-01-02", "2000-01-03"]
///     .iter().map(|s| s.to_string()).collect();
/// let pr = ForcingSeries::from_dates(Variable::Precipitation, &dates, &[0.0, 5.2, 1.1]).unwrap();
/// let pet = ForcingSeries::from_dates(Variable::Pet, &dates[1..], &[3.0, 3.1]).unwrap();
///
/// // Alinea al período común (2000-01-02 .. 2000-01-03).
/// let set = ForcingSet::align(vec![pr, pet]).unwrap();
/// assert_eq!(set.len(), 2);
/// let csv = set.to_csv();
/// assert!(csv.starts_with("date,pr,pet\n2000-01-02,5.2,3\n"));
/// ```
#[derive(Debug, Clone)]
pub struct ForcingSet {
    start_day: i64,
    columns: Vec<(Variable, Vec<f64>)>,
}

impl ForcingSet {
    /// Alinea series al período común (intersección de ejes temporales).
    ///
    /// # Errors
    ///
    /// - [`DownscaleError::InvalidParameter`] si no hay series o hay
    ///   variables repetidas.
    /// - [`DownscaleError::NoOverlap`] si la intersección es vacía.
    pub fn align(series: Vec<ForcingSeries>) -> Result<Self> {
        if series.is_empty() {
            return Err(DownscaleError::InvalidParameter {
                name: "series",
                value: 0.0,
                expected: ">= 1 serie",
            });
        }
        for (i, s) in series.iter().enumerate() {
            if series[..i].iter().any(|t| t.variable() == s.variable()) {
                return Err(DownscaleError::InvalidParameter {
                    name: "series",
                    value: i as f64,
                    expected: "variables sin repetir",
                });
            }
        }
        let start = series
            .iter()
            .map(ForcingSeries::start_day)
            .max()
            .expect("no vacío");
        let end = series
            .iter()
            .map(ForcingSeries::end_day)
            .min()
            .expect("no vacío");
        if end < start {
            return Err(DownscaleError::NoOverlap);
        }
        let len = (end - start + 1) as usize;
        let columns = series
            .into_iter()
            .map(|s| {
                let offset = (start - s.start_day) as usize;
                (s.variable, s.values[offset..offset + len].to_vec())
            })
            .collect();
        Ok(Self {
            start_day: start,
            columns,
        })
    }

    /// Número de días del período común.
    #[must_use]
    pub fn len(&self) -> usize {
        self.columns.first().map_or(0, |(_, v)| v.len())
    }

    /// `true` si no hay días.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Primer día del período común (época Unix).
    #[must_use]
    pub fn start_day(&self) -> i64 {
        self.start_day
    }

    /// Variables presentes, en orden de columna.
    #[must_use]
    pub fn variables(&self) -> Vec<Variable> {
        self.columns.iter().map(|(v, _)| *v).collect()
    }

    /// Valores de una variable, si está presente.
    #[must_use]
    pub fn series(&self, variable: Variable) -> Option<&[f64]> {
        self.columns
            .iter()
            .find(|(v, _)| *v == variable)
            .map(|(_, vals)| vals.as_slice())
    }

    /// CSV canónico del contrato: header `date,<col>...` y una fila por
    /// día. Es el formato que consume el lector de forzantes de rainflow.
    #[must_use]
    pub fn to_csv(&self) -> String {
        let mut out = String::from("date");
        for (v, _) in &self.columns {
            out.push(',');
            out.push_str(v.column_name());
        }
        out.push('\n');
        for i in 0..self.len() {
            out.push_str(&format_date(self.start_day + i as i64));
            for (_, vals) in &self.columns {
                let _ = write!(out, ",{}", vals[i]);
            }
            out.push('\n');
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dates(items: &[&str]) -> Vec<String> {
        items.iter().map(ToString::to_string).collect()
    }

    #[test]
    fn epoch_day_roundtrip_and_known_values() {
        assert_eq!(epoch_day_from_civil(1970, 1, 1), 0);
        assert_eq!(epoch_day_from_civil(2000, 3, 1), 11017);
        assert_eq!(civil_from_epoch_day(0), (1970, 1, 1));
        // Bisiestos: 2000-02-29 existe, 1900-02-29 no (siglo no bisiesto).
        assert_eq!(
            epoch_day_from_civil(2000, 2, 29) + 1,
            epoch_day_from_civil(2000, 3, 1)
        );
        for z in [-150_000, -1, 0, 1, 10_958, 20_000, 150_000] {
            let (y, m, d) = civil_from_epoch_day(z);
            assert_eq!(epoch_day_from_civil(y, m, d), z);
        }
    }

    #[test]
    fn parse_date_rejects_nonexistent_dates() {
        assert!(parse_date("2023-02-28").is_ok());
        assert!(parse_date("2023-02-29").is_err());
        assert!(parse_date("2024-02-29").is_ok());
        assert!(parse_date("2023-13-01").is_err());
        assert!(parse_date("garbage").is_err());
        assert_eq!(format_date(parse_date("1950-01-01").unwrap()), "1950-01-01");
    }

    #[test]
    fn series_detects_gaps_and_disorder() {
        let v = [1.0, 2.0, 3.0];
        let gap = ForcingSeries::from_dates(
            Variable::Precipitation,
            &dates(&["2000-01-01", "2000-01-02", "2000-01-04"]),
            &v,
        )
        .unwrap_err();
        assert_eq!(
            gap,
            DownscaleError::NonContiguous {
                name: "pr",
                date: "2000-01-02".into(),
                gap_days: 2
            }
        );
        assert!(
            ForcingSeries::from_dates(
                Variable::Precipitation,
                &dates(&["2000-01-02", "2000-01-01", "2000-01-03"]),
                &v,
            )
            .is_err()
        );
        // Cruce de año y de bisiesto, contiguo.
        assert!(
            ForcingSeries::from_dates(
                Variable::Precipitation,
                &dates(&["2000-02-28", "2000-02-29", "2000-03-01"]),
                &v,
            )
            .is_ok()
        );
    }

    #[test]
    fn align_intersects_periods() {
        let pr = ForcingSeries::from_dates(
            Variable::Precipitation,
            &dates(&["2000-01-01", "2000-01-02", "2000-01-03", "2000-01-04"]),
            &[1.0, 2.0, 3.0, 4.0],
        )
        .unwrap();
        let pet = ForcingSeries::from_dates(
            Variable::Pet,
            &dates(&["2000-01-03", "2000-01-04", "2000-01-05"]),
            &[30.0, 40.0, 50.0],
        )
        .unwrap();
        let set = ForcingSet::align(vec![pr, pet]).unwrap();
        assert_eq!(set.len(), 2);
        assert_eq!(format_date(set.start_day()), "2000-01-03");
        assert_eq!(set.series(Variable::Precipitation).unwrap(), &[3.0, 4.0]);
        assert_eq!(set.series(Variable::Pet).unwrap(), &[30.0, 40.0]);
    }

    #[test]
    fn align_rejects_disjoint_and_duplicates() {
        let a = ForcingSeries::from_dates(
            Variable::Precipitation,
            &dates(&["2000-01-01", "2000-01-02"]),
            &[1.0, 2.0],
        )
        .unwrap();
        let b = ForcingSeries::from_dates(
            Variable::Pet,
            &dates(&["2010-01-01", "2010-01-02"]),
            &[1.0, 2.0],
        )
        .unwrap();
        assert_eq!(
            ForcingSet::align(vec![a.clone(), b]).unwrap_err(),
            DownscaleError::NoOverlap
        );
        assert!(ForcingSet::align(vec![a.clone(), a]).is_err());
    }

    #[test]
    fn csv_matches_rainflow_contract() {
        let d = dates(&["2000-01-01", "2000-01-02"]);
        let pr = ForcingSeries::from_dates(Variable::Precipitation, &d, &[0.0, 5.5]).unwrap();
        let pet = ForcingSeries::from_dates(Variable::Pet, &d, &[3.0, 3.25]).unwrap();
        let t = ForcingSeries::from_dates(Variable::TemperatureMean, &d, &[14.5, 15.0]).unwrap();
        let set = ForcingSet::align(vec![pr, pet, t]).unwrap();
        assert_eq!(
            set.to_csv(),
            "date,pr,pet,tmean\n2000-01-01,0,3,14.5\n2000-01-02,5.5,3.25,15\n"
        );
    }
}
