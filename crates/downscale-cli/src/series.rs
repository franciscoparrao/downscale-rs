//! Lectura de series temporales `fecha,valor` desde CSV y emparejamiento
//! por fecha.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result, bail};

/// Serie temporal diaria: fechas ISO (strings) y valores.
#[derive(Debug, Clone)]
pub struct Series {
    pub dates: Vec<String>,
    pub values: Vec<f64>,
}

impl Series {
    /// Lee un CSV con columnas `fecha,valor`.
    ///
    /// - Se salta la primera línea si no parsea como dato (header).
    /// - Líneas vacías o sin coma se ignoran.
    /// - Valores no numéricos o centinela (`-9999`, `NA`) se descartan.
    pub fn read_csv(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("no se pudo leer {}", path.display()))?;
        let mut dates = Vec::new();
        let mut values = Vec::new();
        for line in text.lines() {
            let Some((date, value)) = line.split_once(',') else {
                continue;
            };
            let date = date.trim();
            let value = value.trim();
            // Solo filas cuyo primer campo parece fecha ISO (YYYY-MM-DD...).
            if date.len() < 8 || !date.as_bytes()[0].is_ascii_digit() || !date.contains('-') {
                continue;
            }
            let Ok(v) = value.parse::<f64>() else {
                continue; // header, NA, vacío
            };
            if v == -9999.0 || !v.is_finite() {
                continue; // centinela de dato faltante (convención DGA/CR2)
            }
            dates.push(date.to_string());
            values.push(v);
        }
        if dates.is_empty() {
            bail!("{}: ninguna fila 'fecha,valor' válida", path.display());
        }
        Ok(Self { dates, values })
    }
}

/// Empareja dos series por fecha exacta (intersección, orden de `a`).
pub fn pair_by_date(a: &Series, b: &Series) -> Result<(Vec<String>, Vec<f64>, Vec<f64>)> {
    let index: HashMap<&str, f64> = b
        .dates
        .iter()
        .map(String::as_str)
        .zip(b.values.iter().copied())
        .collect();
    let mut dates = Vec::new();
    let mut va = Vec::new();
    let mut vb = Vec::new();
    for (d, &v) in a.dates.iter().zip(&a.values) {
        if let Some(&w) = index.get(d.as_str()) {
            dates.push(d.clone());
            va.push(v);
            vb.push(w);
        }
    }
    if dates.is_empty() {
        bail!("las series no comparten ninguna fecha");
    }
    Ok((dates, va, vb))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pairs_intersection_in_order() {
        let a = Series {
            dates: vec![
                "2000-01-01".into(),
                "2000-01-02".into(),
                "2000-01-03".into(),
            ],
            values: vec![1.0, 2.0, 3.0],
        };
        let b = Series {
            dates: vec!["2000-01-03".into(), "2000-01-01".into()],
            values: vec![30.0, 10.0],
        };
        let (dates, va, vb) = pair_by_date(&a, &b).unwrap();
        assert_eq!(dates, vec!["2000-01-01", "2000-01-03"]);
        assert_eq!(va, vec![1.0, 3.0]);
        assert_eq!(vb, vec![10.0, 30.0]);
    }

    #[test]
    fn read_csv_skips_headers_and_sentinels() {
        let dir = std::env::temp_dir().join("downscale_cli_test_read");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("s.csv");
        std::fs::write(
            &p,
            "lat,lon\n-33,70\n\ndate,pr\n2000-01-01,1.5\n2000-01-02,-9999\n2000-01-03,NA\n2000-01-04,0\n",
        )
        .unwrap();
        let s = Series::read_csv(&p).unwrap();
        assert_eq!(s.dates, vec!["2000-01-01", "2000-01-04"]);
        assert_eq!(s.values, vec![1.5, 0.0]);
    }
}
