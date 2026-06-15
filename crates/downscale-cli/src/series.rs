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

/// Matriz de predictores leída de un CSV `date,col1,col2,...`.
#[derive(Debug, Clone)]
pub struct Matrix {
    pub dates: Vec<String>,
    /// Nombres de las columnas de predictores (sin `date`).
    pub columns: Vec<String>,
    /// Datos aplanados fila-por-día; `n_features = columns.len()`.
    pub data: Vec<f64>,
}

impl Matrix {
    /// Número de predictores por día.
    pub fn n_features(&self) -> usize {
        self.columns.len()
    }

    /// Lee un CSV con primera columna `date` (ISO) y el resto numéricas.
    ///
    /// Salta filas con cualquier valor no numérico o centinela (`-9999`,
    /// `NA`) — las observaciones con huecos deben resolverse aguas arriba.
    pub fn read_csv(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("no se pudo leer {}", path.display()))?;
        let mut lines = text.lines();
        let header = lines
            .next()
            .with_context(|| format!("{}: archivo vacío", path.display()))?;
        let cols: Vec<&str> = header.split(',').map(str::trim).collect();
        if cols.len() < 2 || !cols[0].eq_ignore_ascii_case("date") {
            bail!(
                "{}: se esperaba cabecera 'date,<col>,...'; se encontró {:?}",
                path.display(),
                header
            );
        }
        let columns: Vec<String> = cols[1..].iter().map(|s| (*s).to_string()).collect();
        let n_features = columns.len();

        let mut dates = Vec::new();
        let mut data = Vec::new();
        for line in lines {
            let fields: Vec<&str> = line.split(',').map(str::trim).collect();
            if fields.len() != n_features + 1 {
                continue;
            }
            let date = fields[0];
            if date.len() < 8 || !date.as_bytes()[0].is_ascii_digit() || !date.contains('-') {
                continue;
            }
            let mut row = Vec::with_capacity(n_features);
            let mut ok = true;
            for f in &fields[1..] {
                match f.parse::<f64>() {
                    Ok(v) if v != -9999.0 && v.is_finite() => row.push(v),
                    _ => {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                dates.push(date.to_string());
                data.extend(row);
            }
        }
        if dates.is_empty() {
            bail!("{}: ninguna fila de predictores válida", path.display());
        }
        Ok(Self {
            dates,
            columns,
            data,
        })
    }

    /// Vector de predictores del día `i`.
    fn row(&self, i: usize) -> &[f64] {
        let n = self.n_features();
        &self.data[i * n..(i + 1) * n]
    }
}

/// Parea una matriz de predictores con una serie de observaciones por fecha
/// (intersección, en el orden de la matriz). Devuelve `(predictores
/// aplanados, n_features, obs)`.
pub fn pair_matrix_series(m: &Matrix, obs: &Series) -> Result<(Vec<f64>, usize, Vec<f64>)> {
    let index: HashMap<&str, f64> = obs
        .dates
        .iter()
        .map(String::as_str)
        .zip(obs.values.iter().copied())
        .collect();
    let n = m.n_features();
    let mut data = Vec::new();
    let mut y = Vec::new();
    for (i, d) in m.dates.iter().enumerate() {
        if let Some(&v) = index.get(d.as_str()) {
            data.extend_from_slice(m.row(i));
            y.push(v);
        }
    }
    if y.is_empty() {
        bail!("los predictores y las observaciones no comparten fechas");
    }
    Ok((data, n, y))
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

    #[test]
    fn matrix_reads_and_pairs() {
        let dir = std::env::temp_dir().join("downscale_cli_test_matrix");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("pred.csv");
        std::fs::write(
            &p,
            "date,pmsl,rh\n2000-01-01,1012.0,60\n2000-01-02,1008.0,-9999\n2000-01-03,1015.0,55\n",
        )
        .unwrap();
        let m = Matrix::read_csv(&p).unwrap();
        assert_eq!(m.columns, vec!["pmsl", "rh"]);
        assert_eq!(m.n_features(), 2);
        // La fila con -9999 se descarta.
        assert_eq!(m.dates, vec!["2000-01-01", "2000-01-03"]);
        assert_eq!(m.data, vec![1012.0, 60.0, 1015.0, 55.0]);

        let obs = Series {
            dates: vec!["2000-01-03".into(), "2000-01-01".into()],
            values: vec![3.0, 1.0],
        };
        let (data, n, y) = pair_matrix_series(&m, &obs).unwrap();
        assert_eq!(n, 2);
        assert_eq!(data, vec![1012.0, 60.0, 1015.0, 55.0]);
        assert_eq!(y, vec![1.0, 3.0]);
    }
}
