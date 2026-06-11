# Interfaz de forzantes downscale-rs → rainflow

Fecha: 2026-06-11. Define el contrato entre la salida de downscale-rs y la
entrada de motores hidrológicos (rainflow hoy; snowmelt-rs/nowcast después).
Implementación: `downscale_core::forcing` + subcomando `downscale forcing`.

## Contrato de datos

**Formato**: CSV ancho, una fila por día.

```
date,pr,pet[,tmean]
1950-01-01,0,7.35,22.4
```

| Columna | Variable | Unidad | Obligatoria |
|---|---|---|---|
| `date` | fecha ISO `YYYY-MM-DD` | — | sí |
| `pr` | precipitación | mm/día | sí |
| `pet` | evapotranspiración potencial | mm/día | sí |
| `tmean` | temperatura media | °C | no (la usa HBV) |

Los nombres están dentro de los alias que el lector de rainflow
(`rainflow-cli/src/forcing.rs`) reconoce (`pr` ∈ PRECIP_NAMES,
`pet` ∈ PET_NAMES, `tmean` ∈ TEMP_NAMES).

**Invariantes que downscale-rs garantiza al escribir** (y rainflow exige
al leer):

1. **Eje diario contiguo**: sin huecos ni fechas desordenadas. Validado
   por `ForcingSeries::from_dates` (error `NonContiguous` con la fecha y
   el tamaño del salto).
2. **Sin valores faltantes**: NaN/centinelas rechazados (`NonFinite`).
   Los huecos de observaciones deben resolverse aguas arriba (la serie
   corregida proviene del modelo/reanálisis, que es continuo).
3. **Fechas reales del calendario**: `2023-02-29` se rechaza
   (`InvalidDate`), ida y vuelta por aritmética civil de Hinnant.
4. **Período común**: si las variables cubren rangos distintos,
   `ForcingSet::align` recorta a la intersección (error `NoOverlap` si es
   vacía).

## Uso

```bash
# 1. Corregir el sesgo de la variable del modelo/reanálisis
downscale correct --obs obs_pr.csv --model era5_pr.csv --kind mult \
  -o pr_corrected.csv

# 2. Ensamblar las forzantes para rainflow
downscale forcing --pr pr_corrected.csv --pet era5_pet.csv \
  --temp era5_tmean.csv -o forcing.csv

# 3. Correr el modelo hidrológico
rainflow run --forcing forcing.csv --x1 350 --x2 0 --x3 90 --x4 1.7 \
  --output qsim.csv
```

## Verificación end-to-end (2026-06-11)

Cadena completa probada con Quinta Normal: ERA5 pr corregido por EQM
(calibrado contra la estación CR2 330020) + ERA5 PET/tmean → `downscale
forcing` (24.905 días contiguos 1950–2018, alineación automática) →
`rainflow run` GR4J: 24.905 pasos simulados sin errores, `qsim` escrito.

## API Rust (para acople directo, sin CSV)

`downscale_core::forcing` expone `Variable`, `ForcingSeries` (eje validado)
y `ForcingSet` (alineación + `to_csv()` canónico). Cuando rainflow quiera
consumir forzantes en memoria (PyO3 o crate a crate), el tipo de
intercambio es `ForcingSet`; el CSV es la serialización de referencia.

## Pendiente / v0.2

- PET corregida (hoy se pasa cruda; Hargreaves desde tmin/tmax corregidas
  sería lo propio).
- Multi-sitio (un set por estación/celda) para lo semi-distribuido de
  rainflow v0.2.
- Metadatos de procedencia (método de corrección, período de calibración)
  como comentario de cabecera opcional.
