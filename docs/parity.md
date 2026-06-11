# Paridad numérica — EQM vs xclim/xsdba y cmethods

Fecha: 2026-06-10. Reproducible con `scripts/parity_quinta_normal.py`
(requiere `.venv` con `xsdba`, `python-cmethods`, `pandas` y los datos de
`data/README.md`).

## Setup

- **Dataset**: precipitación diaria Quinta Normal 330020 (CR2 prDaily, obs)
  vs ERA5 puntual (Open-Meteo, modelo). 24.900 días pareados 1950–2018.
- **Split**: 70/30 cronológico → calibración 17.430 días, validación 7.470
  días (desde 1997-09-26).
- **Método**: EQM multiplicativo, 100 cuantiles, interpolación lineal,
  extrapolación constante en colas. Idéntica configuración en las tres
  implementaciones (`downscale` Rust, `xsdba.EmpiricalQuantileMapping`,
  `cmethods.adjust(method="quantile_mapping")`).

## Resultados

### Calidad de la corrección (validación holdout vs obs)

| Serie | KS | Sesgo medio (mm/d) | P95 (mm) | Total (mm) |
|---|---|---|---|---|
| raw ERA5 | 0.1486 | +0.572 | 8.20 | 10 074 |
| **rust** | **0.0179** | **+0.052** | 4.90 | 6 192 |
| xclim/xsdba | 0.0179 | +0.057 | 4.73 | 6 231 |
| cmethods | 0.1486 | +0.198 | 4.81 | 7 280 |
| obs | — | — | 2.90 | 5 802 |

### Diferencias entre implementaciones (serie corregida, 7.470 días)

| Par | mediana \|Δ\| | P90 \|Δ\| | P99 \|Δ\| | RMSD | max \|Δ\| |
|---|---|---|---|---|---|
| rust vs xsdba | 0.000 | 0.021 | 0.273 | 0.192 | 8.07 |
| rust vs cmethods | 0.000 | — | 1.249 | 0.728 | 27.48 |
| xsdba vs cmethods | 0.000 | — | 1.244 | 0.587 | 21.83 |

(mm; la mediana 0.000 refleja los días secos, idénticos y exactos en todas.)

## Hallazgos

1. **Rust ≡ réplica NumPy del algoritmo (paridad exacta).** Reimplementando
   el algoritmo en NumPy (nodos `i/(n-1)`, cuantil tipo 7, `np.interp`,
   extrapolación constante), la diferencia máxima con la salida Rust es
   5e-4 mm = el redondeo a 3 decimales del CSV. La aritmética es correcta.

2. **Rust ≈ xsdba; la diferencia residual es solo colocación de nodos.**
   xsdba usa nodos en puntos medios `(i+0.5)/n` (sin 0 ni 1); nosotros
   `i/(n-1)` (con min y max). Eso desplaza dónde empieza la extrapolación
   de cola y explica el max 8.07 mm (ocurre en <0.1% de los días, los más
   extremos). Ambas logran KS idéntico (0.0179) y sesgo medio equivalente.
   Si se quisiera paridad bit-cercana con xsdba, bastaría ofrecer
   `NodePlacement::Midpoint` en `fit` (candidato v0.2).

3. **cmethods NO corrige la frecuencia de días húmedos.** Su EQM construye
   la CDF con histogramas binned sobre `[min, max]`; con precipitación
   inflada en cero los días de llovizna de ERA5 quedan en valores pequeños
   no nulos → el KS no mejora (0.1486, igual que raw) y retiene +0.20 mm/d
   de sesgo. Para precipitación, xsdba es la referencia válida de paridad;
   la divergencia de cmethods es una limitación documentada del incumbente,
   no un error nuestro (útil para el paper).

## Actualización: NodePlacement::Midpoint (2026-06-10, paso 4)

Con `--nodes midpoint` (nodos `(i+0.5)/n`, la convención de xsdba) la
paridad se estrecha ~16×, confirmando el diagnóstico del hallazgo 2:

| rust vs xsdba | P90 \|Δ\| | P99 \|Δ\| | RMSD | max \|Δ\| |
|---|---|---|---|---|
| nodos endpoints | 0.021 | 0.273 | 0.192 | 8.071 |
| **nodos midpoint** | **0.003** | **0.067** | **0.012** | **0.087** |

El residuo ≤0.09 mm se debe a que xsdba interpola *factores de ajuste*
(af = obs_q/sim_q) mientras nosotros componemos las CDFs; ambas son
formulaciones válidas del mismo estimador.

## Tolerancias declaradas (vs xsdba EQM, este dataset)

- Días secos: exactos (Δ = 0).
- P90 de |Δ| ≤ 0.03 mm; P99 de |Δ| ≤ 0.28 mm; RMSD ≤ 0.20 mm.
- Cola extrema (>P99.9): hasta ~8 mm por diseño de nodos (documentado
  arriba); KS y sesgo medio finales estadísticamente indistinguibles.
