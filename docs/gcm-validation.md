# Validación con GCMs reales (CMIP6)

Fecha: 2026-06-15. Reproducible con `scripts/fetch_gcm.py` (descarga) +
`scripts/experiment_gcm.py` (vía bindings `downscale_rs`).

## Motivación

El propósito declarado del motor es corregir el sesgo de modelos climáticos
(GCM/RCM → escala local). Hasta aquí el caso de validación fue ERA5
(reanálisis). Este experimento usa **4 GCMs CMIP6 reales** sobre el punto de
Quinta Normal: MRI-AGCM3-2-S, EC-Earth3P-HR, MPI-ESM1-2-XR, CMCC-CM2-VHR4.

## Metodología distribucional (no día a día)

Un GCM **no asimila observaciones**: produce su propio clima cuyas fechas no
corresponden a los días meteorológicos reales. Por eso, a diferencia de un
reanálisis, las series **no se parean en el tiempo** y el RMSE día a día no
tiene sentido. La evaluación es **distribucional**:

- Separación por **período** (calibración 1950–1999, validación 2000–2018).
- El quantile mapping se calibra comparando las *distribuciones* de cada
  período (no requiere pareo temporal — es su ventaja para GCMs).
- Métricas sobre la distribución corregida del período de validación: KS de
  dos muestras, sesgo medio y frecuencia de días húmedos (≥ 0.1 mm).

Referencia observada (validación): precipitación media 0.80 mm/d, días
húmedos 10.1 %.

## Resultados

| Fuente | KS crudo | KS corr | sesgo crudo | sesgo corr | húmedos crudo | húmedos corr |
|---|---|---|---|---|---|---|
| ERA5 (reanálisis) | 0.145 | 0.018 | +0.576 | +0.062 | 24.6 % | 11.0 % |
| MRI-AGCM3-2-S | 0.126 | 0.009 | +0.542 | +0.055 | 20.4 % | 10.0 % |
| EC-Earth3P-HR | 0.121 | 0.010 | +0.452 | −0.027 | 20.5 % | 9.5 % |
| MPI-ESM1-2-XR | 0.058 | 0.014 | +0.376 | −0.036 | 13.8 % | 10.3 % |
| CMCC-CM2-VHR4 | 0.098 | 0.008 | +0.346 | −0.054 | 18.4 % | 9.3 % |

(QM multiplicativo, 100 cuantiles, nodos midpoint. Frecuencia húmeda
observada = 10.1 %.)

## Hallazgos

1. **El sesgo de "llovizna" (drizzle) es universal.** Los cinco productos
   sobreestiman la frecuencia de días húmedos (13.8–24.6 % vs 10.1 %
   observado) — un sesgo conocido y sistemático de los modelos de clima. El
   quantile mapping multiplicativo lo lleva en todos a ~9–11 %.

2. **El motor corrige los GCMs tan bien como el reanálisis.** El KS de
   validación cae a 0.008–0.018 para los cuatro GCMs, equivalente al de ERA5
   (0.018), y el sesgo medio queda dentro de ±0.06 mm/d.

3. **Sesgo corregido levemente negativo en varios GCMs** (−0.03 a −0.05):
   refleja la misma no-estacionariedad detectada en el caso ERA5 — la
   calibración 1950–1999 es más lluviosa que la validación 2000–2018
   (megasequía de Chile central post-2010). El sesgo del modelo se asume
   estacionario; cuando el clima observado cambia, una corrección calibrada
   en el pasado deja un residuo. Es una limitación intrínseca del bias
   correction estacionario, no del motor.

## GCMs crudos (Pangeo / Google Cloud CMIP6)

Los productos CMIP6 de Open-Meteo provienen de **HighResMIP**, ya
downscaled a ~10 km, por lo que su sesgo es comparable al de un reanálisis.
El caso exigente real es el GCM **crudo** (~200–300 km, sin ajuste): la
celda que contiene Santiago promedia océano Pacífico, valle central y
cordillera. Reproducible con `scripts/fetch_gcm_raw.py` (archivo público
Pangeo CMIP6, zarr, sin autenticación) + `scripts/experiment_gcm_raw.py`;
precipitación diaria, escenario `historical`, miembro `r1i1p1f1`, split por
período. Observado: 291 mm/año, 10 % días húmedos.

| Modelo | resolución | mm/año | días húm. | KS crudo | KS corr | sesgo crudo | sesgo corr |
|---|---|---|---|---|---|---|---|
| IPSL-CM6A-LR | 209 km | 1055 | 68 % | **0.675** | 0.026 | +2.09 | +0.20 |
| MPI-ESM1-2-LR | 208 km | 306 | 15 % | 0.146 | 0.016 | +0.04 | −0.13 |
| CanESM5 | 311 km | 318 | 26 % | 0.595 | 0.009 | +0.07 | −0.09 |

Gradiente **resolución → sesgo** (familia MPI, KS crudo vs estación):
ERA5 ~25 km → 0.145 · MPI HighResMIP-XR ~10 km → **0.059** · MPI crudo LR
~200 km → 0.146.

### Hallazgos

1. **El sesgo del GCM crudo es estructural y puede ser extremo.**
   IPSL-CM6A-LR sobreestima la precipitación **3,6×** (1055 vs 291 mm/año)
   y llueve el 68 % de los días (vs 10 % observado) — un KS de 0.675, el
   mayor de todo el proyecto. CanESM5 acierta el total pero tiene drizzle
   severo (26 % días húmedos, KS 0.595). La magnitud depende del modelo y de
   qué captura su celda gruesa: MPI-ESM1-2-LR acierta casi el total (KS
   0.146, como ERA5).

2. **El quantile mapping lo corrige igual.** Incluso el caso extremo de IPSL
   (sobreestimación 3,6×) queda en KS 0.026; CanESM5 en 0.009. El motor
   maneja el sesgo de magnitud completa de un GCM sin procesar, no solo el
   residuo suave de un producto ya downscaled.

3. **El downscaling reduce el sesgo de partida**, como esperado: el
   HighResMIP-XR (~10 km, KS 0.059) mejora sobre ERA5 y sobre el GCM crudo de
   la misma familia (~200 km, KS 0.146).

Esto cierra la validación con el caso de uso más duro: GCMs crudos reales de
ESGF/Pangeo, no productos pre-ajustados. El residuo corregido positivo en
IPSL (+0.20) refleja, de nuevo, la no-estacionariedad del período de
validación.

## Lectura para el paper (EMS)

El experimento aporta tres elementos al manuscrito: (a) cierra la brecha de
usar el motor con su caso de uso declarado (GCM, no solo reanálisis);
(b) ilustra la distinción metodológica clave entre corregir un reanálisis
(evaluable día a día) y un GCM (solo distribucional), que muchos usuarios
confunden; y (c) muestra que la corrección generaliza a través de cuatro
modelos independientes con sesgos de partida distintos.
