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

## Limitación honesta

Los productos CMIP6 de Open-Meteo provienen de **HighResMIP** y ya están
estadísticamente downscaled a ~10 km. Por eso su sesgo residual es
comparable al de ERA5 y no mayor, como sería el de un GCM **crudo** (~100–250
km, sin ajuste previo). El experimento demuestra el flujo GCM→local con datos
CMIP6 reales y la metodología distribucional correcta, pero no representa el
sesgo de magnitud completa de un GCM sin procesar. Para ese caso se
requeriría descargar GCMs crudos (p. ej. desde ESGF), pendiente para una
validación más exigente.

## Lectura para el paper (EMS)

El experimento aporta tres elementos al manuscrito: (a) cierra la brecha de
usar el motor con su caso de uso declarado (GCM, no solo reanálisis);
(b) ilustra la distinción metodológica clave entre corregir un reanálisis
(evaluable día a día) y un GCM (solo distribucional), que muchos usuarios
confunden; y (c) muestra que la corrección generaliza a través de cuatro
modelos independientes con sesgos de partida distintos.
