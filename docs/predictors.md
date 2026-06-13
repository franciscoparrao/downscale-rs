# Experimento: predictores sinópticos en downscaling por análogos/regresión

Fecha: 2026-06-13. Reproducible con `scripts/fetch_synoptic_predictors.py`
(descarga) + `scripts/experiment_predictors.py` (vía bindings `downscale_rs`).

## Motivación

Con un solo predictor (la propia precipitación ERA5), los métodos de
downscaling por análogos y regresión no superan al quantile mapping
(ver `docs/parity.md` y el ejemplo `quinta_normal_downscaling`). La
pregunta del experimento: ¿el panorama cambia con predictores **sinópticos
de gran escala** (perfect-prognosis), que es el escenario para el que esos
métodos fueron diseñados?

## Diseño

- **Predictando**: precipitación diaria observada, estación DMC Quinta
  Normal 330020. **Modelo**: ERA5 puntual (Open-Meteo). 24 900 días
  pareados 1950–2018, split temporal 70/30 (17 430 / 7 470).
- **Predictores sinópticos** (media diaria, ERA5): presión a nivel del mar,
  humedad relativa 2 m, punto de rocío, nubosidad, déficit de presión de
  vapor, viento (u, v) y temperatura 2 m, más armónicos del día del año.
  El viento se descompone a componentes u/v antes de promediar.
- **Conjuntos**: `pr` (solo precipitación), `pr+season`, `synoptic`
  (estado atmosférico SIN precipitación), `all` (synoptic + pr).

## Resultados (período de validación, vs observado)

| método | predictores | RMSE (mm) | KS | sesgo (mm) |
|---|---|---|---|---|
| raw ERA5 | — | 3.734 | 0.1486 | +0.572 |
| EQM | pr | 3.069 | 0.0179 | +0.052 |
| QDM | pr | 3.068 | 0.0179 | +0.044 |
| análogos k=10 | pr | 3.607 | 0.1033 | −0.053 |
| análogos k=10 | synoptic | 3.525 | 0.2823 | −0.038 |
| **análogos k=10** | **all** | **3.143** | 0.2756 | −0.025 |
| análogos k=1 | pr | 4.420 | 0.0268 | −0.074 |
| **análogos k=1** | **synoptic** | 4.660 | **0.0058** | −0.004 |
| análogos k=1 | all | 4.269 | 0.0070 | −0.015 |
| regresión | synoptic | 3.862 | 0.5116 | +0.251 |
| **regresión** | **all** | **3.031** | 0.5431 | +0.157 |

(R² de calibración de la regresión: pr 0.30, synoptic 0.17, all 0.32.)

## Hallazgos

1. **Los predictores sinópticos mejoran el RMSE de los análogos.** Para
   k=10, enriquecer el predictor reduce el RMSE monótonamente:
   3.607 (`pr`) → 3.525 (`synoptic`) → **3.143 (`all`)**, una reducción de
   ~13 %. Con un solo predictor los análogos perdían claramente; con el
   estado atmosférico completo se acercan al quantile mapping.

2. **Análogos k=1 con predictores sinópticos da el mejor KS de toda la
   tabla** (0.0058), mejor incluso que EQM/QDM (0.0179). Al remuestrear
   valores observados reales seleccionados por análogos atmosféricos
   informativos, la distribución corregida es casi indistinguible de la
   observada — a costa de mayor error día a día (RMSE 4.66, el ruido
   intrínseco del remuestreo de un solo análogo).

3. **El trade-off RMSE↔distribución se confirma y se vuelve accionable.**
   La regresión minimiza el RMSE (3.031, el mejor) pero colapsa la varianza
   (KS 0.54, regresión a la media); análogos k=1 preserva la distribución
   (KS 0.006) sacrificando RMSE; k=10 con `all` equilibra (RMSE 3.14,
   KS 0.28). La elección de método depende de si el uso aguas abajo
   prioriza el acierto puntual (forzar un modelo lluvia-escorrentía día a
   día) o la fidelidad distribucional (frecuencia de extremos, percentiles).

## Lectura para el paper (EMS)

El experimento responde la objeción natural de un revisor a la tabla
univariada: con predictores apropiados, análogos y regresión cumplen su
rol. Más relevante, ilustra que **ningún método domina en todas las
métricas** — el motor entrega EQM/QDM (corrección distribucional robusta
con un predictor), regresión (mínimo error puntual) y análogos (remuestreo
fiel a la distribución, escalable a multivariado), y la elección es del
usuario según el criterio de evaluación. Es exactamente el tipo de
comparación reproducible que distingue una herramienta de modelado.
