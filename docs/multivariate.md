# Corrección multivariada: MBCn en datos reales

Fecha: 2026-06-16. Reproducible con `scripts/experiment_mbcn.py` (vía
bindings `downscale_rs`). Datos: precipitación + temperatura media diaria de
la estación DMC Quinta Normal 330020 (CR2) vs ERA5 puntual, 17 793 días
pareados 1969–2018 (la temperatura de la estación empieza en 1969), split
temporal 70/30.

## Qué demuestra

La corrección de sesgo univariada (quantile mapping por variable) corrige
las distribuciones marginales pero **no toca la estructura de dependencia
entre variables**. En clima mediterráneo la precipitación y la temperatura
correlacionan negativamente — los días lluviosos son más fríos — y un modelo
puede reproducir mal esa correlación. MBCn (Cannon 2018) corrige marginales
**y** dependencia.

ERA5 es reanálisis, sincronizado día a día con la realidad, así que la
correlación pr–temp día a día es significativa tanto en las observaciones
como en el modelo (a diferencia de un GCM, ver `docs/gcm-validation.md`).

## Resultados (período de validación)

### Correlación pr–temp

| | corr(pr, temp) | dist. a obs |
|---|---|---|
| observado | −0.158 | — |
| ERA5 crudo | −0.269 | 0.111 |
| QM univariado | −0.214 | 0.056 |
| **MBCn** | **−0.175** | **0.017** |

### KS marginal vs observado

| variable | crudo | QM univariado | MBCn |
|---|---|---|---|
| pr | 0.144 | 0.009 | 0.019 |
| temp | 0.053 | 0.030 | 0.050 |

## Hallazgos

1. **MBCn recupera la correlación observada; el univariado no.** ERA5
   sobreestima la anticorrelación pr–temp (−0.269 vs −0.158 observado). MBCn
   la corrige a −0.175 (a 0.017 de lo observado), mientras el QM univariado
   queda en −0.214 — más cerca del modelo que de la realidad. El QM por
   variable es una transformación monótona que altera algo la correlación de
   Pearson, pero no la dirige hacia las observaciones; esa es justamente la
   función de MBCn.

2. **Ambos corrigen bien las marginales.** El QM univariado las clava (KS pr
   0.009); MBCn las deja casi igual de bien (KS pr 0.019, temp 0.050) — el
   pequeño costo del reordenamiento iterativo. La precipitación es donde más
   mejora la corrección (KS 0.144 → ~0.01–0.02); la temperatura de ERA5 ya
   parte con buena distribución (KS crudo 0.053).

3. **El trade-off es claro:** si solo importan las distribuciones marginales,
   el QM univariado basta y es más simple; si el uso aguas abajo depende de
   la covariación entre variables (índices de aridez, balance hídrico,
   confort térmico, riesgo de incendio que combina sequía y calor), MBCn es
   la herramienta correcta.

## Lectura para el paper (EMS)

Cierra el conjunto de métodos con un caso multivariado real y honesto: el
motor entrega tanto la corrección univariada (rápida, marginales exactas)
como la multivariada (MBCn, dependencia corregida), y la elección depende de
si la aplicación aguas abajo es sensible a la covariación. Es el tipo de
distinción metodológica accionable que un software de modelado debe ofrecer
y documentar.
