# Cadena multi-cuenca: downscale-rs → rainflow

Fecha: 2026-06-19. Reproducible con `scripts/experiment_basins.py` (vía
bindings `downscale_rs` + binario `rainflow`). Cuencas CAMELS-CL
near-natural (Alvarez-Garreton et al. 2018), forzante de precipitación de
referencia CR2MET, caudal observado, período 1979–2016.

## Pregunta

El modelado hidrológico operacional en Chile usa CR2MET (precipitación
grillada regional, ~5 km) como forzante. ¿Qué pasa si se fuerza el modelo
con ERA5 cruda (global, gratis, pero sesgada)? ¿Y cuánto recupera
downscale-rs corrigiendo ERA5 contra CR2MET? Es la cadena de la familia de
motores en acción: corrección de sesgo → modelo lluvia-escorrentía.

Dos cuencas pluviales near-natural (GR4J puro): Río Itata en Cholguán
(8123001, 860 km²) y Río Perquilauquén en San Manuel (7330001, 502 km²).
Para cada una: ERA5 puntual en el aforo, corregida a CR2MET con quantile
mapping multiplicativo; tres forzantes (CR2MET / ERA5 cruda / ERA5
corregida) calibran GR4J con el split-sample de Klemeš.

## Resultados

### (A) Sesgo de la precipitación de cuenca vs CR2MET (KS)

| cuenca | ERA5 cruda | ERA5 corregida |
|---|---|---|
| Itata | 0.073 | **0.006** |
| Perquilauquén | 0.062 | **0.010** |

### (B) Impacto en GR4J (KGE de validación, split-sample)

| cuenca | CR2MET | ERA5 cruda | ERA5 corregida |
|---|---|---|---|
| Itata | 0.767 | 0.827 | 0.811 |
| Perquilauquén | 0.842 | 0.813 | 0.851 |

## Hallazgos

1. **La cadena multi-cuenca funciona end-to-end.** downscale-rs corrige la
   forzante de cada cuenca y rainflow calibra GR4J sobre ella — el eslabón
   corrección → hidrología del ecosistema, sobre cuencas reales.

2. **La corrección mejora claramente la forzante.** El sesgo de
   precipitación cae ~10× (KS 0.07 → 0.006, 0.06 → 0.01). El sesgo de
   partida de ERA5 es modesto aquí porque son cuencas lluviosas del
   centro-sur, donde ERA5 ya reproduce bien la distribución (consistente con
   `docs/multistation.md`).

3. **El impacto en el caudal es matizado — y por una razón conocida.** En
   Perquilauquén la corrección ayuda (KGE 0.813 → 0.851, superando incluso a
   CR2MET); en Itata ERA5 cruda ya da buen KGE. Un modelo conceptual
   calibrado contra el caudal **absorbe el sesgo de volumen de la forzante**
   vía sus parámetros (x1, la capacidad del reservorio de producción) —
   compensación de parámetros / equifinalidad. Por eso el KGE de un modelo
   calibrado es relativamente insensible al sesgo sistemático de la lluvia.

## Lectura para el paper (EMS)

El experimento demuestra la integración con valor real y, sobre todo, evita
la conclusión simplista de que "el bias correction siempre mejora el
caudal". Lo correcto es más útil: cuando hay caudal para calibrar, el modelo
compensa el sesgo de la forzante; el bias correction importa donde **no** se
puede calibrar — escenarios climáticos futuros (GCM, sin caudal futuro),
cuencas sin aforo (transferencia de parámetros), o el análisis de la propia
precipitación corregida. Saber distinguir esos regímenes es justo lo que la
cadena de herramientas permite.

> Alcance: dos cuencas pluviales donde el sesgo de ERA5 es moderado; no es
> una evaluación exhaustiva del impacto hidrológico, sino una demostración de
> la cadena y del fenómeno de compensación. Las cuencas nivales (HBV +
> temperatura) y un set mayor quedan como extensión.
