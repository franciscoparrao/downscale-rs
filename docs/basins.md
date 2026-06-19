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

Cuatro cuencas near-natural: dos **pluviales** del centro-sur (GR4J puro,
Itata 8123001 y Perquilauquén 7330001) y dos **nivales** del Norte Chico
andino (HBV con rutina de nieve, Las Ramadas 4511002 y Choapa 4703002, elev.
~3100 m). Para cada una: ERA5 puntual en el aforo, corregida a CR2MET con
quantile mapping (precipitación multiplicativa; en las nivales también
temperatura, aditiva — clave para la nieve); tres forzantes (CR2MET / ERA5
cruda / ERA5 corregida) calibran el modelo con el split-sample de Klemeš.

## Resultados

### (A) Sesgo de la precipitación de cuenca vs CR2MET (KS)

| cuenca | modelo | ERA5 cruda | ERA5 corregida |
|---|---|---|---|
| Itata | GR4J | 0.073 | **0.006** |
| Perquilauquén | GR4J | 0.062 | **0.010** |
| Las Ramadas | HBV | 0.140 | **0.008** |
| Choapa | HBV | 0.266 | **0.008** |

### (B) Impacto en el modelo (KGE de validación, split-sample)

| cuenca | modelo | CR2MET | ERA5 cruda | ERA5 corregida |
|---|---|---|---|---|
| Itata | GR4J | 0.767 | 0.827 | 0.811 |
| Perquilauquén | GR4J | 0.842 | 0.813 | 0.851 |
| Las Ramadas | HBV | 0.439 | 0.342 | 0.393 |
| Choapa | HBV | 0.592 | 0.722 | 0.587 |

## Hallazgos

1. **La cadena multi-cuenca funciona end-to-end, en dos modelos.**
   downscale-rs corrige la forzante de cada cuenca y rainflow calibra GR4J
   (pluvial) o HBV con nieve (nival, donde también se corrige la
   temperatura) — el eslabón corrección → hidrología del ecosistema, sobre
   cuencas reales.

2. **La corrección mejora claramente la forzante, y el sesgo de partida
   sigue el gradiente climático.** El KS de precipitación cae a ~0.008 en
   las cuatro. El sesgo crudo es mayor en las nivales del Norte Chico
   semiárido (0.14 y 0.27) que en las pluviales del centro-sur (0.06–0.07)
   — el mismo gradiente árido-húmedo de `docs/multistation.md`.

3. **El impacto en el caudal es matizado — por una razón conocida.** En
   Perquilauquén y Las Ramadas la corrección ayuda (KGE 0.813 → 0.851,
   0.342 → 0.393); en Itata y Choapa la ERA5 cruda ya da buen KGE. Un modelo
   conceptual calibrado contra el caudal **absorbe el sesgo de volumen de la
   forzante** vía sus parámetros (x1 en GR4J; los 12 de HBV aún más) —
   compensación de parámetros / equifinalidad. Por eso el KGE de un modelo
   calibrado es relativamente insensible al sesgo sistemático de la forzante;
   el efecto es más errático en las nivales, que son intrínsecamente más
   difíciles (KGE ~0.4–0.6 vs ~0.8 pluvial).

## Lectura para el paper (EMS)

El experimento demuestra la integración con valor real y, sobre todo, evita
la conclusión simplista de que "el bias correction siempre mejora el
caudal". Lo correcto es más útil: cuando hay caudal para calibrar, el modelo
compensa el sesgo de la forzante; el bias correction importa donde **no** se
puede calibrar — escenarios climáticos futuros (GCM, sin caudal futuro),
cuencas sin aforo (transferencia de parámetros), o el análisis de la propia
precipitación corregida. Saber distinguir esos regímenes es justo lo que la
cadena de herramientas permite.

> Alcance: cuatro cuencas (dos pluviales, dos nivales); no es una evaluación
> exhaustiva del impacto hidrológico, sino una demostración de la cadena
> sobre ambos modelos y del fenómeno de compensación. Un set mayor de cuencas
> BNA queda como extensión natural (la cadena ya lo soporta vía `areal` +
> `rainflow batch`).
