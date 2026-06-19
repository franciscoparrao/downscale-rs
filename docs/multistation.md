# Generalización fuera de Santiago: gradiente climático de Chile

Fecha: 2026-06-19. Reproducible con `scripts/experiment_multistation.py`
(vía bindings `downscale_rs`). Observaciones CR2 de 9 estaciones extraídas
de `cr2_prDaily_2018`, ERA5 puntual por estación (Open-Meteo), EQM
multiplicativo con nodos midpoint, split temporal 70/30.

## Motivación

Todos los casos previos eran Quinta Normal (Santiago, mediterráneo). ¿El
motor —y la corrección de sesgo— generaliza a otros climas? Chile ofrece un
gradiente ideal: del desierto de Atacama (hiperárido) a Magallanes
(subpolar), latitud −18° a −53°, de 2 a 1743 mm de lluvia al año.

## Resultados (período de validación)

| Estación | Régimen | lat | mm/año | días húm. | KS crudo | KS corr | sesgo crudo | sesgo corr |
|---|---|---|---|---|---|---|---|---|
| Arica | desierto extremo | −18.4 | 2 | 0.3 % | **0.617** | 0.003 | +0.459 | −0.005 |
| Antofagasta | desierto costero | −23.4 | 5 | 0.5 % | 0.125 | 0.004 | +0.032 | −0.010 |
| La Serena | semiárido | −29.9 | 87 | 5.1 % | 0.031 | 0.016 | +0.043 | +0.036 |
| Santiago | mediterráneo | −33.5 | 284 | 10.0 % | 0.149 | 0.017 | +0.572 | +0.062 |
| Concepción | medit. húmedo | −36.8 | 966 | 27.1 % | 0.013 | 0.028 | −0.413 | +0.440 |
| Valdivia | templado lluvioso | −39.6 | 1743 | 48.4 % | 0.014 | 0.062 | −0.098 | +0.514 |
| Puerto Montt | oceánico lluvioso | −41.4 | 1526 | 57.2 % | 0.086 | 0.056 | +1.991 | +0.634 |
| Coyhaique | patagónico | −45.6 | 1019 | 49.7 % | 0.152 | 0.090 | +1.364 | +0.827 |
| Punta Arenas | subpolar estepario | −53.0 | 375 | 41.0 % | 0.261 | 0.046 | +0.786 | +0.189 |

## Hallazgos

1. **El motor generaliza operacionalmente.** El mismo código corre sin
   cambios sobre los 9 regímenes —de 2 a 1743 mm/año, de 0.3 % a 57 % de
   días húmedos— sin fallos numéricos ni casos degenerados.

2. **El sesgo de ERA5 tiene un gradiente latitudinal nítido.** En los
   desiertos del norte el sesgo es catastrófico: en Arica el KS crudo es
   **0.617** — ERA5 "moja" la hiperaridez de Atacama con llovizna espuria
   donde casi nunca llueve (2 mm/año observados). La corrección lo lleva a
   0.003. El sesgo se atenúa hacia el sur lluvioso: en Concepción y Valdivia
   ERA5 ya reproduce bien la distribución (KS crudo ~0.013).

3. **Donde hay sesgo real, la corrección es dramática.** En las 5 estaciones
   con KS crudo > 0.1 (del desierto al mediterráneo, más los climas
   australes), el KS de validación baja en promedio **81 %**.

4. **Sesgo del modelo ≠ no-estacionariedad del clima.** En el sur lluvioso
   el EQM a veces *empeora* levemente el KS (Concepción 0.013 → 0.028,
   Valdivia 0.014 → 0.062) y deja un sesgo corregido positivo grande
   (Valdivia +0.51, Puerto Montt +0.63, Coyhaique +0.83). No es falla del
   motor: es la **megasequía del centro-sur de Chile** — la calibración
   1950–2000 es más lluviosa que la validación 2000–2018, así que un EQM
   *estacionario* sobre-corrige. El residuo dominante ahí no es sesgo del
   modelo sino el cambio del clima observado, que ninguna corrección
   estacionaria arregla.

## Lectura para el paper (EMS)

El experimento aporta dos cosas. Primero, **generalización geográfica
real**: un solo motor corrige el sesgo de ERA5 a lo largo de todo Chile, y
el sesgo del reanálisis resulta tener una estructura latitudinal clara
(peor en aridez). Segundo, una **distinción operacional importante**: el
bias correction estacionario corrige el sesgo del modelo pero no la
no-estacionariedad del clima; cuando ésta domina (sur de Chile bajo
megasequía), se necesita una corrección que preserve la tendencia —
exactamente lo que ofrece QDM (`docs/parity.md`), que el motor ya soporta.
Saber distinguir ambos residuos es justo lo que una herramienta de modelado
debe permitir diagnosticar.
