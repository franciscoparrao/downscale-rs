# Rendimiento

Fecha: 2026-06-13. Reproducible con `cargo bench -p downscale-core`
(microbenchmarks) y `scripts/benchmark_vs_python.py` (cross-language).
Toolchain rustc 1.94.1, perfil release (`lto = true`, `codegen-units = 1`).

## Microbenchmarks (criterion)

Tamaño del caso Quinta Normal: 17 430 días de calibración, 7 470 de
validación, 100 cuantiles. Mediana de la estimación de criterion.

| Operación | Tiempo |
|---|---|
| EQM — calibrar (`fit`) | 434 µs |
| EQM — aplicar (`apply`) | 165 µs |
| QDM — aplicar | 264 µs |
| QM paramétrico gamma — calibrar (MLE) | 330 µs |
| QM paramétrico gamma — aplicar (inversa gamma) | 3.55 ms |
| Análogos k=10 — predecir (4 predictores) | 283 ms |
| Regresión OLS — calibrar | 875 µs |

La aplicación del QM paramétrico es ~20× más lenta que la del EQM por la
inversa de la gamma incompleta (Newton/Halley por valor). Los análogos usan
un k-d tree (búsqueda k-NN exacta, O(log n) promedio por consulta): bajaron
de 1.14 s a 283 ms (~4×) en este microbenchmark, que es el peor caso —
predictores uniformes en 4-D, donde la poda del árbol es menos efectiva.
Con predictores reales correlacionados o `k` menor el speedup es mayor. El
test `kdtree_knn_matches_bruteforce` verifica que da el mismo conjunto de
análogos que la fuerza bruta.

## Comparación cross-language

Corrección completa (calibrar + aplicar) sobre el dataset Quinta Normal,
100 cuantiles, multiplicativo, mediana de 50 repeticiones. Para xclim/xsdba
se fuerza el cómputo (`.values`) porque su evaluación es perezosa (dask).

| Método | Implementación | Tiempo | Speedup |
|---|---|---|---|
| EQM | **downscale-rs** | **1.49 ms** | 1.0× (ref) |
| EQM | cmethods | 16.8 ms | 11× |
| EQM | xsdba (xclim) | 233 ms | 156× |
| QDM | **downscale-rs** | **1.93 ms** | 1.0× (ref) |
| QDM | xsdba (xclim) | 240 ms | 125× |

### Interpretación honesta

- **11× vs cmethods** es la comparación algoritmo-a-algoritmo: cmethods es
  NumPy puro sin framework, igual que downscale-rs opera sobre arrays
  planos. Es el speedup atribuible al lenguaje y al diseño.
- **125–156× vs xsdba** incluye el costo del stack xarray/dask que xclim
  impone (construcción de `DataArray`, grupos, scheduler perezoso). Es el
  speedup real que ve un usuario que adopta el flujo xclim, pero parte de
  la diferencia es overhead de framework, no de algoritmo.

Ambos números son reales y relevantes según con qué se compare. El mensaje
del motor: una corrección que toma cientos de milisegundos en el ecosistema
Python ocurre en ~1–2 ms aquí, sin cadena de dependencias y con un binario
único — lo que habilita corridas masivas (muchas estaciones/celdas, grandes
ensambles de GCM) que de otro modo serían costosas.

## Metodología

- Se mide solo el cómputo (datos ya en memoria); no incluye lectura de CSV
  ni serialización. Para xsdba los `DataArray` se construyen una vez fuera
  del lazo cronometrado.
- 3 repeticiones de calentamiento antes de medir, para excluir el costo de
  importación/JIT y el primer toque de caché.
- Mismo número de cuantiles (100), mismo tipo (multiplicativo), mismo
  dataset real en las tres implementaciones.
