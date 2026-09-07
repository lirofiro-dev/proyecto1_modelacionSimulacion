# Persona 5 — Comparación de métodos y análisis Monte Carlo

## Diseño experimental

Se compararon el método Polar y Aceptación–Rechazo usando 10,000 valores normales estándar por método, la misma cantidad de observaciones y semillas reproducibles. El tiempo reportado es la mediana de 7 ejecuciones. El error conjunto se definió como $\sqrt{(\bar x-0)^2+(S-1)^2}$.

## Comparación numérica

| Método | n | Media | S | Error | Tiempo (s) | KS p-valor |
|---|---:|---:|---:|---:|---:|---:|
| Polar | 10000 | 0.004343 | 1.000123 | 0.004345 | 0.002891 | 0.641234 |
| Aceptacion-Rechazo | 10000 | -0.004907 | 0.991740 | 0.009608 | 0.003515 | 0.949235 |

La media teórica es 0 y la desviación estándar teórica es 1. En la prueba KS individual, un p-valor mayor que 0.05 significa que la muestra no aporta evidencia suficiente para rechazar la distribución normal estándar; no significa que la distribución quede demostrada.

## Pruebas entre métodos

| Prueba | Estadístico | p-valor | Decisión con α = 0.05 |
|---|---:|---:|---|
| t de Welch: igualdad de medias | 0.656768 | 0.511337 | no se rechaza la hipótesis nula |
| Levene: igualdad de varianzas | 0.914144 | 0.339029 | no se rechaza la hipótesis nula |
| KS de dos muestras: igualdad de distribuciones | 0.008600 | 0.853388 | no se rechaza la hipótesis nula |

## Convergencia de P(éxito)

Para estudiar únicamente el efecto del generador normal, se ejecutó el mismo modelo dos veces: una usando Polar para las variables normales V y C, y otra usando Aceptación–Rechazo. En ambos escenarios se conservaron los parámetros de `simulacion.py`; los intervalos mostrados son intervalos de Wilson del 95 %.

| Método | n final | Éxitos | P(éxito) | IC 95 % | Error estándar |
|---|---:|---:|---:|---:|---:|
| Polar | 10000 | 9957 | 0.995700 | [0.994213, 0.996806] | 0.000654 |
| Aceptacion-Rechazo | 10000 | 9951 | 0.995100 | [0.993528, 0.996291] | 0.000698 |

La estimación se estabiliza conforme aumenta n y el intervalo de confianza se estrecha, que es el comportamiento esperado del error Monte Carlo, proporcional a $1/\sqrt{n}$.

## Conclusiones

El método con menor tiempo mediano fue **Polar** y el menor error conjunto en esta ejecución correspondió a **Polar**. Ambos resultados dependen de la muestra concreta, por lo que la decisión principal debe apoyarse conjuntamente en ajuste estadístico, tiempo y simplicidad del algoritmo.

Polar obtuvo media 0.004343 y S = 1.000123; Aceptación–Rechazo obtuvo media -0.004907 y S = 0.991740. Las pruebas anteriores permiten verificar si las diferencias observadas son estadísticamente distinguibles al nivel de 5 %.

## Revisión de consistencia del proyecto

La integración principal de V, C, W y N es coherente: V y C son normales, W es exponencial, N se genera mediante un proceso de Poisson, y las cuatro variables alimentan el índice R. La comparación usa directamente las funciones entregadas por las personas 1 y 2 y los parámetros del simulador de la persona 4.

Existe una diferencia de parámetros que el grupo debe resolver antes de entregar: `transformada_inversa_viento.ipynb` documenta media del viento de 12 m/s y tiempo medio entre anomalías de 300 s, mientras `simulacion.py` usa media del viento de 8 m/s y un proceso de anomalías con λ = 0.5 durante una misión de una unidad. Este análisis adopta `simulacion.py` como fuente de parámetros para no mezclar modelos dentro del mismo experimento.

También debe mantenerse la misma unidad temporal para λ y `TIEMPO_MISION`. Si la misión se expresa en segundos, λ debe expresarse en anomalías por segundo; si se usa otra unidad, ambos valores deben convertirse conjuntamente.

Con los pesos, penalizaciones y umbral actuales, la probabilidad estimada de éxito queda alrededor de 99.5 %. El cálculo es coherente con el código, pero el grupo debe confirmar si desea un escenario tan favorable. Si no fue intencional, deben calibrarse el umbral o las penalizaciones antes de interpretar el modelo; no conviene cambiarlos después de observar los resultados sin justificar el nuevo criterio.

## Gráficas

- `graficas/comparacion_distribuciones.png`
- `graficas/comparacion_tiempos.png`
- `graficas/convergencia_probabilidad_exito.png`
