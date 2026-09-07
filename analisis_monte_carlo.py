import csv
import math
import os
import random
import statistics
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from metodo_aceptacion_rechazo import generar_normales as generar_normales_ar
from metodo_polar import generar_normales as generar_normales_polar
from metodo_inversa import generar_viento
from proceso_poisson import generar_cantidad_anomalias
from simulacion import (
    DESVIACION_COMBUSTIBLE,
    DESVIACION_VELOCIDAD,
    LAMBDA_ANOMALIAS,
    LAMBDA_VIENTO,
    MEDIA_COMBUSTIBLE,
    MEDIA_VELOCIDAD,
    TIEMPO_MISION,
    simular_lanzamiento,
)


SEMILLA = 42
N_COMPARACION = 10000
REPETICIONES_TIEMPO = 7
TAMANOS_CONVERGENCIA = (100, 500, 1000, 2500, 5000, 10000)


def generar_normales_metodo(metodo, n, semilla):
    if metodo == "Polar":
        random.seed(semilla)
        return generar_normales_polar(n)
    if metodo == "Aceptacion-Rechazo":
        muestra, estadisticas = generar_normales_ar(n, semilla=semilla)
        return muestra
    raise ValueError("Método desconocido")


def medir_tiempo(metodo, n, semilla, repeticiones=REPETICIONES_TIEMPO):
    tiempos = []
    for repeticion in range(repeticiones):
        inicio = time.perf_counter()
        generar_normales_metodo(metodo, n, semilla + repeticion)
        tiempos.append(time.perf_counter() - inicio)
    return statistics.median(tiempos)


def estadisticas_muestra(metodo, muestra, tiempo_segundos):
    n = len(muestra)
    media = statistics.fmean(muestra)
    desviacion = statistics.stdev(muestra)
    error_media = abs(media)
    error_desviacion = abs(desviacion - 1.0)
    error_conjunto = math.hypot(error_media, error_desviacion)
    ks = stats.kstest(muestra, "norm")
    return {
        "metodo": metodo,
        "n": n,
        "media": media,
        "desviacion": desviacion,
        "error_media": error_media,
        "error_desviacion": error_desviacion,
        "error_conjunto": error_conjunto,
        "tiempo_segundos": tiempo_segundos,
        "ks_estadistico": ks.statistic,
        "ks_p_valor": ks.pvalue,
    }


def comparar_generadores(n=N_COMPARACION, semilla=SEMILLA):
    muestras = {}
    resumen = []
    for indice, metodo in enumerate(("Polar", "Aceptacion-Rechazo")):
        semilla_metodo = semilla + indice * 1000
        muestra = generar_normales_metodo(metodo, n, semilla_metodo)
        tiempo = medir_tiempo(metodo, n, semilla_metodo)
        muestras[metodo] = muestra
        resumen.append(estadisticas_muestra(metodo, muestra, tiempo))

    polar = muestras["Polar"]
    aceptacion = muestras["Aceptacion-Rechazo"]
    prueba_media = stats.ttest_ind(polar, aceptacion, equal_var=False)
    prueba_varianza = stats.levene(polar, aceptacion, center="median")
    prueba_distribucion = stats.ks_2samp(polar, aceptacion)
    pruebas = [
        {
            "prueba": "t de Welch: igualdad de medias",
            "estadistico": prueba_media.statistic,
            "p_valor": prueba_media.pvalue,
        },
        {
            "prueba": "Levene: igualdad de varianzas",
            "estadistico": prueba_varianza.statistic,
            "p_valor": prueba_varianza.pvalue,
        },
        {
            "prueba": "KS de dos muestras: igualdad de distribuciones",
            "estadistico": prueba_distribucion.statistic,
            "p_valor": prueba_distribucion.pvalue,
        },
    ]
    return muestras, resumen, pruebas


def generar_variables_modelo(metodo, n, semilla):
    normales_velocidad = generar_normales_metodo(metodo, n, semilla)
    normales_combustible = generar_normales_metodo(metodo, n, semilla + 1)
    velocidades = [
        MEDIA_VELOCIDAD + DESVIACION_VELOCIDAD * z
        for z in normales_velocidad
    ]
    combustibles = [
        MEDIA_COMBUSTIBLE + DESVIACION_COMBUSTIBLE * z
        for z in normales_combustible
    ]
    random.seed(semilla + 2)
    vientos = generar_viento(n, media=1.0 / LAMBDA_VIENTO)
    random.seed(semilla + 3)
    anomalias = generar_cantidad_anomalias(
        n, TIEMPO_MISION, LAMBDA_ANOMALIAS
    )
    return velocidades, combustibles, vientos, anomalias


def intervalo_wilson(exitos, n, z=1.96):
    proporcion = exitos / n
    denominador = 1.0 + z**2 / n
    centro = (proporcion + z**2 / (2.0 * n)) / denominador
    margen = z * math.sqrt(
        proporcion * (1.0 - proporcion) / n + z**2 / (4.0 * n**2)
    ) / denominador
    return centro - margen, centro + margen


def analizar_convergencia(tamanos=TAMANOS_CONVERGENCIA, semilla=SEMILLA):
    maximo = max(tamanos)
    filas = []
    for indice, metodo in enumerate(("Polar", "Aceptacion-Rechazo")):
        variables = generar_variables_modelo(
            metodo, maximo, semilla + indice * 1000
        )
        exitos_acumulados = []
        total_exitos = 0
        for valores in zip(*variables):
            total_exitos += int(simular_lanzamiento(*valores)["exito"])
            exitos_acumulados.append(total_exitos)
        for n in tamanos:
            exitos = exitos_acumulados[n - 1]
            inferior, superior = intervalo_wilson(exitos, n)
            filas.append(
                {
                    "metodo": metodo,
                    "n": n,
                    "exitos": exitos,
                    "probabilidad_exito": exitos / n,
                    "ic95_inferior": inferior,
                    "ic95_superior": superior,
                    "error_estandar": math.sqrt(
                        (exitos / n) * (1.0 - exitos / n) / n
                    ),
                }
            )
    return filas


def guardar_csv(filas, ruta):
    os.makedirs(os.path.dirname(ruta), exist_ok=True)
    with open(ruta, "w", newline="", encoding="utf-8") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=filas[0].keys())
        escritor.writeheader()
        escritor.writerows(filas)


def grafica_distribuciones(muestras, ruta):
    figura, eje = plt.subplots(figsize=(9, 5))
    eje.hist(
        muestras["Polar"], bins=50, density=True, alpha=0.5,
        label="Polar", color="#3568A8"
    )
    eje.hist(
        muestras["Aceptacion-Rechazo"], bins=50, density=True, alpha=0.5,
        label="Aceptación–Rechazo", color="#E18B3A"
    )
    x = [(-4.0 + 8.0 * i / 500) for i in range(501)]
    y = [math.exp(-(valor**2) / 2.0) / math.sqrt(2.0 * math.pi) for valor in x]
    eje.plot(x, y, color="black", linewidth=1.6, label="Normal estándar teórica")
    eje.set_xlabel("Valor generado")
    eje.set_ylabel("Densidad")
    eje.set_title("Comparación de generadores normales")
    eje.legend()
    figura.tight_layout()
    figura.savefig(ruta, dpi=180)
    plt.close(figura)


def grafica_tiempos(resumen, ruta):
    figura, eje = plt.subplots(figsize=(7, 4.5))
    metodos = [fila["metodo"] for fila in resumen]
    tiempos = [fila["tiempo_segundos"] * 1000.0 for fila in resumen]
    barras = eje.bar(metodos, tiempos, color=["#3568A8", "#E18B3A"])
    eje.bar_label(barras, fmt="%.2f ms", padding=3)
    eje.set_ylabel("Mediana del tiempo para 10 000 valores (ms)")
    eje.set_title("Costo computacional de los métodos")
    figura.tight_layout()
    figura.savefig(ruta, dpi=180)
    plt.close(figura)


def grafica_convergencia(filas, ruta):
    figura, eje = plt.subplots(figsize=(9, 5))
    colores = {"Polar": "#3568A8", "Aceptacion-Rechazo": "#E18B3A"}
    etiquetas = {"Polar": "Polar", "Aceptacion-Rechazo": "Aceptación–Rechazo"}
    for metodo in colores:
        datos = [fila for fila in filas if fila["metodo"] == metodo]
        x = [fila["n"] for fila in datos]
        y = [fila["probabilidad_exito"] for fila in datos]
        inferior = [fila["ic95_inferior"] for fila in datos]
        superior = [fila["ic95_superior"] for fila in datos]
        eje.plot(x, y, marker="o", color=colores[metodo], label=etiquetas[metodo])
        eje.fill_between(x, inferior, superior, color=colores[metodo], alpha=0.15)
    eje.set_xscale("log")
    eje.set_xlabel("Número de simulaciones")
    eje.set_ylabel("Probabilidad estimada de éxito")
    eje.set_title("Convergencia de la estimación Monte Carlo")
    eje.grid(alpha=0.25)
    eje.legend()
    figura.tight_layout()
    figura.savefig(ruta, dpi=180)
    plt.close(figura)


def conclusion_prueba(p_valor, alfa=0.05):
    if p_valor < alfa:
        return "se rechaza la hipótesis nula"
    return "no se rechaza la hipótesis nula"


def crear_informe(resumen, pruebas, convergencia, ruta):
    por_metodo = {fila["metodo"]: fila for fila in resumen}
    final = {
        fila["metodo"]: fila
        for fila in convergencia
        if fila["n"] == max(TAMANOS_CONVERGENCIA)
    }
    rapido = min(resumen, key=lambda fila: fila["tiempo_segundos"])
    preciso = min(resumen, key=lambda fila: fila["error_conjunto"])
    lineas = [
        "# Persona 5 — Comparación de métodos y análisis Monte Carlo",
        "",
        "## Diseño experimental",
        "",
        f"Se compararon el método Polar y Aceptación–Rechazo usando {N_COMPARACION:,} valores normales estándar por método, la misma cantidad de observaciones y semillas reproducibles. El tiempo reportado es la mediana de {REPETICIONES_TIEMPO} ejecuciones. El error conjunto se definió como $\\sqrt{{(\\bar x-0)^2+(S-1)^2}}$.",
        "",
        "## Comparación numérica",
        "",
        "| Método | n | Media | S | Error | Tiempo (s) | KS p-valor |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for fila in resumen:
        lineas.append(
            f"| {fila['metodo']} | {fila['n']} | {fila['media']:.6f} | "
            f"{fila['desviacion']:.6f} | {fila['error_conjunto']:.6f} | "
            f"{fila['tiempo_segundos']:.6f} | {fila['ks_p_valor']:.6f} |"
        )
    lineas.extend([
        "",
        "La media teórica es 0 y la desviación estándar teórica es 1. En la prueba KS individual, un p-valor mayor que 0.05 significa que la muestra no aporta evidencia suficiente para rechazar la distribución normal estándar; no significa que la distribución quede demostrada.",
        "",
        "## Pruebas entre métodos",
        "",
        "| Prueba | Estadístico | p-valor | Decisión con α = 0.05 |",
        "|---|---:|---:|---|",
    ])
    for fila in pruebas:
        lineas.append(
            f"| {fila['prueba']} | {fila['estadistico']:.6f} | "
            f"{fila['p_valor']:.6f} | {conclusion_prueba(fila['p_valor'])} |"
        )
    lineas.extend([
        "",
        "## Convergencia de P(éxito)",
        "",
        "Para estudiar únicamente el efecto del generador normal, se ejecutó el mismo modelo dos veces: una usando Polar para las variables normales V y C, y otra usando Aceptación–Rechazo. En ambos escenarios se conservaron los parámetros de `simulacion.py`; los intervalos mostrados son intervalos de Wilson del 95 %.",
        "",
        "| Método | n final | Éxitos | P(éxito) | IC 95 % | Error estándar |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for metodo in ("Polar", "Aceptacion-Rechazo"):
        fila = final[metodo]
        lineas.append(
            f"| {metodo} | {fila['n']} | {fila['exitos']} | "
            f"{fila['probabilidad_exito']:.6f} | "
            f"[{fila['ic95_inferior']:.6f}, {fila['ic95_superior']:.6f}] | "
            f"{fila['error_estandar']:.6f} |"
        )
    lineas.extend([
        "",
        "La estimación se estabiliza conforme aumenta n y el intervalo de confianza se estrecha, que es el comportamiento esperado del error Monte Carlo, proporcional a $1/\\sqrt{n}$.",
        "",
        "## Conclusiones",
        "",
        f"El método con menor tiempo mediano fue **{rapido['metodo']}** y el menor error conjunto en esta ejecución correspondió a **{preciso['metodo']}**. Ambos resultados dependen de la muestra concreta, por lo que la decisión principal debe apoyarse conjuntamente en ajuste estadístico, tiempo y simplicidad del algoritmo.",
        "",
        f"Polar obtuvo media {por_metodo['Polar']['media']:.6f} y S = {por_metodo['Polar']['desviacion']:.6f}; Aceptación–Rechazo obtuvo media {por_metodo['Aceptacion-Rechazo']['media']:.6f} y S = {por_metodo['Aceptacion-Rechazo']['desviacion']:.6f}. Las pruebas anteriores permiten verificar si las diferencias observadas son estadísticamente distinguibles al nivel de 5 %.",
        "",
        "## Revisión de consistencia del proyecto",
        "",
        "La integración principal de V, C, W y N es coherente: V y C son normales, W es exponencial, N se genera mediante un proceso de Poisson, y las cuatro variables alimentan el índice R. La comparación usa directamente las funciones entregadas por las personas 1 y 2 y los parámetros del simulador de la persona 4.",
        "",
        "Existe una diferencia de parámetros que el grupo debe resolver antes de entregar: `transformada_inversa_viento.ipynb` documenta media del viento de 12 m/s y tiempo medio entre anomalías de 300 s, mientras `simulacion.py` usa media del viento de 8 m/s y un proceso de anomalías con λ = 0.5 durante una misión de una unidad. Este análisis adopta `simulacion.py` como fuente de parámetros para no mezclar modelos dentro del mismo experimento.",
        "",
        "También debe mantenerse la misma unidad temporal para λ y `TIEMPO_MISION`. Si la misión se expresa en segundos, λ debe expresarse en anomalías por segundo; si se usa otra unidad, ambos valores deben convertirse conjuntamente.",
        "",
        "Con los pesos, penalizaciones y umbral actuales, la probabilidad estimada de éxito queda alrededor de 99.5 %. El cálculo es coherente con el código, pero el grupo debe confirmar si desea un escenario tan favorable. Si no fue intencional, deben calibrarse el umbral o las penalizaciones antes de interpretar el modelo; no conviene cambiarlos después de observar los resultados sin justificar el nuevo criterio.",
        "",
        "## Gráficas",
        "",
        "- `graficas/comparacion_distribuciones.png`",
        "- `graficas/comparacion_tiempos.png`",
        "- `graficas/convergencia_probabilidad_exito.png`",
    ])
    with open(ruta, "w", encoding="utf-8") as archivo:
        archivo.write("\n".join(lineas) + "\n")


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    ruta_resultados = os.path.join(base, "resultados")
    ruta_graficas = os.path.join(base, "graficas")
    os.makedirs(ruta_resultados, exist_ok=True)
    os.makedirs(ruta_graficas, exist_ok=True)

    muestras, resumen, pruebas = comparar_generadores()
    convergencia = analizar_convergencia()

    guardar_csv(resumen, os.path.join(ruta_resultados, "comparacion_metodos.csv"))
    guardar_csv(pruebas, os.path.join(ruta_resultados, "pruebas_estadisticas.csv"))
    guardar_csv(convergencia, os.path.join(ruta_resultados, "convergencia_exito.csv"))
    grafica_distribuciones(
        muestras, os.path.join(ruta_graficas, "comparacion_distribuciones.png")
    )
    grafica_tiempos(
        resumen, os.path.join(ruta_graficas, "comparacion_tiempos.png")
    )
    grafica_convergencia(
        convergencia,
        os.path.join(ruta_graficas, "convergencia_probabilidad_exito.png"),
    )
    crear_informe(
        resumen, pruebas, convergencia, os.path.join(base, "informe_persona_5.md")
    )

    print("Análisis de la persona 5 completado")
    for fila in resumen:
        print(
            f"{fila['metodo']}: media={fila['media']:.5f}, "
            f"S={fila['desviacion']:.5f}, "
            f"error={fila['error_conjunto']:.5f}, "
            f"tiempo={fila['tiempo_segundos']:.5f} s"
        )


if __name__ == "__main__":
    main()
