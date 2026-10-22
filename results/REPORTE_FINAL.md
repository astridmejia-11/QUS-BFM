# Reporte final: el rol del beamformer en la estimación del ACS

Astrid Mejía, LIM, PUCP. Asesores: Christian Soto y Arelí Sánchez. Fecha: 22 de octubre de 2026.

Este reporte responde a los seis puntos planteados por los asesores. Todo número sale de un script del repositorio [QUS-BFM](https://github.com/astridmejia-11/QUS-BFM). Las predicciones se escribieron antes de cada corrida y están, con fecha y hora, en `results/PREDICCIONES.md`; las decisiones no fijadas de antemano están en `results/SUPUESTOS.md`. El ACS se expresa en dB/cm/MHz.

## Base de validación

Antes de responder a los asesores se fijó un resultado de referencia y se depuró el entorno.

- **Reproducción de Merino.** `phantoms.m` corre sin cambios de método sobre `targets_edu` (CIRS 049, referencia 544). Con SWIFT, el fondo queda entre 0.535 y 0.562 frente a 0.55; las inclusiones T6, T7 y T8 dan 0.874, 0.979 y 0.830 frente a 0.97, 0.95 y 0.95 (`merino/phantoms_repro.m`).
- **Estimador propio frente al de Merino.** Con los parámetros de Merino, el estimador propio reproduce su tabla con diferencia nula. Las siete diferencias de implementación (banda, ventana espectral, NFFT, muB, muC, segunda pasada y métrica), probadas de a una o todas juntas, mueven el fondo menos de 0.02. El sesgo de la data propia no se origina en el estimador (`scripts/own_qus_edu.m`).
- **"Zoom" del modo B propio.** No era un error del beamformer: el laboratorio y `das_focused` usan rejillas distintas, y al graficar por índice de pixel una imagen parece ampliada. En coordenadas físicas, la escala axial coincide (razón de diámetros 1.004). Quedan dos diferencias del lado del laboratorio: el origen axial está corrido 1.50 mm (no se restan el retardo del pulso ni el de la lente) y el eje lateral está comprimido 2.7 %, porque se convirtió con la frecuencia del pulso, 7.8125 MHz, en lugar de 7.6 MHz. Ninguna de las dos altera la distancia entre ventanas `L` (`scripts/zoom_diag.m`, figura `results/zoom/bmode_physical.png`).
- **Limpieza (fase 0).** El script que comparaba las dos referencias de `DataTargets2809` (`qus_new.m`) cargaba siempre las referencias reconstruidas a 1540 m/s, aunque ya existían las reconstruidas con la velocidad de cada phantom. Se corrigió junto con la banda activa (7 a 9 MHz en vez de 4 a 9) y el emparejamiento de preSets de `load_vsx_rf` (cuatro convenciones de nombre, con prueba).

## Punto 1. Variar los parámetros de SWIFT, en particular el umbral de ΔBSC

**Qué se hizo.**
- Malla de calibración sobre `targets_edu`: muB en 10^{2.5 a 4.5}, muC en 10^{0 a 3} y umbral en {6, 8, 10} (fase 2, punto 5).
- Barrido del umbral 10, 8 y 6 con los demás parámetros congelados, en la data antigua y en `DataTargets2809` (fase 3.4).

**Resultados.**

| Dataset | Rango de variación al cambiar el umbral de 10 a 6 |
|---|---|
| `targets_edu`, calibración (objetivo) | 0.0756 (8), 0.0762 (6), 0.0775 (10) |
| Data antigua, PW y FOC | ≤ 0.025 en inclusiones, ≤ 0.015 en el fondo |
| `DataTargets2809`, propio y lab | ≤ 0.012 en inclusiones, ≤ 0.008 en el fondo |

En la malla, muC es lo que más mueve el resultado: con muB = 10^4, pasar de 10^1.5 a 10^2 cambia el objetivo de 0.088 a 0.076, más que cualquier cambio de umbral.

**Figura principal:** `results/figs/fase2/edu_T7_calibracion_frozen.png` (mapas de ACS y de ΔBSC con el umbral congelado).

**Predicción frente a resultado.** Se predijo que el umbral óptimo sería 8 o 10 y que su efecto sería secundario (< 0.02 en el fondo). Se cumplió en ambos casos.

**Conclusión.** El umbral de ΔBSC afecta poco el ACS en estos phantoms. Bajarlo no corrige el sesgo de las inclusiones, porque ese sesgo no proviene de bloques contaminados por retrodispersión. El parámetro relevante es la regularización de la retrodispersión (muC).

## Punto 2. Procesar el bloque nuevo con las configuraciones de la presentación y con las nuevas

**Qué se hizo.** Se procesó `DataTargets2809` (8 targets, referencias 261 y 544) con el estimador de Merino en tres configuraciones:
1. data del laboratorio con los parámetros originales de Merino (configuración de la tarea 03);
2. reconstrucción propia con la velocidad de cada phantom y los parámetros de Merino (fase 1);
3. reconstrucción propia con los parámetros calibrados (fase 3.2).

**Resultados** (SWIFT, media de 8 targets):

| Configuración | Fondo, ref. 544 | Fondo, ref. 261 | Constante 261 − 544, fondo |
|---|---|---|---|
| Lab, parámetros de Merino | 0.497 | 0.596 | +0.099 |
| Propio con su c, parámetros de Merino | 0.531 | 0.625 | +0.093 |
| Propio con su c, parámetros calibrados | **0.540** | 0.616 | +0.076 |

Con la tercera configuración, el sesgo medio del fondo con la 544 es −1.8 %. Las inclusiones siguen con sesgo: con la 544, el error relativo medio es 0.150; T6 y T8 dan 0.78 frente a 0.97 y 0.95, y T4 da 0.87 frente a 0.81.

**Figura principal:** `results/figs/fase3/dt2809_own_frozen_cut8_T7_F_ref544.png`.

**Predicción frente a resultado.** Se predijo una constante cercana a +0.09 y un fondo con la 544 cerca de 0.53. Se cumplió (+0.076 y 0.540).

**Conclusión.** Con cada phantom reconstruido a su velocidad y los parámetros calibrados, el fondo del CIRS 049 en el bloque nuevo coincide con el certificado usando la referencia 544. La referencia 261 sigue dando valores más altos, en unos 0.08.

## Punto 3. Comparar el beamformer del laboratorio con el propio en SWIFT

**Qué se hizo.** Se usó el mismo estimador con los mismos parámetros calibrados; cada beamformer usa referencias reconstruidas por sí mismo (fase 3.3). Hay una diferencia de procedimiento: el laboratorio reconstruyó todas sus referencias a 1540 m/s, mientras que las propias usan 1509 para la 261 y 1539 para la 544.

**Resultados** (laboratorio menos propio, media ± desviación sobre 8 targets):

| | Inclusión, 261 | Inclusión, 544 | Fondo, 261 | Fondo, 544 |
|---|---|---|---|---|
| Lab − propio | +0.012 ± 0.025 | +0.003 ± 0.021 | −0.111 ± 0.027 | −0.020 ± 0.026 |
| Desviación individual típica | 0.086 | 0.074 | 0.095 | 0.073 |

Con la 544 (velocidad cercana a 1540 en ambos procedimientos), los dos beamformers difieren menos que la dispersión de una sola estimación. Este resultado es coherente con el anterior (0.026 contra 0.092 de dispersión).

**Figura principal:** `results/figs/fase3/dt2809_lab_frozen_cut8_T7_F_ref544.png`, frente a la propia citada en el punto 2.

**Predicción frente a resultado.**
- Con la 544 se predijo una diferencia < 0.03, y se cumplió.
- Con la 261 se predijo que el laboratorio daría +0.05, y no se cumplió: en la inclusión la diferencia es nula y en el fondo el laboratorio da 0.111 menos.

**Conclusión.** La implementación del beamformer no cambia el ACS más allá de la dispersión cuando la referencia se reconstruye de forma equivalente. Lo que sí lo cambia es reconstruir una referencia con una velocidad que no es la suya: con la 261 a 1540, el fondo se desplaza 0.11.

## Punto 4. Cada referencia con su propia velocidad del sonido y beamformeado homogéneo

**Qué se hizo.**
- Reconstrucción propia desde `ChannelData`: la 261 a 1509, la 544 a 1539 y las muestras a 1542, todas en la misma rejilla axial. Luego se aplicó el estimador de Merino con cada referencia (fase 1).
- Prueba de profundidad: ACS del fondo por franjas de z.
- Comparación campo por campo de los 15 preSets.

**Resultados.** Constante 261 − 544 en el fondo, con TV (el método menos regularizado), media ± desviación sobre 8 targets:

| Reconstrucción | Fondo | Inclusión | z 0.5–1.5 cm | z 1.5–2.5 cm | z 2.5–3.5 cm |
|---|---|---|---|---|---|
| Laboratorio, todo a 1540 | +0.161 ± 0.007 | +0.243 ± 0.013 | +0.111 | +0.332 | +0.275 |
| Propia, cada una con su c | **+0.109 ± 0.002** | **+0.118 ± 0.003** | **+0.034** | **+0.169** | **+0.164** |

Los preSets de muestras y referencias son idénticos en frecuencia, muestreo, profundidades, foco de transmisión (98 λ), apertura activa (23 elementos), voltaje (20 V), pulso y TGC.

**Figura principal:** `results/figs/fase1/dt2809_own_merino_T7_F_ref261.png` y `..._ref544.png`.

**Predicción frente a resultado.**
- Se predijo que la constante bajaría a un valor entre 0.08 y 0.12 sin desaparecer, y que dependería de la profundidad. Ambas cosas se cumplieron.
- Se predijo un máximo cerca del foco (2 cm). No se cumplió: la constante es pequeña por encima de 1.5 cm y casi constante por debajo.

**Conclusión.**
- Reconstruir cada referencia con su velocidad reduce la constante en un tercio y la iguala entre inclusión y fondo. La velocidad del sonido es, por lo tanto, una causa parcial de la discrepancia.
- La parte restante no proviene del montaje (los preSets son iguales) y cambia con la profundidad, lo que apunta a un efecto de haz en transmisión. El foco se fijó con retardos calculados para 1540 m/s y no puede corregirse en el beamformer.

## Punto 5. Parámetros de SWIFT calibrados con `targets_edu` y su uso en la data antigua (PW contra FOC)

**Qué se hizo.**
- Calibración con partición fijada antes de mirar resultados: T1, T3, T5 y T7 para calibrar; T2, T4, T6 y T8 para probar.
- Objetivo fijado de antemano: media del valor absoluto del sesgo relativo en inclusión y fondo.
- El control de la calibración reproduce exactamente la tabla de Merino.
- La combinación congelada se evaluó una sola vez en los targets de prueba. Luego se aplicó sin cambios a la data antigua (CIRS 049, T1, T3, T5, T7 y T8), con referencia Wisconsin y cada phantom reconstruido a su velocidad.

**Calibración** (congelado: muB = 10^4, muC = 10^2, umbral 8):

| Conjunto | Congelado | Parámetros de Merino |
|---|---|---|
| Calibración | 0.076 | 0.097 |
| Prueba | **0.073** | **0.071** |

**Data antigua.** La banda 2.5 a 7.5 MHz cae fuera del espectro útil de este montaje, por lo que se usó la banda a −20 dB de la referencia: 3.5 a 9.9 MHz en PW y 3.6 a 10.0 MHz en FOC.

| Target | Cert. | PW inclusión | FOC inclusión | PW − FOC | PW fondo | FOC fondo |
|---|---|---|---|---|---|---|
| T1 | 0.52 | 0.641 ± 0.126 | 0.682 ± 0.128 | −0.040 | 0.470 | 0.437 |
| T3 | 0.74 | 0.484 ± 0.126 | 0.404 ± 0.123 | +0.080 | 0.445 | 0.425 |
| T5 | 0.75 | 0.513 ± 0.128 | 0.551 ± 0.087 | −0.038 | 0.443 | 0.416 |
| T7 | 0.95 | 0.747 ± 0.135 | 0.831 ± 0.112 | −0.084 | 0.464 | 0.421 |
| T8 | 0.95 | 0.593 ± 0.120 | 0.614 ± 0.097 | −0.021 | 0.434 | 0.400 |

**Figura principal:** `results/figs/fase3/old_PW_T7_cut8.png` y `results/figs/fase3/old_FOC_T7_cut8.png`.

**Predicción frente a resultado.**
- Calibración: se predijo un objetivo entre 0.05 y 0.08 y una prueba entre 0.07 y 0.11, sin gran diferencia con Merino. Se cumplió; no hay sobreajuste, pero tampoco mejora en los targets de prueba.
- Data antigua: se predijo |PW − FOC| menor que la dispersión individual, y se cumplió en los cinco targets. No se previó que ambos esquemas subestimaran de forma marcada: fondo entre 0.40 y 0.47 frente a 0.55, y T3 entre 0.40 y 0.48 frente a 0.74.

**Conclusión.**
- Los parámetros calibrados rinden igual que los de Merino en datos que no se usaron para elegirlos, así que los de Merino ya eran adecuados para este phantom.
- En la data antigua, el esquema de transmisión no cambia el ACS más allá de la dispersión. Ambos esquemas comparten un sesgo negativo que no se explica por el estimador, los parámetros ni la velocidad de reconstrucción, y que sí desaparece en el bloque nuevo con la referencia 544. El candidato es la referencia Wisconsin o la adquisición antigua.

## Punto 6. Repositorio de GitHub

**Qué se hizo.**
- Repositorio [astridmejia-11/QUS-BFM](https://github.com/astridmejia-11/QUS-BFM), ordenado como los de Merino y Esteban Aviles: utilidades en `src/`, adaptaciones de Merino en `merino/`, scripts de experimentos en `scripts/`, código del laboratorio sin modificar en `lab/`, briefs en `tasks/` y resultados ignorados en `results/`.
- `README.md` en español con la configuración, el origen de los datos y cómo correr cada experimento.
- Un commit y un push al cerrar cada fase.
- El `.gitignore` es una lista blanca: el historial no contiene ningún archivo de datos ni imágenes, y el objeto más grande pesa 18 KB.
- Tras mover los archivos se comprobó que `phantoms_repro.m` y la corrida de la data antigua reproducen sus tablas con diferencia nula.

**Conclusión.** El repositorio permite reproducir cada número de este reporte a partir de los datos crudos y de `config.m`.

## Bloque `ASTRID` (fase 4)

`targets_edu` proviene de un Ultrasonix SonixTouch con sonda L14-5/38 y beamformado interno; `ASTRID` es Verasonics con L11-5v, datos de canal y foco a 30 mm. Solo coinciden la frecuencia central y el muestreo; el TGC y el foco de `targets_edu` no están registrados. Por la regla de usar el mismo montaje en muestra y referencia, **`ASTRID` no se procesó** con la referencia de `targets_edu`. Se predijo que los montajes no coincidirían, y así fue.

## Qué cambia respecto a lo ya presentado

| Afirmación anterior | Estado | Fundamento |
|---|---|---|
| El CIRS 049 podría haber envejecido y su ACS ya no ser la certificada | **Se retira** | El mismo phantom da el fondo correcto con el código de Merino en `targets_edu` (sesgo de −2.7 % a +2.2 %) y en el bloque nuevo con la referencia 544 (−1.8 %). |
| La velocidad del sonido estaba descartada como causa de la constante entre referencias | **Se corrige** | La corrida que la descartaba nunca usó las referencias reconstruidas con su velocidad (error de caché en `qus_new.m`). Con cada referencia a su c, la constante TV baja de +0.161 a +0.109. La velocidad es una causa parcial. |
| El umbral del detector se mantiene en 10 | **Se confirma con matiz** | El umbral mueve el ACS menos de 0.025 en los tres datasets. La calibración eligió 8, pero la diferencia con 10 es irrelevante. Se confirma que bajarlo no corrige la contaminación. |
| El beamformer propio difiere poco del laboratorio (0.026 frente a 0.092 de dispersión) | **Se confirma** | Con los parámetros calibrados y la referencia 544: +0.003 en la inclusión y −0.020 en el fondo, frente a una dispersión de 0.07. |
| El modo B propio tiene un "zoom" por error del beamformer | **Se retira** | Era un efecto de graficar por índice de pixel. Las diferencias reales de ejes son del laboratorio (origen 1.50 mm, eje lateral 2.7 %). |
| El sesgo de −17 % en el fondo es del pipeline propio | **Se corrige** | No proviene del estimador ni del beamformer. En el bloque nuevo desaparece al reconstruir con la c correcta y calibrar SWIFT con la referencia 544. En la data antigua persiste con el estimador de Merino. |
| PW y FOC dan el mismo ACS | **Se confirma** | \|PW − FOC\| es menor que la dispersión individual en los cinco targets, con los parámetros calibrados y cada phantom a su c. |

## Preguntas abiertas

1. **¿Por qué la referencia 261 da unos 0.08 a 0.11 más que la 544 aun reconstruida con su velocidad, y solo por debajo de 1.5 cm?**
   - Experimento: adquirir la 261 y la 544 con onda plana, o con retardos de transmisión calculados a 1509 m/s, y repetir la prueba de profundidad. Si la constante desaparece, el origen es el foco de transmisión.
   - Además, medir el ACS de la 261 por sustitución en banco, para verificar el certificado.
2. **¿De dónde viene el sesgo negativo de la data antigua (fondo entre −15 % y −27 %) con la referencia Wisconsin?**
   - Experimento: adquirir la referencia 544 con el montaje antiguo (L11-4v a 7.24 MHz) y repetir la comparación PW frente a FOC.
   - Comparar además los preSets de muestras y de Wisconsin (TGC, voltaje) como se hizo con el bloque nuevo.
3. **¿Por qué se subestiman las inclusiones de ACS alta (T6 y T8: 0.78 a 0.82 frente a 0.95 a 0.97) en todos los datasets, incluido `targets_edu`?**
   - Experimento: barrer el tamaño de bloque y muB sobre una simulación con inclusiones de ACS conocido, para separar la pérdida de resolución por regularización de un sesgo físico.
4. **¿Por qué el fondo del laboratorio con la 261 queda 0.11 por debajo del propio, mientras que con la 544 coinciden?**
   - Experimento: reconstruir la 261 propia a 1540 m/s (como hizo el laboratorio). Si la diferencia desaparece, queda atribuida por completo a la velocidad de reconstrucción de la referencia.
5. **¿Puede estimarse el ACS en `ASTRID`?**
   - Experimento: adquirir una referencia homogénea (544) con el mismo montaje de `ASTRID` (L11-5v, foco 30 mm, mismo TGC y voltaje).
6. **¿El eje lateral del laboratorio usa la frecuencia del pulso por diseño?**
   - Experimento: confirmar con el laboratorio el código que genera `x`. Si se corrige, repetir la comparación de beamformers con ejes iguales.
