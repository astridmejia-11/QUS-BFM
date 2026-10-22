# Tarea 03: aislar en qué difiere el estimador propio del de Merino

## Contexto

Resultado de la tarea 01: con `targets_edu` y `phantoms.m`, el fondo del CIRS 049 (serial 211270316, el "ID316") sale con sesgo entre −2.7 % y +2.2 %. Con el pipeline propio, el mismo phantom daba alrededor de −17 %.

Consecuencias:
- **La hipótesis de envejecimiento del CIRS queda descartada.** Mismo phantom, otro pipeline, valor correcto. El sesgo está en el procesamiento propio o en la referencia usada.
- Queda abierta la **constante de 0.158** entre las referencias 261 y 544 en `DataTargets2809`.

Resultado de la tarea 02: el zoom era de visualización; la escala axial propia es correcta. No hace falta tocar el beamformer.

Esta tarea separa beamformer de estimador: todo se hace sobre data **ya beamformeada**.

## Parte A: el pipeline de Merino sobre `DataTargets2809`

Pregunta: ¿la discrepancia de 0.158 entre referencias aparece también con el estimador de Merino?

1. Copiar `merino/phantoms_repro.m` a `merino/phantoms_dt2809.m`.
2. Entrada: `DataTargets2809/Samples/BeamformedData` y `Reference/BeamformedData/261` y `/544`, la versión del laboratorio.
3. Cambios **necesarios y solo esos**, documentados en la cabecera:
   - rutas;
   - nombres de variables si difieren (`RF`, `x`, `z`, `fs`);
   - región `x_inf`, `x_sup`: esta data va de −1.85 a +1.85 cm, no de 0 a 3.8;
   - centro de la inclusión `c1x`, `c1z` medido sobre el modo B de cada target, no copiado;
   - ACS de la referencia: 0.54 para 261, 0.53 para 544;
   - procesar los ocho targets, no solo los tres últimos.
4. Correr cada target dos veces, una por referencia, y reportar la diferencia por target.

Predicción antes de correr, a escribir por Code: si la constante sigue en torno a 0.15, el problema está en la data o los certificados; si baja hacia 0.01, está en el estimador propio.

## Parte B: el estimador propio sobre `targets_edu`

Pregunta: ¿qué diferencia de implementación produce el sesgo?

1. Correr el estimador propio (el de `qus_new` / `qus_all`) sobre `targets_edu` **con todos los parámetros de Merino**: banda 2.5–7.5 MHz, ventana Hamming, `NFFT = 2^(nextpow2(nz/2)+2)`, muB 10^3.5, muC según target, segunda pasada con `w`, métrica de media en máscaras circulares.
2. Debe reproducir la tabla de la tarea 01 dentro de ±0.01. Si no la reproduce, detenerse y reportar: hay una diferencia de implementación no listada.
3. Una vez reproducida, cambiar **un solo factor a la vez** al valor propio y registrar el efecto en fondo y en cada inclusión:

| # | Factor | Merino | Propio |
|---|---|---|---|
| 1 | Banda | 2.5–7.5 MHz | 4–9 MHz |
| 2 | Ventana espectral | Hamming | Tukey 0.25 (`getSpectrum`) |
| 3 | NFFT | `nextpow2(nz/2)+2` | `nextpow2(nz/2)+1` |
| 4 | muB | 10^3.5 | 10^3 |
| 5 | muC | 10^2 / 10^1 | 10^1 |
| 6 | Segunda pasada | `w` | `wExt` |
| 7 | Métrica | media, círculos | mediana, rectángulo |

4. Tabla de salida: factor, fondo, T6, T7, T8, cambio respecto a la línea base de Merino.

## Criterio

El factor que mueve el fondo del orden de −0.1 dB/cm/MHz es el candidato principal al −17 % propio. Si ninguno lo hace solo, probar combinaciones de los dos más grandes.

## Qué no hacer

- No elegir la combinación final por cercanía al certificado. El objetivo es entender qué causa la diferencia, no minimizarla.
- No modificar los archivos originales de Merino ni del laboratorio.

## Entregable

Las dos tablas (Parte A y Parte B), la predicción escrita antes de cada corrida, y una conclusión de dos o tres líneas para Astrid en español.