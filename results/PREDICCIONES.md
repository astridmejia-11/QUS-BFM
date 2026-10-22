# Predicciones y resultados (tarea 04)

Cada entrada se escribe antes de correr; el resultado se agrega después, al lado.

## Fase 0

### 2026-10-21 · prueba de `load_vsx_rf` con `PLANE/RFdata_26_T8_..._15-53-50.mat`
- **Predicción:** antes del cambio no habría emparejado `...preSave.mat` (solo buscaba `preSet` y `_preSet`) y habría caído al respaldo con el warning `assumedPreSet`. Con el cambio, usa su propio `preSave` sin warning.
- **Resultado:** usa `RFdata_26_T8_-August-2026_15-53-50preSave.mat`, sin warning. PASS (`scripts/test_load_presave.m`). La variable dentro del `preSave` se llama `preSet`, así que no hubo que tocar la lectura.

## Fase 1: referencias con su propia velocidad

### 2026-10-21 18:40 · estimador de Merino sobre reconstrucciones propias (261 a 1509, 544 a 1539, muestras a 1542)
- **Predicción:** la constante 261 − 544 baja, pero **no desaparece**. Con TV, en el fondo, pasa de +0.16 a algo entre +0.08 y +0.12.
  - Reconstruir la 261 a 1509 corrige la escala axial (~2 %, poco efecto en L) y el foco en **recepción**.
  - El foco en **transmisión** quedó fijado en la adquisición con retardos calculados para 1540; en un medio de 1509 el haz converge a otra profundidad, y eso no se corrige en el beamformer.
  - Por eso espero que la constante **dependa de z**, con la mayor diferencia cerca del foco (≈2 cm) y no constante entre franjas.
  - El sesgo negativo del fondo con la 544 (0.42–0.57) casi no debería cambiar: la 544 (1539) y las muestras (1542) ya estaban cerca de 1540.

## Fase 2: calibración de SWIFT con `targets_edu`

### 2026-10-21 19:01 · partición y objetivo, fijados antes de mirar resultados
- **Calibración:** T1, T3, T5, T7. **Prueba:** T2, T4, T6, T8. Referencia 544 (4 archivos, frame 1).
- **Malla:** muB ∈ 10^{2.5, 3, 3.5, 4, 4.5}; muC ∈ 10^{0, 0.5, 1, 1.5, 2, 2.5, 3}; ratioCutOff ∈ {6, 8, 10}. Todo lo demás como `phantoms.m`.
- **Objetivo:** media del |sesgo relativo| en inclusión y fondo de los 4 targets de calibración (8 valores). Desempate por menor desviación estándar media.
- **Centros:** medidos sobre el modo B de cada target (`merino/centres_edu.m`), no los de Merino.
- **Control:** con los parámetros y centros de Merino, el pipeline debe reproducir exactamente la tabla de la tarea 01.
- **Predicción:**
  - el objetivo de calibración queda entre 0.05 y 0.08;
  - el error de prueba sale algo mayor (0.07–0.11), porque el óptimo se ajusta a 4 targets;
  - los parámetros originales de Merino, sobre los mismos targets de prueba, quedan cerca del congelado (diferencia < 0.02), porque Merino ya los eligió sobre este phantom;
  - el umbral óptimo es 8 o 10; con 6 se penalizan bordes reales y empeora la inclusión.
- **Resultado** (`merino/phase1_summary.m`; constante 261 − 544, media ± desviación sobre 8 targets):

  | Juego | TV fondo | TV inclusión | SWIFT fondo | SWIFT inclusión |
  |---|---|---|---|---|
  | lab (todo a 1540) | +0.161 ± 0.007 | +0.243 ± 0.013 | +0.099 ± 0.049 | +0.105 ± 0.018 |
  | propio (cada uno con su c) | **+0.109 ± 0.002** | **+0.118 ± 0.003** | +0.093 ± 0.010 | +0.118 ± 0.030 |

  Prueba de profundidad, TV, fondo, juego propio: 0.5–1.5 cm **+0.034**; 1.5–2.5 cm **+0.169**; 2.5–3.5 cm **+0.164**. En el juego lab: +0.111, +0.332 y +0.275.
  - **Se cumple:** la constante baja (de 0.161 a 0.109, dentro del rango predicho) sin desaparecer, y depende de z.
  - **No se cumple del todo:** no tiene un máximo cerca del foco. Es chica arriba de 1.5 cm y casi constante debajo.
  - **Sesgo del fondo con la 544:** pasa de 0.471 a 0.492 con TV y de 0.497 a 0.531 con SWIFT. Mejora más de lo previsto.
  - **PreSets** (`scripts/compare_presets_dt2809.m`): los 15 son idénticos en fc, fs, profundidades, foco de transmisión (98 λ), elementos activos (23), voltaje (20 V), pulso y TGC. La diferencia no viene del montaje.
- **Resultado** (`merino/calib_edu.m`; control contra la tarea 01: exacto):
  - **Congelado:** muB = 10^4, muC = 10^2, ratioCutOff = 8. Objetivo de calibración 0.0756 (desempate innecesario; el siguiente es 0.0762 con umbral 6).
  - **Error medio |sesgo relativo|:**

    | Conjunto | Congelado | Merino original |
    |---|---|---|
    | Calibración (T1, T3, T5, T7) | 0.076 | 0.097 |
    | Prueba (T2, T4, T6, T8) | **0.073** | **0.071** |

  - **Se cumple:** objetivo en el rango predicho (0.05–0.08); prueba en el rango (0.07–0.11) y sin sobreajuste; diferencia con Merino en prueba < 0.02; umbral óptimo 8.
  - **Lectura:** la calibración mejora solo los targets con que se calibró. En los de prueba, los parámetros de Merino rinden igual. Las inclusiones siguen con sesgo (T6 0.82 y T8 0.81 frente a 0.97 y 0.95); el fondo queda entre 0.50 y 0.59.

## Fase 3: parámetros congelados

### 2026-10-21 19:32 · predicciones antes de correr
- **3.1, data antigua PW contra FOC:** la diferencia PW − FOC en cada inclusión queda por debajo de la desviación individual (~0.08) en al menos 4 de 5 targets. Ningún esquema es sistemáticamente mejor.
- **3.2, `DataTargets2809` con el beamformer propio:** la constante 261 − 544 con SWIFT sigue en torno a +0.09; el fondo con la 544 queda cerca de 0.53.
- **3.3, lab contra propio:** con la 544 difieren poco (< 0.03), porque la 544 está cerca de 1540. Con la 261, el lab da valores más altos que el propio (≈ +0.05), por la referencia reconstruida a 1540 en vez de 1509.
- **3.4, umbral 10/8/6:** cambiar el umbral mueve las medias menos de 0.02 en el fondo; en las inclusiones hasta ≈ 0.05 en los targets con ΔBSC grande (T1, T4, T7, T8). Se confirma que el umbral es secundario.
- **Resultado 3.1** (`merino/old_pwfoc.m`, umbral congelado 8; la banda de la fase 2 no cabe en el espectro de este montaje: −20 dB = 3.5–9.9 MHz en PW y 3.6–10.0 en FOC, que fue la que se usó):

  | Target | Cert. | PW inclusión | FOC inclusión | PW − FOC | PW fondo | FOC fondo |
  |---|---|---|---|---|---|---|
  | T1 | 0.52 | 0.641 ± 0.126 | 0.682 ± 0.128 | −0.040 | 0.470 | 0.437 |
  | T3 | 0.74 | 0.484 ± 0.126 | 0.404 ± 0.123 | +0.080 | 0.445 | 0.425 |
  | T5 | 0.75 | 0.513 ± 0.128 | 0.551 ± 0.087 | −0.038 | 0.443 | 0.416 |
  | T7 | 0.95 | 0.747 ± 0.135 | 0.831 ± 0.112 | −0.084 | 0.464 | 0.421 |
  | T8 | 0.95 | 0.593 ± 0.120 | 0.614 ± 0.097 | −0.021 | 0.434 | 0.400 |

  - **Se cumple:** |PW − FOC| < desviación individual en **los 5** targets; ningún esquema es sistemáticamente mejor en la inclusión. En el fondo, PW queda siempre un poco más alto (+0.02 a +0.04).
  - **No previsto:** los dos esquemas subestiman mucho. El fondo da 0.40–0.47 contra 0.55 (−15 % a −27 %), e inclusiones como T3 dan 0.40–0.48 contra 0.74. El sesgo de la data antigua sigue con el estimador de Merino, los parámetros calibrados y cada phantom a su c.
  - **3.4 en data antigua:** los umbrales 10, 8 y 6 mueven las medias ≤ 0.025 en la inclusión y ≤ 0.015 en el fondo (PW − FOC estable). El umbral es secundario.

## Fase 4: bloque `ASTRID`

### 2026-10-22 · comparación de montajes (`scripts/inspect_astrid_preset.m` y metadata XML de `targets_edu`)
- **Predicción:** no coinciden. `targets_edu` no trae preSet de Verasonics, así que probablemente sea de otro sistema.
- **Resultado:**

  | | `targets_edu` | `ASTRID` |
  |---|---|---|
  | Sistema | Ultrasonix SonixTouch (SW 6.1.1, `study.xml`) | Verasonics Vantage |
  | Sonda | L14-5/38, 128 elementos, pitch 0.305 mm | L11-5v, 128 elementos, pitch 0.300 mm |
  | Frecuencia | fcentre 9.5 MHz (`fc` en el `.mat`: 10 MHz), fs 40 MHz | `Trans.frequency` 10 MHz, pulso 10 MHz, fs 40 MHz |
  | Formación de haz | interna del Sonix, entrega RF ya beamformeado | datos de canal, foco 30 mm (194.8 λ), 33 elementos activos |
  | TGC, foco, voltaje | no registrados en la metadata exportada | TGC conocido, 20 V |

  - **Se cumple la predicción.** Solo coinciden fc y fs; difieren el sistema, la sonda, el beamformer y la cadena de recepción. El TGC y el foco de `targets_edu` no se conocen.
  - Por la regla 2 de CLAUDE.md (muestra y referencia con el mismo montaje), **`ASTRID` no se procesa** con la referencia 544 de `targets_edu`. Hace falta una referencia adquirida con el mismo montaje de `ASTRID`.

### 2026-10-22 · resultados 3.2–3.4 (`merino/phase3_dt_summary.m`, SWIFT muB 10^4, muC 10^2)
- **3.2, beamformer propio (cada phantom con su c), umbral 8**, media sobre 8 targets:

  | Referencia | Inclusión | Fondo | |sesgo rel.| inclusión | |sesgo rel.| fondo |
  |---|---|---|---|---|
  | 261 | 0.799 | 0.616 | 0.133 | 0.120 |
  | 544 | 0.719 | **0.540** | 0.150 | **0.044** |

  Constante 261 − 544: **+0.076** en el fondo y **+0.080** en la inclusión, casi igual en los 8 targets (fondo entre +0.054 y +0.090).
  - **Se cumple:** constante cercana a +0.09 y fondo con la 544 cerca de 0.53. Con la 544, el fondo queda a −1.8 % del certificado.
  - Las inclusiones siguen con sesgo: subestimadas las de ACS alta (T6 0.78, T8 0.78 frente a 0.97 y 0.95) y sobrestimada T4 (0.87 frente a 0.81).
- **3.3, lab menos propio, mismo estimador y parámetros:**

  | | Inclusión 261 | Inclusión 544 | Fondo 261 | Fondo 544 |
  |---|---|---|---|---|
  | Media ± desv. (8 targets) | +0.012 ± 0.025 | +0.003 ± 0.021 | **−0.111 ± 0.027** | −0.020 ± 0.026 |
  | Desviación individual típica | 0.086 | 0.074 | 0.095 | 0.073 |

  - **Se cumple** con la 544: lab y propio difieren mucho menos que la dispersión individual.
  - **No se cumple** con la 261: predije que el lab daría valores más altos (+0.05), y en la inclusión la diferencia es nula mientras que en el fondo el lab da **0.111 menos**. La diferencia de procedimiento (la 261 del lab reconstruida a 1540 y no a 1509) pesa en el fondo, que abarca más rango de profundidad que la inclusión, y no en la inclusión. Con el lab, la constante 261 − 544 cambia de signo entre regiones (inclusión +0.089, fondo −0.015).
- **3.4, umbral 10/8/6:** en `DataTargets2809`, los cambios son ≤ 0.012 en la inclusión y ≤ 0.008 en el fondo; en la data antigua, ≤ 0.025 y ≤ 0.015. **Se cumple:** el umbral es secundario. El 8 elegido en la fase 2 queda, con margen mínimo, como el mejor o casi el mejor en todos los casos.
