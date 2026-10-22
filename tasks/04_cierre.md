# Tarea 04: cerrar lo pedido por los asesores, de principio a fin

## Cómo trabajar esta tarea

- Trabaja **de forma autónoma** hasta terminar todas las fases. Astrid no quiere reportes intermedios.
- **No te detengas a preguntar** salvo ante una decisión que no se pueda deshacer (borrar datos, reescribir historia de git, sobrescribir resultados ya presentados) o si falta un archivo imprescindible. En cualquier otro caso, toma la opción más razonable, anótala en `results/SUPUESTOS.md` y sigue.
- Antes de cada corrida, escribe la predicción en `results/PREDICCIONES.md` con fecha y hora. Después, el resultado al lado.
- Al terminar cada fase, haz un commit con un mensaje claro. Al final, push.
- El único entregable para Astrid es `results/REPORTE_FINAL.md` (fase 7). Todo lo demás es soporte.
- Siguen vigentes todas las reglas de `CLAUDE.md`, con la precisión de la fase 2 sobre la regla 1.

## Lo que pidieron los asesores (fuente de esta tarea)

1. Variar los parámetros de SWIFT, sobre todo el umbral de ΔBSC, para que detecte los bordes.
2. Procesar el bloque nuevo con las configuraciones de la presentación y con las nuevas.
3. Comparar la data beamformeada del laboratorio con el beamformer propio en SWIFT.
4. Verificar el beamformeado para que cada referencia use su propia velocidad del sonido; homogeneizar el beamformeado y repetir los resultados con los parámetros corregidos.
5. Con `targets_edu` y el repositorio de Merino, obtener parámetros de SWIFT que den resultados cercanos al certificado, y usar esos mismos parámetros sobre la data antigua para la comparación PW contra FOC.
6. Un repositorio de GitHub con los avances, con estructura guiada por los de Merino o Esteban Aviles.

Los puntos 1 a 3 ya tienen resultados; se rehacen con lo que salga de las fases 1 y 2.

---

## Fase 0: limpieza (sin esto, lo demás no es confiable)

1. `qus_new.m`: el bucle principal debe cargar `bf_new_ref<tag>_c<cRef>.mat`. Restaurar la banda a la de la configuración vigente (no 7–9 MHz).
2. `qus_bfcomp.m`: corregir la rama de respaldo que usa `Dp.params.*`.
3. `load_vsx_rf.m`: comprobar que empareja las cuatro convenciones de preSet (`preSet`, `_preSet`, `preSave`, `_preSave`). Prueba: cargar `PLANE/RFdata_26_T8_-August-2026_15-53-50.mat` y confirmar que usa su propio `...preSave.mat` sin warning.
4. `CLAUDE.md` línea 17: dice que datos y repo están en lugares distintos; hoy todo está en `C:\Users\LENOVO\amejia`. Corregir.
5. Commit.

## Fase 1: referencias con su propia velocidad (punto 4)

Pregunta: ¿la constante 261 − 544 ≈ +0.16 y el sesgo del fondo en `DataTargets2809` vienen de haber reconstruido la 261 (1509 m/s) a 1540?

1. Con `das_focused`, desde `ChannelData`, reconstruir: 261 a 1509, 544 a 1539, muestras a la c del fondo del CIRS 049 (1542). Misma rejilla axial para todas (`dz` fijo), parámetros de beamformer de la configuración vigente.
2. Correr el estimador de Merino (`merino/phantoms_dt2809.m`, con entrada cambiada a estas reconstrucciones) con 261 y con 544 por separado.
3. Reportar por target: fondo e inclusión con cada referencia, y la constante 261 − 544 con TV y SWIFT.
4. **Prueba de profundidad:** estimar la ACS del fondo por franjas de z (por ejemplo 0.5–1.5, 1.5–2.5, 2.5–3.5 cm) con cada referencia. Si la constante cambia con z, el origen es un efecto de haz (foco, difracción), no un offset de certificado.
5. Si la constante no baja: comparar los preSets de 261, 544 y muestras campo por campo (TGC, foco, apertura, voltaje, fecha de adquisición) y dejarlo documentado.

Criterio: la regla 3 se aplica pase lo que pase. Los resultados de las fases siguientes usan estas reconstrucciones.

## Fase 2: calibrar SWIFT con `targets_edu` (punto 5)

Precisión sobre la regla 1: los asesores piden elegir parámetros por cercanía al certificado. Eso es legítimo **solo si la calibración y la evaluación usan datos distintos**. Por eso:

1. **Partición fijada antes de mirar resultados:** calibración = T1, T3, T5, T7 de `targets_edu`; prueba = T2, T4, T6, T8. Escribirla en `PREDICCIONES.md` antes de correr.
2. Estimador de Merino (`phantoms.m` adaptado a los 8 targets, centros de inclusión medidos sobre el modo B de cada uno). Referencia 544.
3. Malla: muB ∈ 10^{2.5, 3, 3.5, 4, 4.5}, muC ∈ 10^{0, 0.5, 1, 1.5, 2, 2.5, 3}, ratioCutOff ∈ {6, 8, 10}. Banda, ventana y demás como Merino.
4. Objetivo, fijado antes: media del |sesgo relativo| en inclusiones y fondo de los targets de calibración. Desempate por menor desviación estándar.
5. Congelar la mejor combinación. Evaluarla **una sola vez** en los targets de prueba y reportar ambos conjuntos. Si el error de prueba es mucho mayor que el de calibración, decirlo: hay sobreajuste.
6. Comparar contra los parámetros originales de Merino en los mismos targets de prueba.

## Fase 3: aplicar los parámetros congelados (puntos 2, 3 y 5)

Con los parámetros de la fase 2, sin cambiarlos:

1. **Data antigua, PW contra FOC:** CIRS 049, T1, T3, T5, T7, T8, referencia Wisconsin de cada esquema. Banda: si la de la fase 2 cae fuera del espectro útil de este montaje, usar la banda a −20 dB del espectro de la referencia y documentarlo en `SUPUESTOS.md`. Tabla PW, FOC, diferencia, desviaciones y certificado.
2. **`DataTargets2809` con el beamformer propio** (reconstrucciones de la fase 1), 8 targets, con 261 y con 544.
3. **Comparación de beamformers** sobre `DataTargets2809`: laboratorio contra propio, mismo estimador, cada uno con referencias reconstruidas por sí mismo. Notar que las referencias del laboratorio están a 1540 m/s; reportarlo como diferencia de procedimiento.
4. Barrido del umbral (10, 8, 6) con los demás parámetros congelados, en la data antigua y en `DataTargets2809`, para confirmar o corregir la conclusión previa.

## Fase 4: bloque `ASTRID`

`targets_edu` usa fc 10 MHz y fs 40 MHz, igual que `ASTRID` (L11-5v). Comparar sus preSets (transductor, foco, TGC, apertura). Si coinciden en lo que afecta al cociente espectral, procesar `ASTRID` con la referencia 544 de `targets_edu` y los parámetros congelados. Si no coinciden, documentar qué difiere y no procesar.

## Fase 5: figuras

Para cada resultado de las fases 1 a 4: modo B con región de análisis y de medida, mapa de ACS, mapa de ΔBSC. Coordenadas físicas en todos los ejes. Guardar en `results/figs/<fase>/` con nombres descriptivos.

## Fase 6: repositorio de GitHub (punto 6)

1. Revisar la estructura de `SebastianMerino/SWIFT_UFFC2024` y de los repos de `EstebanAviles` en GitHub, y ordenar el repo de forma parecida: `src/`, `scripts/`, `merino/`, `lab/`, `tasks/`, `results/` (ignorado).
2. Mover archivos solo si después todo sigue corriendo; probar al menos `phantoms_repro.m` y una corrida de la fase 3 tras mover.
3. `README.md` en español: objetivo, cómo configurar `config.m`, cómo correr cada experimento, de dónde salen los datos (sin subirlos), y un resumen de resultados con enlace al reporte.
4. `.gitignore` con `*.mat`, imágenes, `.xlsx` y `results/`. Verificar que ningún archivo de datos quede en el historial.
5. Commit y push.

## Fase 7: reporte final

`results/REPORTE_FINAL.md`, en español, prosa académica sin guiones largos, conciso. Una sección por cada punto de los asesores, con:
- qué se hizo, tabla de resultados y figura principal (ruta);
- predicción frente a resultado;
- conclusión en dos o tres líneas.

Además, una sección **"Qué cambia respecto a lo ya presentado"**, con cada afirmación anterior que se confirma, se corrige o se retira. Como mínimo:
- la hipótesis de envejecimiento del CIRS (la tarea 03 la debilitó);
- la afirmación de que la velocidad del sonido estaba descartada como causa de la constante (la corrida nunca usó las referencias con su c);
- la elección del umbral en 10.

Y una sección final **"Preguntas abiertas"**, cada una con el experimento que la resolvería.

Último commit y push con el reporte.
