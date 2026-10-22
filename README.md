# qus-bfm: el rol del beamformer en la estimación de atenuación (QUS)

Proyecto de Astrid Mejía, Laboratorio de Imágenes Médicas (LIM), PUCP. Asesores: Christian Soto y Arelí Sánchez.

**Pregunta:** ¿la elección del beamformer afecta la estimación del coeficiente de atenuación (ACS) con el método de referencia y SWIFT?

El estimador de referencia es el de Sebastián Merino ([SWIFT_UFFC2024](https://github.com/SebastianMerino/SWIFT_UFFC2024), `phantoms.m`). Todo resultado propio se valida primero contra su reproducción.

## Estructura

```
config.m        DATA_ROOT y MERINO_DIR (único lugar con rutas absolutas)
qus_paths.m     path para los scripts de scripts/legacy
src/            funciones propias: load_vsx_rf, das_focused, das_planewave,
                merino_qus (estimador de Merino como función), get_spectrum_opt,
                measure_inclusion, fit_inclusion_circle, plot_qus_case
merino/         reproducción y adaptaciones de phantoms.m (calibración, data antigua,
                DataTargets2809) y sus resúmenes
scripts/        reconstrucción, diagnósticos, pruebas y figuras; legacy/ = scripts
                anteriores a la tarea 04 (qus_*, fig_*)
lab/            código del laboratorio sin modificar (SWIFT_Simple, SetUp Verasonics)
tasks/          briefs de cada tarea
results/        salidas (ignorado por git, salvo los .md de reporte)
```

Se sigue el patrón de los repos de referencia: utilidades en una carpeta (`AttenuationUtils`/`UltrasoundUtils` en Merino, `utils/` en Esteban Aviles) y scripts por experimento.

## Configuración

1. Clonar este repo y, aparte, el de Merino:
   `git clone https://github.com/SebastianMerino/SWIFT_UFFC2024`
2. Editar `config.m`:
   ```matlab
   DATA_ROOT  = 'C:\ruta\a\los\datos';               % carpetas de adquisición
   MERINO_DIR = fullfile(DATA_ROOT, 'SWIFT_UFFC2024'); % o donde se clonó
   ```
3. MATLAB R2026a (se probó ahí) con Signal Processing e Image Processing Toolbox.
4. Comprobar el path: `which spectra -all` debe resolver al repo de Merino cuando se corre algo de `merino/`. `lab/SWIFT_Simple` tiene copias con el mismo nombre; los scripts nuevos no lo agregan al path.

## Datos (no se suben)

Los datos viven en `DATA_ROOT` y nunca entran al repo (`.gitignore` en lista blanca: solo código y documentos).

| Carpeta | Contenido |
|---|---|
| `targets_edu/` | CIRS 049 (ID316), Ultrasonix L14-5/38, ya beamformeado; `data/` T1–T8, `referencia/544/` |
| `DataTargets2809/` | CIRS 049, Verasonics L11-4v 7.6 MHz, focalizado; `ChannelData` y `BeamformedData` del lab; referencias 261 y 544 |
| `PLANE/`, `TARGETS/` | data antigua, CIRS 049 T1 T3 T5 T7 T8, onda plana y focalizado |
| `Referencia Plane/`, `Ref Targets ATS/` | Wisconsin, referencia de PW y de FOC (data antigua) |
| `ASTRID/` | L11-5v 10 MHz, sin referencia |

## Cómo correr cada experimento

Desde la raíz del repo, en MATLAB:

| Experimento | Comando |
|---|---|
| Reproducir Merino (tarea 01) | `run('merino/phantoms_repro.m')` |
| Prueba del estimador como función | `run('scripts/test_merino_qus.m')` |
| Diagnóstico del zoom (tarea 02) | `run('scripts/zoom_diag.m')` |
| Estimador propio frente a Merino (tarea 03 B) | `run('scripts/own_qus_edu.m')` |
| Reconstrucción con la c de cada phantom | `run('scripts/bf_own_dt2809.m')`, `run('scripts/bf_own_old.m')` |
| Estimador de Merino en DataTargets2809 | `setenv('QUS_BFSET','own'); setenv('QUS_SWIFT','[4 2 8]'); setenv('QUS_TAG','frozen_cut8'); run('merino/phantoms_dt2809.m')` (`QUS_BFSET='lab'` usa la data del lab; `QUS_SWIFT` vacío = parámetros de Merino) |
| Centros de inclusión | `run('merino/centres_dt2809.m')`, `run('merino/centres_edu.m')` |
| Constante entre referencias y profundidad (fase 1) | `run('merino/phase1_summary.m')` |
| Calibración de SWIFT (fase 2) | `run('merino/calib_edu.m')` |
| Data antigua PW contra FOC (fase 3) | `run('merino/old_pwfoc.m')` |
| Resumen DataTargets2809 congelado (fase 3) | `run('merino/phase3_dt_summary.m')` |
| Figuras | `run('scripts/make_figures.m')` |
| Scripts anteriores | `qus_paths; run('scripts/legacy/qus_new.m')` |

## Resultados

El reporte completo está en [`results/REPORTE_FINAL.md`](results/REPORTE_FINAL.md). Predicciones y resultados de cada corrida: [`results/PREDICCIONES.md`](results/PREDICCIONES.md); decisiones tomadas: [`results/SUPUESTOS.md`](results/SUPUESTOS.md).

En breve:
- El estimador de Merino se reproduce exactamente y el estimador propio lo iguala con sus parámetros. Las diferencias de implementación mueven el fondo menos de 0.02 dB/cm/MHz.
- El beamformer (lab o propio, PW o FOC) cambia la ACS menos que la dispersión de una sola estimación.
- La constante entre las referencias 261 y 544 baja de +0.16 a +0.11 al reconstruir cada una con su velocidad, pero no desaparece y depende de la profundidad.
- La calibración de SWIFT en `targets_edu` (muB 10^4, muC 10^2, umbral 8) no mejora a los parámetros de Merino en los targets de prueba.
