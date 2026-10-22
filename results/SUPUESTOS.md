# Supuestos y decisiones (tarea 04)

Cada decisión que no estaba fijada en el brief, con su motivo.

## Fase 0
- **Repositorio git.** `amejia` no era un repositorio. Se inicializó en la propia carpeta (que también es `DATA_ROOT`) con un `.gitignore` en lista blanca. Primer commit = estado inicial sin las correcciones, para que el historial muestre cada una.
- **Remoto y push.** Al inicio no había remoto y `gh` no tenía sesión, así que las fases 0, 1, 2 y 6 se commitearon en local. El 2026-10-22 Astrid creó el remoto `https://github.com/astridmejia-11/QUS-BFM.git`. Se conectó como `origin`, se integró su commit inicial con `--allow-unrelated-histories` (conservando nuestro README) y se hizo push normal, sin force y sin `gh`. Antes se verificó que ningún `.mat`, imagen o archivo de más de 50 MB estuviera en el historial (máximo 18 KB). Desde ahí, push al cerrar cada fase. El clon auxiliar `amejia/QUS-BFM` se borró con autorización.
- **`qus_new.m`, banda.** "Configuración vigente" = 4 a 9 MHz (la que CLAUDE.md lista como propia). La salida se renombró a `qus_new_f49.mat` / `qus_new_f49_table.csv` para no sobrescribir `qus_new_f79.mat` ni `qus_new_table.csv`, que ya se presentaron.
- **`qus_bfcomp.m`, rama de respaldo.** Se reemplazó `Dp.params.fs`, `.startDepth`, `.xe` por los campos que `load_vsx_rf` sí devuelve (`Dp.fs`, `Dp.startDepth_m`, `Dp.elemX`). La rama no se ejecuta con la data actual, porque los archivos del lab traen `x` y `z`.

## Fase 1
- **Parámetros del beamformer propio** = los de `qus_new` (F/2, Hanning, z de 3 a 60 mm, `dz` = 1540/(2·31.25 MHz) fijo para todas las reconstrucciones). Muestras a 1542 m/s (c del fondo del CIRS 049).
- **Frame.** Los archivos de canal de `DataTargets2809` no traen `RcvLastFrame`; `load_vsx_rf` usa el frame 1 (el único en el buffer).
- **Centros en el juego propio.** Se parte del centro del lab llevado a los ejes propios (z − 1.50 mm, x × 0.300/0.2918, tarea 02) y se ajusta el círculo sobre la imagen propia. Regla fijada antes de usar los centros: si el ajuste mueve el centro más de 2 mm, se rechaza y se usa el punto de partida. Solo afectó a T3 (isoecoico).
- **Entrada del estimador.** `merino/phantoms_dt2809.m` elige la entrada con la variable de entorno `QUS_BFSET` (`lab`/`own`), porque `setup` hace `clear`. Con `lab` reproduce exactamente la tabla de la tarea 03.
- **Fondo por franjas.** Se usan puntos de la rejilla ACS a ≥ 1.05 cm del centro (radio de fondo de Merino). Franjas de 0.5–1.5, 1.5–2.5 y 2.5–3.5 cm.
- **Fechas de los preSets.** La fecha del archivo es de copia, no de adquisición; no sirve para ordenar las adquisiciones.

## Fase 2
- **Centros en `targets_edu`.** Ajuste de círculo sobre el frame 1. Para targets con contraste lateral < 3 dB (T3 y T6, isoecoicos), regla fijada antes de calibrar: reajustar desde la mediana de los centros con contraste; si se mueve más de 2 mm, usar el centro publicado por Merino (T6: 1.8, 1.9 cm) o, si no existe, la mediana (T3: 1.914, 1.981 cm). Ambos cayeron en el respaldo.
- **muC de Merino para T1–T5** (comparación de la fase 2 y juego `lab`): su caso 2 (10^1), igual que en `phantoms_dt2809`.

## Fase 3
- **Banda en la data antigua.** Regla: si 2.5–7.5 MHz no cabe dentro de la banda a −20 dB del espectro medio de la referencia (1–3 cm de profundidad, todas las columnas), se usa la banda a −20 dB. Se aplicó: PW 3.53–9.85 MHz, FOC 3.61–10.01 MHz. Como en Merino, la ventana de 8×12 λ se calcula con λ = 1540/centro de la banda, así que también cambia el tamaño de bloque.
- **Reconstrucción de la data antigua** (`scripts/bf_own_old.m`): parámetros de `qus_all` (F/2, Hanning, z de 3 a 58 mm, `dz` = 1515/(2·27.78 MHz)). Muestras a 1542 y Wisconsin a 1535.7 (regla 3; los cachés de `qus_all` estaban a 1540). Frame = `RcvLastFrame` de cada archivo. Las referencias Wisconsin comparten un preSet por carpeta (emparejamiento de respaldo de `load_vsx_rf`).
- **Referencia Wisconsin:** ACS 0.5246 con exponente 1.0019 (reporte UW).
- **Centros en la data antigua:** misma regla que en `targets_edu`. T3 PW y T7 FOC cayeron en el reajuste desde la mediana; ambos se aceptaron (movimiento ≤ 2 mm).
- **Corridas de `DataTargets2809` con parámetros congelados:** `QUS_SWIFT = [4 2 umbral]` aplica muB = 10^4 y muC = 10^2 a SWIFT en todos los targets. TV y SWTV conservan los valores de Merino.
- **MATLAB -batch que no cierra.** Algunas corridas de `phantoms_dt2809` escribieron sus resultados y el proceso no terminó, lo que dejó la cadena colgada una noche. Se relanzó con un vigilante: si `maps.mat` existe y el proceso no sale en 2 minutos, se termina. No afecta los resultados, porque `maps.mat` es lo último que se escribe.

## Fase 5
- **Figuras.** 166 en `results/figs/<fase>/`: modo B con la región de análisis (rectángulo verde) y las máscaras de medida (círculo continuo r = 0.85 cm, inclusión; discontinuo r = 1.05 cm, borde interior del fondo), mapa de ACS SWIFT y mapa de ΔBSC de la primera pasada, todo en cm. No se suben al repo (`*.png` ignorado); se regeneran con `scripts/make_figures.m`. Se fuerza el tema claro de MATLAB porque las sesiones `-batch` usan el oscuro (texto gris sobre blanco).
- **Fase 4 sin figura.** `ASTRID` no se procesó.
