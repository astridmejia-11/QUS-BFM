# Tarea 02: diagnosticar el "zoom" del beamformer propio

## Síntoma

Sobre la misma data cruda de `DataTargets2809`, el modo B reconstruido con `das_focused` se ve ampliado respecto al reconstruido por el laboratorio. Los asesores lo atribuyen a un error del beamformer propio.

## Por qué importa

Si la escala axial está mal, la distancia entre ventanas `L` también lo está, y la ACS hereda el error. Puede contribuir a la discrepancia de 0.158 entre referencias.

## Hipótesis, en orden de probabilidad

1. **Diferencia de rejilla, no de imagen.** El laboratorio conserva el Nz de la data de canal (2560 muestras desde startDepth). `das_focused` reconstruye solo `zLim = [3 60] mm` con `dz` fijado a mano. Si las dos imágenes se grafican por índice de pixel en vez de por coordenada en mm, una parece ampliada sin que haya error. Descartar esto primero.
2. **Origen temporal.** En este bloque `startDepth = 1.0 λ`; en la data antigua era 0. Si `load_vsx_rf` no suma el startDepth al eje temporal, toda la imagen queda corrida en profundidad.
3. **Escala axial.** `dz` fijado con 1515 o 1540 y `fs` de otro dataset (27.78 en vez de 31.25 MHz). Un `fs` equivocado escala toda la profundidad en la razón 31.25/27.78 = 1.125, un 12.5 %, que se vería como zoom.
4. **Términos de t0.** `TW.peak` (1.8696 λ) y `lensCorrection` (7.2965 λ) distintos a los de la data antigua; comprobar que se leen del preSet de este bloque.

## Procedimiento

Usar un solo archivo, por ejemplo `Samples/ChannelData/T7_F.mat`, con su versión de `Samples/BeamformedData`.

1. Cargar la versión del laboratorio y listar sus variables (`whos -file`). Anotar tamaño de la matriz, y si trae ejes `x`, `z` y en qué unidades.
2. Reconstruir con `das_focused` y anotar tamaño, `x(1)`, `x(end)`, `z(1)`, `z(end)`, `dz` y el `fs` que usó.
3. Graficar las dos **en coordenadas físicas** (`imagesc(x, z, ...)`) lado a lado, con los mismos límites de ejes.
4. Medir en ambas la posición del centro de la inclusión y su diámetro, en mm. Hacerlo con un perfil axial y uno lateral por el centro, no a ojo.
5. Comparar:
   - si coinciden en mm, el zoom era de visualización (hipótesis 1) y no hay error;
   - si el diámetro difiere en una razón constante, es escala (hipótesis 3);
   - si el diámetro coincide pero el centro está corrido, es origen (hipótesis 2 o 4).

## Datos de control

- Para la data antigua ya hay validación con los hilos del ATS 539: intercepto −0.25 mm y velocidades implicadas 1515/1516 m/s. Ese resultado es del montaje antiguo y **no** valida este bloque, que tiene otro `fs` y otro startDepth.
- Si hay duda sobre la escala, verificar `dz = c/(2*fs)` con el `fs` de **este** preSet.

## Entregable

Una tabla corta con las dos reconstrucciones (tamaño, ejes, centro y diámetro de la inclusión), la hipótesis confirmada y, si hay error, el cambio mínimo en `load_vsx_rf` o `das_focused` que lo corrige. Después de corregir, repetir la comparación y confirmar que coinciden.

## Qué no hacer

- No cambiar la velocidad del sonido para "hacer coincidir" las imágenes.
- No tocar el código del laboratorio.