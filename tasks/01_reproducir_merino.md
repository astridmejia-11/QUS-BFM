# Tarea 01: reproducir `phantoms.m` de Merino con `targets_edu`

## Objetivo

Obtener, con el código de Merino **sin modificar su método**, valores de ACS cercanos al certificado del CIRS 049. Si se logra, el estimador queda validado de forma independiente al beamformer propio, y ese resultado pasa a ser la referencia de todo lo demás.

## Por qué primero

Hasta ahora se depuraba sin un resultado correcto conocido. Este es el primero: misma data, mismo código, resultado publicado. Si no se reproduce, el problema está en el entorno (rutas, versiones, archivos) y hay que resolverlo antes de tocar nada propio.

## Paso 1: inspeccionar `targets_edu`

No asumir estructura. Listar todo y mirar las variables de un archivo de muestra y uno de referencia.

```matlab
run('config.m');
root = fullfile(DATA_ROOT, 'targets_edu');

d = dir(fullfile(root, '**', '*.mat'));
fprintf('%d archivos\n', numel(d));
for i = 1:numel(d)
    fprintf('[%s]\n', erase(fullfile(d(i).folder, d(i).name), [root filesep]));
end

w = whos('-file', fullfile(d(1).folder, d(1).name));
for i = 1:numel(w)
    fprintf('  %-10s %-8s %s\n', w(i).name, w(i).class, mat2str(w(i).size));
end
```

Reportar a Astrid:
- qué carpeta tiene las muestras y cuál las referencias;
- si los archivos contienen `RF`, `x`, `z`, `fs` (lo que espera `phantoms.m`);
- cuántos archivos de muestra hay y en qué orden los devuelve `dir`.

**Ojo:** `phantoms.m` toma `targetFiles(end-2:end)`, los tres últimos archivos de la carpeta. El orden alfabético decide qué targets se procesan, y el ground truth `[0.97, 0.95, 0.95, 0.55]` asume que son T6, T7 y T8. Verificar que coincida.

## Paso 2: preparar la copia

1. Copiar `phantoms.m` a `merino/phantoms_repro.m`.
2. Cambiar **solo** estas tres líneas:

```matlab
dataDir    = fullfile(DATA_ROOT, 'targets_edu', '<carpeta de muestras>');
refDir     = fullfile(DATA_ROOT, 'targets_edu', '<carpeta de referencias>');
resultsDir = fullfile(pwd, 'results', 'merino');
```

3. Las rutas en `phantoms.m` usan `[dataDir,'\*.mat']`. Funcionan en Windows; no cambiarlas.
4. Al inicio, antes de `setup`, agregar el repositorio de Merino al path y cargar la configuración:

```matlab
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(genpath('C:\Users\LENOVO\repos\SWIFT_UFFC2024'));
```

5. Documentar los cambios en un comentario al inicio del archivo.

No cambiar banda, pesos, ventana, región ni métrica.

## Paso 3: dependencias

El script necesita, del repo de Merino: `setup`, `spectra`, `AlterOpti_ADMM`, `AlterOptiAdmmAnisWeighted`, `optimAdmmTvTikhonov`, `optimAdmmWeightedTvTikhonov`, `save_all_figures_to_directory`.

Antes de correr, comprobar con `which <nombre>` que cada una resuelve **al repo de Merino** y no a una copia del laboratorio en `lab/` con el mismo nombre pero otra versión. Si hay conflicto de nombres, reportarlo antes de seguir.

## Paso 4: correr y reportar

Predicción, escrita antes de correr: los valores de SWIFT deberían quedar cerca de `[0.97, 0.95, 0.95]` en las inclusiones y de 0.55 en el fondo.

Correr. La tabla se guarda en `results/merino/phantoms.xlsx`. Reportar a Astrid, para cada target y método (TV, SWTV, SWIFT): `meanInc`, `stdInc`, `meanBack`, `stdBack`, `biasInc`, `biasBack`.

## Criterio de éxito

- Sesgo de SWIFT en inclusiones y fondo claramente menor que el de la data propia (que rondaba −17 % en el fondo).
- Si no se acerca al certificado, **no ajustar parámetros**. Revisar primero: orden de archivos, carpeta de referencias correcta, versión de las funciones.

## Qué no hacer

- No reemplazar funciones de Merino por las del laboratorio ni al revés.
- No cambiar parámetros para "mejorar" el resultado.
- No subir `targets_edu` ni resultados al repositorio.