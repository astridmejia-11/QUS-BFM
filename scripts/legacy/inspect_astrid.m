% INSPECT_ASTRID  What is in the new acquisition block, and can the existing
%                 pipeline process it unchanged?
%
%   The reference-phantom method is only valid when sample and reference
%   share every acquisition setting. A TGC mismatch in particular enters the
%   spectral ratio and is indistinguishable from attenuation. So before
%   running anything, this checks the new files against the settings the
%   existing pipeline was built around.
%
%   Run from the amejia folder:  inspect_astrid
%
%   No local functions, so it also runs if executed section by section.
%   Plain ASCII on purpose.

clear; clc;

NEWDIR = fullfile(pwd, 'DataTargets2809');

% ---- what the existing pipeline assumes ---------------------------------
REF_fc   = 7.24e6;
REF_fs   = 27.7778e6;
REF_tgcA = [200 375 490 560 630 700 765 830];   % ATS acquisitions
REF_tgcB = [233 436 570 651 733 814 889 965];   % CIRS and thyroid

files = dir(fullfile(NEWDIR, '**', '*.mat'));
files = files(~[files.isdir]);
if isempty(files), error('Nothing found under %s', NEWDIR); end

isPre = false(numel(files),1); isDat = false(numel(files),1);
for k = 1:numel(files)
    try
        w = whos('-file', fullfile(files(k).folder, files(k).name));
        nm = {w.name};
    catch
        nm = {};
    end
    isPre(k) = any(strcmp(nm,'preSet'));
    isDat(k) = any(strcmp(nm,'RcvData'));
end
fprintf('%d data, %d preSet, %d unrecognised\n\n', ...
        nnz(isDat), nnz(isPre), nnz(~isDat & ~isPre));

dataIdx = find(isDat);
preList = arrayfun(@(f) string(fullfile(f.folder,f.name)), files(isPre));

for i = 1:numel(dataIdx)
    k = dataIdx(i);
    [~, base] = fileparts(files(k).name);
    fprintf('==================== %s ====================\n', files(k).name);

    % ---- pair with its preSet: three naming conventions coexist ---------
    cand = string(fullfile(files(k).folder, ...
           {[base 'preSet.mat'], [base '_preSet.mat'], [base '_preSave.mat']}));
    hit = cand(isfile(cand));
    if isempty(hit)
        m = preList(contains(preList, base));
        if isempty(m)
            fprintf('  NO preSet FOUND -- cannot process\n\n');
            continue
        end
        hit = m(1);
    end
    S = load(char(hit(1)), 'preSet');
    p = S.preSet;
    fprintf('  preSet: %s\n', char(erase(hit(1), [pwd filesep])));

    % ---- transducer and sampling ----------------------------------------
    fc = p.Trans.frequency*1e6;
    fs = p.Receive(1).decimSampleRate*1e6;
    if abs(fc-REF_fc)/REF_fc < 1e-3, sfc = '(igual)'; else, sfc = '(DIFIERE)'; end
    if abs(fs-REF_fs)/REF_fs < 1e-3, sfs = '(igual)'; else, sfs = '(DIFIERE)'; end
    fprintf('  transducer %s | %d elementos\n', p.Trans.name, size(p.Trans.ElementPos,1));
    fprintf('  fc %.3f MHz %s | fs %.4f MHz %s\n', fc/1e6, sfc, fs/1e6, sfs);

    % ---- transmit scheme: this is what decides which reference is valid --
    lambdaMM = 1540/fc*1e3;
    nTX = numel(p.TX);
    foc = 0;
    if isfield(p.TX(1),'focus'), foc = p.TX(1).focus; end
    ang = 0;
    if isfield(p.TX(1),'Steer'), ang = p.TX(1).Steer(1); end
    if foc == 0 || isinf(foc)
        fprintf('  TRANSMISION: onda plana | %d disparos | steer %.2f deg\n', ...
                nTX, rad2deg(ang));
    else
        fprintf('  TRANSMISION: focalizada | %d disparos | foco %.1f wl = %.1f mm\n', ...
                nTX, foc, foc*lambdaMM);
    end

    % ---- depth window and t0 terms --------------------------------------
    nS = p.Receive(1).endSample - p.Receive(1).startSample + 1;
    fprintf('  muestras %d | startDepth %.1f wl | endDepth %.1f wl = %.1f mm\n', ...
            nS, p.Receive(1).startDepth, p.Receive(1).endDepth, ...
            p.Receive(1).endDepth*lambdaMM);
    lc = 0;
    if isfield(p.Trans,'lensCorrection'), lc = p.Trans.lensCorrection; end
    fprintf('  TW.peak %.4f wl | lensCorrection %.4f wl\n', p.TW(1).peak, lc);

    % ---- TGC: the setting that must match the reference ------------------
    g = [];
    if isfield(p,'TGCparam')
        raw = p.TGCparam;
        if isnumeric(raw)
            g = raw(:)';
        elseif iscell(raw)
            fprintf('  TGC: cell con %d entradas\n', numel(raw));
            if isnumeric(raw{1}), g = raw{1}(:)'; end
        elseif isstruct(raw)
            fn = fieldnames(raw);
            fprintf('  TGC: struct con campos %s\n', strjoin(fn', ', '));
            for q = 1:numel(fn)
                v = raw(1).(fn{q});
                if isnumeric(v) && numel(v) >= 8
                    fprintf('    %s = %s\n', fn{q}, mat2str(v(:)'));
                    if isempty(g), g = v(:)'; end
                end
            end
        end
    else
        fprintf('  TGC: campo ausente\n');
    end

    if ~isempty(g)
        fprintf('  TGC %s\n', mat2str(g));
        if isequal(g, REF_tgcA)
            fprintf('  -> coincide con curva A (set ATS)\n');
        elseif isequal(g, REF_tgcB)
            fprintf('  -> coincide con curva B (CIRS/tiroides): referencias existentes usables\n');
        else
            fprintf('  -> CURVA TGC NUEVA: necesita referencia adquirida con esta misma curva\n');
            if numel(g) == numel(REF_tgcB)
                fprintf('     razon respecto a curva B: %s\n', mat2str(round(g./REF_tgcB,3)));
            end
        end
    end
    fprintf('\n');
end

fprintf(['\nQUE MIRAR EN LA SALIDA\n' ...
 '  1. Todos los archivos focalizados al mismo foco? El nombre dice 30 mm,\n' ...
 '     pero lo que lee el beamformer es el preSet.\n' ...
 '  2. La TGC es la misma curva de las adquisiciones CIRS existentes? Si no,\n' ...
 '     ninguna referencia actual es valida para este bloque.\n' ...
 '  3. Hay alguna referencia Wisconsin en esta carpeta? Ocho targets sin\n' ...
 '     referencia significa que el bloque no se puede procesar solo.\n']);