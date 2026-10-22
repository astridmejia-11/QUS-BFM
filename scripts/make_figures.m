% make_figures.m - Task 04, phase 5: one figure per result of phases 1-3
% (B-mode with analysis and measurement regions, ACS map, Delta-BSC map,
% physical coordinates). Output: results/figs/<phase>/<name>.png
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(genpath(fullfile(MERINO_DIR, 'AttenuationUtils')));
addpath(fullfile(DATA_ROOT, 'src'));
figRoot = fullfile(DATA_ROOT, 'results', 'figs');
R0 = fullfile(DATA_ROOT, 'results');

% ---- phase 1 and 3 (DataTargets2809): every maps.mat under merino_dt2809* --
runs = {'fase1', 'merino_dt2809', 'lab_merino'; 'fase1', 'merino_dt2809_own', 'own_merino'};
d = dir(fullfile(R0, 'merino_dt2809*', '**', 'maps.mat'));
for k = 1:numel(d)
    rel = erase(d(k).folder, [R0 filesep]);
    if any(strcmp(rel, runs(:,2))), continue; end
    parts = strsplit(rel, filesep);
    grp = 'lab'; if contains(parts{1}, '_own'), grp = 'own'; end
    runs(end+1,:) = {'fase3', rel, [grp '_' strjoin(parts(2:end), '_')]}; %#ok<SAGROW>
end
for r = 1:size(runs,1)
    f = fullfile(R0, runs{r,2}, 'maps.mat');
    if ~isfile(f), continue; end
    M = load(f); M = M.Maps;
    for i = 1:numel(M)
        [~, nm] = fileparts(M(i).file);
        C = struct('Bmode', M(i).Bmode, 'x', M(i).x, 'z', M(i).z, 'x_ACS', M(i).x_ACS, ...
            'z_ACS', M(i).z_ACS, 'acs', M(i).SWIFT, 'bsc', M(i).bscSWIFT, ...
            'centre', [M(i).c1x M(i).c1z]);
        ttl = sprintf('DataTargets2809 %s, %s, ref %s, SWIFT muB 10^%.1f muC 10^%.1f cut %g', ...
            runs{r,3}, nm, M(i).ref, log10(M(i).muBswift), log10(M(i).muCswift), M(i).ratioCutOff);
        plot_qus_case(C, ttl, fullfile(figRoot, runs{r,1}, ...
            sprintf('dt2809_%s_%s_ref%s.png', runs{r,3}, nm, M(i).ref)));
    end
    fprintf('%s: %d figures\n', runs{r,2}, numel(M));
end

% ---- phase 2 (targets_edu, frozen parameters, all eight targets) ----------
F = load(fullfile(R0, 'calib_edu', 'calib_edu.mat'), 'best', 'TEST');
cen = readtable(fullfile(R0, 'calib_edu', 'centres.csv'));
root = fullfile(DATA_ROOT, 'targets_edu');
dr = dir(fullfile(root, 'referencia', '544', '*.mat'));
for j = 1:numel(dr)
    S = load(fullfile(dr(j).folder, dr(j).name), 'RF');
    if j == 1, ref = zeros([size(S.RF,1) size(S.RF,2) numel(dr)]); end
    ref(:,:,j) = S.RF(:,:,1);
end
for t = 1:8
    S = load(fullfile(root, 'data', sprintf('T%d.mat', t)), 'RF', 'x', 'z', 'fs');
    o = struct('lmuB', F.best.lmuB, 'lmuC', F.best.lmuC, 'cut', F.best.cut, ...
        'centre', [cen.c1x(t) cen.c1z(t)]);
    R = merino_qus(S.RF(:,:,1), ref, S.x, S.z, S.fs, o);
    grp = 'calibracion'; if ismember(t, F.TEST), grp = 'prueba'; end
    C = struct('Bmode', R.Bmode, 'x', R.x, 'z', R.z, 'x_ACS', R.x_ACS, 'z_ACS', R.z_ACS, ...
        'acs', R.SWIFT, 'bsc', R.bsc, 'centre', o.centre);
    plot_qus_case(C, sprintf('targets_edu T%d (%s), SWIFT congelado: muB 10^%.1f muC 10^%.1f cut %d', ...
        t, grp, F.best.lmuB, F.best.lmuC, F.best.cut), ...
        fullfile(figRoot, 'fase2', sprintf('edu_T%d_%s_frozen.png', t, grp)));
end
fprintf('fase2: 8 figures\n');

% ---- phase 3 (old data, PW vs FOC) ----------------------------------------
f = fullfile(R0, 'old_pwfoc', 'old_pwfoc.mat');
if isfile(f)
    O = load(f);
    for i = 1:numel(O.res)
        q = O.res(i); R = q.map;
        C = struct('Bmode', R.Bmode, 'x', R.x, 'z', R.z, 'x_ACS', R.x_ACS, 'z_ACS', R.z_ACS, ...
            'acs', R.SWIFT, 'bsc', R.bsc, 'centre', [q.c1x q.c1z]);
        plot_qus_case(C, sprintf('Data antigua %s %s, cut %d, banda %.1f-%.1f MHz', ...
            q.scheme, q.target, q.cut, q.bandL/1e6, q.bandH/1e6), ...
            fullfile(figRoot, 'fase3', sprintf('old_%s_%s_cut%d.png', q.scheme, q.target, q.cut)));
    end
    fprintf('fase3 old: %d figures\n', numel(O.res));
end
