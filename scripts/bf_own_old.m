% bf_own_old.m - Task 04, phase 3.1: reconstruct the old data (L11-4v,
% 7.24 MHz, fs 27.78 MHz) with each phantom at its own speed of sound
% (rule 3): CIRS 049 samples at 1542, Wisconsin references at 1535.7.
% PW with das_planewave, FOC with das_focused. Beamformer settings of
% qus_all (F/2, Hanning, z 3-58 mm, fixed dz = 1515/(2*fs)).
% Output: results/bf_old/<PW|FOC>/<samples|ref>/<name>.mat with rf, x, z, fs.
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(fullfile(DATA_ROOT, 'src'));
outDir = fullfile(DATA_ROOT, 'results', 'bf_old');
bfOpt = {'zLim', [3e-3 58e-3], 'dz', 1515/(2*27.7778e6), 'FNumber', 2, ...
         'Apod', 'hanning', 'Verbose', false};
dataOnly = @(d) d(~contains({d.name}, 'preSet') & ~contains({d.name}, 'preSave'));

sch = struct('name', {'PW', 'FOC'}, 'bf', {@das_planewave, @das_focused}, ...
    'samples', {'PLANE', 'TARGETS'}, 'ref', {'Referencia Plane', 'Ref Targets ATS'});
jobs = {};
for s = 1:2
    d = dataOnly(dir(fullfile(DATA_ROOT, sch(s).samples, '*.mat')));
    for k = 1:numel(d)
        jobs(end+1,:) = {s, fullfile(d(k).folder, d(k).name), ...
            fullfile(outDir, sch(s).name, 'samples', d(k).name), 1542}; %#ok<AGROW>
    end
    d = dataOnly(dir(fullfile(DATA_ROOT, sch(s).ref, '*.mat')));
    for k = 1:numel(d)
        jobs(end+1,:) = {s, fullfile(d(k).folder, d(k).name), ...
            fullfile(outDir, sch(s).name, 'ref', d(k).name), 1535.7}; %#ok<AGROW>
    end
end
for j = 1:size(jobs,1)
    od = fileparts(jobs{j,3});
    if ~isfolder(od), mkdir(od); end
    [~, nm] = fileparts(jobs{j,3});
    if isfile(jobs{j,3}), fprintf('[%2d/%d] exists %s\n', j, size(jobs,1), nm); continue; end
    D = load_vsx_rf(jobs{j,2}, 'c', jobs{j,4}, 'Verbose', false);
    B = sch(jobs{j,1}).bf(D, bfOpt{:});
    rf = B.rf; x = B.x(:).'; z = B.z(:).'; fs = B.params.fs;
    c = jobs{j,4}; frame = D.frame; preSetFile = D.preSetFile; srcFile = jobs{j,2}; %#ok<NASGU>
    save(jobs{j,3}, 'rf', 'x', 'z', 'fs', 'c', 'frame', 'preSetFile', 'srcFile');
    fprintf('[%2d/%d] %-4s %-50s c %.1f frame %d %dx%d\n', j, size(jobs,1), ...
        sch(jobs{j,1}).name, nm, c, frame, size(rf,1), size(rf,2));
    clear D B
end
