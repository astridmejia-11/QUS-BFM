% bf_own_dt2809.m - Task 04, phase 1: own reconstruction of DataTargets2809
% from channel data, each phantom at the speed of sound of its certificate
% (rule 3): reference 261 at 1509, reference 544 at 1539, CIRS 049 samples at
% 1542 (background). Same axial grid for all (fixed dz), beamformer settings
% of the current configuration (qus_new): F/2, Hanning, z 3-60 mm.
%
% Output mirrors the lab layout and variable names (rf, x, z, fs) so the
% Merino-based scripts read it unchanged:
%   results/bf_own/Samples/<name>.mat
%   results/bf_own/Reference/<261|544>/<name>.mat
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(fullfile(DATA_ROOT, 'src'));
base   = fullfile(DATA_ROOT, 'DataTargets2809');
outDir = fullfile(DATA_ROOT, 'results', 'bf_own');

FS = 31.25e6;
bfOpt = {'zLim', [3e-3 60e-3], 'dz', 1540/(2*FS), 'FNumber', 2, 'Apod', 'hanning', ...
         'Verbose', false};

jobs = {};   % {channel file, output file, c}
s = dir(fullfile(base, 'Samples', 'ChannelData', '*.mat'));
s = s(~contains({s.name}, 'preSet') & ~contains({s.name}, 'preSave'));
for k = 1:numel(s)
    jobs(end+1,:) = {fullfile(s(k).folder, s(k).name), ...
        fullfile(outDir, 'Samples', s(k).name), 1542}; %#ok<AGROW>
end
refs = {'261', 1509; '544', 1539};
for r = 1:size(refs,1)
    d = dir(fullfile(base, 'Reference', 'ChannelData', refs{r,1}, '*.mat'));
    d = d(~contains({d.name}, 'preSet') & ~contains({d.name}, 'preSave'));
    for k = 1:numel(d)
        jobs(end+1,:) = {fullfile(d(k).folder, d(k).name), ...
            fullfile(outDir, 'Reference', refs{r,1}, d(k).name), refs{r,2}}; %#ok<AGROW>
    end
end

for j = 1:size(jobs,1)
    [od, nm] = fileparts(jobs{j,2});
    if ~isfolder(od), mkdir(od); end
    if isfile(jobs{j,2}), fprintf('[%2d/%d] %s exists\n', j, size(jobs,1), nm); continue; end
    D = load_vsx_rf(jobs{j,1}, 'c', jobs{j,3}, 'Verbose', false);
    B = das_focused(D, bfOpt{:});
    rf = B.rf; x = B.x(:).'; z = B.z(:).'; fs = B.params.fs;
    c = jobs{j,3}; frame = D.frame; srcFile = jobs{j,1}; %#ok<NASGU>
    save(jobs{j,2}, 'rf', 'x', 'z', 'fs', 'c', 'frame', 'srcFile');
    fprintf('[%2d/%d] %-10s c %d  frame %d  %dx%d\n', j, size(jobs,1), nm, c, frame, size(rf,1), size(rf,2));
    clear D B
end
