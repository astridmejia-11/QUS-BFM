% zoom_inspect.m - Task 02 step 1: what is inside the lab file and the preSet
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
base = fullfile(DATA_ROOT, 'DataTargets2809', 'Samples');
fLab = fullfile(base, 'BeamformedData', 'T7_F.mat');
fPre = fullfile(base, 'ChannelData', 'T7_F_preSet.mat');
fCh  = fullfile(base, 'ChannelData', 'T7_F.mat');

for f = {fLab, fCh, fPre}
    fprintf('\n--- %s\n', f{1});
    w = whos('-file', f{1});
    for i = 1:numel(w)
        fprintf('  %-18s %-8s %s\n', w(i).name, w(i).class, mat2str(w(i).size));
    end
end

S = load(fLab);
fn = fieldnames(S);
for i = 1:numel(fn)
    v = S.(fn{i});
    if isnumeric(v) && isvector(v) && numel(v) > 1
        fprintf('  %s: %.6g .. %.6g, step %.6g, n %d\n', fn{i}, v(1), v(end), v(2)-v(1), numel(v));
    elseif isnumeric(v) && isscalar(v)
        fprintf('  %s = %.6g\n', fn{i}, v);
    elseif isstruct(v)
        fprintf('  %s: struct with fields %s\n', fn{i}, strjoin(fieldnames(v), ', '));
    end
end

P = load(fPre, 'preSet'); P = P.preSet;
R = P.Receive(1);
fprintf('\npreSet: Trans.frequency %.4g MHz, numelements %d, units %s\n', ...
    P.Trans.frequency, P.Trans.numelements, char(string(P.Trans.units)));
fprintf('  Receive(1): startDepth %.4g, endDepth %.4g, startSample %d, endSample %d\n', ...
    R.startDepth, R.endDepth, R.startSample, R.endSample);
if isfield(R, 'decimSampleRate'), fprintf('  decimSampleRate %.6g\n', R.decimSampleRate); end
if isfield(R, 'samplesPerWave'), fprintf('  samplesPerWave %.6g\n', R.samplesPerWave); end
if isfield(R, 'sampleMode'), fprintf('  sampleMode %s\n', char(R.sampleMode)); end
fprintf('  TW.peak %.4g, lensCorrection %.4g\n', P.TW(1).peak, P.Trans.lensCorrection);
if isfield(P, 'P'), disp(P.P); end
fprintf('  TX(1).focus %.4g, Origin %s; TX(64).Origin %s\n', P.TX(1).focus, ...
    mat2str(P.TX(1).Origin), mat2str(P.TX(64).Origin));
fprintf('  ElementPos x: %.4g .. %.4g\n', P.Trans.ElementPos(1,1), P.Trans.ElementPos(end,1));
fprintf('  numel(TX) %d, numel(Receive) %d\n', numel(P.TX), numel(P.Receive));
