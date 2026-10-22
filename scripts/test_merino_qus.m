% test_merino_qus.m - src/merino_qus must reproduce Task 01 (phantoms.m on
% targets_edu, T6 T7 T8, Merino's parameters and centres) for TV and SWIFT.
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(genpath(fullfile(MERINO_DIR, 'AttenuationUtils')));
addpath(fullfile(DATA_ROOT, 'src'));
root = fullfile(DATA_ROOT, 'targets_edu');
d = dir(fullfile(root, 'referencia', '544', '*.mat'));
for j = 1:numel(d)
    S = load(fullfile(d(j).folder, d(j).name), 'RF');
    if j == 1, ref = zeros([size(S.RF,1) size(S.RF,2) numel(d)]); end
    ref(:,:,j) = S.RF(:,:,1);
end
T01 = readtable(fullfile(DATA_ROOT, 'results', 'merino', 'phantoms.xlsx'));
cases = {6, 2, [1.8 1.9], [3.5 3.5]; 7, 1, [1.95 1.95], [3 3]; 8, 1, [1.85 1.9], [3.5 3.5]};
for k = 1:3
    S = load(fullfile(root, 'data', sprintf('T%d.mat', cases{k,1})), 'RF', 'x', 'z', 'fs');
    o = struct('lmuC', cases{k,2}, 'centre', cases{k,3}, 'lmuTV', cases{k,4});
    R = merino_qus(S.RF(:,:,1), ref, S.x, S.z, S.fs, o);
    tv = T01(strcmp(T01.method,'TV'),:); sw = T01(strcmp(T01.method,'SWIFT'),:);
    e = [R.incTV-tv.meanInc(k), R.backTV-tv.meanBack(k), R.incSWIFT-sw.meanInc(k), R.backSWIFT-sw.meanBack(k)];
    fprintf('T%d: TV inc %.4f back %.4f | SWIFT inc %.4f back %.4f | max |diff| %.2g\n', ...
        cases{k,1}, R.incTV, R.backTV, R.incSWIFT, R.backSWIFT, max(abs(e)));
    assert(max(abs(e)) < 1e-6);
end
fprintf('PASS: merino_qus reproduces Task 01 (TV and SWIFT)\n');
