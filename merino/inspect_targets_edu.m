% inspect_targets_edu.m - list targets_edu files and their variables (Task 01, step 1)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
root = fullfile(DATA_ROOT, 'targets_edu');

d = dir(fullfile(root, '**', '*.mat'));
fprintf('%d files\n', numel(d));
for i = 1:numel(d)
    fprintf('[%s]\n', erase(fullfile(d(i).folder, d(i).name), [root filesep]));
end

% Order returned by dir for the sample folder, as phantoms.m sees it
t = dir(fullfile(root, 'data', '*.mat'));
fprintf('\nSample order: %s\n', strjoin({t.name}, ', '));
fprintf('targetFiles(end-2:end): %s\n', strjoin({t(end-2:end).name}, ', '));

files = {fullfile(root,'data','T6.mat'), fullfile(root,'data','T8.mat'), ...
         fullfile(root,'referencia','544','16-19-52.mat')};
for k = 1:numel(files)
    fprintf('\n--- %s\n', files{k});
    w = whos('-file', files{k});
    for i = 1:numel(w)
        fprintf('  %-12s %-8s %s\n', w(i).name, w(i).class, mat2str(w(i).size));
    end
    S = load(files{k}, 'x', 'z', 'fs');
    if isfield(S,'x'), fprintf('  x: %.4g .. %.4g, dx %.4g\n', S.x(1), S.x(end), S.x(2)-S.x(1)); end
    if isfield(S,'z'), fprintf('  z: %.4g .. %.4g, dz %.4g\n', S.z(1), S.z(end), S.z(2)-S.z(1)); end
    if isfield(S,'fs'), fprintf('  fs: %.6g\n', S.fs); end
end

% Same grid check across all sample and reference files
fprintf('\nGrid per file (size RF, x range, z range, fs):\n');
for i = 1:numel(d)
    f = fullfile(d(i).folder, d(i).name);
    S = load(f, 'x', 'z', 'fs');
    m = matfile(f); sz = size(m, 'RF');
    fprintf('  %-28s RF %s  x[%.4f %.4f] z[%.4f %.4f] fs %.4g\n', d(i).name, ...
        mat2str(sz), S.x(1), S.x(end), S.z(1), S.z(end), S.fs);
end
