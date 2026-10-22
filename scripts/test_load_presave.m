% test_load_presave.m - Phase 0 check: load_vsx_rf pairs preSave files
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(fullfile(DATA_ROOT, 'src'));
f = fullfile(DATA_ROOT, 'PLANE', 'RFdata_26_T8_-August-2026_15-53-50.mat');
w = whos('-file', strrep(f, '.mat', 'preSave.mat'));
fprintf('variables in preSave: %s\n', strjoin({w.name}, ', '));
lastwarn('');
D = load_vsx_rf(f, 'Verbose', true);
[msg, id] = lastwarn;
fprintf('preSet used: %s\n', D.preSetFile);
fprintf('warning: "%s" (%s)\n', msg, id);
assert(endsWith(D.preSetFile, 'RFdata_26_T8_-August-2026_15-53-50preSave.mat'));
assert(~strcmp(id, 'load_vsx_rf:assumedPreSet'));
fprintf('PASS\n');
