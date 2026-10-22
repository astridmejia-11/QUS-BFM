% check_deps.m - verify that phantoms_repro dependencies resolve to the Merino repo
here = fileparts(mfilename('fullpath'));
cd(fullfile(here, '..', 'SWIFT_UFFC2024'));
setup,
here = fileparts(mfilename('fullpath'));
cd(fullfile(here, '..'));
run(fullfile(pwd, 'config.m'));
names = {'setup','spectra','AlterOpti_ADMM','AlterOptiAdmmAnisWeighted', ...
    'optimAdmmTvTikhonov','optimAdmmWeightedTvTikhonov','save_all_figures_to_directory'};
ok = true;
for i = 1:numel(names)
    all = which(names{i}, '-all');
    fprintf('%-30s -> %s\n', names{i}, all{1});
    for k = 2:numel(all), fprintf('%30s    shadowed: %s\n', '', all{k}); end
    ok = ok && startsWith(all{1}, MERINO_DIR);
end
fprintf('All resolve to Merino repo: %d\n', ok);
fprintf('MATLAB %s\n', version);
