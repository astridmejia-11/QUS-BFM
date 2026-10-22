% qus_paths.m - path for the legacy scripts in scripts/legacy (qus_*, fig_*).
% They expect load_vsx_rf, das_* and the lab SWIFT_Simple code on the path
% and are run from the project root, e.g.
%   qus_paths; run('scripts/legacy/qus_new.m')
% The new scripts (merino/, scripts/) set their own path and do not need it.
% Do NOT combine with the Merino repo on the path: SWIFT_Simple holds
% same-name copies of spectra, AlterOpti_ADMM, optimAdmm*.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'lab', 'SWIFT_Simple'));
clear root
