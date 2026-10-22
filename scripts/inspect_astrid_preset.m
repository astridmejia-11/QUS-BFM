% inspect_astrid_preset.m - Task 04, phase 4: acquisition settings of the
% ASTRID block (Verasonics) to compare against targets_edu (Ultrasonix).
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
f = fullfile(DATA_ROOT, 'ASTRID', 'ASTRID', 'T1_L_Raw_F30mm_preSet.mat');
S = load(f, 'preSet'); P = S.preSet;
fprintf('Trans.name %s, frequency %.3f MHz, elements %d\n', char(string(P.Trans.name)), ...
    P.Trans.frequency, P.Trans.numelements);
pitch = median(diff(P.Trans.ElementPos(:,1)));
fprintf('pitch %.4f (units %s)\n', pitch, char(string(P.Trans.units)));
R = P.Receive(1);
fprintf('fs %.3f MHz, startDepth %.2f, endDepth %.2f wl\n', R.decimSampleRate, R.startDepth, R.endDepth);
fprintf('TX(1).focus %.2f wl, active elements %d, TW %s\n', P.TX(1).focus, nnz(P.TX(1).Apod), ...
    mat2str(P.TW(1).Parameters));
if isfield(P,'TPC'), fprintf('hv %.1f V\n', P.TPC(1).hv); end
if isfield(P,'TGCparam'), fprintf('TGC %s\n', jsonencode(P.TGCparam.TGC.CntrlPts)); end
fprintf('numel(TX) %d\n', numel(P.TX));
