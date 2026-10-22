% compare_presets_dt2809.m - Task 04, phase 1, step 5: compare the preSets of
% references 261, 544 and the samples field by field (TGC, focus, aperture,
% voltage, sampling, acquisition date) and report what differs.
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
base = fullfile(DATA_ROOT, 'DataTargets2809');
outDir = fullfile(DATA_ROOT, 'results', 'phase1');
if ~isfolder(outDir), mkdir(outDir); end

pre = [dir(fullfile(base, 'Samples', 'ChannelData', '*preSet.mat')); ...
       dir(fullfile(base, 'Reference', 'ChannelData', '*', '*preSet.mat'))];
rows = {};
for k = 1:numel(pre)
    pf = fullfile(pre(k).folder, pre(k).name);
    S = load(pf, 'preSet'); P = S.preSet;
    df = strrep(pf, '_preSet.mat', '.mat');
    if isfile(df), dd = dir(df); dt = datestr(dd.datenum, 'yyyy-mm-dd HH:MM'); else, dt = '(no data)'; end
    [~, grp] = fileparts(pre(k).folder);
    if strcmp(grp, 'ChannelData'), grp = 'sample'; end
    tgc = tgcString(P);
    tx = P.TX(64);
    nAct = nnz(tx.Apod);
    hv = NaN; if isfield(P, 'TPC') && isfield(P.TPC, 'hv'), hv = P.TPC(1).hv; end
    R = P.Receive(1);
    tw = ''; if isfield(P, 'TW') && isfield(P.TW, 'Parameters'), tw = mat2str(P.TW(1).Parameters); end
    sos = NaN; if isfield(P, 'P') && isfield(P.P, 'speedOfSound'), sos = P.P.speedOfSound; end
    rows(end+1,:) = {string(grp), string(pre(k).name), string(dt), P.Trans.frequency, ...
        R.decimSampleRate, R.startDepth, R.endDepth, tx.focus, nAct, hv, string(tw), ...
        string(tgc), sos}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'group','preSet','dataDate','fc','fs', ...
    'startDepth','endDepth','txFocus_wl','txActiveEl','hv','TWparams','TGC','speedOfSound'});
T = sortrows(T, {'group','preSet'});
disp(T(:, {'group','preSet','dataDate','fc','fs','startDepth','endDepth','txFocus_wl','txActiveEl','hv','speedOfSound'}));
fprintf('\nTW parameters per file:\n');
for k = 1:height(T), fprintf('  %-22s %s\n', T.preSet(k), T.TWparams(k)); end
fprintf('\nTGC per file:\n');
for k = 1:height(T), fprintf('  %-22s %s\n', T.preSet(k), T.TGC(k)); end
writetable(T, fullfile(outDir, 'presets_dt2809.csv'));

% fields that differ between any two files
fprintf('\nColumns that are not identical across all files:\n');
for v = T.Properties.VariableNames(4:end)
    col = T.(v{1});
    if isstring(col), u = unique(col); else, u = unique(col(~isnan(col))); end
    if numel(u) > 1, fprintf('  %s: %d distinct values\n', v{1}, numel(u)); end
end

function s = tgcString(P)
s = '(none)';
if isfield(P, 'TGC') && ~isempty(P.TGC)
    g = P.TGC(1);
    if isfield(g, 'CntrlPts'), s = ['CntrlPts ' mat2str(g.CntrlPts)]; end
    if isfield(g, 'rangeMax'), s = [s ' rangeMax ' num2str(g.rangeMax)]; end
elseif isfield(P, 'TGCparam') && ~isempty(P.TGCparam)
    s = ['TGCparam ' jsonencode(P.TGCparam)];
end
end
