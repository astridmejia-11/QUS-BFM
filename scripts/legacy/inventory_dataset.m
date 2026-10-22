function T = inventory_dataset(rootDir, outCsv)
%INVENTORY_DATASET  Summarise a LIM-PUCP Verasonics dataset.
%
%   T = inventory_dataset(rootDir)
%   T = inventory_dataset(rootDir, 'inventory.csv')
%
%   This dataset stores each acquisition as TWO files:
%       <name>.mat                          -> RcvData, RcvLastFrame
%       <name>preSet.mat or <name>_preSet.mat -> preSet.{Trans,TX,Receive,...}
%   Both naming conventions occur in the same dataset, so both are tried.
%   Several data files may legitimately share one preSet (reference
%   acquisitions: identical settings, repeated over independent speckle).
%
%   Returns ONE ROW PER DATA FILE with the settings read from its matched
%   preSet, converted to millimetres where that aids interpretation.
%   RcvData is never loaded.
%
%   Rationale: the reference-phantom method is only valid when a sample and
%   its reference share every acquisition setting -- TGC above all, since a
%   TGC mismatch enters the spectral ratio and is indistinguishable from the
%   attenuation you are trying to measure.

if nargin < 1 || isempty(rootDir), rootDir = pwd; end
if nargin < 2, outCsv = ''; end

wState  = warning('off', 'MATLAB:load:variableNotFound');
cleaner = onCleanup(@() warning(wState)); %#ok<NASGU>

files = dir(fullfile(rootDir, '**', '*.mat'));
files = files(~[files.isdir]);
if isempty(files), T = table(); return; end
fprintf('Scanning %d .mat files under %s ...\n', numel(files), rootDir);

% ---- 1. classify each file ----------------------------------------------
isPre  = false(numel(files), 1);
isData = false(numel(files), 1);
for k = 1:numel(files)
    try
        w  = whos('-file', fullfile(files(k).folder, files(k).name));
        nm = {w.name};
    catch
        nm = {};
    end
    isPre(k)  = any(strcmp(nm, 'preSet'));
    isData(k) = any(strcmp(nm, 'RcvData'));
end
fprintf('  %d data file(s), %d preSet file(s), %d unrecognised\n', ...
    nnz(isData), nnz(isPre), nnz(~isData & ~isPre));

dataIdx = find(isData);
preIdx  = find(isPre);
preFull   = arrayfun(@(f) string(fullfile(f.folder, f.name)), files(preIdx));
preFolder = arrayfun(@(f) string(f.folder), files(preIdx));
preTime   = arrayfun(@(f) f.datenum,        files(preIdx));

rows = repmat(defaultRow(), max(numel(dataIdx), 1), 1);
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:numel(dataIdx)
    k = dataIdx(i);
    r = defaultRow();
    r.folder = string(erase(files(k).folder, rootDir));
    r.file   = string(files(k).name);
    r.fileMB = round(files(k).bytes/1e6, 1);

    % ---- 2. pair with a preSet ------------------------------------------
    % Two naming conventions coexist: "<base>preSet.mat" and
    % "<base>_preSet.mat". Try both before falling back to a guess.
    [~, base] = fileparts(files(k).name);
    cand = [string(fullfile(files(k).folder, [base 'preSet.mat'])), ...
            string(fullfile(files(k).folder, [base '_preSet.mat']))];
    j = [];
    for cc = 1:numel(cand)
        j = find(preFull == cand(cc), 1);
        if ~isempty(j), break; end
    end

    if ~isempty(j)
        r.pairing = "exact";
    else
        same = find(preFolder == string(files(k).folder));
        if isempty(same)
            r.status = "no preSet in this folder";
            rows(i) = r; continue
        end
        before = same(preTime(same) <= files(k).datenum);
        if ~isempty(before)
            [~, m] = max(preTime(before));            j = before(m);
        else
            [~, m] = min(abs(preTime(same) - files(k).datenum)); j = same(m);
        end
        r.pairing = "assumed";
    end
    r.preSetFile = string(files(preIdx(j)).name);

    % ---- 3. read the setup (cached: many files share one preSet) --------
    key = char(preFull(j));
    try
        if cache.isKey(key)
            m = cache(key);
        else
            S = load(key, 'preSet');
            if ~isfield(S, 'preSet')
                r.status = "file has no preSet variable"; rows(i) = r; continue
            end
            m = readSetup(S.preSet);
            cache(key) = m;
        end
        fn = fieldnames(m);
        for q = 1:numel(fn), r.(fn{q}) = m.(fn{q}); end
        r.status = "ok";
    catch ME
        r.status = "ERROR: " + string(ME.message);
    end

    rows(i) = r;
end

T = struct2table(rows);
T = sortrows(T, {'folder', 'file'});

if ~isempty(outCsv)
    writetable(T, outCsv);
    fprintf('Wrote %s (%d rows)\n', outCsv, height(T));
end

report(T);
end

% =========================================================================
function m = readSetup(p)
m = struct();
lambda_mm = NaN; pitch_mm = NaN;

% Speed of sound assumed at acquisition. Without a Resource struct it cannot
% be read back, so it is an ASSUMPTION until confirmed with the lab. Every
% wavelength-to-millimetre conversion below depends on it.
if isfield(p, 'P') && isfield(p.P, 'speedOfSound')
    m.c_setup = p.P.speedOfSound;  m.cSource = "preSet.P";
else
    m.c_setup = 1540;              m.cSource = "ASSUMED 1540";
end

if isfield(p, 'Trans')
    m.probe  = string(p.Trans.name);
    m.fc_MHz = p.Trans.frequency;
    m.nElem  = p.Trans.numelements;
    if isfield(p.Trans, 'units'), m.transUnits = string(p.Trans.units); end
    if isfield(p.Trans, 'lensCorrection'), m.lensCorr_wl = p.Trans.lensCorrection; end
    lambda_mm = m.c_setup / (m.fc_MHz * 1e6) * 1e3;
    % Element pitch, needed for the f-number.
    if isfield(p.Trans, 'ElementPos') && size(p.Trans.ElementPos, 1) > 1
        d = median(abs(diff(p.Trans.ElementPos(:, 1))));
        if isfield(p.Trans, 'units') && strcmpi(char(p.Trans.units), 'mm')
            pitch_mm = d;
        else
            pitch_mm = d * lambda_mm;
        end
    end
end
m.lambda_mm = lambda_mm;
m.pitch_mm  = pitch_mm;

if isfield(p, 'Receive') && ~isempty(p.Receive)
    R = p.Receive(1);
    m.sampleMode = string(R.sampleMode);
    m.spw        = samplesPerWave(R);
    if isfield(R, 'decimSampleRate') && ~isempty(R.decimSampleRate)
        m.fs_MHz = R.decimSampleRate;
    else
        m.fs_MHz = m.fc_MHz * m.spw;
    end
    m.startDepth_mm = R.startDepth * lambda_mm;
    m.endDepth_mm   = R.endDepth   * lambda_mm;
    m.nSamples      = R.endSample - R.startSample + 1;
    m.nRecvTotal    = numel(p.Receive);
    if isfield(R, 'framenum'), m.nFrames      = max([p.Receive.framenum]); end
    if isfield(R, 'acqNum'),   m.nAcqPerFrame = max([p.Receive.acqNum]);   end
end

if isfield(p, 'TX') && ~isempty(p.TX)
    m.nTX = numel(p.TX);
    foc = [p.TX.focus];
    if isfield(p.TX, 'Steer'), st = arrayfun(@(t) t.Steer(1), p.TX);
    else,                      st = zeros(size(foc));
    end
    % Use the MAXIMUM active aperture, not the mean: scanlines near the array
    % edge are truncated, so the mean understates the intended aperture and
    % gives a meaningless f-number.
    apod = arrayfun(@(t) nnz(t.Apod), p.TX);
    m.apMaxElem = max(apod);
    m.apMinElem = min(apod);
    m.steerDeg  = max(abs(st)) * 180/pi;
    m.focus_mm  = mean(foc) * lambda_mm;

    nOrigins = numel(unique(round(arrayfun(@(t) t.Origin(1), p.TX), 4)));
    if     all(foc == 0),  m.txType = "planewave";
    elseif all(foc <  0),  m.txType = "diverging";
    elseif nOrigins  > 1,  m.txType = "focused-scanlines";
    else,                  m.txType = "focused";
    end

    if m.focus_mm > 0 && isfinite(pitch_mm)
        m.fnumber = m.focus_mm / (m.apMaxElem * pitch_mm);
    else
        m.fnumber = NaN;   % undefined for an unfocused transmit
    end
end

if isfield(p, 'TW') && ~isempty(p.TW) && isfield(p.TW, 'peak')
    m.twPeak_wl = p.TW(1).peak;
end

[m.tgc, m.tgcInfo] = tgcString(p);
m.hv_V = NaN;
if isfield(p, 'TPC') && ~isempty(p.TPC) && isfield(p.TPC, 'hv')
    m.hv_V = p.TPC(1).hv;
end
end

% =========================================================================
function report(T)
bad = T.status ~= "ok";
if any(bad)
    fprintf('\n--- %d data file(s) could not be characterised ---\n', nnz(bad));
    for i = find(bad)'
        fprintf('  %s\\%s  ->  %s\n', T.folder(i), T.file(i), T.status(i));
    end
end
Tok = T(~bad, :);
if isempty(Tok), fprintf('\nNothing readable.\n'); return; end

fprintf('\n--- distinct acquisition configurations ---\n');
fprintf('(a sample and its reference MUST share every field below)\n\n');
key = strcat(sstr(Tok.txType), ...
             " | focus=",  sstr(round(Tok.focus_mm, 1)), "mm", ...
             " | ap=",     sstr(Tok.apMaxElem), "el", ...
             " | f/",      sstr(round(Tok.fnumber, 2)), ...
             " | fs=",     sstr(round(Tok.fs_MHz, 2)), "MHz", ...
             " | depth=",  sstr(round(Tok.startDepth_mm, 1)), "-", ...
                           sstr(round(Tok.endDepth_mm, 1)), "mm", ...
             " | TGC=",    sstr(Tok.tgc), ...
             " | hv=",     sstr(Tok.hv_V), "V");
[u, ~, ic] = unique(key);
for i = 1:numel(u)
    sel  = ic == i;
    fprintf('CONFIG %d  (%d file(s))\n  %s\n', i, nnz(sel), u(i));
    fold = unique(Tok.folder(sel));
    for q = 1:numel(fold)
        fprintf('     %2d in %s\n', nnz(sel & Tok.folder == fold(q)), fold(q));
    end
    fprintf('\n');
end

nAssumed = nnz(Tok.pairing == "assumed");
if nAssumed > 0
    fprintf(['NOTE: %d file(s) had no exactly-named preSet; paired with the\n' ...
             '      nearest one in the same folder. Verify these.\n'], nAssumed);
end
if all(Tok.tgc == "?" | Tok.tgc == "")
    fprintf(['\nWARNING: TGC could not be read from any preSet. Its raw shape is\n' ...
             '         in column tgcInfo. A TGC mismatch between sample and\n' ...
             '         reference is indistinguishable from tissue attenuation,\n' ...
             '         so this must be resolved before estimating ACS.\n']);
    fprintf('         tgcInfo: %s\n', unique(Tok.tgcInfo));
end
end

% =========================================================================
function r = defaultRow()
r = struct( ...
    'folder', "", 'file', "", 'fileMB', NaN, 'preSetFile', "", 'pairing', "", ...
    'probe', "", 'fc_MHz', NaN, 'nElem', NaN, 'transUnits', "", ...
    'lensCorr_wl', NaN, 'lambda_mm', NaN, 'pitch_mm', NaN, ...
    'sampleMode', "", 'spw', NaN, 'fs_MHz', NaN, ...
    'startDepth_mm', NaN, 'endDepth_mm', NaN, 'nSamples', NaN, ...
    'nRecvTotal', NaN, 'nFrames', NaN, 'nAcqPerFrame', NaN, ...
    'nTX', NaN, 'txType', "", 'focus_mm', NaN, 'fnumber', NaN, ...
    'apMaxElem', NaN, 'apMinElem', NaN, 'steerDeg', NaN, ...
    'twPeak_wl', NaN, 'c_setup', NaN, 'cSource', "", ...
    'tgc', "", 'tgcInfo', "", 'hv_V', NaN, 'status', "");
end

% =========================================================================
function [s, info] = tgcString(p)
%TGCSTRING  Extract the TGC control points, whatever shape they are stored in.
s = ""; info = "absent";
v = [];
if     isfield(p, 'TGCparam'), v = p.TGCparam; info = describeVar(v);
elseif isfield(p, 'TGC'),      v = p.TGC;      info = describeVar(v);
end
if isempty(v), return; end

if isnumeric(v)
    s = string(strtrim(num2str(double(v(:)'), '%g ')));
elseif isstruct(v)
    fn = fieldnames(v);
    for i = 1:numel(fn)
        w = v(1).(fn{i});
        if isnumeric(w) && numel(w) >= 2 && numel(w) <= 16
            s = string(fn{i}) + "=" + ...
                string(strtrim(num2str(double(w(:)'), '%g ')));
            return
        end
    end
elseif iscell(v) && ~isempty(v) && isnumeric(v{1})
    s = string(strtrim(num2str(double(v{1}(:)'), '%g ')));
end
end

% =========================================================================
function s = describeVar(v)
s = string(class(v)) + " " + string(mat2str(size(v)));
if isstruct(v)
    s = s + " {" + strjoin(fieldnames(v)', ",") + "}";
end
end

% =========================================================================
function spw = samplesPerWave(R)
if isfield(R, 'samplesPerWave') && ~isempty(R.samplesPerWave)
    spw = R.samplesPerWave; return
end
switch upper(char(R.sampleMode))
    case 'NS200BW',  spw = 4;
    case 'NS200BWI', spw = 2;
    case 'BS100BW',  spw = 2;
    case 'BS67BW',   spw = 4/3;
    case 'BS50BW',   spw = 1;
    otherwise,       spw = NaN;
end
end

% =========================================================================
function s = sstr(x)
%SSTR  Display-safe conversion: NaN/empty -> "?" instead of <missing>, which
%      would poison any strcat it takes part in.
s = string(x);
s(ismissing(s)) = "?";
s(s == "")      = "?";
end