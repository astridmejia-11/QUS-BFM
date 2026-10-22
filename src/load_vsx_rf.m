function D = load_vsx_rf(dataFile, varargin)
%LOAD_VSX_RF  Load one LIM-PUCP acquisition into a beamformer-agnostic struct.
%
%   D = load_vsx_rf(dataFile)
%   D = load_vsx_rf(dataFile, 'Frame', 5, 'c', 1542)
%   D = load_vsx_rf(dataFile, 'PreSet', 'path\to\...preSet.mat')
%
%   Name-value options
%     'PreSet'  : explicit preSet path. Default: same base name with
%                 "preSet" or "_preSet" appended; if neither exists, the
%                 nearest preSet in the same folder (a warning is raised).
%     'Frame'   : frame index. Default: RcvLastFrame, i.e. the last frame the
%                 system actually wrote. The receive buffer is CIRCULAR, so
%                 frame 1 is not necessarily valid or recent -- it may hold a
%                 stale or partially overwritten acquisition.
%     'c'       : sound speed for reconstruction [m/s]. Default: the value
%                 assumed at acquisition.
%     'Events'  : acquisition indices within the frame. Default: all.
%     'Verbose' : print a summary. Default true.
%
%   The returned struct D is the ONLY thing the beamformer should see.
%   Everything Verasonics-specific -- the split data/preSet layout, channel
%   mapping, sample slicing, wavelength units, t0 offsets -- is resolved here
%   once. A focused acquisition and a plane-wave acquisition come out with
%   identical field names in SI units, so one DAS and one QUS pipeline can
%   consume both. That is what makes "focused vs plane wave" a fair
%   comparison: only the contents of D change, never the code.
%
%   Fields
%     D.rf        [nSamples x nElements x nEvents] double, raw ADC counts
%     D.t         [nSamples x 1] receive time axis [s], t = 0 at emission
%     D.fs D.fc   sampling and centre frequency [Hz]
%     D.c D.cAcq  reconstruction and acquisition sound speed [m/s]
%     D.elemX/Z   [nElements x 1] element positions [m]
%     D.tx(e)     .delay [s/element] .apod .origin [m] .focus [m] .steer [rad]
%                 .type ("focused" | "planewave")
%     D.tgcRaw D.hv D.satFraction D.maxAbs D.setup

% ---- options ------------------------------------------------------------
p = inputParser;
p.addParameter('PreSet',  '');
p.addParameter('Frame',   []);
p.addParameter('c',       []);
p.addParameter('Events',  []);
p.addParameter('Verbose', true);
p.parse(varargin{:});
opt = p.Results;

if ~isfile(dataFile)
    error('load_vsx_rf:noFile', 'Data file not found: %s', dataFile);
end

% ---- 1. locate the preSet ----------------------------------------------
preFile = opt.PreSet;
if isempty(preFile)
    [fdir, base] = fileparts(dataFile);
    % Four naming conventions coexist in the old data: preSet, _preSet,
    % preSave and _preSave.
    cand = {fullfile(fdir, [base 'preSet.mat']), ...
            fullfile(fdir, [base '_preSet.mat']), ...
            fullfile(fdir, [base 'preSave.mat']), ...
            fullfile(fdir, [base '_preSave.mat'])};
    hit = find(cellfun(@isfile, cand), 1);
    if ~isempty(hit)
        preFile = cand{hit};
    else
        % Reference folders legitimately hold one preSet for many
        % acquisitions. Fall back to the nearest one saved at or before this
        % file, and say so: an assumed pairing is not a verified one.
        pl = [dir(fullfile(fdir, '*preSet.mat')); dir(fullfile(fdir, '*preSave.mat'))];
        if isempty(pl)
            error('load_vsx_rf:noPreSet', 'No preSet found in %s', fdir);
        end
        dinfo  = dir(dataFile);
        before = pl([pl.datenum] <= dinfo.datenum);
        if ~isempty(before)
            [~, m] = max([before.datenum]);  pick = before(m);
        else
            [~, m] = min(abs([pl.datenum] - dinfo.datenum));  pick = pl(m);
        end
        preFile = fullfile(pick.folder, pick.name);
        warning('load_vsx_rf:assumedPreSet', ...
            'No exactly-named preSet; using %s', pick.name);
    end
end

Sp = load(preFile, 'preSet');
if ~isfield(Sp, 'preSet')
    error('load_vsx_rf:badPreSet', '%s has no preSet variable', preFile);
end
P = Sp.preSet;
for f = {'Trans', 'Receive', 'TX'}
    if ~isfield(P, f{1})
        error('load_vsx_rf:missing', 'preSet lacks %s', f{1});
    end
end

% ---- 2. frame layout ----------------------------------------------------
% There is no Resource struct in this dataset, so the frame/acquisition
% layout is recovered from the Receive array itself.
nFrames = 1; nAcq = numel(P.Receive);
if isfield(P.Receive, 'framenum'), nFrames = max([P.Receive.framenum]); end
if isfield(P.Receive, 'acqNum'),   nAcq    = max([P.Receive.acqNum]);   end

frame = opt.Frame;
if isempty(frame)
    Sl = load(dataFile, 'RcvLastFrame');
    if isfield(Sl, 'RcvLastFrame') && ~isempty(Sl.RcvLastFrame)
        frame = double(Sl.RcvLastFrame);
    else
        frame = nFrames;
    end
end
frame = max(1, min(frame, nFrames));

% ---- 3. read only the frame we need -------------------------------------
% These files are large (up to ~1.6 GB uncompressed). If RcvData is stored as
% a plain array in a v7.3 file, matfile lets us pull a single frame instead of
% the whole buffer.
mf   = matfile(dataFile);
info = whos(mf, 'RcvData');
if isempty(info)
    error('load_vsx_rf:noRcvData', '%s has no RcvData', dataFile);
end
if strcmp(info.class, 'cell')
    Sd  = load(dataFile, 'RcvData');
    buf = Sd.RcvData;
    if iscell(buf), buf = buf{1}; end
    frame = max(1, min(frame, size(buf, 3)));
    frameData = buf(:, :, frame);
    clear Sd buf
else
    frame = max(1, min(frame, info.size(3)));
    frameData = mf.RcvData(:, :, frame);
end

% ---- 4. element <-> channel mapping -------------------------------------
% RcvData columns are SYSTEM CHANNELS, not elements. Trans.Connector maps
% element index -> channel index. A wrong mapping scrambles the aperture and
% produces a defocused image you would otherwise blame on t0.
if isfield(P.Trans, 'Connector') && ~isempty(P.Trans.Connector)
    chMap = double(P.Trans.Connector(:)).';
else
    chMap = 1:P.Trans.numelements;
end

% ---- 5. sampling --------------------------------------------------------
fc = P.Trans.frequency * 1e6;
R1 = P.Receive(1);
if isfield(R1, 'decimSampleRate') && ~isempty(R1.decimSampleRate)
    fs = R1.decimSampleRate * 1e6;
else
    fs = fc * samplesPerWave(R1);
end

% ---- 6. the two roles of the sound speed --------------------------------
% cAcq is the speed that Verasonics used to define the "wavelength" in which
% every geometric quantity in the setup is expressed. It is a UNIT, frozen at
% acquisition time. c is the speed you assume for time-of-flight when
% reconstructing. Conflating them is the classic source of a systematic depth
% error that looks like a t0 problem.
if isfield(P, 'P') && isfield(P.P, 'speedOfSound')
    cAcq = P.P.speedOfSound;
else
    cAcq = 1540;   % assumption: this dataset does not record it
end
if ~isempty(opt.c), c = opt.c; else, c = cAcq; end
lambdaAcq = cAcq / fc;

% ---- 7. slice the RF ----------------------------------------------------
ev = opt.Events;
if isempty(ev), ev = 1:nAcq; end
nS = R1.endSample - R1.startSample + 1;
rf = zeros(nS, numel(chMap), numel(ev));
for i = 1:numel(ev)
    rIdx = (frame - 1)*nAcq + ev(i);
    if rIdx > numel(P.Receive)
        error('load_vsx_rf:recvIndex', ...
            'Receive index %d exceeds %d.', rIdx, numel(P.Receive));
    end
    idx = P.Receive(rIdx).startSample : P.Receive(rIdx).endSample;
    rf(:, :, i) = double(frameData(idx, chMap));
end
clear frameData

% ---- 8. saturation ------------------------------------------------------
% A clipped channel is a hard nonlinearity: it injects energy at frequencies
% that were never backscattered. Invisible in B-mode, fatal for a spectral
% ACS/BSC fit. Check this before trusting any estimate.
maxAbs      = max(abs(rf(:)));
satFraction = nnz(abs(rf) >= 0.99*maxAbs) / numel(rf);

% ---- 9. receive time axis ----------------------------------------------
% Sample 1 corresponds to the round-trip time of Receive.startDepth, given in
% wavelengths of Trans.frequency:
%       t_start = 2*startDepth/fc        (independent of c, by construction)
% Two fixed lags are removed so that t = 0 is the acoustic emission:
%       TW.peak          lag from pulse start to its peak            [wl]
%       2*lensCorrection two-way delay through lens / matching layer [wl]
% VALIDATE THIS with a point target before trusting any depth.
t_start = 2 * R1.startDepth / fc;
lag = 0;
if isfield(P, 'TW') && ~isempty(P.TW) && isfield(P.TW, 'peak')
    lag = lag + P.TW(1).peak / fc;
end
if isfield(P.Trans, 'lensCorrection') && ~isempty(P.Trans.lensCorrection)
    lag = lag + 2 * P.Trans.lensCorrection / fc;
end
t = t_start - lag + (0:nS-1).' / fs;

% ---- 10. element positions in metres ------------------------------------
EP = P.Trans.ElementPos;
if isfield(P.Trans, 'units') && strcmpi(char(P.Trans.units), 'mm')
    geoScale = 1e-3;
else
    geoScale = lambdaAcq;
end

% ---- 11. transmit description -------------------------------------------
% ASSUMPTION: acquisition event i uses TX(i). Holds for the stock scanline
% and flash scripts. Confirm in your SetUp*.m before relying on it.
tx = struct('delay', {}, 'apod', {}, 'origin', {}, 'focus', {}, ...
            'steer', {}, 'type', {});
for i = 1:numel(ev)
    Ti = P.TX(min(ev(i), numel(P.TX)));
    tx(i).delay  = Ti.Delay(:).' / fc;    % wavelengths -> seconds, exactly
    tx(i).apod   = Ti.Apod(:).';
    tx(i).origin = Ti.Origin(:).' * geoScale;
    tx(i).focus  = Ti.focus * geoScale;   % 0 = plane wave
    if isfield(Ti, 'Steer'), tx(i).steer = Ti.Steer(1); else, tx(i).steer = 0; end
    if Ti.focus == 0, tx(i).type = "planewave"; else, tx(i).type = "focused"; end
end

% ---- 12. assemble -------------------------------------------------------
D = struct();
D.dataFile     = string(dataFile);
D.preSetFile   = string(preFile);
D.rf           = rf;
D.t            = t;
D.fs           = fs;
D.fc           = fc;
D.c            = c;
D.cAcq         = cAcq;
D.cIsAssumed   = ~(isfield(P, 'P') && isfield(P.P, 'speedOfSound'));
D.lambdaAcq    = lambdaAcq;
D.elemX        = EP(:, 1) * geoScale;
D.elemZ        = EP(:, 3) * geoScale;
D.pitch        = median(abs(diff(EP(:, 1)))) * geoScale;
D.tx           = tx;
D.frame        = frame;
D.nFrames      = nFrames;
D.events       = ev(:).';
D.maxAbs       = maxAbs;
D.satFraction  = satFraction;
D.startDepth_m = R1.startDepth * lambdaAcq;
D.endDepth_m   = R1.endDepth   * lambdaAcq;
D.hv           = NaN;
if isfield(P, 'TPC') && ~isempty(P.TPC) && isfield(P.TPC, 'hv')
    D.hv = P.TPC(1).hv;
end
% Keep the TGC untouched: it must be compared byte-for-byte between a sample
% and its reference, so no interpretation happens here.
if     isfield(P, 'TGCparam'), D.tgcRaw = P.TGCparam;
elseif isfield(P, 'TGC'),      D.tgcRaw = P.TGC;
else,                          D.tgcRaw = [];
end
% Small setup structs kept for traceability. RcvData is never kept.
D.setup = struct('Trans', P.Trans, 'Receive', P.Receive(1), 'TX', P.TX);
if isfield(P, 'TW'), D.setup.TW = P.TW; end
if isfield(P, 'P'),  D.setup.P  = P.P;  end

if opt.Verbose
    [~, dn, de] = fileparts(dataFile);
    fprintf('%s%s\n', dn, de);
    fprintf('  rf      : %d samples x %d elements x %d events  (%s)\n', ...
        size(rf,1), size(rf,2), size(rf,3), tx(1).type);
    fprintf('  frame   : %d of %d\n', frame, nFrames);
    fprintf('  fs / fc : %.2f / %.2f MHz\n', fs/1e6, fc/1e6);
    fprintf('  c       : %.0f m/s recon, %.0f m/s at acquisition%s\n', ...
        c, cAcq, ternary(D.cIsAssumed, ' (ASSUMED)', ''));
    fprintf('  depth   : %.1f - %.1f mm   pitch %.3f mm\n', ...
        D.startDepth_m*1e3, D.endDepth_m*1e3, D.pitch*1e3);
    fprintf('  t range : %.2f - %.2f us\n', t(1)*1e6, t(end)*1e6);
    fprintf('  maxAbs  : %g   at >=99%% of max: %.3f %%\n', ...
        maxAbs, 100*satFraction);
    if satFraction > 1e-4
        fprintf('  >> possible clipping: inspect before any spectral estimate\n');
    end
end
end

% =========================================================================
function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
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