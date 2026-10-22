% QUS_ALL  Single final run: every sample, one frozen parameter set.
%
%   Every number reported in the presentation comes out of THIS script, so
%   the whole table is reproducible from one command. Nothing is tuned per
%   case. The only thing that changes between cases is which file is loaded,
%   which reference set it is divided by, and where the analysis region sits.
%
%   Requires on the path: the SWIFT folder, load_vsx_rf, das_planewave,
%   das_focused.
%
%   Run from the amejia folder:  qus_all
%
%   Outputs:
%       qus_all_results.mat   maps + stats for every case
%       qus_all_table.csv     the table for the report
%
%   Plain ASCII on purpose.

clear; clc;

% ======================= FROZEN PARAMETERS ===============================
% These are the conclusions of the study. They do not change per case.
P.fNumber   = 2;            % receive f-number   (module 7: comparable to Merino)
P.apod      = 'hanning';    % receive apodisation
P.muB       = 10^3;         % SWIFT, ACS  term   (module 5: joint sweep)
P.muC       = 10^1;         % SWIFT, BSC  term
P.blockWL   = [8 12];       % [lateral axial] block size, in wavelengths
P.overlap   = 0.8;
P.freqL     = 4e6;
P.freqH     = 9e6;
P.acsRef    = 0.5246;       % Wisconsin, UW report 07/09/2025
P.gammaRef  = 1.0019;       % attenuation exponent, same report
P.ratioCut  = 8;           % SWIFT detector threshold, LEFT UNTOUCHED on purpose
P.reject    = 0.1;
P.extension = 3;
P.zBF       = [3e-3 58e-3]; % depth reconstructed by the beamformer
P.tol       = 1e-3;

ROOT   = pwd;
OUTMAT = fullfile(ROOT, 'qus_all_cut8.mat');
OUTCSV = fullfile(ROOT, 'qus_all_cut8.csv');

NptodB = log10(exp(1))*20;
% Common axial grid. dz defaults to c/(2*fs), so a different speed of sound
% would give a different number of rows and the blocks of sample and
% reference would not sit at the same depths. Fixing it decouples the grid
% from c: the speed still enters where it should, in the delays.
% Value taken from the LOWEST c in use (1515), so no case is undersampled.
P.dz = 1515 / (2*27.7778e6);     % 27.27 um
% =========================================================================


% ============================== CASES ====================================
% EDIT THIS BLOCK. File names must match what is actually on disk.
%
%   scheme  : 'PW' or 'FOC'      -> picks the beamformer and the reference set
%   sample  : file name inside its folder, or '' for the self-reference case
%   c       : speed of sound used for reconstruction [m/s]
%   ana     : analysis region [xInf xSup zInf zSup] in cm
%   region  : where the reported value is measured
%               'all' every block in the analysis region
%               'bg'  |x| > xBg
%               'tg'  |x| < xTg and zTg(1) < z < zTg(2)
%   cert    : reference value to compare against, NaN if none exists

k = 0;

% --- module 8: bias localisation -----------------------------------------
k=k+1; C(k) = mkcase('WIS-self', 'PW', '', 1535.7, ...
        [-1.8 1.8 0.6 4.5], 'all', 0.5246, 'UW report');

% Clean window, verified free of wires and cysts in BOTH acquisitions:
% the wires sit at x = -10 mm with one more at +10 mm, and the cysts of
% the "circle" acquisition are below z = 3 cm.
ATSBOX = [-0.5 0.5 1.2 3.0];

k=k+1; C(k) = mkcase('ATS-circle-1515', 'PW', ...
        'PLANE ATS/RFdata_26_circle-August-2026_16-11-53.mat', 1515, ...
        [-1.8 1.8 0.6 4.5], 'box', 0.50, 'ATS datasheet', ATSBOX);

k=k+1; C(k) = mkcase('ATS-wires-1515', 'PW', ...
        'PLANE ATS/RFdata_26_puntitos-August-2026_16-10-26.mat', 1515, ...
        [-1.8 1.8 0.6 4.5], 'box', 0.50, 'ATS datasheet', ATSBOX);

k=k+1; C(k) = mkcase('ATS-circle-1450', 'PW', ...
    'PLANE ATS/RFdata_26_circle-August-2026_16-11-53.mat', 1450, ...
    [-1.8 1.8 0.6 4.5], 'box', 0.50, 'ATS datasheet', ATSBOX);

k=k+1; C(k) = mkcase('ATS-wires-1450', 'PW', ...
    'PLANE ATS/RFdata_26_puntitos-August-2026_16-10-26.mat', 1450, ...
    [-1.8 1.8 0.6 4.5], 'box', 0.50, 'ATS datasheet', ATSBOX);

k=k+1; C(k) = mkcase('CIRS-bg', 'PW', ...
        'PLANE/RFdata_26_T7_-August-2026_15-53-50.mat', 1542, ...
        [-1.8 1.8 0.6 4.5], 'bg', 0.55, 'CIRS certificate');

% --- module 6: targets, both transmission schemes ------------------------
pw = {'PLANE/RFdata_26_T1_-August-2026_15-53-50.mat'
      'PLANE/RFdata_26_T3_-August-2026_15-53-50.mat'
      'PLANE/RFdata_26_T5_-August-2026_15-53-50.mat'
      'PLANE/RFdata_26_T7_-August-2026_15-53-50.mat'
      'PLANE/RFdata_26_T8_-August-2026_15-53-50.mat'};

foc = {'TARGETS/RFdata_26_T1_F20_fn3-August-2026_15-53-50.mat'
       'TARGETS/RFdata_26_T3_F20_fn3-August-2026_15-47-36.mat'
       'TARGETS/RFdata_26_T5_F20_fn3-August-2026_15-49-37.mat'
       'TARGETS/RFdata_26_T7_F20_fn3-August-2026_15-53-50.mat'
       'TARGETS/RFdata_26_T8_F20_fn3-August-2026_15-53-50.mat'};

tg   = {'T1','T3','T5','T7','T8'};
cert = [0.52 0.74 0.75 0.95 0.95];
for i = 1:numel(tg)
    k=k+1; C(k) = mkcase(['CIRS-' tg{i} '-PW'], 'PW', pw{i}, 1542, ...
            [-1.8 1.8 0.6 4.5], 'tg', cert(i), 'CIRS certificate');
    k=k+1; C(k) = mkcase(['CIRS-' tg{i} '-FOC'], 'FOC', foc{i}, 1542, ...
            [-1.8 1.8 0.6 4.5], 'tg', cert(i), 'CIRS certificate');
end

% --- module 9: in vivo ---------------------------------------------------
% THE ONLY CASES WITH A RESTRICTED ANALYSIS REGION. Total variation
% propagates error from anechoic structures (trachea, vessels) into the
% parenchyma, and the weight detector fires on 17-22 % of blocks here
% against 0.1 % in homogeneous phantom. ROIs drawn on the B-mode by anatomy.
k=k+1; C(k) = mkcase('THY-right', 'FOC', 'Thyroid/Focused/RFdata_26-August-2026_17-12-57.mat', 1540, [0.0 0.7 1.3 2.2], 'all', 1.34, 'Rouyer 2015');

k=k+1; C(k) = mkcase('THY-left',  'FOC', 'Thyroid/Focused/RFdata_26-August-2026_17-15-15.mat', 1540, [0.3 0.9 1.4 2.4], 'all', 1.34, 'Rouyer 2015');
% =========================================================================


% ------------------------- reference sets --------------------------------
REF.PW.dir  = fullfile(ROOT, 'Referencia Plane');
REF.PW.bf   = @das_planewave;
REF.FOC.dir = fullfile(ROOT, 'Ref Targets ATS');
REF.FOC.bf  = @das_focused;

dataOnly = @(d) d(~contains({d.name}, 'preSet') & ~contains({d.name}, 'preSave'));


% ===================== beamform the references once ======================
% Cached: this is the slow part and it is identical for every case that
% shares a transmission scheme.
for s = {'PW','FOC'}
    sch   = s{1};
    cache = fullfile(ROOT, ['bf_cache_' sch '.mat']);
    if isfile(cache)
        fprintf('Cached references, %s\n', sch);
        continue
    end
    rl = dataOnly(dir(fullfile(REF.(sch).dir, '*.mat')));
    fprintf('Beamforming %d references, %s...\n', numel(rl), sch);
    Bref = cell(1, numel(rl));
    for j = 1:numel(rl)
        fprintf('  [%d/%d] %s\n', j, numel(rl), rl(j).name);
        Dj = load_vsx_rf(fullfile(rl(j).folder, rl(j).name), ...
                         'c', 1540, 'Verbose', false);
        Bref{j} = REF.(sch).bf(Dj, 'zLim', P.zBF, 'dz', P.dz, ...
                       'FNumber', P.fNumber, 'Apod', P.apod, 'Verbose', false);
        clear Dj
    end
    save(cache, 'Bref', 'rl', '-v7.3');
end


% ============================= main loop =================================
R = struct([]);

for ic = 1:numel(C)
    cs = C(ic);
    fprintf('\n===== [%d/%d] %s =====\n', ic, numel(C), cs.name);

    L0   = load(fullfile(ROOT, ['bf_cache_' cs.scheme '.mat']));
    Bref = L0.Bref;
    bf   = REF.(cs.scheme).bf;

    % ---- the sample -----------------------------------------------------
    if isempty(cs.sample)
        % Self-reference: acquisition 1 as sample, the rest as references.
        % Ground truth is exactly zero difference, so any deviation is the
        % estimator's own.
        Bs   = Bref{1};
        Bref = Bref(2:end);
        fprintf('Self-reference: 1 sample, %d references\n', numel(Bref));
    else
        f = fullfile(ROOT, cs.sample);
        if ~isfile(f), warning('Missing: %s  -- skipped', cs.sample); continue; end
        Ds = load_vsx_rf(f, 'c', cs.c, 'Verbose', false);
        Bs = bf(Ds, 'zLim', P.zBF, 'dz', P.dz, ...
        'FNumber', P.fNumber, 'Apod', P.apod, 'Verbose', false);
        clear Ds
    end

    % ---- shapes the lab code expects ------------------------------------
    sam1 = Bs.rf;
    xBm  = Bs.x(:)' * 100;      % row,    cm
    zBm  = Bs.z(:)   * 100;     % column, cm
    fs   = Bs.params.fs;
    
    if ~isequal(size(Bref{1}.rf), size(sam1))
        error('Grid mismatch: sample %dx%d, reference %dx%d. Same dz?', ...
            size(sam1,1), size(sam1,2), size(Bref{1}.rf,1), size(Bref{1}.rf,2));
    end

    rfRef = zeros([size(sam1) numel(Bref)]);
    for j = 1:numel(Bref), rfRef(:,:,j) = Bref{j}.rf; end

    wl = cs.c / ((P.freqL + P.freqH)/2);
    bp.xInf = cs.ana(1); bp.xSup = cs.ana(2);
    bp.zInf = cs.ana(3); bp.zSup = cs.ana(4);
    bp.overlap   = P.overlap;
    bp.blocksize = P.blockWL * wl;

    [Sp, Sd, x_ACS, z_ACS, f] = getSpectrum(sam1,  xBm, zBm, fs, bp);
    [SpR,SdR,~,     ~,     ~] = getSpectrum(rfRef, xBm, zBm, fs, bp);

    % ---- axial separation L, CORRECTED ----------------------------------
    % The lab code derived L from the rounded block step, which introduced a
    % 2.9 % scale error. Here it comes straight from the block size in
    % samples: nz samples, half of them between window centres.
    dz_cm = zBm(2) - zBm(1);
    nz    = round(bp.blocksize(2)*100 / dz_cm);
    L     = (nz/2) * dz_cm;

    att_ref_map(1,1,:) = P.acsRef * f.^P.gammaRef / NptodB;

    sld  = log(Sp)  - log(Sd);
    sldR = log(SpR) - log(SdR);
    comp = sldR - 4*L*att_ref_map;

    rng = f > P.freqL/1e6 & f < P.freqH/1e6;
    b   = sld(:,:,rng) - comp(:,:,rng);
    ufr = f(rng);

    [m,n,p] = size(b);
    A1   = kron(4*L*ufr, speye(m*n));
    A2   = kron(ones(size(ufr)), speye(m*n));
    mask = ones(m*n*p,1);

    % ---- SWIFT, two passes ----------------------------------------------
    [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),P.muB,P.muC,m,n,P.tol,mask);
    bscMap = reshape(Cn*NptodB,m,n);

    w    = (1-P.reject)*(abs(bscMap) < P.ratioCut) + P.reject;
    wExt = movmin(w, P.extension);
    pctPenal = 100*mean(wExt(:) < 1);     % audit: how much the detector fires

    W   = repmat(wExt,[1 1 p]);
    W   = spdiags(W(:),0,m*n*p,m*n*p);
    [Bn2,Cn2] = optimAdmmWeightedTvTikhonov(W*A1, W*A2, W*b(:), ...
                    P.muB, P.muC, m, n, P.tol, mask, wExt);

    acsMap = reshape(Bn2*NptodB,m,n);
    bscMap2= reshape(Cn2*NptodB,m,n);

    % ---- measurement region ---------------------------------------------
    [Xa, Za] = meshgrid(x_ACS, z_ACS);
    switch cs.region
        case 'all', sel = true(size(Xa));
        case 'bg',  sel = abs(Xa) > 1.2;
        case 'tg',  sel = abs(Xa) < 0.6 & Za > 1.4 & Za < 2.4;
        case 'box'
            sel = Xa > cs.box(1) & Xa < cs.box(2) & ...
                Za > cs.box(3) & Za < cs.box(4);
    end
    v = acsMap(sel); v = v(isfinite(v));

    % ---- stats ----------------------------------------------------------
    st.median = median(v);
    st.std    = std(v);
    st.p10    = prctile(v,10);
    st.p90    = prctile(v,90);
    st.n      = numel(v);
    st.bias   = 100*(st.median - cs.cert)/cs.cert;

    fprintf('ACS %.3f  sd %.3f  p10-p90 [%.3f %.3f]  cert %.2f  bias %+.1f %%\n', ...
            st.median, st.std, st.p10, st.p90, cs.cert, st.bias);
    fprintf('blocks %d x %d   penalised %.1f %%   L = %.4f cm\n', ...
            m, n, pctPenal, L);

    % ---- store: maps too, they feed the module-8 figure ------------------
    bMode = 20*log10(abs(hilbert(sam1))); bMode = bMode - max(bMode(:));

    R(ic).name     = cs.name;      R(ic).scheme = cs.scheme;
    R(ic).cert     = cs.cert;      R(ic).certSrc= cs.certSrc;
    R(ic).region   = cs.region;    R(ic).ana    = cs.ana;
    R(ic).stats    = st;           R(ic).pctPenal = pctPenal;
    R(ic).acsMap   = acsMap;       R(ic).bscMap = bscMap2;
    R(ic).x_ACS    = x_ACS;        R(ic).z_ACS  = z_ACS;
    R(ic).bMode    = bMode;        R(ic).xBm    = xBm;  R(ic).zBm = zBm;
    R(ic).selMask  = sel;

    clear Bs Bref rfRef sld sldR b A1 A2 W
end

save(OUTMAT, 'R', 'P', '-v7.3');


% ============================== table ====================================
fid = fopen(OUTCSV, 'w');
fprintf(fid, 'case,scheme,region,ACS,sd,p10,p90,nblocks,reference,source,bias_pct,penalised_pct\n');
fprintf('\n\n%-16s %6s %7s %7s %7s %7s %8s %9s\n', ...
        'case','scheme','ACS','sd','p10','p90','cert','bias %');
for i = 1:numel(R)
    if isempty(R(i).name), continue; end
    s = R(i).stats;
    fprintf('%-16s %6s %7.3f %7.3f %7.3f %7.3f %8.2f %+8.1f\n', ...
            R(i).name, R(i).scheme, s.median, s.std, s.p10, s.p90, ...
            R(i).cert, s.bias);
    fprintf(fid, '%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%s,%.2f,%.2f\n', ...
            R(i).name, R(i).scheme, R(i).region, s.median, s.std, ...
            s.p10, s.p90, s.n, R(i).cert, R(i).certSrc, s.bias, R(i).pctPenal);
end
fclose(fid);
fprintf('\nWritten: %s\n        %s\n', OUTMAT, OUTCSV);


% ---------------------------- helper -------------------------------------
function c = mkcase(name, scheme, sample, cSound, ana, region, cert, src, box)
    c.name    = name;
    c.scheme  = scheme;
    c.sample  = sample;
    c.c       = cSound;
    c.ana     = ana;
    c.region  = region;
    c.cert    = cert;
    c.certSrc = src;
    if nargin < 9, box = []; end
    c.box = box;                 % [xMin xMax zMin zMax] in cm, for 'box'
end
