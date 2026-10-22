% QUS_SWEEP  Runs SWIFT and RSLD over all five CIRS targets, both schemes.
%
%   Produces the comparison table requested by the advisors: certified ACS
%   next to the estimated ACS for every target, in both transmission schemes.
%
%   Requires on the path: the SWIFT folder, load_vsx_rf, das_planewave,
%   das_focused. Run from the amejia folder:  qus_sweep
%
%   Beamformed references are cached per scheme. The first run of each scheme
%   is slow; later runs reuse the cache.
%
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT    = pwd;
SCHEMES = {'PW', 'FOC'};
C_SOUND = 1542;                   % CIRS background, certificate
Z_BF    = [3e-3 58e-3];

% Targets and their certified ACS (CIRS Model 049, serial 211270316)
TARGETS  = {'T1',  'T3',  'T5',  'T7',  'T8' };
ACS_CERT = [0.52,  0.74,  0.75,  0.95,  0.95 ];
ACS_BG   = 0.55;

% Wisconsin reference phantom, acoustic characterization report 07/09/2025.
% Specific attenuation fit: 0.5246 * f^1.0019, R2 = 0.9993.
ACS_REF   = 0.5246;
GAMMA_REF = 1.0019;

% Analysis region, in cm
X_INF = -1.8; X_SUP = 1.8;
Z_INF =  1.0; Z_SUP =  3.5;
X_TG  = 0.6;  Z_TG  = [1.3 2.7];   % target region
X_BG  = 1.2;                       % background region: |x| greater than this

FREQ_L = 4e6; FREQ_H = 9e6;
OVERLAP  = 0.8;
BLOCK_WL = [8 12];

MU_RSLD = 10^3;
MU_B    = 10^3.5;
MU_C    = 10^2;
RATIO_CUT = 10; REJECT = 0.1; EXTENSION = 3;
% =========================================================================

R = struct('scheme',{},'target',{},'cert',{}, ...
           'rsldBg',{},'rsldTg',{},'swiftBg',{},'swiftTg',{}, ...
           'swiftBgIqr',{},'swiftTgIqr',{});
maps = {};

for iS = 1:numel(SCHEMES)
    scheme = SCHEMES{iS};
    switch scheme
        case 'PW'
            sampleDir = fullfile(ROOT,'PLANE');
            refDir    = fullfile(ROOT,'Referencia Plane');
            bf        = @(D,varargin) das_planewave(D,varargin{:});
            cacheFile = fullfile(ROOT,'bf_cache_PW.mat');
        case 'FOC'
            sampleDir = fullfile(ROOT,'TARGETS');
            refDir    = fullfile(ROOT,'Ref Targets ATS');
            bf        = @(D,varargin) das_focused(D,varargin{:});
            cacheFile = fullfile(ROOT,'bf_cache_FOC.mat');
    end

    % ---- references -----------------------------------------------------
    dataOnly = @(d) d(~contains({d.name},'preSet'));
    if isfile(cacheFile)
        fprintf('\n[%s] loading cached references\n', scheme);
        L0 = load(cacheFile,'Bref'); Bref = L0.Bref;
    else
        rfl = dataOnly(dir(fullfile(refDir,'*.mat')));
        fprintf('\n[%s] beamforming %d references (slow, once)\n', ...
            scheme, numel(rfl));
        Bref = cell(1,numel(rfl));
        for k = 1:numel(rfl)
            fprintf('  [%d/%d]\n', k, numel(rfl));
            Dk = load_vsx_rf(fullfile(rfl(k).folder,rfl(k).name), ...
                             'c',C_SOUND,'Verbose',false);
            Bref{k} = bf(Dk,'zLim',Z_BF,'Verbose',false);
            clear Dk
        end
        save(cacheFile,'Bref','-v7.3');
    end

    sf = dataOnly(dir(fullfile(sampleDir,'*.mat')));

    % ---- targets --------------------------------------------------------
    for iT = 1:numel(TARGETS)
        tg  = TARGETS{iT};
        hit = find(contains({sf.name},['_' tg '_']) | ...
                   contains({sf.name},['_' tg '-']), 1);
        if isempty(hit)
            fprintf('  [%s] %s NOT FOUND, skipped\n', scheme, tg);
            continue
        end

        fprintf('\n--- %s / %s ---\n', scheme, tg);
        Ds = load_vsx_rf(fullfile(sf(hit).folder,sf(hit).name), ...
                         'c',C_SOUND,'Verbose',false);
        Bs = bf(Ds,'zLim',Z_BF,'Verbose',false);
        clear Ds

        out = estimateOne(Bs, Bref, C_SOUND, ...
            X_INF, X_SUP, Z_INF, Z_SUP, FREQ_L, FREQ_H, ...
            ACS_REF, GAMMA_REF, OVERLAP, BLOCK_WL, ...
            MU_RSLD, MU_B, MU_C, RATIO_CUT, REJECT, EXTENSION);

        % ---- regions ----------------------------------------------------
        [Xa,Za] = meshgrid(out.x_ACS, out.z_ACS);
        mBG = abs(Xa) > X_BG;
        mTG = abs(Xa) < X_TG & Za > Z_TG(1) & Za < Z_TG(2);
        mBG2 = abs(Xa) < 0.6 & Za > 3.2 & Za < 4.3;
        med = @(v) median(v(isfinite(v)));
        spr = @(v) iqr(v(isfinite(v)));

        r = struct();
        r.scheme = scheme;  r.target = tg;  r.cert = ACS_CERT(iT);
        r.rsldBg  = med(out.BRSLD(mBG));   r.rsldTg  = med(out.BRSLD(mTG));
        r.swiftBg = med(out.BSWIFT(mBG));  r.swiftTg = med(out.BSWIFT(mTG));
        r.swiftBgIqr = spr(out.BSWIFT(mBG));
        r.swiftTgIqr = spr(out.BSWIFT(mTG));
        R(end+1) = r; %#ok<SAGROW>

        maps{end+1} = struct('scheme',scheme,'target',tg, ...
            'BSWIFT',out.BSWIFT,'BRSLD',out.BRSLD,'CSWIFT',out.CSWIFT, ...
            'x_ACS',out.x_ACS,'z_ACS',out.z_ACS, ...
            'bMode',out.bMode,'xBm',out.xBm,'zBm',out.zBm); %#ok<SAGROW>

        fprintf('  SWIFT  bg %.3f  target %.3f (cert %.2f)\n', ...
            r.swiftBg, r.swiftTg, r.cert);
    end
end

save(fullfile(ROOT,'qus_sweep_results.mat'),'R','maps','-v7.3');

% ---- comparison table ---------------------------------------------------
fprintf('\n\n================ ACS COMPARISON TABLE ================\n');
fprintf('Reference: Wisconsin phantom, %.4f dB/cm/MHz (UW report 07/09/2025)\n', ...
    ACS_REF);
fprintf('Background certified %.2f dB/cm/MHz\n\n', ACS_BG);

fprintf('%-4s %-7s %7s %8s %8s %8s %9s %9s\n', ...
    'sch','target','cert','SWIFT','IQR','NBias','contrast','cert ctr');
for k = 1:numel(R)
    ctr     = R(k).swiftTg - R(k).swiftBg;
    certCtr = R(k).cert - ACS_BG;
    fprintf('%-4s %-7s %7.2f %8.3f %8.3f %7.1f%% %9.3f %9.3f\n', ...
        R(k).scheme, R(k).target, R(k).cert, R(k).swiftTg, ...
        R(k).swiftTgIqr, 100*(R(k).swiftTg-R(k).cert)/R(k).cert, ctr, certCtr);
end

fprintf('\nBackground estimate per acquisition (certified %.2f):\n', ACS_BG);
for k = 1:numel(R)
    fprintf('  %-4s %-4s  SWIFT %.3f (%+.1f%%)   RSLD %.3f (%+.1f%%)\n', ...
        R(k).scheme, R(k).target, ...
        R(k).swiftBg, 100*(R(k).swiftBg-ACS_BG)/ACS_BG, ...
        R(k).rsldBg,  100*(R(k).rsldBg -ACS_BG)/ACS_BG);
end

% ---- scheme comparison --------------------------------------------------
% This is the project question: does the transmission scheme change the
% estimate? Compare the same target across schemes.
fprintf('\n--- estimate by scheme ---\n');
for iT = 1:numel(TARGETS)
    sel = strcmp({R.target}, TARGETS{iT});
    if nnz(sel) == 2
        v = [R(sel).swiftTg];
        s = {R(sel).scheme};
        fprintf('  %-4s  %s %.3f   %s %.3f   diff %+.3f\n', ...
            TARGETS{iT}, s{1}, v(1), s{2}, v(2), v(2)-v(1));
    end
end

% ---- figure -------------------------------------------------------------
nM = numel(maps);
if nM > 0
    figure('Position',[40 40 240*min(nM,5) 480]);
    tiledlayout(2, min(nM,5), 'TileSpacing','compact', 'Padding','compact');
    for k = 1:min(nM,10)
        nexttile
        imagesc(maps{k}.x_ACS, maps{k}.z_ACS, maps{k}.BSWIFT, [0 1.5]);
        axis image; colormap turbo;
        title(sprintf('%s %s', maps{k}.scheme, maps{k}.target));
        if k == 1, ylabel('Axial [cm]'); end
        xlabel('Lateral [cm]');
    end
    cb = colorbar; cb.Label.String = 'ACS [dB/cm/MHz]';
end


% #########################################################################
%                            LOCAL FUNCTION
% #########################################################################

function out = estimateOne(Bs, Bref, c0, xInf, xSup, zInf, zSup, ...
    freqL, freqH, acsRef, gammaRef, overlap, blockWl, ...
    muRsld, muB, muC, ratioCut, reject, extension)
%ESTIMATEONE  One sample through the lab's RSLD and SWIFT estimators.
%
%   The lab code expects xBm as a ROW vector in cm (it indexes x(1,...)) and
%   zBm as a COLUMN vector in cm (later fits use [z_ACS ones(m,1)]).

sam1 = Bs.rf;
xBm  = Bs.x(:)' * 100;
zBm  = Bs.z(:)   * 100;
fs   = Bs.params.fs;

rfRef = zeros(size(sam1,1), size(sam1,2), numel(Bref));
for k = 1:numel(Bref)
    rfRef(:,:,k) = Bref{k}.rf;
end

wl = c0 / ((freqL + freqH)/2);
NptodB = log10(exp(1))*20;

blockParams.xInf = xInf;  blockParams.xSup = xSup;
blockParams.zInf = zInf;  blockParams.zSup = zSup;
blockParams.overlap   = overlap;
blockParams.blocksize = blockWl * wl;

% Two spectra per block, proximal and distal. Their log difference carries
% the local attenuation.
[Sp,   Sd,   x_ACS, z_ACS, f] = getSpectrum(sam1,  xBm, zBm, fs, blockParams);
[SpRef,SdRef,~,     ~,     ~] = getSpectrum(rfRef, xBm, zBm, fs, blockParams);

L = (z_ACS(2) - z_ACS(1))/(1 - overlap)/2;     % cm
att_ref_map(1,1,:) = acsRef * f.^gammaRef / NptodB;   % Np/cm

sld    = log(Sp)    - log(Sd);
sldRef = log(SpRef) - log(SdRef);
compensation = sldRef - 4*L*att_ref_map;

range = f > freqL/1e6 & f < freqH/1e6;
b   = sld(:,:,range) - compensation(:,:,range);
ufr = f(range);

[m,n,p] = size(b);
A1   = kron( 4*L*ufr , speye(m*n) );
A2   = kron( ones(size(ufr)) , speye(m*n) );
mask = ones(m,n,p);
tol  = 1e-3;

[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muRsld,muRsld,m,n,tol,mask(:));
BRSLD = reshape(Bn*NptodB,m,n);

% SWIFT: first pass gives a BSC map; blocks with extreme BSC are where the
% uniform-backscatter assumption fails, typically at target boundaries. Those
% are down-weighted and the system is solved again.
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);
w      = (1-reject)*(abs(bscMap) < ratioCut) + reject;
wExt   = movmin(w, extension);

W   = repmat(wExt,[1 1 p]);
W   = spdiags(W(:),0,m*n*p,m*n*p);
[Bn2,Cn2] = optimAdmmWeightedTvTikhonov(W*A1, W*A2, W*b(:), ...
    muB, muC, m, n, tol, mask(:), wExt);

out.BRSLD  = BRSLD;
out.BSWIFT = reshape(Bn2*NptodB,m,n);
out.CSWIFT = reshape(Cn2*NptodB,m,n);
out.x_ACS  = x_ACS;
out.z_ACS  = z_ACS;
out.xBm    = xBm;
out.zBm    = zBm;
bMode = 20*log10(abs(hilbert(sam1)));
out.bMode = bMode - max(bMode(:));
end
