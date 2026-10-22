% QUS_BFSWEEP  Does the beamformer change the ACS estimate?
%
%   The project question, as stated by the advisors: evaluate the differences
%   between BEAMFORMERS, independently of the acquisition type. Receive
%   f-number and receive apodisation are beamformer parameters: they are
%   chosen at reconstruction time and have nothing to do with how the data
%   were transmitted. Sweeping them over one fixed raw acquisition isolates
%   the beamformer as the only variable.
%
%   Hypothesis: the reference phantom method cancels beam effects as long as
%   they are identical in sample and reference. If that holds, the B-mode
%   should change a lot across the sweep while the ACS stays put. Either
%   outcome is a result.
%
%   CRITICAL: sample and reference are re-beamformed together for every
%   configuration. Comparing a sample at f/1 against a reference at f/2 does
%   not measure tissue, it measures the mismatch.
%
%   Run from the amejia folder:  qus_bfsweep
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT    = pwd;
SCHEME  = 'PW';
TARGET  = 'T7';
C_SOUND = 1542;
Z_BF    = [3e-3 58e-3];

% The sweep: beamformer parameters only
FNUMS = [1.0 1.5 2.0 3.0];
APODS = {'hanning', 'boxcar'};

ACS_REF   = 0.5246;
GAMMA_REF = 1.0019;
ACS_BG    = 0.55;
ACS_TG    = 0.95;

X_INF = -1.8; X_SUP = 1.8;
Z_INF =  1.0; Z_SUP =  4.5;
X_TG  = 0.6;  Z_TG  = [1.3 2.7];
X_BG  = 1.2;

FREQ_L = 4e6; FREQ_H = 9e6;
OVERLAP  = 0.8;
BLOCK_WL = [8 12];
MU_B = 10^3.5; MU_C = 10^2;
RATIO_CUT = 10; REJECT = 0.1; EXTENSION = 3;

NREF = 4;      % references used per configuration; lower = faster sweep
% =========================================================================

switch upper(SCHEME)
    case 'PW'
        sampleDir = fullfile(ROOT,'PLANE');
        refDir    = fullfile(ROOT,'Referencia Plane');
        bf        = @(D,varargin) das_planewave(D,varargin{:});
    case 'FOC'
        sampleDir = fullfile(ROOT,'TARGETS');
        refDir    = fullfile(ROOT,'Ref Targets ATS');
        bf        = @(D,varargin) das_focused(D,varargin{:});
end

dataOnly = @(d) d(~contains({d.name},'preSet'));
sf  = dataOnly(dir(fullfile(sampleDir,'*.mat')));
hit = find(contains({sf.name},['_' TARGET '_']) | ...
           contains({sf.name},['_' TARGET '-']), 1);
rfl = dataOnly(dir(fullfile(refDir,'*.mat')));
rfl = rfl(1:min(NREF,numel(rfl)));

% ---- load the RAW data ONCE --------------------------------------------
% This is what makes the comparison clean: the same samples go through every
% beamformer configuration. Nothing upstream changes.
fprintf('Loading raw data (once)...\n');
Ds = load_vsx_rf(fullfile(sf(hit).folder,sf(hit).name),'c',C_SOUND,'Verbose',false);
Dr = cell(1,numel(rfl));
for k = 1:numel(rfl)
    Dr{k} = load_vsx_rf(fullfile(rfl(k).folder,rfl(k).name), ...
                        'c',C_SOUND,'Verbose',false);
end
fprintf('Sample: %s\nReferences: %d\n\n', sf(hit).name, numel(Dr));

% ---- sweep --------------------------------------------------------------
nCfg = numel(FNUMS)*numel(APODS);
Rs = struct('fnum',{},'apod',{},'acsBg',{},'iqrBg',{},'acsTg',{}, ...
            'iqrTg',{},'contrast',{},'res_mm',{},'bMode',{},'BSWIFT',{}, ...
            'x_ACS',{},'z_ACS',{},'xBm',{},'zBm',{});

for ia = 1:numel(APODS)
    for jf = 1:numel(FNUMS)
        fn = FNUMS(jf); ap = APODS{ia};
        fprintf('--- f/%.1f, %s ---\n', fn, ap);

        Bs = bf(Ds, 'zLim',Z_BF, 'FNumber',fn, 'Apod',ap, 'Verbose',false);
        Bref = cell(1,numel(Dr));
        for k = 1:numel(Dr)
            Bref{k} = bf(Dr{k}, 'zLim',Z_BF, 'FNumber',fn, 'Apod',ap, ...
                         'Verbose',false);
        end

        out = estimateOne(Bs, Bref, C_SOUND, X_INF, X_SUP, Z_INF, Z_SUP, ...
            FREQ_L, FREQ_H, ACS_REF, GAMMA_REF, OVERLAP, BLOCK_WL, ...
            MU_B, MU_C, RATIO_CUT, REJECT, EXTENSION);

        [Xa,Za] = meshgrid(out.x_ACS, out.z_ACS);
        mBG = abs(Xa) > X_BG;
        mTG = abs(Xa) < X_TG & Za > Z_TG(1) & Za < Z_TG(2);
        med = @(v) median(v(isfinite(v)));
        spr = @(v) std(v(isfinite(v)));

        r.fnum = fn; r.apod = string(ap);
        r.acsBg = med(out.BSWIFT(mBG)); r.iqrBg = spr(out.BSWIFT(mBG));
        r.acsTg = med(out.BSWIFT(mTG)); r.iqrTg = spr(out.BSWIFT(mTG));
        r.contrast = r.acsTg - r.acsBg;

        % Lateral resolution proxy: width of the lateral autocorrelation of
        % the envelope at -6 dB, measured in background. It is the number that
        % should move a lot across the sweep while the ACS does not.
        r.res_mm = lateralRes(Bs);

        r.bMode = out.bMode; r.BSWIFT = out.BSWIFT;
        r.x_ACS = out.x_ACS; r.z_ACS = out.z_ACS;
        r.xBm = out.xBm; r.zBm = out.zBm;
        Rs(end+1) = r; %#ok<SAGROW>

        fprintf('   res %.2f mm | bg %.3f | target %.3f | contrast %.3f\n\n', ...
            r.res_mm, r.acsBg, r.acsTg, r.contrast);
    end
end

save(fullfile(ROOT,'qus_bfsweep_results.mat'),'Rs','-v7.3');

% ---- table --------------------------------------------------------------
fprintf('\n============ BEAMFORMER SWEEP: %s %s ============\n', SCHEME, TARGET);
fprintf('Same raw data throughout. Only the beamformer changes.\n\n');
fprintf('%-9s %6s %8s %9s %8s %9s %10s\n', ...
    'apod','f/#','res[mm]','bg ACS','IQR','target','contrast');
for k = 1:numel(Rs)
    fprintf('%-9s %6.1f %8.2f %9.3f %8.3f %9.3f %10.3f\n', ...
        Rs(k).apod, Rs(k).fnum, Rs(k).res_mm, Rs(k).acsBg, Rs(k).iqrBg, ...
        Rs(k).acsTg, Rs(k).contrast);
end

bg = [Rs.acsBg]; ct = [Rs.contrast]; rs = [Rs.res_mm];
fprintf('\nAcross the whole sweep:\n');
fprintf('  lateral resolution : %.2f to %.2f mm  (factor %.1f)\n', ...
    min(rs), max(rs), max(rs)/min(rs));
fprintf('  background ACS     : %.3f to %.3f     (spread %.3f, %.1f%%)\n', ...
    min(bg), max(bg), max(bg)-min(bg), 100*(max(bg)-min(bg))/mean(bg));
fprintf('  contrast           : %.3f to %.3f     (spread %.3f)\n', ...
    min(ct), max(ct), max(ct)-min(ct));
fprintf(['\nIf the resolution changes by a large factor while the ACS barely\n' ...
         'moves, the reference phantom method is insensitive to the receive\n' ...
         'beamformer, which is the expected behaviour and worth stating.\n']);

% ---- figures ------------------------------------------------------------
figure('Position',[40 40 1150 420]);
subplot(1,2,1);
mk = {'o-','s--'};
for ia = 1:numel(APODS)
    sel = [Rs.apod] == string(APODS{ia});
    errorbar(FNUMS, [Rs(sel).acsBg], [Rs(sel).iqrBg]/2, mk{ia}, 'LineWidth',1.4);
    hold on;
end
yline(ACS_BG,'k-','certified 0.55','LineWidth',1.2);
hold off; grid on; ylim([0 0.9]);
xlabel('Receive f-number'); ylabel('Background ACS [dB/cm/MHz]');
legend([APODS, {'certified'}], 'Location','southeast');
title('ACS against beamformer settings');

subplot(1,2,2);
for ia = 1:numel(APODS)
    sel = [Rs.apod] == string(APODS{ia});
    plot(FNUMS, [Rs(sel).res_mm], mk{ia}, 'LineWidth',1.4); hold on;
end
hold off; grid on;
xlabel('Receive f-number'); ylabel('Lateral resolution [mm]');
legend(APODS, 'Location','northwest');
title('Lateral resolution against beamformer settings');

figure('Position',[40 40 260*numel(Rs) 300]);
tiledlayout(1,numel(Rs),'TileSpacing','compact','Padding','compact');
for k = 1:numel(Rs)
    nexttile
    imagesc(Rs(k).x_ACS, Rs(k).z_ACS, Rs(k).BSWIFT, [0 1.5]);
    axis image; colormap turbo;
    title(sprintf('f/%.1f %s', Rs(k).fnum, Rs(k).apod), 'FontSize',9);
    xlabel('x [cm]');
end


% #########################################################################
function res_mm = lateralRes(Bs)
%LATERALRES  Lateral autocorrelation width of the envelope at -6 dB, in mm.
%   A proxy for lateral resolution: in fully developed speckle the
%   autocorrelation width tracks the width of the point spread function.
env = Bs.env(:, :);
zSel = Bs.z > 15e-3 & Bs.z < 40e-3;
E = env(zSel, :);
E = E - mean(E, 2);

maxLag = 20;
c = zeros(1, maxLag+1);
for L = 0:maxLag
    a = E(:, 1:end-L);
    b = E(:, 1+L:end);
    c(L+1) = mean(a(:).*b(:));
end
c = c / c(1);

half = find(c < 10^(-6/20), 1);
if isempty(half)
    res_mm = NaN;
else
    % linear interpolation between the two lags straddling the threshold
    L0 = half - 2; L1 = half - 1;
    t  = (c(L0+1) - 10^(-6/20)) / (c(L0+1) - c(L1+1));
    res_mm = 2 * (L0 + t) * (Bs.x(2)-Bs.x(1)) * 1e3;
end
end

% #########################################################################
function out = estimateOne(Bs, Bref, c0, xInf, xSup, zInf, zSup, ...
    freqL, freqH, acsRef, gammaRef, overlap, blockWl, ...
    muB, muC, ratioCut, reject, extension)

sam1 = Bs.rf;
xBm  = Bs.x(:)' * 100;
zBm  = Bs.z(:)   * 100;
fs   = Bs.params.fs;
rfRef = zeros(size(sam1,1), size(sam1,2), numel(Bref));
for k = 1:numel(Bref), rfRef(:,:,k) = Bref{k}.rf; end

wl = c0 / ((freqL + freqH)/2);
NptodB = log10(exp(1))*20;
bp.xInf = xInf; bp.xSup = xSup; bp.zInf = zInf; bp.zSup = zSup;
bp.overlap = overlap; bp.blocksize = blockWl * wl;

evalc('[Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam1,xBm,zBm,fs,bp);');
evalc('[SpRef,SdRef,~,~,~]   = getSpectrum(rfRef,xBm,zBm,fs,bp);');

% Window separation taken from the window geometry, not from the block step,
% which is rounded to whole samples and leaks a few percent of scale error
% straight into the attenuation.
dz_cm = zBm(2) - zBm(1);
nz    = 2*round(bp.blocksize(2)/(dz_cm/100)/2);
L     = (nz/2) * dz_cm;

att_ref_map(1,1,:) = acsRef * f.^gammaRef / NptodB;
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

[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);
w    = (1-reject)*(abs(bscMap) < ratioCut) + reject;
wExt = movmin(w, extension);
W    = repmat(wExt,[1 1 p]);
W    = spdiags(W(:),0,m*n*p,m*n*p);
[Bn2,~] = optimAdmmWeightedTvTikhonov(W*A1, W*A2, W*b(:), ...
    muB, muC, m, n, tol, mask(:), wExt);

out.BSWIFT = reshape(Bn2*NptodB,m,n);
out.x_ACS = x_ACS; out.z_ACS = z_ACS;
out.xBm = xBm;     out.zBm = zBm;
bMode = 20*log10(abs(hilbert(sam1)));
out.bMode = bMode - max(bMode(:));
end
