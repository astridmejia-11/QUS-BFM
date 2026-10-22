% QUS_MUB  Does the regularisation compress the dynamic range?
%
%   T1 certifies 0.52 and T8 certifies 0.95, a spread of 0.43. At muB = 1000
%   they measure 0.684 and 0.655: a spread of 0.03. Two candidate causes:
%
%     (a) total-variation regularisation pulls neighbouring values together,
%         so raising muB should shrink the spread and lowering it should open
%         the spread back up;
%
%     (b) backscatter contamination. T1 is +11.0 dB and T8 is -9.0 dB, so the
%         hyperechoic target is pushed up and the hypoechoic one down. That
%         has nothing to do with muB and would not move with it.
%
%   This script sweeps muB with everything else frozen. Spectra are computed
%   once per target, since they do not depend on the weights; only the solve
%   repeats.
%
%   Run from the amejia folder, after qus_all has built bf_cache_PW.mat.
%
%   Plain ASCII on purpose.

clear; clc;

% ============================= SWEEP =====================================
MUB_LIST = [1e2 3.16e2 1e3 3.16e3 1e4];   % muB values to test
MUC      = 1e1;                           % frozen at the final value
% =========================================================================

% ---- frozen parameters, identical to qus_all ----------------------------
P.fNumber = 2;    P.apod = 'hanning';
P.blockWL = [8 12];  P.overlap = 0.8;
P.freqL = 4e6;    P.freqH = 9e6;
P.acsRef = 0.5246;   P.gammaRef = 1.0019;
P.ratioCut = 10;  P.reject = 0.1;  P.extension = 3;
P.zBF = [3e-3 58e-3];   P.dz = 1515/(2*27.7778e6);   P.tol = 1e-3;
C_SOUND = 1542;

ANA = [-1.8 1.8 0.6 4.5];        % analysis region [xInf xSup zInf zSup], cm
TG  = struct('x', 0.6, 'z', [1.4 2.4]);   % measurement region inside target

ROOT   = pwd;
NptodB = log10(exp(1))*20;

CASES = { 'T1', 'PLANE/RFdata_26_T1_-August-2026_15-53-50.mat', 0.52, +11.0
          'T8', 'PLANE/RFdata_26_T8_-August-2026_15-53-50.mat', 0.95,  -9.0 };

% ---- references ---------------------------------------------------------
L0   = load(fullfile(ROOT,'bf_cache_PW.mat'), 'Bref');
Bref = L0.Bref;
fprintf('References: %d\n', numel(Bref));

res = nan(size(CASES,1), numel(MUB_LIST));

for ic = 1:size(CASES,1)
    name = CASES{ic,1};
    fprintf('\n===== %s =====\n', name);

    Ds = load_vsx_rf(fullfile(ROOT, CASES{ic,2}), 'c', C_SOUND, 'Verbose', false);
    Bs = das_planewave(Ds, 'zLim', P.zBF, 'dz', P.dz, ...
                       'FNumber', P.fNumber, 'Apod', P.apod, 'Verbose', false);
    clear Ds

    sam1 = Bs.rf;
    xBm  = Bs.x(:)' * 100;
    zBm  = Bs.z(:)   * 100;
    fs   = Bs.params.fs;

    rfRef = zeros([size(sam1) numel(Bref)]);
    for j = 1:numel(Bref), rfRef(:,:,j) = Bref{j}.rf; end

    wl = C_SOUND / ((P.freqL + P.freqH)/2);
    bp.xInf = ANA(1); bp.xSup = ANA(2); bp.zInf = ANA(3); bp.zSup = ANA(4);
    bp.overlap = P.overlap;  bp.blocksize = P.blockWL * wl;

    % ---- spectra: computed ONCE, they do not depend on the weights ------
    [Sp, Sd, x_ACS, z_ACS, f] = getSpectrum(sam1,  xBm, zBm, fs, bp);
    [SpR,SdR,~,     ~,     ~] = getSpectrum(rfRef, xBm, zBm, fs, bp);

    dz_cm = zBm(2) - zBm(1);
    nz    = round(bp.blocksize(2)*100 / dz_cm);
    L     = (nz/2) * dz_cm;

    att_ref_map(1,1,:) = P.acsRef * f.^P.gammaRef / NptodB;
    comp = (log(SpR) - log(SdR)) - 4*L*att_ref_map;
    sld  =  log(Sp)  - log(Sd);

    rng = f > P.freqL/1e6 & f < P.freqH/1e6;
    b   = sld(:,:,rng) - comp(:,:,rng);
    ufr = f(rng);

    [m,n,p] = size(b);
    A1   = kron(4*L*ufr, speye(m*n));
    A2   = kron(ones(size(ufr)), speye(m*n));
    mask = ones(m*n*p,1);

    [Xa, Za] = meshgrid(x_ACS, z_ACS);
    sel = abs(Xa) < TG.x & Za > TG.z(1) & Za < TG.z(2);

    % ---- sweep ----------------------------------------------------------
    fprintf('%10s %8s %8s %8s\n', 'muB', 'ACS', 'sd', 'penal%');
    for im = 1:numel(MUB_LIST)
        muB = MUB_LIST(im);

        [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,MUC,m,n,P.tol,mask);
        bscMap = reshape(Cn*NptodB,m,n);

        w    = (1-P.reject)*(abs(bscMap) < P.ratioCut) + P.reject;
        wExt = movmin(w, P.extension);
        pct  = 100*mean(wExt(:) < 1);

        W = repmat(wExt,[1 1 p]);
        W = spdiags(W(:),0,m*n*p,m*n*p);
        [Bn2,~] = optimAdmmWeightedTvTikhonov(W*A1, W*A2, W*b(:), ...
                      muB, MUC, m, n, P.tol, mask, wExt);
        acsMap = reshape(Bn2*NptodB,m,n);

        v = acsMap(sel); v = v(isfinite(v));
        res(ic,im) = median(v);
        fprintf('%10.0f %8.3f %8.3f %8.1f\n', muB, median(v), std(v), pct);
    end

    clear Bs rfRef sld b A1 A2 W
end

% ============================ verdict ====================================
certSpread = CASES{2,3} - CASES{1,3};       % 0.95 - 0.52 = 0.43

fprintf('\n\n%10s %8s %8s %10s %10s\n', ...
        'muB', 'T1', 'T8', 'spread', '% of cert');
for im = 1:numel(MUB_LIST)
    sp = res(2,im) - res(1,im);
    fprintf('%10.0f %8.3f %8.3f %10.3f %9.1f%%\n', ...
            MUB_LIST(im), res(1,im), res(2,im), sp, 100*sp/certSpread);
end
fprintf('\nCertified spread: %.3f (T8 %.2f - T1 %.2f)\n', ...
        certSpread, CASES{2,3}, CASES{1,3});

fprintf(['\nHOW TO READ THIS\n' ...
  '  Spread widens as muB drops -> the regularisation was compressing the\n' ...
  '  range, and the compression is a smoothing artefact.\n' ...
  '  Spread stays flat          -> muB is not the cause. The remaining\n' ...
  '  candidate is backscatter contamination: T1 is +11.0 dB and T8 is\n' ...
  '  -9.0 dB, and the bias follows that sign.\n']);

figure;
plot(MUB_LIST, res(1,:), 'o-', 'LineWidth', 1.5); hold on;
plot(MUB_LIST, res(2,:), 's-', 'LineWidth', 1.5);
yline(CASES{1,3}, '--'); yline(CASES{2,3}, '--');
set(gca,'XScale','log'); grid on;
xlabel('\mu_B'); ylabel('ACS [dB/cm/MHz]');
legend('T1 (+11.0 dB)','T8 (-9.0 dB)','certified','Location','best');
title('Dynamic range vs regularisation weight');
