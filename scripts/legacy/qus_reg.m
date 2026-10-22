% QUS_REG  Regularisation parameters of SWIFT: weight map and mu sweep.
%
%   Two things the advisors asked to verify:
%
%   1. THE WEIGHT MAP. SWIFT runs a first pass, looks at the resulting BSC
%      map, and down-weights the blocks whose |BSC| exceeds ratioCutOff,
%      because those are where the uniform-backscatter assumption fails.
%      With ratioCutOff = 10 dB and a target whose certified contrast is
%      +10.1 dB, the threshold sits exactly on the real contrast: the weight
%      map may be discarding valid blocks. This has to be looked at, not
%      assumed.
%
%   2. THE mu VALUES. muB = 10^3.5 and muC = 10^2 were taken from the lab's
%      simulation script without adjusting them to real phantom data. muC in
%      particular is suspect: the BSC map comes out nearly flat when the
%      certificate says +10.1 dB for T7.
%
%   Run from the amejia folder:  qus_reg
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT    = pwd;
TARGET  = 'T7';
C_SOUND = 1542;
Z_BF    = [3e-3 58e-3];

ACS_REF = 0.5246; GAMMA_REF = 1.0019;
ACS_BG  = 0.55;   ACS_TG = 0.95;   DBSC_TG = 10.1;   % certificado CIRS

X_INF = -1.8; X_SUP = 1.8;
Z_INF =  0.6; Z_SUP =  5.0;        % region ampliada
X_TG  = 0.6;  Z_TG  = [1.3 2.7];
X_BG  = 1.2;

FREQ_L = 4e6; FREQ_H = 9e6;
OVERLAP = 0.8; BLOCK_WL = [8 12];
FNUM = 2.0; APOD = 'hanning';

MU_B_REF = 10^3.5; MU_C_REF = 10^2;
MU_B_LIST = 10.^[2.5 3.0 3.5];
MU_C_LIST = 10.^[-1 0 1 2];

RATIO_CUT = 10; REJECT = 0.1; EXTENSION = 3;
NREF = 4;
% =========================================================================

dataOnly = @(d) d(~contains({d.name},'preSet') & ~contains({d.name},'preSave'));
sf  = dataOnly(dir(fullfile(ROOT,'PLANE','*.mat')));
hit = find(contains({sf.name},['_' TARGET '_']) | contains({sf.name},['_' TARGET '-']),1);
rfl = dataOnly(dir(fullfile(ROOT,'Referencia Plane','*.mat')));
rfl = rfl(1:min(NREF,numel(rfl)));

fprintf('Beamformando...\n');
Ds = load_vsx_rf(fullfile(sf(hit).folder,sf(hit).name),'c',C_SOUND,'Verbose',false);
Bs = das_planewave(Ds,'zLim',Z_BF,'FNumber',FNUM,'Apod',APOD,'Verbose',false);
clear Ds
Bref = cell(1,numel(rfl));
for k = 1:numel(rfl)
    Dk = load_vsx_rf(fullfile(rfl(k).folder,rfl(k).name),'c',C_SOUND,'Verbose',false);
    Bref{k} = das_planewave(Dk,'zLim',Z_BF,'FNumber',FNUM,'Apod',APOD,'Verbose',false);
    clear Dk
end

% ---- build the linear system ONCE ---------------------------------------
% The system does not depend on mu: only the solver does. Building it once
% makes the sweep cheap and guarantees every configuration sees identical data.
[S, x_ACS, z_ACS, bMode, xBm, zBm] = buildSystem(Bs, Bref, C_SOUND, ...
    X_INF, X_SUP, Z_INF, Z_SUP, FREQ_L, FREQ_H, ACS_REF, GAMMA_REF, ...
    OVERLAP, BLOCK_WL);
fprintf('Sistema: %d x %d bloques, %d frecuencias\n\n', S.m, S.n, S.p);

[Xa, Za] = meshgrid(x_ACS, z_ACS);
mBG = abs(Xa) > X_BG;
mTG = abs(Xa) < X_TG & Za > Z_TG(1) & Za < Z_TG(2);
stat = @(v) deal(mean(v(isfinite(v))), std(v(isfinite(v))));

% =========================================================================
% PART 1: the weight map
% =========================================================================
fprintf('=== 1. MAPA DE PESOS (muB %.0f, muC %.0f, cutoff %g dB) ===\n', ...
    MU_B_REF, MU_C_REF, RATIO_CUT);

[~, Cn] = optimAdmmTvTikhonov(S.A1,S.A2,S.b,MU_B_REF,MU_C_REF,S.m,S.n,S.tol,S.mask);
bsc1 = reshape(Cn*S.NptodB, S.m, S.n);
w    = (1-REJECT)*(abs(bsc1) < RATIO_CUT) + REJECT;
wExt = movmin(w, EXTENSION);

fracRej = mean(wExt(:) < 1);
fracRejTG = mean(wExt(mTG) < 1);
fprintf('  bloques castigados: %.1f%% del total, %.1f%% dentro del target\n', ...
    100*fracRej, 100*fracRejTG);
fprintf('  rango del BSC de la 1a pasada: %.1f a %.1f dB\n', ...
    min(bsc1(:)), max(bsc1(:)));
fprintf('  contraste certificado del target: %+.1f dB\n\n', DBSC_TG);

figure('Position',[40 40 1250 380]);
subplot(1,3,1);
imagesc(xBm, zBm, bMode, [-60 0]); axis image; colormap(gca,gray); colorbar;
hold on;
plot([X_INF X_SUP X_SUP X_INF X_INF],[Z_INF Z_INF Z_SUP Z_SUP Z_INF],'y--','LineWidth',1.2);
plot([-X_TG X_TG X_TG -X_TG -X_TG],[Z_TG(1) Z_TG(1) Z_TG(2) Z_TG(2) Z_TG(1)],'c-','LineWidth',1.4);
plot([-X_SUP -X_BG -X_BG -X_SUP -X_SUP],[Z_INF Z_INF Z_SUP Z_SUP Z_INF],'g-','LineWidth',1.4);
plot([X_BG X_SUP X_SUP X_BG X_BG],[Z_INF Z_INF Z_SUP Z_SUP Z_INF],'g-','LineWidth',1.4);
hold off;
xlabel('Lateral [cm]'); ylabel('Axial [cm]');
title('B-mode: amarillo = analisis, cyan = target, verde = fondo');

subplot(1,3,2);
imagesc(x_ACS, z_ACS, bsc1, [-15 15]); axis image; colormap(gca,turbo); colorbar;
xlabel('Lateral [cm]'); title(sprintf('BSC 1a pasada [dB], cert %+.1f', DBSC_TG));

subplot(1,3,3);
imagesc(x_ACS, z_ACS, wExt, [0 1]); axis image; colormap(gca,parula); colorbar;
xlabel('Lateral [cm]');
title(sprintf('Mapa de pesos: %.0f%% castigado', 100*fracRej));

% =========================================================================
% PART 2: mu sweep
% =========================================================================
fprintf('\n=== BARRIDO CONJUNTO muB x muC ===\n');
fprintf('%8s %8s %8s %7s %8s %7s %9s %9s\n', ...
    'muB','muC','fondo','std','target','std','contr ACS','contr BSC');

best = []; bestScore = Inf; RES = [];
for i = 1:numel(MU_B_LIST)
    for j = 1:numel(MU_C_LIST)
        [acs,bsc] = solveSwift(S, MU_B_LIST(i), MU_C_LIST(j), ...
            RATIO_CUT, REJECT, EXTENSION);
        [mb,sb] = stat(acs(mBG));
        [mt,st] = stat(acs(mTG));
        cb = mean(bsc(mTG),'omitnan') - mean(bsc(mBG),'omitnan');

        fprintf('%8.0f %8.2f %8.3f %7.3f %8.3f %7.3f %9.3f %9.2f\n', ...
            MU_B_LIST(i), MU_C_LIST(j), mb, sb, mt, st, mt-mb, cb);
        RES(end+1,:) = [MU_B_LIST(i) MU_C_LIST(j) mb sb mt st cb]; %#ok<SAGROW>
        % error combinado contra los dos valores certificados
        score = abs(mb-ACS_BG)/ACS_BG + abs(mt-ACS_TG)/ACS_TG;
        if score < bestScore
            bestScore = score;
            best = [MU_B_LIST(i) MU_C_LIST(j) mb mt];
        end
    end
end

fprintf('\nMenor error combinado: muB = %.0f, muC = %.2f\n', best(1), best(2));
fprintf('  fondo %.3f (cert %.2f), target %.3f (cert %.2f)\n', ...
    best(3), ACS_BG, best(4), ACS_TG);

figure('Position',[40 40 1050 400]);
subplot(1,2,1);
for i = 1:numel(MU_B_LIST)
    sel = RES(:,1) == MU_B_LIST(i);
    errorbar(log10(RES(sel,2)), RES(sel,3), RES(sel,4), 'o-','LineWidth',1.4);
    hold on;
end
yline(ACS_BG,'k-','fondo cert 0.55','LineWidth',1.2);
hold off; grid on; xlabel('log_{10}(muC)'); ylabel('ACS fondo');
legend(arrayfun(@(v) sprintf('muB=10^{%.1f}',log10(v)), MU_B_LIST, 'uni',0));
title('Fondo');

subplot(1,2,2);
for i = 1:numel(MU_B_LIST)
    sel = RES(:,1) == MU_B_LIST(i);
    errorbar(log10(RES(sel,2)), RES(sel,5), RES(sel,6), 'o-','LineWidth',1.4);
    hold on;
end
yline(ACS_TG,'k-','target cert 0.95','LineWidth',1.2);
hold off; grid on; xlabel('log_{10}(muC)'); ylabel('ACS target');
title('Target');

% #########################################################################
function [S, x_ACS, z_ACS, bMode, xBm, zBm] = buildSystem(Bs, Bref, c0, ...
    xInf, xSup, zInf, zSup, freqL, freqH, acsRef, gammaRef, overlap, blockWl)

sam1 = Bs.rf;  xBm = Bs.x(:)'*100;  zBm = Bs.z(:)*100;  fs = Bs.params.fs;
rfRef = zeros(size(sam1,1), size(sam1,2), numel(Bref));
for k = 1:numel(Bref), rfRef(:,:,k) = Bref{k}.rf; end

wl = c0/((freqL+freqH)/2);
S.NptodB = log10(exp(1))*20;
bp.xInf=xInf; bp.xSup=xSup; bp.zInf=zInf; bp.zSup=zSup;
bp.overlap=overlap; bp.blocksize=blockWl*wl;

evalc('[Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam1,xBm,zBm,fs,bp);');
evalc('[SpRef,SdRef,~,~,~]   = getSpectrum(rfRef,xBm,zBm,fs,bp);');

dz_cm = zBm(2)-zBm(1);
nz = 2*round(bp.blocksize(2)/(dz_cm/100)/2);
L  = (nz/2)*dz_cm;

att_ref_map(1,1,:) = acsRef*f.^gammaRef/S.NptodB;
sld = log(Sp)-log(Sd);  sldRef = log(SpRef)-log(SdRef);
compensation = sldRef - 4*L*att_ref_map;
range = f>freqL/1e6 & f<freqH/1e6;
b = sld(:,:,range) - compensation(:,:,range);
ufr = f(range);

[S.m,S.n,S.p] = size(b);
S.A1 = kron(4*L*ufr, speye(S.m*S.n));
S.A2 = kron(ones(size(ufr)), speye(S.m*S.n));
S.b  = b(:);
S.mask = ones(S.m*S.n*S.p,1);
S.tol = 1e-3;

bMode = 20*log10(abs(hilbert(sam1)));  bMode = bMode - max(bMode(:));
end

% #########################################################################
function [acs, bsc] = solveSwift(S, muB, muC, ratioCut, reject, extension)
[~,Cn] = optimAdmmTvTikhonov(S.A1,S.A2,S.b,muB,muC,S.m,S.n,S.tol,S.mask);
bscMap = reshape(Cn*S.NptodB,S.m,S.n);
w    = (1-reject)*(abs(bscMap) < ratioCut) + reject;
wExt = movmin(w, extension);
W = repmat(wExt,[1 1 S.p]);
W = spdiags(W(:),0,S.m*S.n*S.p,S.m*S.n*S.p);
[Bn2,Cn2] = optimAdmmWeightedTvTikhonov(W*S.A1, W*S.A2, W*S.b, ...
    muB, muC, S.m, S.n, S.tol, S.mask, wExt);
acs = reshape(Bn2*S.NptodB,S.m,S.n);
bsc = reshape(Cn2*S.NptodB,S.m,S.n);
end
