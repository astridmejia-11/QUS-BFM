% QUS_DIAG  Where does the systematic ACS underestimation come from?
%
%   Hypothesis: the analysis reaches into a corner of the frequency-depth
%   plane where there is no signal left, only electronic noise. The method
%   measures attenuation as ln(Sp) - ln(Sd), so a noise floor that props up
%   the distal window makes the difference smaller and the attenuation LOWER.
%   That matches the observation: the bias is always negative, never positive.
%
%   Two tests:
%     (1) power spectrum against depth, for sample and reference. The noise
%         floor shows up as a region that stops decaying.
%     (2) a sweep of the upper band limit and the maximum depth. If the
%         hypothesis holds, the bias shrinks as both are pulled back.
%
%   Run from the amejia folder:  qus_diag
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT    = pwd;
TARGET  = 'T7';
C_SOUND = 1542;
Z_BF    = [3e-3 58e-3];

ACS_REF   = 0.5246;               % Wisconsin, UW report 07/09/2025
GAMMA_REF = 1.0019;
ACS_BG    = 0.55;                 % CIRS background certificate

X_INF = -1.8; X_SUP = 1.8;
Z_INF = 1.0;
X_BG  = 1.2;                      % background region: |x| greater than this

FREQ_L   = 4e6;
OVERLAP  = 0.8;
BLOCK_WL = [8 12];
MU_B = 10^3.5; MU_C = 10^2;
RATIO_CUT = 10; REJECT = 0.1; EXTENSION = 3;

% Sweep grid
FREQ_H_LIST = [6e6 7e6 8e6 9e6];
Z_SUP_LIST  = [2.5 3.0 3.5 4.5];
% =========================================================================

% ---- beamform -----------------------------------------------------------
dataOnly = @(d) d(~contains({d.name},'preSet'));
sf  = dataOnly(dir(fullfile(ROOT,'PLANE','*.mat')));
hit = find(contains({sf.name},['_' TARGET '_']) | ...
           contains({sf.name},['_' TARGET '-']), 1);
Ds = load_vsx_rf(fullfile(sf(hit).folder,sf(hit).name),'c',C_SOUND,'Verbose',false);
Bs = das_planewave(Ds,'zLim',Z_BF,'Verbose',false);
clear Ds

cacheFile = fullfile(ROOT,'bf_cache_PW.mat');
if isfile(cacheFile)
    L0 = load(cacheFile,'Bref'); Bref = L0.Bref;
else
    rfl = dataOnly(dir(fullfile(ROOT,'Referencia Plane','*.mat')));
    fprintf('Beamforming %d references (once)...\n', numel(rfl));
    Bref = cell(1,numel(rfl));
    for k = 1:numel(rfl)
        fprintf('  [%d/%d]\n', k, numel(rfl));
        Dk = load_vsx_rf(fullfile(rfl(k).folder,rfl(k).name), ...
                         'c',C_SOUND,'Verbose',false);
        Bref{k} = das_planewave(Dk,'zLim',Z_BF,'Verbose',false);
        clear Dk
    end
    save(cacheFile,'Bref','-v7.3');
end

sam1 = Bs.rf;
xBm  = Bs.x(:)' * 100;
zBm  = Bs.z(:)   * 100;
fs   = Bs.params.fs;
rfRef = zeros(size(sam1,1), size(sam1,2), numel(Bref));
for k = 1:numel(Bref), rfRef(:,:,k) = Bref{k}.rf; end

NptodB = log10(exp(1))*20;

% ---- TEST 1: spectrum against depth ------------------------------------
wl = C_SOUND / ((FREQ_L + 9e6)/2);
bp.xInf = X_INF; bp.xSup = X_SUP; bp.zInf = Z_INF; bp.zSup = 4.5;
bp.overlap = OVERLAP; bp.blocksize = BLOCK_WL * wl;

fprintf('--- spectra by depth ---\n');
evalc('[Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam1,xBm,zBm,fs,bp);');
evalc('[SpRef,SdRef,~,~,~]   = getSpectrum(rfRef,xBm,zBm,fs,bp);');

specSam = 10*log10(squeeze(mean(Sp/2 + Sd/2, 2)));
specRef = 10*log10(squeeze(mean(SpRef/2 + SdRef/2, 2)));
specSam = specSam - max(specSam(:));      % global normalisation, NOT per row:
specRef = specRef - max(specRef(:));      % we want to see the decay itself

figure('Position',[60 60 1100 420]);
subplot(1,2,1);
imagesc(f, z_ACS, specSam, [-60 0]); colorbar; colormap turbo;
hold on; xline(FREQ_L/1e6,'w--','LineWidth',2); xline(9,'w--','LineWidth',2); hold off;
xlim([0 12]); xlabel('f [MHz]'); ylabel('z [cm]');
title('Sample power spectrum by depth [dB]');

subplot(1,2,2);
imagesc(f, z_ACS, specRef, [-60 0]); colorbar; colormap turbo;
hold on; xline(FREQ_L/1e6,'w--','LineWidth',2); xline(9,'w--','LineWidth',2); hold off;
xlim([0 12]); xlabel('f [MHz]'); ylabel('z [cm]');
title('Reference power spectrum by depth [dB]');

% Where does the spectrum stop decaying with depth? That is the noise floor.
figure('Position',[60 60 700 420]);
plot(f, specSam(1,:),   'LineWidth',1.5); hold on;
plot(f, specSam(round(end/2),:), 'LineWidth',1.5);
plot(f, specSam(end,:), 'LineWidth',1.5);
xline(FREQ_L/1e6,'k--'); xline(9,'k--');
hold off; grid on; xlim([0 12]); ylim([-70 0]);
xlabel('f [MHz]'); ylabel('Power [dB]');
legend({sprintf('z = %.1f cm', z_ACS(1)), ...
        sprintf('z = %.1f cm', z_ACS(round(end/2))), ...
        sprintf('z = %.1f cm', z_ACS(end))}, 'Location','southwest');
title('If the deep curve flattens instead of falling, that is the noise floor');

% ---- TEST 2: sweep of band and depth -----------------------------------
fprintf('\n--- background ACS vs band and depth (certified %.2f) ---\n', ACS_BG);
fprintf('%8s', 'zSup\fH');
for jf = 1:numel(FREQ_H_LIST)
    fprintf('%9.0fMHz', FREQ_H_LIST(jf)/1e6);
end
fprintf('\n');

res = nan(numel(Z_SUP_LIST), numel(FREQ_H_LIST));
for iz = 1:numel(Z_SUP_LIST)
    fprintf('%7.1fcm', Z_SUP_LIST(iz));
    for jf = 1:numel(FREQ_H_LIST)
        try
            v = runOne(sam1, rfRef, xBm, zBm, fs, C_SOUND, ...
                X_INF, X_SUP, Z_INF, Z_SUP_LIST(iz), ...
                FREQ_L, FREQ_H_LIST(jf), ACS_REF, GAMMA_REF, ...
                OVERLAP, BLOCK_WL, MU_B, MU_C, ...
                RATIO_CUT, REJECT, EXTENSION, X_BG);
            res(iz,jf) = v;
            fprintf('%12.3f', v);
        catch
            fprintf('%12s', '--');
        end
    end
    fprintf('\n');
end

figure;
imagesc(FREQ_H_LIST/1e6, Z_SUP_LIST, res, [0.3 0.7]);
colorbar; colormap turbo;
xlabel('Upper band limit [MHz]'); ylabel('Max depth [cm]');
title(sprintf('Background ACS (certified %.2f)', ACS_BG));
set(gca,'YDir','normal');

fprintf(['\nIf the estimate rises toward %.2f as the band and the depth are\n' ...
         'pulled back, the bias comes from noise-dominated data, not from\n' ...
         'the method or the phantom.\n'], ACS_BG);


% #########################################################################
function acsBg = runOne(sam1, rfRef, xBm, zBm, fs, c0, ...
    xInf, xSup, zInf, zSup, freqL, freqH, acsRef, gammaRef, ...
    overlap, blockWl, muB, muC, ratioCut, reject, extension, xBG)

wl = c0 / ((freqL + freqH)/2);
NptodB = log10(exp(1))*20;

bp.xInf = xInf; bp.xSup = xSup; bp.zInf = zInf; bp.zSup = zSup;
bp.overlap = overlap; bp.blocksize = blockWl * wl;

evalc('[Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam1,xBm,zBm,fs,bp);');
evalc('[SpRef,SdRef,~,~,~]   = getSpectrum(rfRef,xBm,zBm,fs,bp);');

% Axial separation between the proximal and distal windows, in cm.
% Taken directly from the window geometry rather than from the block spacing:
% the block step is rounded to an integer number of samples, and that rounding
% leaks straight into the attenuation as a few percent of scale error.
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
BSWIFT = reshape(Bn2*NptodB,m,n);

[Xa,~] = meshgrid(x_ACS, z_ACS);
v = BSWIFT(abs(Xa) > xBG);
acsBg = median(v(isfinite(v)));
end
