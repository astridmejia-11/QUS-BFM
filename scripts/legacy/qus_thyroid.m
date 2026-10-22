% QUS_THYROID  ACS estimation on the in vivo thyroid acquisitions.
%
%   The two thyroid acquisitions share configuration with TARGETS and with
%   the focused Wisconsin references: focused, focus at 20 mm, f/2.9, and the
%   same TGC curve. That was established in the metadata audit, and it is what
%   makes the same reference usable here.
%
%   WHAT IS DIFFERENT ABOUT THIS CASE
%   There is no ground truth. Live tissue has no certificate, so NBias cannot
%   be computed. The phantoms validate the pipeline; the thyroid shows that
%   the validated pipeline transfers to real tissue. Two different roles.
%
%   Reported instead: the estimated value, its spread, and whether it falls
%   inside the range published for healthy thyroid (roughly 0.4 to 0.9
%   dB/cm/MHz depending on the study and the frequency band).
%
%   This is also the case where the uniform-backscatter assumption clearly
%   fails: there are lobe, muscle, trachea and vessels in the field. The map
%   is therefore more informative than any single global number, and it is
%   where a regularised estimator earns its place.
%
%   Run from the amejia folder:  qus_thyroid
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT    = pwd;
C_SOUND = 1542;                   % soft tissue, standard value
Z_BF    = [3e-3 58e-3];

ACS_REF   = 0.5246;               % Wisconsin, UW report 07/09/2025
GAMMA_REF = 1.0019;

% Region de parenquima homogeneo, una por adquisicion.
% Evitar musculos infrahioideos por delante, traquea y vasos.
% El prior de variacion total propaga los errores de las zonas anecoicas
% al parenquima, asi que la region debe contener solo glandula.
ROI = [ 0.0  0.7  1.3  2.2;
    0.3  0.9  1.4  2.4 ];
FREQ_L = 4e6;
FREQ_H = 9e6;
OVERLAP  = 0.8;
BLOCK_WL = [8 12];
F_NUM = 2.0; 
APOD = 'hanning';
DZ_COMMON = 1515 / (2*27.7778e6);     % 27.27 um, same as qus_all
MU_B = 1000;
MU_C = 10;
RATIO_CUT = 10; REJECT = 0.1; EXTENSION = 3;

LIT_RANGE = [1.19 1.49];            % published range, healthy thyroid
% =========================================================================

dataOnly = @(d) d(~contains({d.name},'preSet') & ~contains({d.name},'preSave'));
whos Z_BF F_NUM APOD DZ_COMMON
bfFOC = @(D) das_focused(D, 'zLim', Z_BF, 'dz', DZ_COMMON, 'FNumber', F_NUM, 'Apod', APOD, 'Verbose', false);

% ---- references ---------------------------------------------------------
% The thyroid was acquired with the focused configuration, so it needs the
% focused references, not the plane wave ones.
cacheFile = fullfile(ROOT,'bf_cache_FOC.mat');
if isfile(cacheFile)
    fprintf('Loading cached focused references\n');
    L0 = load(cacheFile,'Bref'); Bref = L0.Bref;
else
    rfl = dataOnly(dir(fullfile(ROOT,'Ref Targets ATS','*.mat')));
    fprintf('Beamforming %d focused references (slow, once)...\n', numel(rfl));
    Bref = cell(1,numel(rfl));
    for k = 1:numel(rfl)
        fprintf('  [%d/%d]\n', k, numel(rfl));
        Dk = load_vsx_rf(fullfile(rfl(k).folder,rfl(k).name), ...
                         'c',C_SOUND,'Verbose',false);
        Bref{k} = bfFOC(Dk); clear Dk
    end
    save(cacheFile,'Bref','-v7.3');
end

% ---- thyroid acquisitions ----------------------------------------------
tf = dataOnly(dir(fullfile(ROOT,'Thyroid','Focused','*.mat')));
if isempty(tf)
    tf = dataOnly(dir(fullfile(ROOT,'Thyroid','**','*.mat')));
end
if isempty(tf)
    error('No thyroid files found under %s', fullfile(ROOT,'Thyroid'));
end
fprintf('Found %d thyroid acquisition(s)\n\n', numel(tf));

res = struct('name',{},'med',{},'std',{},'p10',{},'p90',{});
figure('Position',[40 40 640*numel(tf) 420]);

for k = 1:numel(tf)
    fprintf('--- %s ---\n', tf(k).name);
    Dk = load_vsx_rf(fullfile(tf(k).folder,tf(k).name),'c',C_SOUND,'Verbose',false);
    Bk = bfFOC(Dk); clear Dk
   xInf = ROI(k,1); xSup = ROI(k,2);
    zInf = ROI(k,3); zSup = ROI(k,4);

    A = estimateOne(Bk, Bref, C_SOUND, xInf, xSup, zInf, zSup, ...
        FREQ_L, FREQ_H, ACS_REF, GAMMA_REF, OVERLAP, BLOCK_WL, ...
        MU_B, MU_C, RATIO_CUT, REJECT, EXTENSION);

    v = A.BSWIFT(:); v = v(isfinite(v));
    r.name = string(tf(k).name);
    r.med = median(v); r.std = std(v);
    r.p10 = prctile(v,10); r.p90 = prctile(v,90);
    res(end+1) = r; %#ok<SAGROW>

    fprintf('  ACS median %.3f, IQR %.3f, p10-p90 %.3f-%.3f dB/cm/MHz\n\n', ...
        r.med, r.std, r.p10, r.p90);

    subplot(1, 2*numel(tf), 2*k-1);
    imagesc(A.xBm, A.zBm, A.bMode, [-60 0]);
    axis image; colormap(gca,gray); colorbar;
    xlabel('Lateral [cm]'); ylabel('Axial [cm]');
    title(sprintf('B-mode, thyroid %d', k));
    hold on;
    plot([xInf xSup xSup xInf xInf], [zInf zInf zSup zSup zInf], ...
     'y--', 'LineWidth', 1.5);
    hold off;

    subplot(1, 2*numel(tf), 2*k);
    imagesc(A.x_ACS, A.z_ACS, A.BSWIFT, [0 3]);
    axis image; colormap(gca,turbo); colorbar;
    xlabel('Lateral [cm]');
    title(sprintf('SWIFT ACS, thyroid %d', k));
end

% ---- summary ------------------------------------------------------------
fprintf('\n============ THYROID, IN VIVO ============\n');
fprintf('No certificate exists for live tissue, so no NBias is computed.\n');
fprintf('Published range for healthy thyroid: %.1f to %.1f dB/cm/MHz\n\n', ...
    LIT_RANGE(1), LIT_RANGE(2));
fprintf('%-40s %8s %8s %14s\n','acquisition','median','STD','p10-p90');
for k = 1:numel(res)
    nm = char(res(k).name); if numel(nm) > 38, nm = nm(1:38); end
    fprintf('%-40s %8.3f %8.3f %6.2f-%.2f\n', nm, res(k).med, res(k).std, ...
        res(k).p10, res(k).p90);
end

if numel(res) == 2
    fprintf('\nDifference between the two acquisitions: %.3f dB/cm/MHz\n', ...
        abs(res(1).med - res(2).med));
    fprintf(['Both are the same subject, so this difference is a measure of\n' ...
             'repeatability under probe repositioning.\n']);
end

allMed = [res.med];
inRange = allMed >= LIT_RANGE(1) & allMed <= LIT_RANGE(2);
fprintf('\n%d of %d acquisitions fall inside the published range.\n', ...
    nnz(inRange), numel(allMed));

fprintf('\n--- parametros usados ---\n');
fprintf('  muB %.0f, muC %.0f, ratioCutOff %g, reject %.1f, extension %d\n', ...
    MU_B, MU_C, RATIO_CUT, REJECT, EXTENSION);
fprintf('  ventana %g x %g long. de onda, solapamiento %.0f%%\n', ...
    BLOCK_WL(1), BLOCK_WL(2), 100*OVERLAP);
fprintf('  banda %.1f - %.1f MHz, f/%.1f %s\n', ...
    FREQ_L/1e6, FREQ_H/1e6, F_NUM, APOD);
fprintf('  referencia Wisconsin %.4f dB/cm/MHz, c = %d m/s\n', ACS_REF, C_SOUND);
for k = 1:size(ROI,1)
    fprintf('  ROI %d: x %.1f a %.1f cm, z %.1f a %.1f cm\n', k, ROI(k,:));
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
