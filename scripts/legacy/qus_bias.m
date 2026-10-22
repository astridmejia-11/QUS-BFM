% QUS_BIAS  Where does the systematic underestimation come from?
%
%   Four hypotheses have already been ruled out experimentally:
%     - the provisional reference constant (the UW report gives 0.5246,
%       essentially the value that was assumed)
%     - a truncated-aperture edge effect (centre and edge agree)
%     - a noise floor (spectra at 1.1, 2.7 and 4.3 cm lie on top of each other)
%     - the extent of the analysis region (flat between 3.0 and 4.5 cm)
%
%   Two remain, and these two experiments separate them.
%
%   A. REFERENCE AGAINST ITSELF
%      One Wisconsin acquisition as the sample, the rest as the reference.
%      Same material, so the correct answer is exactly its certified value,
%      0.5246. Any departure is produced by the estimation chain itself,
%      since there is no real difference between sample and reference.
%
%   B. ATS AS A THIRD SAMPLE
%      A different phantom, with its own certificate (0.5 dB/cm/MHz +/-10%),
%      acquired in the same configuration as the Wisconsin references.
%      If it lands near 0.5, the pipeline measures correctly and the problem
%      is specific to the CIRS. If it is also 20% low, the problem is shared.
%
%   Run from the amejia folder:  qus_bias
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT    = pwd;
C_SOUND = 1542;
Z_BF    = [3e-3 58e-3];

ACS_REF_CERT = 0.5246;            % Wisconsin, UW report 07/09/2025
GAMMA_REF    = 1.0019;
ACS_ATS_CERT = 0.50;              % ATS 539 datasheet, +/-10%
ACS_CIRS_BG  = 0.55;              % CIRS certificate

X_INF = -1.8; X_SUP = 1.8;
Z_INF =  1.0; Z_SUP =  4.5;
X_BG  = 1.2;

FREQ_L = 4e6; FREQ_H = 9e6;
OVERLAP  = 0.8;
BLOCK_WL = [8 12];
FNUM = 2.0; APOD = 'hanning';     % kept fixed: the sweep showed it does not matter
MU_B = 10^3.5; MU_C = 10^2;
RATIO_CUT = 10; REJECT = 0.1; EXTENSION = 3;
% =========================================================================

dataOnly = @(d) d(~contains({d.name},'preSet') & ~contains({d.name},'preSave'));
bfPW = @(D) das_planewave(D,'zLim',Z_BF,'FNumber',FNUM,'Apod',APOD,'Verbose',false);

% ---- beamform every Wisconsin reference --------------------------------
rfl = dataOnly(dir(fullfile(ROOT,'Referencia Plane','*.mat')));
fprintf('Beamforming %d Wisconsin acquisitions...\n', numel(rfl));
Bw = cell(1,numel(rfl));
for k = 1:numel(rfl)
    Dk = load_vsx_rf(fullfile(rfl(k).folder,rfl(k).name), ...
                     'c',C_SOUND,'Verbose',false);
    Bw{k} = bfPW(Dk);
    clear Dk
end

% =========================================================================
% EXPERIMENT A: reference against itself
% =========================================================================
fprintf('\n=========== A. REFERENCE AGAINST ITSELF ===========\n');
fprintf('Correct answer: %.4f dB/cm/MHz for every case.\n\n', ACS_REF_CERT);
fprintf('%-10s %10s %10s %10s\n','sample as','SWIFT','error','NBias');

vA = nan(1,numel(Bw));
for k = 1:numel(Bw)
    others = Bw(setdiff(1:numel(Bw), k));
    A = estimateOne(Bw{k}, others, C_SOUND, X_INF, X_SUP, Z_INF, Z_SUP, ...
        FREQ_L, FREQ_H, ACS_REF_CERT, GAMMA_REF, OVERLAP, BLOCK_WL, ...
        MU_B, MU_C, RATIO_CUT, REJECT, EXTENSION);
    [Xa,~] = meshgrid(A.x_ACS, A.z_ACS);
    v = A.BSWIFT(abs(Xa) > X_BG);
    vA(k) = median(v(isfinite(v)));
    fprintf('W%-9d %10.3f %10.3f %9.1f%%\n', k, vA(k), ...
        vA(k)-ACS_REF_CERT, 100*(vA(k)-ACS_REF_CERT)/ACS_REF_CERT);
end
fprintf('\nmedian %.3f, spread %.3f\n', median(vA), max(vA)-min(vA));

% =========================================================================
% EXPERIMENT B: ATS as a third sample
% =========================================================================
fprintf('\n=========== B. ATS AS A THIRD SAMPLE ===========\n');
atsDir = fullfile(ROOT,'PLANE ATS');
af = dataOnly(dir(fullfile(atsDir,'*.mat')));
if isempty(af)
    fprintf('No files in %s, experiment B skipped.\n', atsDir);
    vB = [];
else
    fprintf('Certificate: %.2f dB/cm/MHz (+/-10%%)\n\n', ACS_ATS_CERT);
    fprintf('%-28s %10s %10s\n','file','SWIFT','NBias');
    vB = nan(1,numel(af));
    for k = 1:numel(af)
        Dk = load_vsx_rf(fullfile(af(k).folder,af(k).name), ...
                         'c',C_SOUND,'Verbose',false);
        Bk = bfPW(Dk); clear Dk
        A = estimateOne(Bk, Bw, C_SOUND, X_INF, X_SUP, Z_INF, Z_SUP, ...
            FREQ_L, FREQ_H, ACS_REF_CERT, GAMMA_REF, OVERLAP, BLOCK_WL, ...
            MU_B, MU_C, RATIO_CUT, REJECT, EXTENSION);
        % The ATS "circle" and "puntitos" acquisitions contain wires and
        % cysts, so take the median over the whole map rather than the edges:
        % the median is robust to those localised structures.
        v = A.BSWIFT(:);
        vB(k) = median(v(isfinite(v)));
        nm = af(k).name; if numel(nm) > 26, nm = nm(1:26); end
        fprintf('%-28s %10.3f %9.1f%%\n', nm, vB(k), ...
            100*(vB(k)-ACS_ATS_CERT)/ACS_ATS_CERT);
    end
end

% =========================================================================
% VERDICT
% =========================================================================
fprintf('\n=================== SUMMARY ===================\n');
fprintf('%-34s %8s %8s %8s\n','case','cert','measured','NBias');
fprintf('%-34s %8.3f %8.3f %7.1f%%\n', 'A. Wisconsin against itself', ...
    ACS_REF_CERT, median(vA), 100*(median(vA)-ACS_REF_CERT)/ACS_REF_CERT);
if ~isempty(vB)
    fprintf('%-34s %8.3f %8.3f %7.1f%%\n', 'B. ATS 539 vs Wisconsin', ...
        ACS_ATS_CERT, median(vB), 100*(median(vB)-ACS_ATS_CERT)/ACS_ATS_CERT);
end
fprintf('%-34s %8.3f %8.3f %7.1f%%\n', 'C. CIRS background vs Wisconsin', ...
    ACS_CIRS_BG, 0.44, 100*(0.44-ACS_CIRS_BG)/ACS_CIRS_BG);
fprintf('   (C is the value already obtained in qus_sweep)\n');

fprintf(['\nHow to read this:\n' ...
  '  A biased  -> the estimation chain itself is biased, since sample and\n' ...
  '               reference are the same material.\n' ...
  '  A clean, B clean, C biased -> the pipeline measures correctly and the\n' ...
  '               discrepancy is specific to the CIRS phantom.\n' ...
  '  A clean, B and C biased    -> something common to both samples, most\n' ...
  '               likely a mismatch against the Wisconsin reference.\n']);

figure;
bar([median(vA), ACS_REF_CERT; ...
     median([vB nan],'omitnan'), ACS_ATS_CERT; ...
     0.44, ACS_CIRS_BG]);
set(gca,'XTickLabel',{'Wisconsin vs itself','ATS 539','CIRS background'});
ylabel('ACS [dB/cm/MHz]'); legend({'measured','certified'},'Location','northwest');
grid on; title('Measured against certified, three materials');


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

% Window separation from the window geometry, not from the block step, which
% is rounded to whole samples and leaks a few percent of scale error into the
% attenuation.
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
end
