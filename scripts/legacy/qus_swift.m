% QUS_SWIFT  Runs the lab's RSLD and SWIFT estimators on our own beamformed data.
%
%   This is an ADAPTER, not a reimplementation. The estimation is done by the
%   lab's own functions (getSpectrum, AlterOpti_ADMM, optimAdmmTvTikhonov,
%   optimAdmmWeightedTvTikhonov). All this script does is hand them our data
%   in the units and shapes they expect:
%
%       sam1   RF image                         [nz x nx]
%       rfRef  reference RF images stacked      [nz x nx x nRef]
%       xBm    lateral axis in cm, ROW vector
%       zBm    axial axis in cm, COLUMN vector
%       fs     sampling frequency in Hz
%
%   Requires on the path: the SWIFT folder files, plus load_vsx_rf,
%   das_planewave and das_focused.
%
%   Run from the amejia folder:  qus_swift
%
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT      = pwd;
SCHEME    = 'PW';                 % 'PW' or 'FOC'
TARGET    = 'T7';
C_SOUND   = 1542;                 % CIRS background, certificate
Z_BF      = [3e-3 58e-3];         % depth reconstructed by the beamformer

% Analysis region, in cm (this is what the lab code expects)
X_INF = -1.8;  X_SUP = 1.8;
Z_INF =  1.0;  Z_SUP =  4.5;

FREQ_L = 4e6;  FREQ_H = 9e6;      % analysis band
ACS_REF = 0.53;                   % reference phantom ACS, provisional
GAMMA_REF = 1.00;                 % reference attenuation exponent
OVERLAP = 0.8;
BLOCK_WL = [8 12];                % [lateral axial] in wavelengths

% Regularisation weights. These are the lab's defaults for simulation; real
% phantom data may need a sweep.
MU_RSLD   = 10^3;
MU_B      = 10^3.5;               % SWIFT, ACS term
MU_C      = 10^2;                 % SWIFT, BSC term
RATIO_CUT = 10; REJECT = 0.1; EXTENSION = 3;

% Ground truth for the comparison table
ACS_BG = 0.55;
ACS_TG = 0.95;                    % change with TARGET
X_TG   = 0.6;                     % target region: |x| < this, in cm
Z_TG   = [1.3 2.7];               % target region depth, in cm
X_BG   = 1.2;                     % background region: |x| > this, in cm
% =========================================================================

switch upper(SCHEME)
    case 'PW'
        sampleDir = fullfile(ROOT, 'PLANE');
        refDir    = fullfile(ROOT, 'Referencia Plane');
        bf        = @(D, varargin) das_planewave(D, varargin{:});
        cacheFile = fullfile(ROOT, 'bf_cache_PW.mat');
    case 'FOC'
        sampleDir = fullfile(ROOT, 'TARGETS');
        refDir    = fullfile(ROOT, 'Ref Targets ATS');
        bf        = @(D, varargin) das_focused(D, varargin{:});
        cacheFile = fullfile(ROOT, 'bf_cache_FOC.mat');
end

% ---- beamform -----------------------------------------------------------
dataOnly = @(d) d(~contains({d.name}, 'preSet'));
sf  = dataOnly(dir(fullfile(sampleDir, '*.mat')));
hit = find(contains({sf.name}, ['_' TARGET '_']) | ...
           contains({sf.name}, ['_' TARGET '-']), 1);
if isempty(hit)
    fprintf('  %s\n', sf.name); error('Target %s not found.', TARGET);
end

fprintf('Sample: %s\n\n', sf(hit).name);
Ds = load_vsx_rf(fullfile(sf(hit).folder, sf(hit).name), 'c', C_SOUND);
Bs = bf(Ds, 'zLim', Z_BF);
clear Ds

if isfile(cacheFile)
    fprintf('\nLoading cached references\n');
    L0 = load(cacheFile, 'Bref'); Bref = L0.Bref;
else
    rfl = dataOnly(dir(fullfile(refDir, '*.mat')));
    fprintf('\nBeamforming %d references...\n', numel(rfl));
    Bref = cell(1, numel(rfl));
    for k = 1:numel(rfl)
        fprintf('  [%d/%d]\n', k, numel(rfl));
        Dk = load_vsx_rf(fullfile(rfl(k).folder, rfl(k).name), ...
                         'c', C_SOUND, 'Verbose', false);
        Bref{k} = bf(Dk, 'zLim', Z_BF, 'Verbose', false);
        clear Dk
    end
    save(cacheFile, 'Bref', '-v7.3');
end

% ---- shape the data the way the lab code expects ------------------------
% xBm must be a ROW vector because getSpectrum indexes it as x(1,...).
% zBm must be a COLUMN vector because the later fits use [z_ACS ones(m,1)].
% Both in centimetres.
sam1 = Bs.rf;
xBm  = Bs.x(:)' * 100;
zBm  = Bs.z(:)   * 100;
fs   = Bs.params.fs;

rfRef = zeros(size(sam1,1), size(sam1,2), numel(Bref));
for k = 1:numel(Bref)
    rfRef(:,:,k) = Bref{k}.rf;
end

c0 = C_SOUND;
wl = c0 / ((FREQ_L + FREQ_H)/2);          % wavelength at the band centre [m]
NptodB = log10(exp(1))*20;                % 8.686

blockParams.xInf = X_INF;  blockParams.xSup = X_SUP;
blockParams.zInf = Z_INF;  blockParams.zSup = Z_SUP;
blockParams.overlap   = OVERLAP;
blockParams.blocksize = BLOCK_WL * wl;    % [lateral axial] in metres

fprintf('\nWavelength at band centre: %.3f mm\n', wl*1e3);
fprintf('Block size: %.2f x %.2f mm\n', blockParams.blocksize*1e3);

addpath(genpath("C:\Users\LENOVO\OneDrive\Desktop\amejia\SWIFT_Simple"))

% ---- spectra ------------------------------------------------------------
% getSpectrum returns TWO spectra per block: proximal and distal windows,
% separated axially. The spectral log difference between them is what carries
% the local attenuation.
[Sp,   Sd,   x_ACS, z_ACS, f] = getSpectrum(sam1,  xBm, zBm, fs, blockParams);
[SpRef,SdRef,~,     ~,     ~] = getSpectrum(rfRef, xBm, zBm, fs, blockParams);

% ---- diffraction compensation ------------------------------------------
% L is the axial separation between the proximal and distal windows, in cm.
L = (z_ACS(2) - z_ACS(1))/(1 - OVERLAP)/2;

% Known attenuation of the reference phantom, in Np/cm, as a function of f.
att_ref_map(1,1,:) = ACS_REF * f.^GAMMA_REF / NptodB;

sld    = log(Sp)    - log(Sd);
sldRef = log(SpRef) - log(SdRef);
compensation = sldRef - 4*L*att_ref_map;

range = f > FREQ_L/1e6 & f < FREQ_H/1e6;
b   = sld(:,:,range) - compensation(:,:,range);
ufr = f(range);

[m,n,p] = size(b);
fprintf('System: %d axial x %d lateral blocks, %d frequency bins\n', m, n, p);

% A1 multiplies the ACS unknowns, A2 the BSC unknowns. Solving both together
% is what makes the estimate identifiable: they are coupled in the same
% equation and cannot be separated one at a time.
A1   = kron( 4*L*ufr , speye(m*n) );
A2   = kron( ones(size(ufr)) , speye(m*n) );
mask = ones(m,n,p);
tol  = 1e-3;

% ---- RSLD ---------------------------------------------------------------
fprintf('\nRunning RSLD...\n'); tic
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),MU_RSLD,MU_RSLD,m,n,tol,mask(:));
BRSLD = reshape(Bn*NptodB,m,n);
toc

% ---- SWIFT --------------------------------------------------------------
% Two passes. The first solves with plain TV to get a BSC map; blocks whose
% BSC is extreme are the ones where the assumption of uniform backscatter
% fails, typically at target boundaries. Those get down-weighted, and the
% system is solved again. That is what the "spatially weighted" in SWIFT means.
fprintf('\nRunning SWIFT, first pass...\n'); tic
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),MU_B,MU_C,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

w    = (1-REJECT)*(abs(bscMap) < RATIO_CUT) + REJECT;
wExt = movmin(w, EXTENSION);

W   = repmat(wExt,[1 1 p]);
W   = spdiags(W(:),0,m*n*p,m*n*p);
bw  = W*b(:);
A1w = W*A1;
A2w = W*A2;

fprintf('SWIFT, second pass...\n');
[Bn2,Cn2] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,MU_B,MU_C,m,n,tol,mask(:),wExt);
BSWIFT = reshape(Bn2*NptodB,m,n);
CSWIFT = reshape(Cn2*NptodB,m,n);
toc

% ---- figures ------------------------------------------------------------
bMode = 20*log10(abs(hilbert(sam1)));
bMode = bMode - max(bMode(:));
attRange = [0 1.5];

figure('Position',[50 50 1400 420]);

subplot(1,4,1);
imagesc(xBm, zBm, bMode, [-60 0]); axis image; colormap(gca,gray); colorbar;
xlabel('Lateral [cm]'); ylabel('Axial [cm]');
title(sprintf('B-mode %s %s', TARGET, SCHEME));

subplot(1,4,2);
imagesc(x_ACS, z_ACS, BRSLD, attRange); axis image;
colormap(gca,turbo); colorbar;
xlabel('Lateral [cm]'); title('RSLD, ACS');

subplot(1,4,3);
imagesc(x_ACS, z_ACS, BSWIFT, attRange); axis image;
colormap(gca,turbo); colorbar;
xlabel('Lateral [cm]'); title('SWIFT, ACS');

subplot(1,4,4);
imagesc(x_ACS, z_ACS, CSWIFT, [-10 10]); axis image;
colormap(gca,turbo); colorbar;
xlabel('Lateral [cm]'); title('SWIFT, BSC [dB]');

% ---- comparison table ---------------------------------------------------
[Xa, Za] = meshgrid(x_ACS, z_ACS);
maskBG = abs(Xa) > X_BG;
maskTG = abs(Xa) < X_TG & Za > Z_TG(1) & Za < Z_TG(2);

fprintf('\n%-18s %8s %8s %8s %8s\n', 'method / region', 'median', 'IQR', 'cert', 'NBias');
names = {'RSLD', 'SWIFT'};
maps  = {BRSLD, BSWIFT};
for k = 1:2
    M = maps{k};
    for r = 1:2
        if r == 1
            v = M(maskBG); cert = ACS_BG; lab = 'background';
        else
            v = M(maskTG); cert = ACS_TG; lab = 'target';
        end
        v = v(isfinite(v));
        fprintf('%-8s %-9s %8.3f %8.3f %8.2f %7.1f%%\n', names{k}, lab, ...
            median(v), iqr(v), cert, 100*(median(v)-cert)/cert);
    end
end

vb = BSWIFT(maskBG); vt = BSWIFT(maskTG);
fprintf('\nSWIFT contrast target - background: %+.3f  (certificate %+.3f)\n', ...
    median(vt(isfinite(vt))) - median(vb(isfinite(vb))), ACS_TG - ACS_BG);
