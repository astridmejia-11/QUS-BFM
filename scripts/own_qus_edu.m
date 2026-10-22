% own_qus_edu.m - Task 03, part B: own estimator (qus_new pipeline) on
% targets_edu. Baseline uses every Merino parameter and must reproduce the
% SWIFT rows of results/merino/phantoms.xlsx within +-0.01. Then one factor
% at a time is switched to the own value.
%
% Pipeline copied from qus_new.m: getSpectrum (here get_spectrum_opt, same
% code with window/NFFT exposed), reference compensation, first pass
% optimAdmmTvTikhonov, weights w/wExt, second pass
% optimAdmmWeightedTvTikhonov, metric. Data: frame 1 of each file, as
% phantoms.m uses (see Task 01).
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(genpath(fullfile(MERINO_DIR, 'AttenuationUtils')));
addpath(fullfile(DATA_ROOT, 'src'));
addpath(fullfile(DATA_ROOT, 'lab', 'SWIFT_Simple'), '-end');   % only for getSpectrum check
outDir = fullfile(DATA_ROOT, 'results', 'own_edu');
if ~isfolder(outDir), mkdir(outDir); end
NptodB = log10(exp(1))*20;

% ---- data ---------------------------------------------------------------
root = fullfile(DATA_ROOT, 'targets_edu');
tNames = {'T6','T7','T8'};
gtInc  = [0.97 0.95 0.95];  gtBack = 0.55;
cxy    = [1.8 1.9; 1.95 1.95; 1.85 1.9];      % Merino's centres, cm
muCm   = [10^2 10^1 10^1];                    % Merino SWIFT muC per target
S = cell(1,3);
for k = 1:3
    L0 = load(fullfile(root, 'data', [tNames{k} '.mat']), 'RF', 'x', 'z', 'fs', 'fc');
    S{k} = L0.RF(:,:,1);
    if k == 1, x = L0.x*100; z = L0.z*100; fs = L0.fs; fc = L0.fc; end
end
d = dir(fullfile(root, 'referencia', '544', '*.mat'));
REF = zeros([size(S{1}) numel(d)]);
for j = 1:numel(d)
    L0 = load(fullfile(d(j).folder, d(j).name), 'RF');
    REF(:,:,j) = L0.RF(:,:,1);
end
clear L0
fprintf('fs %.2f MHz, fc %.2f MHz, %d refs\n', fs/1e6, fc/1e6, numel(d));

% ---- Merino reference values (Task 01) ---------------------------------
T1 = readtable(fullfile(DATA_ROOT, 'results', 'merino', 'phantoms.xlsx'));
T1 = T1(strcmp(T1.method, 'SWIFT'), :);

% ---- check: get_spectrum_opt defaults == getSpectrum -------------------
bp0 = struct('xInf',0.1,'xSup',3.8,'zInf',0.2,'zSup',3.5,'overlap',0.8, ...
             'blocksize',[8 12]*1540/5e6,'verbose',false);
[a1,a2] = getSpectrum(S{1}, x, z, fs, bp0);
[b1,b2] = get_spectrum_opt(S{1}, x, z, fs, bp0);
fprintf('get_spectrum_opt vs getSpectrum, max |diff|: %g %g\n', ...
    max(abs(a1(:)-b1(:))), max(abs(a2(:)-b2(:))));

% ---- configurations -----------------------------------------------------
base = struct('band',[2.5 7.5]*1e6, 'wlFromBand',false, 'window','hamming', ...
    'nfftExtra',2, 'muB',10^3.5, 'muC','merino', 'pass2','w', 'metric','circle');
cfg = {'0 Merino baseline', base};
c = base; c.band = [4 9]*1e6;                cfg(end+1,:) = {'1a Band 4-9 (block fixed)', c};
c.wlFromBand = true;                          cfg(end+1,:) = {'1b Band 4-9 + block from band', c};
c = base; c.window = 'tukey';                 cfg(end+1,:) = {'2 Tukey 0.25', c};
c = base; c.nfftExtra = 1;                    cfg(end+1,:) = {'3 NFFT +1', c};
c = base; c.muB = 10^3;                       cfg(end+1,:) = {'4 muB 10^3', c};
c = base; c.muC = 10^1;                       cfg(end+1,:) = {'5 muC 10^1', c};
c = base; c.pass2 = 'wExt';                   cfg(end+1,:) = {'6 Second pass wExt', c};
c = base; c.metric = 'rect';                  cfg(end+1,:) = {'7 Median, rectangle', c};
c = struct('band',[4 9]*1e6, 'wlFromBand',true, 'window','tukey', 'nfftExtra',1, ...
    'muB',10^3, 'muC',10^1, 'pass2','wExt', 'metric','rect');
cfg(end+1,:) = {'8 All own', c};

nC = size(cfg,1);
inc = nan(nC,3); bck = nan(nC,3);
for ic = 1:nC
    o = cfg{ic,2};
    for k = 1:3
        muC = o.muC; if ischar(muC), muC = muCm(k); end
        [inc(ic,k), bck(ic,k)] = runOne(S{k}, REF, x, z, fs, o, o.muB, muC, cxy(k,:), NptodB);
    end
    fprintf('%-32s back %.3f %.3f %.3f | inc %.3f %.3f %.3f\n', cfg{ic,1}, bck(ic,:), inc(ic,:));
end

% ---- baseline vs Merino -------------------------------------------------
dInc  = inc(1,:) - T1.meanInc.';
dBack = bck(1,:) - T1.meanBack.';
fprintf('\nBaseline - Merino: inc %+.4f %+.4f %+.4f | back %+.4f %+.4f %+.4f\n', dInc, dBack);
reproduced = all(abs([dInc dBack]) <= 0.01);
fprintf('Reproduced within 0.01: %d\n', reproduced);

% ---- table --------------------------------------------------------------
backM = mean(bck, 2);
R = table(string(cfg(:,1)), backM, inc(:,1), inc(:,2), inc(:,3), ...
    backM - backM(1), inc(:,1)-inc(1,1), inc(:,2)-inc(1,2), inc(:,3)-inc(1,3), ...
    'VariableNames', {'factor','back','T6','T7','T8','dBack','dT6','dT7','dT8'});
disp(R);
writetable(R, fullfile(outDir, 'factors.csv'));
save(fullfile(outDir, 'factors.mat'), 'R', 'inc', 'bck', 'cfg', 'reproduced');

% =========================================================================
function [mInc, mBack] = runOne(sam, ref, x, z, fs, o, muB, muC, cxy, NptodB)
c0 = 1540;
fC = 5e6; if o.wlFromBand, fC = mean(o.band); end
wl = c0 / fC;
bp = struct('xInf',0.1,'xSup',3.8,'zInf',0.2,'zSup',3.5,'overlap',0.8, ...
    'blocksize',[8 12]*wl,'window',o.window,'nfftExtra',o.nfftExtra,'verbose',false);
[Sp,Sd,x_ACS,z_ACS,f] = get_spectrum_opt(sam, x, z, fs, bp);
[SpR,SdR]             = get_spectrum_opt(ref, x, z, fs, bp);

% axial half-length as in qus_new
dz_cm = z(2) - z(1);
nz = round(bp.blocksize(2)*100 / dz_cm);
L  = (nz/2) * dz_cm;

att_ref_map(1,1,:) = 0.53 * f / NptodB;
comp = (log(SpR) - log(SdR)) - 4*L*att_ref_map;
sld  = log(Sp) - log(Sd);
rngf = f > o.band(1)/1e6 & f < o.band(2)/1e6;
b    = sld(:,:,rngf) - comp(:,:,rngf);
ufr  = f(rngf);

[m,n,p] = size(b);
A1 = kron(4*L*ufr, speye(m*n));
A2 = kron(ones(size(ufr)), speye(m*n));
mask = ones(m*n*p,1);
tol = 1e-3;

[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask);
bsc1 = reshape(Cn*NptodB,m,n);
w    = (1-0.1)*(abs(bsc1) < 10) + 0.1;
wExt = movmin(w, 3);
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
if strcmp(o.pass2,'w'), wReg = w; else, wReg = wExt; end
Bn = optimAdmmWeightedTvTikhonov(W*A1, W*A2, W*b(:), muB, muC, m, n, tol, mask, wReg);
acs = reshape(Bn*NptodB,m,n);

switch o.metric
    case 'circle'   % Merino: interpolate to pixel grid, mean in circles
        xr = x(0.1 <= x & x <= 3.8);  zr = z(0.2 <= z & z <= 3.5);
        [X,Z] = meshgrid(x_ACS, z_ACS);  [Xq,Zq] = meshgrid(xr, zr);
        A = interp2(X,Z,acs,Xq,Zq);
        rr = (Xq-cxy(1)).^2 + (Zq-cxy(2)).^2;
        mInc  = mean(A(rr <= 0.85^2), 'omitnan');
        mBack = mean(A(rr >= 1.05^2), 'omitnan');
    case 'rect'     % own: median on ACS grid, rectangle in the inclusion
        [Xa,Za] = meshgrid(x_ACS, z_ACS);
        selT = abs(Xa-cxy(1)) < 0.6 & abs(Za-cxy(2)) < 0.5;
        selB = abs(Xa-cxy(1)) > 1.2;
        vT = acs(selT); vB = acs(selB);
        mInc  = median(vT(isfinite(vT)));
        mBack = median(vB(isfinite(vB)));
end
end
