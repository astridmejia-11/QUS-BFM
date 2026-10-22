% calib_edu.m - Task 04, phase 2: calibrate SWIFT on targets_edu with a
% calibration / test split fixed beforehand (results/PREDICCIONES.md):
%   calibration = T1 T3 T5 T7,  test = T2 T4 T6 T8.
%
% Estimator: Merino's phantoms.m pipeline, reproduced step by step (crop,
% 8x12 wavelength blocks at c0 = 1540 and 5 MHz, 80 % overlap, Hamming,
% NFFT = 2^(nextpow2(nz/2)+2), band 2.5-7.5 MHz, reference 544 at 0.53
% averaged over its files, frame 1 everywhere, first pass
% optimAdmmTvTikhonov, w = (1-reject)*(|bsc|<cut)+reject, wExt = movmin,
% second pass optimAdmmWeightedTvTikhonov regularised with w, metric = mean
% over circular masks r <= 0.85 / r >= 1.05 cm after interp2 to the pixel
% grid). Centres measured on each B-mode (merino/centres_edu.m).
%
% Check: Merino's own parameters and hand-picked centres on T6 T7 T8 must
% give the Task 01 SWIFT table.
%
% Grid: muB 10^{2.5..4.5}, muC 10^{0..3 step 0.5}, ratioCutOff {6 8 10}.
% Objective (fixed beforehand): mean |relative bias| over inclusion and
% background of the calibration targets; tie-break by mean std.
% The frozen combination is evaluated ONCE on the test targets, next to
% Merino's original parameters on the same targets.
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(genpath(fullfile(MERINO_DIR, 'AttenuationUtils')));
outDir = fullfile(DATA_ROOT, 'results', 'calib_edu');
NptodB = log10(exp(1))*20;

CAL  = [1 3 5 7];  TEST = [2 4 6 8];
GT   = [0.52 0.55 0.74 0.81 0.75 0.97 0.95 0.95];  GTB = 0.55;
LMUB = 2.5:0.5:4.5;  LMUC = 0:0.5:3;  CUTS = [6 8 10];
% Merino's original SWIFT parameters per target (his cases: T6 -> 1,
% T7 -> 2, T8 -> 3; others take case 2, as in phantoms_dt2809)
merinoLmuC = [1 1 1 1 1 2 1 1];

P = prepare(DATA_ROOT, NptodB);
C = readtable(fullfile(outDir, 'centres.csv'));
cen = [C.c1x C.c1z];

% ---- check against Task 01 ---------------------------------------------
T01 = readtable(fullfile(DATA_ROOT, 'results', 'merino', 'phantoms.xlsx'));
T01 = T01(strcmp(T01.method, 'SWIFT'), :);
hand = [1.8 1.9; 1.95 1.95; 1.85 1.9];
for k = 1:3
    t = 5 + k;
    r = swiftEval(P.t{t}, P, 3.5, merinoLmuC(t), 10, hand(k,:), NptodB);
    fprintf('check T%d: inc %.4f (Task01 %.4f)  back %.4f (Task01 %.4f)\n', t, ...
        r.inc, T01.meanInc(k), r.back, T01.meanBack(k));
    assert(abs(r.inc - T01.meanInc(k)) < 1e-6 && abs(r.back - T01.meanBack(k)) < 1e-6);
end
fprintf('check passed: pipeline reproduces Task 01 exactly\n\n');

% ---- calibration grid, calibration targets only -------------------------
rows = {};
for ib = 1:numel(LMUB)
    for ic = 1:numel(LMUC)
        first = cell(1, 8);
        for t = CAL
            first{t} = firstPass(P.t{t}, P, LMUB(ib), LMUC(ic));
        end
        for iu = 1:numel(CUTS)
            relI = []; relB = []; sds = []; v = [];
            for t = CAL
                r = swiftEval(P.t{t}, P, LMUB(ib), LMUC(ic), CUTS(iu), cen(t,:), NptodB, first{t});
                relI(end+1) = abs(r.inc - GT(t))/GT(t);  %#ok<AGROW>
                relB(end+1) = abs(r.back - GTB)/GTB;     %#ok<AGROW>
                sds(end+1)  = mean([r.stdInc r.stdBack]); %#ok<AGROW>
                v = [v r.inc r.back]; %#ok<AGROW>
            end
            obj = mean([relI relB]);
            rows(end+1,:) = {LMUB(ib), LMUC(ic), CUTS(iu), obj, mean(sds), mean(relI), mean(relB), v}; %#ok<AGROW>
        end
        fprintf('muB 10^%.1f muC 10^%.1f done\n', LMUB(ib), LMUC(ic));
    end
end
G = cell2table(rows(:,1:7), 'VariableNames', ...
    {'lmuB','lmuC','cut','objective','meanStd','relBiasInc','relBiasBack'});
V = cell2mat(rows(:,8));
for j = 1:numel(CAL)
    G.(sprintf('incT%d', CAL(j)))  = V(:, 2*j-1);
    G.(sprintf('backT%d', CAL(j))) = V(:, 2*j);
end
G = sortrows(G, {'objective','meanStd'});
writetable(G, fullfile(outDir, 'grid_calibration.csv'));
best = G(1,:);
fprintf('\nFROZEN: muB 10^%.1f  muC 10^%.1f  cut %d  | objective %.4f  std %.4f\n', ...
    best.lmuB, best.lmuC, best.cut, best.objective, best.meanStd);
disp(G(1:10, 1:7));

% ---- single evaluation: frozen vs Merino, all eight targets --------------
E = table();
for t = 1:8
    rf = swiftEval(P.t{t}, P, best.lmuB, best.lmuC, best.cut, cen(t,:), NptodB);
    rm = swiftEval(P.t{t}, P, 3.5, merinoLmuC(t), 10, cen(t,:), NptodB);
    set = "calibration"; if ismember(t, TEST), set = "test"; end
    E = [E; table(string(sprintf('T%d', t)), set, GT(t), ...
        rf.inc, rf.stdInc, rf.back, rf.stdBack, ...
        rm.inc, rm.stdInc, rm.back, rm.stdBack, ...
        'VariableNames', {'target','set','gt','incFrozen','sdIncFrozen', ...
        'backFrozen','sdBackFrozen','incMerino','sdIncMerino','backMerino','sdBackMerino'})]; %#ok<AGROW>
    M(t) = struct('frozen', rf, 'merino', rm); %#ok<SAGROW>
end
rel = @(e, g) abs(e - g)./g;
for s = ["calibration" "test"]
    sel = E.set == s;
    eF = mean([rel(E.incFrozen(sel), E.gt(sel)); rel(E.backFrozen(sel), GTB)]);
    eM = mean([rel(E.incMerino(sel), E.gt(sel)); rel(E.backMerino(sel), GTB)]);
    fprintf('%-12s mean |rel bias|: frozen %.4f   Merino %.4f\n', s, eF, eM);
end
disp(E);
writetable(E, fullfile(outDir, 'evaluation.csv'));
save(fullfile(outDir, 'calib_edu.mat'), 'G', 'E', 'best', 'M', 'CAL', 'TEST');

% =========================================================================
function P = prepare(DATA_ROOT, NptodB)
% Everything that does not depend on the regularisation: crop, blocks,
% reference spectra (shared grid), sample spectra and the system b, A1, A2.
root = fullfile(DATA_ROOT, 'targets_edu');
freq_L = 2.5e6; freq_H = 7.5e6; freq_C = (freq_L + freq_H)/2;
overlap_pc = 0.8; ratio_zx = 12/8; blocksize = 8;
x_inf = 0.1; x_sup = 3.8; z_inf = 0.2; z_sup = 3.5;

S = load(fullfile(root, 'data', 'T1.mat'), 'x', 'z', 'fs');
dx = S.x(2)-S.x(1); dz = S.z(2)-S.z(1); fs = S.fs;
x = S.x*1e2; z = S.z*1e2;
ind_x = x_inf <= x & x <= x_sup;  ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x); z = z(ind_z);
wl = 1540/freq_C;
wx = round(blocksize*wl*(1-overlap_pc)/dx);  nx = round(blocksize*wl/dx);
x0 = 1:wx:length(x)-nx;  x_ACS = x(1, x0+round(nx/2));  n = length(x0);
wz = round(blocksize*wl*(1-overlap_pc)/dz * ratio_zx);
nz = 2*round(blocksize*wl/dz /2 * ratio_zx);
L = (nz/2)*dz*100;
z0p = 1:wz:length(z)-nz;  z0d = z0p + nz/2;  z_ACS = z(z0p + nz/2);  m = length(z0p);
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;
rang = band > freq_L & band < freq_H;
f = band(rang)*1e-6;  p = length(f);
windowing = hamming(nz/2)*ones(1,nx);

spec = @(sam) blockSpectra(sam, x0, z0p, z0d, nz, nx, windowing, NFFT, rang, m, n, p);

refFiles = dir(fullfile(root, 'referencia', '544', '*.mat'));
Sp = zeros(m,n,p,numel(refFiles)); Sd = Sp;
for i = 1:numel(refFiles)
    R = load(fullfile(refFiles(i).folder, refFiles(i).name), 'RF');
    [Sp(:,:,:,i), Sd(:,:,:,i)] = spec(R.RF(ind_z, ind_x, 1));
end
att_ref_map = repmat(reshape(0.53*f/NptodB, 1, 1, p), m, n);
compensation = (log(mean(Sp,4)) - log(mean(Sd,4))) - 4*L*att_ref_map;

P = struct('x', x, 'z', z, 'x_ACS', x_ACS, 'z_ACS', z_ACS, 'm', m, 'n', n, 'p', p, ...
    'A1', kron(4*L*f, speye(m*n)), 'A2', kron(ones(size(f)), speye(m*n)));
P.t = cell(1, 8);
for t = 1:8
    R = load(fullfile(root, 'data', sprintf('T%d.mat', t)), 'RF');
    [sp, sd] = spec(R.RF(ind_z, ind_x, 1));
    b = (log(sp) - log(sd)) - compensation;
    P.t{t} = struct('b', b(:));
end
end

function [Sp, Sd] = blockSpectra(sam, x0, z0p, z0d, nz, nx, windowing, NFFT, rang, m, n, p)
Sp = zeros(m,n,p); Sd = zeros(m,n,p);
for jj = 1:n
    for ii = 1:m
        xw = x0(jj); zp = z0p(ii); zd = z0d(ii);
        [tp,~] = spectra(sam(zp:zp+nz/2-1, xw:xw+nx-1), windowing, 0, nz/2, NFFT);
        [td,~] = spectra(sam(zd:zd+nz/2-1, xw:xw+nx-1), windowing, 0, nz/2, NFFT);
        Sp(ii,jj,:) = tp(rang); Sd(ii,jj,:) = td(rang);
    end
end
end

function Cn = firstPass(T, P, lmuB, lmuC)
mask = ones(P.m*P.n*P.p, 1);
[~, Cn] = optimAdmmTvTikhonov(P.A1, P.A2, T.b, 10^lmuB, 10^lmuC, P.m, P.n, 1e-3, mask);
end

function r = swiftEval(T, P, lmuB, lmuC, cut, cxy, NptodB, Cn)
if nargin < 8, Cn = firstPass(T, P, lmuB, lmuC); end
m = P.m; n = P.n; p = P.p;
mask = ones(m*n*p, 1);
bscMap = reshape(Cn*NptodB, m, n);
w = (1-0.1)*(abs(bscMap) < cut) + 0.1;
wExt = movmin(w, 3);
W = repmat(wExt, [1 1 p]);
W = spdiags(W(:), 0, m*n*p, m*n*p);
Bn = optimAdmmWeightedTvTikhonov(W*P.A1, W*P.A2, W*T.b, 10^lmuB, 10^lmuC, m, n, 1e-3, mask, w);
acs = reshape(Bn*NptodB, m, n);
[X, Z] = meshgrid(P.x_ACS, P.z_ACS);  [Xq, Zq] = meshgrid(P.x, P.z);
A = interp2(X, Z, acs, Xq, Zq);
rr = (Xq-cxy(1)).^2 + (Zq-cxy(2)).^2;
inc = rr <= 0.85^2;  back = rr >= 1.05^2;
r = struct('inc', mean(A(inc),'omitnan'), 'stdInc', std(A(inc),'omitnan'), ...
    'back', mean(A(back),'omitnan'), 'stdBack', std(A(back),'omitnan'), ...
    'acs', acs, 'bsc', bscMap, 'w', w);
end
