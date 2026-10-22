% zoom_diag.m - Task 02: compare lab and own (das_focused) B-modes of T7_F
% in physical coordinates and measure inclusion centre and diameter with
% axial and lateral profiles.
%
% Own reconstruction uses the same parameters as qus_new.m.
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(fullfile(DATA_ROOT, 'src'));
base   = fullfile(DATA_ROOT, 'DataTargets2809', 'Samples');
outDir = fullfile(DATA_ROOT, 'results', 'zoom');
if ~isfolder(outDir), mkdir(outDir); end

C  = 1540;
FS = 31.25e6;

% ---- lab ---------------------------------------------------------------
L = load(fullfile(base, 'BeamformedData', 'T7_F.mat'));
lab.x   = L.x(:).';  lab.z = L.z(:);
lab.env = abs(hilbert(L.rf));
lab.name = 'Lab';

% ---- own ---------------------------------------------------------------
D = load_vsx_rf(fullfile(base, 'ChannelData', 'T7_F.mat'), 'c', C, 'Verbose', true);
B = das_focused(D, 'zLim', [3e-3 60e-3], 'dz', C/(2*FS), ...
                'FNumber', 2, 'Apod', 'hanning', 'Verbose', true);
own.x   = B.x(:).';  own.z = B.z(:);
own.env = B.env;
own.name = 'Own (das\_focused)';
fprintf('own t(1) = %.4f us -> depth of sample 1 = %.3f mm\n', D.t(1)*1e6, D.t(1)*C/2*1e3);
fprintf('own fs used = %.4f MHz, fc = %.4f MHz, dz = %.3f um\n', ...
    B.params.fs/1e6, B.params.fc/1e6, B.params.dz*1e6);

% ---- measurements ------------------------------------------------------
imgs = {lab, own};
for k = 1:2
    imgs{k} = measureInclusion(imgs{k}, [0 19e-3]);
    I = imgs{k};
    fprintf('\n%s\n  size %d x %d\n  x %.3f .. %.3f mm (dx %.4f)\n  z %.3f .. %.3f mm (dz %.4f)\n', ...
        I.name, size(I.env,1), size(I.env,2), I.x(1)*1e3, I.x(end)*1e3, ...
        (I.x(2)-I.x(1))*1e3, I.z(1)*1e3, I.z(end)*1e3, (I.z(2)-I.z(1))*1e3);
    fprintf('  centre (x, z) = (%.3f, %.3f) mm\n  diameter lateral %.3f mm, axial %.3f mm\n', ...
        I.xc*1e3, I.zc*1e3, I.dLat*1e3, I.dAx*1e3);
end
a = imgs{1}; b = imgs{2};
fprintf('\nOwn - Lab: dxc %.3f mm, dzc %.3f mm; ratio dLat %.4f, dAx %.4f\n', ...
    (b.xc-a.xc)*1e3, (b.zc-a.zc)*1e3, b.dLat/a.dLat, b.dAx/a.dAx);
fprintf('Lab x step / L11-4v pitch (0.300 mm): %.4f;  1540/7.8125 / (1540/7.6) = %.4f\n', ...
    (a.x(2)-a.x(1))/0.3e-3, 7.6/7.8125);
save(fullfile(outDir, 'zoom_diag.mat'), 'imgs');

% ---- figure 1: physical coordinates, same limits -----------------------
xl = [-19 19]; zl = [0 60];
f1 = figure('Units','centimeters','Position',[2 2 22 12]);
for k = 1:2
    I = imgs{k};
    subplot(1,2,k);
    bm = 20*log10(I.env/max(I.env(:)));
    imagesc(I.x*1e3, I.z*1e3, bm, [-60 0]); colormap gray; axis image;
    xlim(xl); ylim(zl); hold on
    rectangle('Position',[I.xc-I.dLat/2, I.zc-I.dAx/2, I.dLat, I.dAx]*1e3, ...
        'Curvature',1,'EdgeColor','r','LineStyle','--');
    plot(I.xc*1e3, I.zc*1e3, 'r+'); hold off
    xlabel('x [mm]'); ylabel('z [mm]'); title(I.name);
end
exportgraphics(f1, fullfile(outDir, 'bmode_physical.png'), 'Resolution', 150);

% ---- figure 2: pixel index, how it looks without axes ------------------
f2 = figure('Units','centimeters','Position',[2 2 22 12]);
for k = 1:2
    I = imgs{k};
    subplot(1,2,k);
    imagesc(20*log10(I.env/max(I.env(:))), [-60 0]); colormap gray;
    xlabel('column'); ylabel('row'); title([I.name ' (pixel index)']);
end
exportgraphics(f2, fullfile(outDir, 'bmode_pixels.png'), 'Resolution', 150);

% ---- figure 3: profiles ------------------------------------------------
f3 = figure('Units','centimeters','Position',[2 2 22 9]);
subplot(1,2,1); hold on
for k = 1:2, plot(imgs{k}.latX*1e3, imgs{k}.latP); end
xlabel('x [mm]'); ylabel('smoothed env [dB]'); title('Lateral profile through centre');
legend(cellfun(@(s) s.name, imgs, 'UniformOutput', false), 'Location','south'); grid on
subplot(1,2,2); hold on
for k = 1:2, plot(imgs{k}.axZ*1e3, imgs{k}.axP); end
xlabel('z [mm]'); ylabel('smoothed env [dB]'); title('Axial profile through centre'); grid on
exportgraphics(f3, fullfile(outDir, 'profiles.png'), 'Resolution', 150);

% =========================================================================
function I = measureInclusion(I, c0)
% Centre and diameter from half-contrast crossings of profiles taken on a
% smoothed log-envelope. Smoothing is defined in mm, so both grids get the
% same physical kernel. Two passes refine the centre.
sig = 1.0e-3;                       % Gaussian sigma [m]
dx = I.x(2)-I.x(1); dz = I.z(2)-I.z(1);
lg = 20*log10(I.env/max(I.env(:)) + eps);
lg = imgaussfilt(lg, [sig/dz, sig/dx]);
band = 2e-3;                        % half-width of the averaging strip [m]
xc = c0(1); zc = c0(2);
for pass = 1:3
    rows = abs(I.z - zc) <= band;
    latP = mean(lg(rows, :), 1);
    [xl, xr] = crossings(I.x, latP, xc);
    xc = (xl + xr)/2;
    cols = abs(I.x - xc) <= band;
    axP = mean(lg(:, cols), 2).';
    [zt, zb] = crossings(I.z.', axP, zc);
    zc = (zt + zb)/2;
end
I.xc = xc; I.zc = zc; I.dLat = xr - xl; I.dAx = zb - zt;
I.latX = I.x; I.latP = latP; I.axZ = I.z.'; I.axP = axP;
end

function [lo, hi] = crossings(u, p, uc)
% Half-level crossings left and right of uc. Inside level: value at uc.
% Outside level: median of the profile 8 to 14 mm away from uc on both sides.
uf = linspace(u(1), u(end), 20*numel(u));
pf = interp1(u, p, uf);
[~, ic] = min(abs(uf - uc));
inL  = median(pf(abs(uf - uc) <= 2e-3));
outM = abs(uf - uc) >= 8e-3 & abs(uf - uc) <= 14e-3;
outL = median(pf(outM));
h = (inL + outL)/2;
s = sign(inL - outL);
iL = find(s*(pf(1:ic) - h) < 0, 1, 'last');
iR = ic - 1 + find(s*(pf(ic:end) - h) < 0, 1, 'first');
lo = interp1(pf(iL:iL+1), uf(iL:iL+1), h);
hi = interp1(pf(iR-1:iR), uf(iR-1:iR), h);
end
