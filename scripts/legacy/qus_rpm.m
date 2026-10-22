% QUS_RPM  Reference Phantom Method, self-contained.
%
%   Everything the estimator needs lives in this one file: the driver at the
%   top and the estimation functions at the bottom as local functions. Only
%   load_vsx_rf, das_planewave and das_focused are external.
%
%   Run from the amejia folder:  qus_rpm
%
%   Written in plain ASCII on purpose, so that no encoding surprise can break
%   the parser.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT      = pwd;
SCHEME    = 'PW';                 % 'PW' (plane wave) or 'FOC' (focused)
TARGET    = 'T1';                 % T1, T3, T5, T7, T8
C_SOUND   = 1542;                 % CIRS background, from the certificate
Z_BF      = [3e-3  58e-3];        % depth reconstructed by the beamformer
Z_ROI     = [10e-3 45e-3];        % depth actually analysed
ACS_REF   = 0.53;                 % provisional until the Wisconsin table arrives
ACS_TRUE  = 0.55;                 % certificate: CIRS background
BAND_MAN  = [4e6 9e6];
BAND_DROP = 14;                   % dB below peak for the automatic band
MAP_SPAN  = 15e-3;                % depth span for the depth-resolved map
LAT_AVG   = 3;                    % lateral blocks pooled per spectrum

% Region analysis. The spectral log difference assumes the BSC term does not
% vary with depth inside the fitted span, so each region must lie ENTIRELY
% within one material. Crossing the target boundary injects the BSC step into
% the slope and biases the ACS with an arbitrary sign.
X_BG      = 12e-3;                % background: |x| greater than this
X_TG      =  6e-3;                % target: |x| smaller than this
Z_TG      = [13e-3 25e-3];        % depth span fully inside the target
ACS_TG    = 0.52;                 % certificate value for this target
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
    otherwise
        error('SCHEME must be PW or FOC');
end

% ---- locate the files ---------------------------------------------------
dataOnly = @(d) d(~contains({d.name}, 'preSet'));

sf  = dataOnly(dir(fullfile(sampleDir, '*.mat')));
hit = find(contains({sf.name}, ['_' TARGET '_']) | ...
           contains({sf.name}, ['_' TARGET '-']), 1);
if isempty(hit)
    fprintf('Files in %s:\n', sampleDir);
    fprintf('  %s\n', sf.name);
    error('Target %s not found.', TARGET);
end
sampleFile = fullfile(sf(hit).folder, sf(hit).name);

rfl = dataOnly(dir(fullfile(refDir, '*.mat')));
if isempty(rfl), error('No reference files in %s', refDir); end
fprintf('Sample    : %s\n', sf(hit).name);
fprintf('References: %d file(s)\n\n', numel(rfl));

% ---- beamform -----------------------------------------------------------
Ds = load_vsx_rf(sampleFile, 'c', C_SOUND);
Bs = bf(Ds, 'zLim', Z_BF);
clear Ds

if isfile(cacheFile)
    fprintf('\nLoading cached references\n');
    L = load(cacheFile, 'Bref');
    Bref = L.Bref;
else
    fprintf('\nBeamforming %d references...\n', numel(rfl));
    Bref = cell(1, numel(rfl));
    for k = 1:numel(rfl)
        fprintf('  [%d/%d] %s\n', k, numel(rfl), rfl(k).name);
        Dk = load_vsx_rf(fullfile(rfl(k).folder, rfl(k).name), ...
                         'c', C_SOUND, 'Verbose', false);
        Bref{k} = bf(Dk, 'zLim', Z_BF, 'Verbose', false);
        clear Dk
    end
    save(cacheFile, 'Bref', '-v7.3');
end

% ---- estimate -----------------------------------------------------------
common = {'AcsRef', ACS_REF, 'LatAvg', LAT_AVG, 'BandDrop', BAND_DROP};

fprintf('\n=== (a) full-region fit, automatic band ===\n');
Aprof = rpm_estimate(Bs, Bref, common{:}, 'zLim', Z_ROI);

fprintf('\n=== (b) short-span fit, depth-resolved map ===\n');
Amap  = rpm_estimate(Bs, Bref, common{:}, 'zLim', Z_ROI, 'FitSpan', MAP_SPAN);

fprintf('\n=== (c) full-region fit, manual band ===\n');
Aman  = rpm_estimate(Bs, Bref, common{:}, 'zLim', Z_ROI, 'Band', BAND_MAN);

% ---- figures ------------------------------------------------------------
figure('Position', [60 60 1300 430]);

subplot(1,3,1);
imagesc(Bs.x*1e3, Bs.z*1e3, 20*log10(Bs.env/max(Bs.env(:))), [-60 0]);
axis image; colormap(gca, gray); colorbar;
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title(sprintf('B-mode %s %s', TARGET, SCHEME));
hold on;
plot([Bs.x(1) Bs.x(end) Bs.x(end) Bs.x(1) Bs.x(1)]*1e3, ...
     [Z_ROI(1) Z_ROI(1) Z_ROI(2) Z_ROI(2) Z_ROI(1)]*1e3, ...
     'y--', 'LineWidth', 1);
hold off;

subplot(1,3,2);
plot(Aprof.x*1e3, Aprof.acs, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
hold on;
plot(Aman.x*1e3, Aman.acs, 's--', 'LineWidth', 1.0, 'MarkerSize', 4);
yline(ACS_TRUE, 'k-', 'background 0.55', 'LineWidth', 1.2);
hold off; grid on; ylim([0 1.4]);
xlabel('Lateral [mm]'); ylabel('ACS [dB/cm/MHz]');
title(sprintf('Profile, fit over %.0f mm', Aprof.params.fitSpan_mm));
legend({sprintf('auto %.1f-%.1f MHz', Aprof.band/1e6), ...
        sprintf('manual %.1f-%.1f MHz', Aman.band/1e6)}, ...
       'Location', 'south', 'FontSize', 8);

subplot(1,3,3);
imagesc(Amap.x*1e3, Amap.z*1e3, Amap.acs, [0.2 1.0]);
axis image; colormap(gca, turbo); colorbar;
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title(sprintf('Map, fit over %.0f mm', Amap.params.fitSpan_mm));

% ---- region analysis ----------------------------------------------------
% Both regions use the SAME fixed band, so any difference between them comes
% from the tissue and not from the band selection.
fprintf('\n=== region analysis (band %.1f-%.1f MHz) ===\n', BAND_MAN/1e6);

Abg = rpm_estimate(Bs, Bref, common{:}, 'Band', BAND_MAN, ...
                   'zLim', Z_ROI, 'xSel', @(x) abs(x) > X_BG, 'Verbose', false);
Atg = rpm_estimate(Bs, Bref, common{:}, 'Band', BAND_MAN, ...
                   'zLim', Z_TG,  'xSel', @(x) abs(x) < X_TG, 'Verbose', false);

% CONTROL: the background lateral region, but over the SAME short depth span
% used for the target. Target and background differ in two things at once
% (lateral position and depth span), so this isolates the span. If the
% control is also inflated, the bias belongs to the short fit span and not to
% the target material.
Act = rpm_estimate(Bs, Bref, common{:}, 'Band', BAND_MAN, ...
                   'zLim', Z_TG,  'xSel', @(x) abs(x) > X_BG, 'Verbose', false);

vbg = Abg.acs(:); vbg = vbg(isfinite(vbg));
vtg = Atg.acs(:); vtg = vtg(isfinite(vtg));
vct = Act.acs(:); vct = vct(isfinite(vct));

fprintf('\n%-14s %8s %8s %8s %8s %8s\n', ...
    'region', 'span', 'median', 'IQR', 'cert', 'NBias');
fprintf('%-14s %6.1fmm %8.3f %8.3f %8.2f %7.1f%%\n', 'background', ...
    Abg.params.fitSpan_mm, median(vbg), iqr(vbg), ACS_TRUE, ...
    100*(median(vbg)-ACS_TRUE)/ACS_TRUE);
fprintf('%-14s %6.1fmm %8.3f %8.3f %8.2f %7.1f%%\n', 'target', ...
    Atg.params.fitSpan_mm, median(vtg), iqr(vtg), ACS_TG, ...
    100*(median(vtg)-ACS_TG)/ACS_TG);
fprintf('%-14s %6.1fmm %8.3f %8.3f %8.2f %7.1f%%   <- control\n', 'bg, short', ...
    Act.params.fitSpan_mm, median(vct), iqr(vct), ACS_TRUE, ...
    100*(median(vct)-ACS_TRUE)/ACS_TRUE);
fprintf(['\nIf "bg, short" is inflated like the target, the bias belongs to\n' ...
         'the short fit span, not to the target material.\n']);

% The contrast is the number that does NOT depend on the reference: the
% provisional AcsRef cancels in the subtraction, so this can be compared with
% the certificate today, without waiting for the Wisconsin table.
fprintf('\ncontrast target - background\n');
fprintf('  measured   %+.3f dB/cm/MHz\n', median(vtg) - median(vbg));
fprintf('  certificate %+.3f  (%.2f - %.2f)\n', ...
    ACS_TG - ACS_TRUE, ACS_TG, ACS_TRUE);

% What the reference phantom would have to be for the background to land on
% its certified value.
fprintf('\nimplied reference ACS: %.3f dB/cm/MHz\n', ...
    ACS_REF + (ACS_TRUE - median(vbg)));

Bs2   = Bref{1};
Brest = Bref(2:end);

spans = linspace(8e-3, 34e-3, 10);
iqrs  = zeros(size(spans));
zc    = 27.5e-3;

for k = 1:numel(spans)
    zk = zc + [-1 1]*spans(k)/2;
    Ak = rpm_estimate(Bs2, Brest, 'AcsRef', 0, 'LatAvg', 1, ...
                      'Band', BAND_MAN, 'zLim', zk, 'Verbose', false);
    v = Ak.acs(:); v = v(isfinite(v));
    iqrs(k) = iqr(v);
    fprintf('span %5.1f mm   n=%4d   mediana %+.3f   IQR %.3f\n', ...
        Ak.params.fitSpan_mm, numel(v), median(v), iqrs(k));
end

ok = spans <= 22e-3;
pw = polyfit(log(spans(ok)), log(iqrs(ok)), 1);
fprintf('Exponente (span <= 22 mm): %.2f\n', pw(1));

figure;
loglog(spans*1e3, iqrs, 'o', 'MarkerSize', 7, 'LineWidth', 1.5); hold on;
loglog(spans*1e3, exp(polyval(pw, log(spans))), '-',  'LineWidth', 1.3);
loglog(spans*1e3, iqrs(end)*(spans/spans(end)).^-1.5, '--', 'LineWidth', 1.3);
hold off; grid on; yline(median(iqrs(spans > 22e-3)), ':', 'piso ~0.12', 'LineWidth', 1.2);
xlim([7 40]); xlabel('Span del ajuste [mm]'); ylabel('IQR de ACS [dB/cm/MHz]');
legend({'medido', sprintf('ajuste L^{%.2f}', pw(1)), 'teoria L^{-3/2}'}, ...
       'Location', 'southwest');



% #########################################################################
%                            LOCAL FUNCTIONS
% #########################################################################

function A = rpm_estimate(Bs, Bref, varargin)
%RPM_ESTIMATE  ACS by the Reference Phantom Method.
%
%   The power spectrum measured in a window at depth z is
%
%       S(f,z) = G(f) * D(f,z) * BSC(f) * exp(-4*alpha(f)*z)
%
%   G is the system response (transducer, pulse, TGC) and D is beam
%   diffraction. Neither is known. Dividing by a reference acquired with
%   IDENTICAL settings cancels both:
%
%       ln[Ss/Sr](f,z) = ln[BSCs/BSCr](f) - 4*[alpha_s(f)-alpha_r(f)]*z
%
%   The first term does not depend on z, so it lands in the intercept and the
%   SLOPE against depth gives the attenuation difference. That slope is taken
%   at every frequency, converted to dB/cm, and fitted against frequency: the
%   second slope is dACS. Adding the reference ACS gives the sample ACS.
%
%   Because the final fit is a slope against frequency, only the reference
%   ACS SLOPE enters, not its absolute level. A provisional AcsRef therefore
%   shifts every result by a known constant, correctable later without
%   recomputing anything.
%
%   FitSpan is the parameter that matters most. The estimate is the slope of
%   a line fitted against depth, so its variance scales as 1/(N*var(z)) and
%   collapses when the span is short. With 12-wavelength windows at 80 %
%   overlap the step is about 0.5 mm, so a 5-window fit spans only 4.5 mm.
%   Over that distance a 0.02 dB/cm/MHz difference moves the log spectral
%   ratio by hundredths, while speckle noise on that ratio is of order 0.5.
%   Fitting over the full region instead raises the span eightfold and the
%   sample count from 5 to about 70, cutting the slope error roughly 20x.
%   The price is depth resolution, which is exactly the trade that
%   regularised methods such as SWIFT exist to avoid.

p = inputParser;
p.addParameter('Band',       'auto');
p.addParameter('BandDrop',   20);
p.addParameter('AcsRef',     0.53);
p.addParameter('WinAxLam',   12);
p.addParameter('WinLatLam',  8);
p.addParameter('Overlap',    0.8);
p.addParameter('FitSpan',    Inf);
p.addParameter('LatAvg',     1);
p.addParameter('RefLatAvg',  true);
p.addParameter('zLim',       []);
p.addParameter('xSel',       []);
p.addParameter('Verbose',    true);
p.parse(varargin{:});
o = p.Results;

if ~iscell(Bref), Bref = {Bref}; end

c   = Bs.params.c;
fc  = Bs.params.fc;
lam = c/fc;
dz  = Bs.z(2) - Bs.z(1);
dx  = Bs.x(2) - Bs.x(1);

for k = 1:numel(Bref)
    if numel(Bref{k}.z) ~= numel(Bs.z) || numel(Bref{k}.x) ~= numel(Bs.x)
        error('rpm:grid', ['Sample and reference are on different grids. ' ...
            'Beamform both with the same options.']);
    end
end

zL = o.zLim;
if isempty(zL), zL = [Bs.z(1), Bs.z(end)]; end
zmask = find(Bs.z >= zL(1) & Bs.z <= zL(2));

winZ  = round(o.WinAxLam  * lam / dz);
winX  = max(1, round(o.WinLatLam * lam / dx));
stepZ = max(1, round(winZ * (1 - o.Overlap)));
stepX = max(1, round(winX * (1 - o.Overlap)));

z0 = zmask(1) : stepZ : (zmask(end) - winZ + 1);
x0 = 1        : stepX : (numel(Bs.x) - winX + 1);

% Restrict to a lateral region. Needed so that the fitted span lies entirely
% inside one material: the BSC term of the log spectral ratio must not vary
% with depth, and it does the moment the region crosses a target boundary.
if ~isempty(o.xSel)
    xc   = Bs.x(x0) + winX*dx/2;
    keep = o.xSel(xc);
    x0   = x0(keep);
    if isempty(x0)
        error('rpm:xSel', 'xSel removed every lateral block.');
    end
end

nzb = numel(z0);
nxb = numel(x0);
if nzb < 3
    error('rpm:tooFewWindows', ...
        'Only %d depth windows fit. Widen zLim or lower Overlap.', nzb);
end

if isinf(o.FitSpan)
    nFit = nzb;
else
    nFit = max(3, round(o.FitSpan / (stepZ*dz)) + 1);
    nFit = min(nFit, nzb);
end
fitSpan_mm = (nFit-1) * stepZ * dz * 1e3 + o.WinAxLam*lam*1e3;

nfft  = 2^nextpow2(winZ) * 2;
f     = (0:nfft-1).' * (Bs.params.fs / nfft);
half  = 1 : floor(nfft/2);
f     = f(half);
taper = hann(winZ);

% ---- power spectra ------------------------------------------------------
Ss = blockSpectra(Bs.rf, z0, x0, winZ, winX, taper, nfft, half);
if o.LatAvg > 1
    % Pool neighbouring lateral blocks: more independent speckle per
    % spectrum, at the cost of lateral resolution.
    Ss = movmean(Ss, o.LatAvg, 3);
end

Sr = zeros(size(Ss));
for k = 1:numel(Bref)
    Sr = Sr + blockSpectra(Bref{k}.rf, z0, x0, winZ, winX, taper, nfft, half);
end
Sr = Sr / numel(Bref);

if o.RefLatAvg
    % The reference phantom is homogeneous, so its expected spectrum depends
    % on depth only. Averaging laterally keeps the diffraction correction
    % intact while removing the speckle variance it would carry into the ratio.
    Sr = repmat(mean(Sr, 3), 1, 1, nxb);
end

% ---- analysis band ------------------------------------------------------
if ischar(o.Band) || isstring(o.Band)
    meanS  = mean(mean(Ss, 3), 2);
    meanS  = meanS / max(meanS);
    inBand = 10*log10(meanS + eps) >= -o.BandDrop;
    [~, ipk] = max(meanS);
    lo = ipk; while lo > 1        && inBand(lo-1), lo = lo - 1; end
    hi = ipk; while hi < numel(f) && inBand(hi+1), hi = hi + 1; end
    band = [f(lo) f(hi)];
    bandSource = sprintf('auto (-%g dB)', o.BandDrop);
else
    band = o.Band(:).';
    bandSource = 'manual';
end
fsel = find(f >= band(1) & f <= band(2));
if numel(fsel) < 4
    error('rpm:band', 'Band contains only %d bins.', numel(fsel));
end
fMHz = f(fsel) / 1e6;

% ---- log spectral ratio -------------------------------------------------
X = log(Ss(fsel,:,:) ./ (Sr(fsel,:,:) + eps));

% ---- fit 1: slope against depth, at every frequency ---------------------
zBlock = (Bs.z(z0) + winZ*dz/2) * 100;      % block centre, cm
nOut   = nzb - nFit + 1;
acs    = nan(nOut, nxb);
r2     = nan(nOut, nxb);
zOut   = zeros(nOut, 1);

fBar = mean(fMHz);
fDev = fMHz - fBar;
fDen = sum(fDev.^2);

for iz = 1:nOut
    kk   = iz : iz + nFit - 1;
    zk   = zBlock(kk);
    zDev = zk(:) - mean(zk);
    zDen = sum(zDev.^2);
    zOut(iz) = mean(zk);

    for ix = 1:nxb
        Xk = squeeze(X(:, kk, ix));
        % slope of ln(ratio) against depth, one value per frequency [1/cm]
        m  = (Xk - mean(Xk, 2)) * zDev / zDen;

        % m = -4*dAlpha with dAlpha in Np/cm  ->  dAlpha[dB/cm] = -8.686*m/4
        dAlpha = -8.686 * m / 4;

        % fit 2: dAlpha against frequency. That slope is dACS [dB/cm/MHz].
        dAcs = sum(fDev .* (dAlpha - mean(dAlpha))) / fDen;
        acs(iz, ix) = dAcs + o.AcsRef;

        pred = mean(dAlpha) + dAcs * fDev;
        ssr  = sum((dAlpha - pred).^2);
        sst  = sum((dAlpha - mean(dAlpha)).^2);
        r2(iz, ix) = 1 - ssr / max(sst, eps);
    end
end

A = struct();
A.acs  = acs;
A.r2   = r2;
A.z    = zOut / 100;
A.x    = Bs.x(x0) + winX*dx/2;
A.band = band;
A.f    = f(fsel);
A.params = struct('acsRef', o.AcsRef, 'bandSource', string(bandSource), ...
                  'winAxLam', o.WinAxLam, 'winLatLam', o.WinLatLam, ...
                  'winZ_samples', winZ, 'winX_lines', winX, ...
                  'overlap', o.Overlap, 'nFit', nFit, ...
                  'fitSpan_mm', fitSpan_mm, 'latAvg', o.LatAvg, ...
                  'refLatAvg', o.RefLatAvg, 'nRef', numel(Bref), ...
                  'lambda_mm', lam*1e3);

if o.Verbose
    fprintf('RPM\n');
    fprintf('  window  : %.1f x %.1f mm (%d samples x %d lines), overlap %.0f%%\n', ...
        winZ*dz*1e3, winX*dx*1e3, winZ, winX, 100*o.Overlap);
    fprintf('  blocks  : %d axial x %d lateral -> output %d x %d\n', ...
        nzb, nxb, nOut, nxb);
    fprintf('  fit span: %.1f mm over %d windows (lateral pooling %d)\n', ...
        fitSpan_mm, nFit, o.LatAvg);
    fprintf('  band    : %.2f - %.2f MHz [%s], %d bins\n', ...
        band(1)/1e6, band(2)/1e6, bandSource, numel(fsel));
    fprintf('  refs    : %d acquisition(s)\n', numel(Bref));
    fprintf('  result  : median %.3f, IQR %.3f dB/cm/MHz, median R2 %.2f\n', ...
        median(acs(:), 'omitnan'), iqr(acs(:)), median(r2(:), 'omitnan'));
end
end

% =========================================================================
function S = blockSpectra(rf, z0, x0, winZ, winX, taper, nfft, half)
%BLOCKSPECTRA  Mean power spectrum of every analysis block. [nf x nzb x nxb]
nzb = numel(z0);
nxb = numel(x0);
S   = zeros(numel(half), nzb, nxb);
for iz = 1:nzb
    rows = z0(iz) : z0(iz) + winZ - 1;
    for ix = 1:nxb
        cols = x0(ix) : x0(ix) + winX - 1;
        seg  = rf(rows, cols) .* taper;   % taper suppresses spectral leakage
        F    = fft(seg, nfft, 1);
        % Average the periodograms of the lines in the block: a single-line
        % spectrum is dominated by its own speckle realisation.
        S(:, iz, ix) = mean(abs(F(half, :)).^2, 2);
    end
end
end