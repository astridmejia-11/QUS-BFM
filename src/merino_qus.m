function R = merino_qus(sam, refStack, x, z, fs, o)
%MERINO_QUS  Merino's phantoms.m estimator (TV and SWIFT) as a function.
%
%   R = merino_qus(sam, refStack, x, z, fs, o)
%   sam      : beamformed RF [nz x nx]
%   refStack : reference RF [nz x nx x nRef], same grid as sam
%   x, z     : axes in METRES (as stored in the .mat files)
%   fs       : sampling frequency [Hz]
%   o        : struct with fields (defaults = phantoms.m)
%     band      [fL fH] Hz                  [2.5e6 7.5e6]
%     region    [xInf xSup zInf zSup] cm    [0.1 3.8 0.2 3.5]
%     refAcs    reference ACS dB/cm/MHz     0.53
%     refGamma  frequency exponent           1
%     lmuB, lmuC  log10 SWIFT weights        3.5, 1
%     cut       ratioCutOff                  10
%     lmuTV     [log10 muB log10 muC] for TV [3.5 3.5]
%     centre    [c1x c1z] cm (inclusion), [] = no masks
%
%   Steps copied from phantoms.m: crop, 8x12 wavelength blocks with
%   wl = 1540/mean(band), 80 % overlap, Hamming, NFFT = 2^(nextpow2(nz/2)+2),
%   reference spectra averaged over files, compensation with 4*L*ACS_ref,
%   TV = AlterOpti_ADMM, SWIFT = optimAdmmTvTikhonov -> w, wExt -> weighted
%   optimAdmmWeightedTvTikhonov regularised with w. Metric: interp2 of the
%   map to the pixel grid, mean over r <= 0.85 cm (inclusion) and
%   r >= 1.05 cm (background). Verified against Task 01 in
%   scripts/test_merino_qus.m.
d = struct('band',[2.5e6 7.5e6], 'region',[0.1 3.8 0.2 3.5], 'refAcs',0.53, ...
    'refGamma',1, 'lmuB',3.5, 'lmuC',1, 'cut',10, 'lmuTV',[3.5 3.5], 'centre',[]);
fn = fieldnames(d);
for i = 1:numel(fn), if ~isfield(o, fn{i}), o.(fn{i}) = d.(fn{i}); end, end
NptodB = log10(exp(1))*20;
blocksize = 8; overlap_pc = 0.8; ratio_zx = 12/8;

dx = x(2)-x(1); dz = z(2)-z(1);
x = x(:).'*1e2; z = z(:).'*1e2;
ind_x = o.region(1) <= x & x <= o.region(2);
ind_z = o.region(3) <= z & z <= o.region(4);
x = x(ind_x); z = z(ind_z);
sam = sam(ind_z, ind_x);
refStack = refStack(ind_z, ind_x, :);

wl = 1540/mean(o.band);
wx = round(blocksize*wl*(1-overlap_pc)/dx);  nx = round(blocksize*wl/dx);
x0 = 1:wx:length(x)-nx;  x_ACS = x(1, x0+round(nx/2));  n = length(x0);
wz = round(blocksize*wl*(1-overlap_pc)/dz * ratio_zx);
nz = 2*round(blocksize*wl/dz /2 * ratio_zx);
L = (nz/2)*dz*100;
z0p = 1:wz:length(z)-nz;  z0d = z0p + nz/2;  z_ACS = z(z0p + nz/2);  m = length(z0p);
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;
rang = band > o.band(1) & band < o.band(2);
f = band(rang)*1e-6;  p = length(f);
windowing = hamming(nz/2)*ones(1,nx);

Sp = zeros(m,n,p,size(refStack,3)); Sd = Sp;
for i = 1:size(refStack,3)
    [Sp(:,:,:,i), Sd(:,:,:,i)] = blockSpectra(refStack(:,:,i));
end
att_ref_map = repmat(reshape(o.refAcs*f.^o.refGamma/NptodB, 1, 1, p), m, n);
compensation = (log(mean(Sp,4)) - log(mean(Sd,4))) - 4*L*att_ref_map;
[sp, sd] = blockSpectra(sam);
b = (log(sp) - log(sd)) - compensation;

A1 = kron(4*L*f, speye(m*n));  A2 = kron(ones(size(f)), speye(m*n));
mask = ones(m*n*p, 1); tol = 1e-3;

% TV (RSLD)
Bn = AlterOpti_ADMM(A1, A2, b(:), 10^o.lmuTV(1), 10^o.lmuTV(2), m, n, tol, mask);
R.TV = reshape(Bn*NptodB, m, n);

% SWIFT
[~, Cn] = optimAdmmTvTikhonov(A1, A2, b(:), 10^o.lmuB, 10^o.lmuC, m, n, tol, mask);
R.bsc = reshape(Cn*NptodB, m, n);
w = (1-0.1)*(abs(R.bsc) < o.cut) + 0.1;
wExt = movmin(w, 3);
W = repmat(wExt, [1 1 p]);  W = spdiags(W(:), 0, m*n*p, m*n*p);
Bn = optimAdmmWeightedTvTikhonov(W*A1, W*A2, W*b(:), 10^o.lmuB, 10^o.lmuC, m, n, tol, mask, w);
R.SWIFT = reshape(Bn*NptodB, m, n);
R.w = w; R.x_ACS = x_ACS; R.z_ACS = z_ACS; R.x = x; R.z = z; R.f = f; R.L = L;
R.pctPenal = 100*mean(wExt(:) < 1);
bm = 20*log10(abs(hilbert(sam))); R.Bmode = bm - max(bm(:));

if ~isempty(o.centre)
    [X, Z] = meshgrid(x_ACS, z_ACS);  [Xq, Zq] = meshgrid(x, z);
    rr = (Xq-o.centre(1)).^2 + (Zq-o.centre(2)).^2;
    inc = rr <= 0.85^2;  back = rr >= 1.05^2;
    for meth = {'TV', 'SWIFT'}
        A = interp2(X, Z, R.(meth{1}), Xq, Zq);
        R.(['inc' meth{1}])     = mean(A(inc), 'omitnan');
        R.(['sdInc' meth{1}])   = std(A(inc), 'omitnan');
        R.(['back' meth{1}])    = mean(A(back), 'omitnan');
        R.(['sdBack' meth{1}])  = std(A(back), 'omitnan');
    end
end

    function [Sp1, Sd1] = blockSpectra(s)
        Sp1 = zeros(m,n,p); Sd1 = zeros(m,n,p);
        for jj = 1:n
            for ii = 1:m
                xw = x0(jj); zp = z0p(ii); zd = z0d(ii);
                [tp,~] = spectra(s(zp:zp+nz/2-1, xw:xw+nx-1), windowing, 0, nz/2, NFFT);
                [td,~] = spectra(s(zd:zd+nz/2-1, xw:xw+nx-1), windowing, 0, nz/2, NFFT);
                Sp1(ii,jj,:) = tp(rang); Sd1(ii,jj,:) = td(rang);
            end
        end
    end
end
