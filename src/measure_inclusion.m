function M = measure_inclusion(env, x, z, c0, varargin)
%MEASURE_INCLUSION  Centre and diameter of a circular inclusion in a B-mode.
%
%   M = measure_inclusion(env, x, z, c0)
%   env : envelope [nz x nx];  x, z : axes [m];  c0 : initial guess [x z] [m]
%
%   Profiles are taken on a log-envelope smoothed with a Gaussian defined in
%   mm, so different grids get the same physical kernel. The centre is the
%   midpoint of the half-contrast crossings of a lateral and an axial
%   profile; three passes refine it. M.contrastLat / M.contrastAx report the
%   inside-minus-outside level [dB]: below ~3 dB the crossings are not
%   reliable and the caller should not trust the result.
p = inputParser;
p.addParameter('Sigma', 1.0e-3);    % smoothing sigma [m]
p.addParameter('Band',  2e-3);      % half-width of the averaging strip [m]
p.parse(varargin{:});
o = p.Results;

x = x(:).'; z = z(:);
dx = x(2)-x(1); dz = z(2)-z(1);
lg = 20*log10(env/max(env(:)) + eps);
lg = imgaussfilt(lg, [o.Sigma/dz, o.Sigma/dx]);
xc = c0(1); zc = c0(2);
M.ok = true;
try
    for pass = 1:3
        rows = abs(z - zc) <= o.Band;
        latP = mean(lg(rows, :), 1);
        [xl, xr, cL] = crossings(x, latP, xc);
        xc = (xl + xr)/2;
        cols = abs(x - xc) <= o.Band;
        axP = mean(lg(:, cols), 2).';
        [zt, zb, cA] = crossings(z.', axP, zc);
        zc = (zt + zb)/2;
    end
catch
    M.ok = false; xl = NaN; xr = NaN; zt = NaN; zb = NaN; cL = NaN; cA = NaN;
    latP = []; axP = [];
end
M.xc = xc; M.zc = zc; M.dLat = xr - xl; M.dAx = zb - zt;
M.contrastLat = cL; M.contrastAx = cA;
M.latX = x; M.latP = latP; M.axZ = z.'; M.axP = axP;
end

function [lo, hi, contrast] = crossings(u, p, uc)
% Half-level crossings left and right of uc. Inside level: median within
% 2 mm of uc. Outside level: median 8 to 14 mm away on both sides.
uf = linspace(u(1), u(end), 20*numel(u));
pf = interp1(u, p, uf);
[~, ic] = min(abs(uf - uc));
inL  = median(pf(abs(uf - uc) <= 2e-3));
outM = abs(uf - uc) >= 8e-3 & abs(uf - uc) <= 14e-3;
outL = median(pf(outM));
contrast = inL - outL;
h = (inL + outL)/2;
s = sign(inL - outL);
iL = find(s*(pf(1:ic) - h) < 0, 1, 'last');
iR = ic - 1 + find(s*(pf(ic:end) - h) < 0, 1, 'first');
lo = interp1(pf(iL:iL+1), uf(iL:iL+1), h);
hi = interp1(pf(iR-1:iR), uf(iR-1:iR), h);
end
