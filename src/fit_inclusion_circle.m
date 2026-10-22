function M = fit_inclusion_circle(env, x, z, c0, sgn, varargin)
%FIT_INCLUSION_CIRCLE  Circle fit to the edge of an inclusion in a B-mode.
%
%   M = fit_inclusion_circle(env, x, z, c0, sgn)
%   env : envelope [nz x nx];  x, z : axes [m];  c0 : initial centre [x z] [m]
%   sgn : +1 hyperechoic inclusion, -1 hypoechoic
%
%   Rays leave the current centre over the upper half and the sides only
%   (the region below an inclusion is shadowed or enhanced by its
%   attenuation, which biases any edge found there). Along each ray the
%   edge is the radius of steepest change of the smoothed log-envelope with
%   the expected sign. A least-squares (Kasa) circle is fitted to the edge
%   points; three passes refine the centre. M.r is a built-in check: it
%   should be close to the nominal inclusion radius for every target.
p = inputParser;
p.addParameter('Sigma',  1.0e-3);           % smoothing sigma [m]
p.addParameter('RRange', [6e-3 14e-3]);     % radial search range [m]
p.addParameter('Angles', -200:5:20);        % degrees; 0 = +x, -90 = up (z down)
p.parse(varargin{:});
o = p.Results;

x = x(:).'; z = z(:);
dx = x(2)-x(1); dz = z(2)-z(1);
lg = 20*log10(env/max(env(:)) + eps);
lg = imgaussfilt(lg, [o.Sigma/dz, o.Sigma/dx]);
[X, Z] = meshgrid(x, z);
F = griddedInterpolant(X.', Z.', lg.', 'linear', 'none');

r  = linspace(o.RRange(1), o.RRange(2), 200);
xc = c0(1); zc = c0(2);
th = deg2rad(o.Angles);
for pass = 1:3
    ex = nan(size(th)); ez = ex;
    for k = 1:numel(th)
        px = xc + r*cos(th(k));  pz = zc + r*sin(th(k));
        v  = F(px, pz);
        if any(isnan(v)), continue; end
        dv = gradient(v, r);
        [~, i] = min(sgn*dv);            % steepest drop going outward (hyper)
        ex(k) = px(i); ez(k) = pz(i);
    end
    ok = ~isnan(ex);
    A = [ex(ok).' ez(ok).' ones(nnz(ok),1)];
    s = A \ (ex(ok).'.^2 + ez(ok).'.^2);
    xc = s(1)/2; zc = s(2)/2;
    R = sqrt(s(3) + xc^2 + zc^2);
end
res = hypot(ex(ok) - xc, ez(ok) - zc) - R;
M.xc = xc; M.zc = zc; M.r = R; M.rmsRes = sqrt(mean(res.^2));
M.nPts = nnz(ok); M.ex = ex; M.ez = ez;
end
