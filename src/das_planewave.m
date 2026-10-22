function B = das_planewave(D, varargin)
%DAS_PLANEWAVE  Delay-and-sum beamforming of plane-wave data on a Cartesian grid.
%
%   B = das_planewave(D)
%   B = das_planewave(D, 'FNumber', 2, 'dx', 0.3e-3)
%
%   D comes from load_vsx_rf. Options:
%     'FNumber' : receive f-number (focus depth / active aperture). Default 2.
%     'Apod'    : 'hanning' (default) or 'boxcar' receive apodisation.
%     'dx'      : lateral pixel spacing [m]. Default: the element pitch, so the
%                 grid matches the 128 scanlines of the focused acquisition and
%                 both schemes can be compared on identical sampling.
%     'dz'      : axial pixel spacing [m]. Default c/(2*fs), the native axial
%                 sample spacing. Do NOT coarsen this: the spectral estimator
%                 needs the full bandwidth.
%     'zLim'    : [zmin zmax] in metres. Default [1e-3, D.endDepth_m].
%     'xLim'    : [xmin xmax] in metres. Default the array extent.
%     'fDemod'  : demodulation frequency [Hz]. Default D.fc.
%
%   Output B:
%     B.x, B.z     grid axes [m]
%     B.iq         [nz x nx] complex analytic beamformed signal
%     B.rf         [nz x nx] real RF  = real(B.iq)   <- feed this to SWIFT
%     B.env        [nz x nx] envelope = abs(B.iq)
%     B.nActive    [nz x nx] elements contributing to each pixel
%     B.params     every setting used, for the report and for reproducibility
%+
%   WHY IQ INSTEAD OF RAW RF INTERPOLATION
%   The delay for a pixel almost never lands on a sample. At 27.8 MHz with a
%   7.24 MHz carrier there are only ~3.8 samples per cycle, so interpolating
%   the raw RF linearly commits a phase error that varies with sub-sample
%   position. Summed over 128 channels those errors destroy coherence. The fix
%   is to strip the carrier first: the complex envelope varies slowly and
%   interpolates cleanly, and the carrier phase is restored analytically with
%   exp(1j*2*pi*f*tau), which has no interpolation error at all.
%   This matters more for QUS than for B-mode. A coherence loss that varies
%   with depth is indistinguishable from attenuation, which is exactly the
%   quantity being estimated.

% ---- options ------------------------------------------------------------
p = inputParser;
p.addParameter('FNumber', 2);
p.addParameter('Apod',    'hanning');
p.addParameter('dx',      []);
p.addParameter('dz',      []);
p.addParameter('zLim',    []);
p.addParameter('xLim',    []);
p.addParameter('fDemod',  []);
p.addParameter('Verbose', true);
p.parse(varargin{:});
o = p.Results;

if ~all(arrayfun(@(s) s.type == "planewave", D.tx))
    error('das_planewave:notPlaneWave', ...
        'This function expects plane-wave transmits. Use das_focused instead.');
end

c  = D.c;
fs = D.fs;
fd = o.fDemod; if isempty(fd), fd = D.fc; end

ex = D.elemX(:);            % element positions [m]
ez = D.elemZ(:);
nEl = numel(ex);

% ---- 1. grid ------------------------------------------------------------
dx = o.dx; if isempty(dx), dx = D.pitch; end
dz = o.dz; if isempty(dz), dz = c/(2*fs); end
zL = o.zLim; if isempty(zL), zL = [1e-3, D.endDepth_m]; end
xL = o.xLim; if isempty(xL), xL = [min(ex), max(ex)]; end

x = (xL(1):dx:xL(2)).';
z = (zL(1):dz:zL(2)).';
nx = numel(x); nz = numel(z);
[X, Z] = meshgrid(x, z);        % both [nz x nx]

% ---- 2. demodulate to baseband -----------------------------------------
% hilbert() builds the analytic signal column by column, so each element's
% channel is transformed independently. Multiplying by exp(-j*2*pi*fd*t)
% then shifts the spectrum down to zero, leaving the slowly varying complex
% envelope that is safe to interpolate.
t  = D.t(:);
nS = numel(t);
nEv = numel(D.tx);

BB = complex(zeros(nS, nEl, nEv));
carrierDown = exp(-1j*2*pi*fd*t);
for e = 1:nEv
    BB(:, :, e) = hilbert(D.rf(:, :, e)) .* carrierDown;
end

% ---- 3. beamform --------------------------------------------------------
iq      = complex(zeros(nz, nx));
nActive = zeros(nz, nx);
t0 = t(1);
dt = 1/fs;

for e = 1:nEv
    % Transmit time of flight for a plane wave steered by theta. The wavefront
    % is a straight line, so the travel time to a point is just the projection
    % of its position onto the propagation direction. For theta = 0 this is
    % simply z/c, which matches a transmit with all delays set to zero.
    th = D.tx(e).steer;
    tauTx = (Z*cos(th) + X*sin(th)) / c;

    for i = 1:nEl
        % Receive time of flight: straight line from pixel back to element i.
        tauRx = hypot(X - ex(i), Z - ez(i)) / c;
        tau   = tauTx + tauRx;

        % Dynamic receive aperture. An element only contributes to a pixel if
        % it lies inside the cone set by the f-number, so the aperture grows
        % with depth and the lateral resolution stays roughly constant.
        halfAp = Z / (2*o.FNumber);
        dxe    = X - ex(i);
        inAp   = abs(dxe) <= halfAp;

        switch lower(o.Apod)
            case 'boxcar'
                w = double(inAp);
            case 'hanning'
                % Smooth taper across the aperture: avoids the discontinuity a
                % boxcar creates as elements enter the aperture with depth.
                w = zeros(nz, nx);
                w(inAp) = 0.5 * (1 + cos(pi * dxe(inAp) ./ halfAp(inAp)));
            otherwise
                error('das_planewave:apod', 'Unknown apodisation: %s', o.Apod);
        end

        % Linear interpolation of the BASEBAND signal, then analytic phase
        % restoration. Done by hand rather than with interp1: same result,
        % roughly an order of magnitude faster over ~40M lookups.
        fidx = (tau - t0)/dt + 1;
        i0   = floor(fidx);
        frac = fidx - i0;
        ok   = (i0 >= 1) & (i0 < nS) & (w > 0);

        if ~any(ok(:)), continue; end

        bb  = BB(:, i, e);
        val = complex(zeros(nz, nx));
        val(ok) = bb(i0(ok)) .* (1 - frac(ok)) + bb(i0(ok) + 1) .* frac(ok);

        % Put the carrier back at the exact delay. No interpolation error here.
        val = val .* exp(1j*2*pi*fd*tau);

        iq      = iq + w .* val;
        nActive = nActive + double(ok);
    end
end

% ---- 4. outputs ---------------------------------------------------------
B = struct();
B.x   = x;
B.z   = z;
B.iq  = iq;
B.rf  = real(iq);     % the analytic signal's real part IS the RF
B.env = abs(iq);
B.nActive = nActive;
B.params = struct('fNumber', o.FNumber, 'apod', string(o.Apod), ...
                  'dx', dx, 'dz', dz, 'c', c, 'fs', fs, 'fc', D.fc, ...
                  'fDemod', fd, 'nEvents', nEv, ...
                  'dataFile', D.dataFile, 'preSetFile', D.preSetFile, ...
                  'frame', D.frame, 'cIsAssumed', D.cIsAssumed);

if o.Verbose
    fprintf('DAS plane wave\n');
    fprintf('  grid    : %d x %d  (dz %.1f um, dx %.1f um)\n', ...
        nz, nx, dz*1e6, dx*1e6);
    fprintf('  extent  : x %.1f..%.1f mm,  z %.1f..%.1f mm\n', ...
        x(1)*1e3, x(end)*1e3, z(1)*1e3, z(end)*1e3);
    fprintf('  f-number: %.2f (%s apodisation)\n', o.FNumber, o.Apod);
    fprintf('  aperture: %d..%d elements per pixel\n', ...
        min(nActive(:)), max(nActive(:)));
    fprintf('  c       : %.0f m/s%s\n', c, ...
        repmat(' (ASSUMED at acquisition)', 1, D.cIsAssumed));
end
end
