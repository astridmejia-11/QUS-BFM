function B = das_focused(D, varargin)
%DAS_FOCUSED  Delay-and-sum beamforming of focused scanline data.
%
%   B = das_focused(D)
%   B = das_focused(D, 'FNumber', 2, 'zLim', [1e-3 60e-3])
%
%   D comes from load_vsx_rf and must hold focused transmits (one focused beam
%   per lateral origin). One A-line is reconstructed per transmit event, at the
%   lateral position of that event's transmit origin, so the output grid has
%   one column per scanline.
%
%   Options (same meaning as das_planewave):
%     'FNumber' receive f-number, default 2
%     'Apod'    'hanning' (default) or 'boxcar'
%     'dz'      axial spacing [m], default c/(2*fs)
%     'zLim'    [zmin zmax] in metres
%     'fDemod'  demodulation frequency [Hz], default D.fc
%
%   THE TRANSMIT MODEL
%   A focused transmit is not a simple wavefront: it converges to the focus and
%   diverges after it. The standard treatment is the VIRTUAL SOURCE model. The
%   focal point behaves as a point source that fires at a single instant t_f;
%   before the focal plane the wave is arriving at it, after the focal plane it
%   is leaving it. So the transmit time of flight to a pixel is
%
%       tau_tx = t_f + sign(z - z_f) * distance(pixel, focus) / c
%
%   The sign is what makes this work: above the focus the pixel is reached
%   BEFORE t_f, below it AFTER.
%
%   t_f is not assumed here, it is measured from the TX struct: for every
%   active element, delay_i + distance(element_i, focus)/c must give the same
%   number, because that is what "focused" means. The spread of those numbers
%   is reported as a diagnostic; a large spread means the delays, the focus, or
%   the assumed sound speed disagree with each other.

p = inputParser;
p.addParameter('FNumber', 2);
p.addParameter('Apod',    'hanning');
p.addParameter('dz',      []);
p.addParameter('zLim',    []);
p.addParameter('fDemod',  []);
p.addParameter('Verbose', true);
p.parse(varargin{:});
o = p.Results;

if ~any(arrayfun(@(s) s.type == "focused", D.tx))
    error('das_focused:notFocused', ...
        'This function expects focused transmits. Use das_planewave instead.');
end

c  = D.c;
fs = D.fs;
fd = o.fDemod; if isempty(fd), fd = D.fc; end
ex = D.elemX(:);
ez = D.elemZ(:);
nEl = numel(ex);

dz = o.dz; if isempty(dz), dz = c/(2*fs); end
zL = o.zLim; if isempty(zL), zL = [1e-3, D.endDepth_m]; end
z  = (zL(1):dz:zL(2)).';
nz = numel(z);

nEv = numel(D.tx);
t  = D.t(:);
nS = numel(t);
t0 = t(1);
dt = 1/fs;
carrierDown = exp(-1j*2*pi*fd*t);

x       = zeros(1, nEv);
iq      = complex(zeros(nz, nEv));
nActive = zeros(nz, nEv);
tfSpread = zeros(1, nEv);

for e = 1:nEv
    T  = D.tx(e);
    th = T.steer;

    % ---- focus geometry -------------------------------------------------
    xf = T.origin(1) + T.focus * sin(th);
    zf = T.focus * cos(th);
    x(e) = T.origin(1);

    % ---- the instant the wave passes through the focus ------------------
    act = find(T.apod ~= 0);
    if isempty(act)
        continue    % no active elements: nothing was transmitted
    end
    tArrive = T.delay(act).' + hypot(ex(act) - xf, ez(act) - zf) / c;
    tF      = mean(tArrive);
    tfSpread(e) = std(tArrive);

    % ---- transmit time of flight to every pixel on this line ------------
    dFocus = hypot(x(e) - xf, z - zf);
    tauTx  = tF + sign(z - zf) .* dFocus / c;

    % ---- receive delays for all elements at once ------------------------
    dxe   = x(e) - ex(:).';                   % [1 x nEl]
    tauRx = hypot(repmat(dxe, nz, 1), z - ez(:).') / c;   % [nz x nEl]
    tau   = tauTx + tauRx;

    % ---- dynamic receive aperture ---------------------------------------
    halfAp = z / (2*o.FNumber);               % [nz x 1]
    inAp   = abs(dxe) <= halfAp;              % [nz x nEl] by broadcasting
    switch lower(o.Apod)
        case 'boxcar'
            w = double(inAp);
        case 'hanning'
            w = zeros(nz, nEl);
            R = repmat(dxe, nz, 1) ./ repmat(halfAp, 1, nEl);
            w(inAp) = 0.5 * (1 + cos(pi * R(inAp)));
        otherwise
            error('das_focused:apod', 'Unknown apodisation: %s', o.Apod);
    end

    % ---- demodulate this event only -------------------------------------
    % Done inside the loop on purpose: pre-computing the baseband for all 128
    % events would need several hundred MB of complex doubles.
    bb = hilbert(D.rf(:, :, e)) .* carrierDown;    % [nS x nEl]

    fidx = (tau - t0)/dt + 1;
    i0   = floor(fidx);
    frac = fidx - i0;
    ok   = (i0 >= 1) & (i0 < nS) & (w > 0);

    val = complex(zeros(nz, nEl));
    if any(ok(:))
        % Column offsets turn the per-element index into a linear index into bb
        colOff = repmat((0:nEl-1)*nS, nz, 1);
        li     = i0 + colOff;
        val(ok) = bb(li(ok)) .* (1 - frac(ok)) + bb(li(ok) + 1) .* frac(ok);
    end
    val = val .* exp(1j*2*pi*fd*tau);

    iq(:, e)      = sum(w .* val, 2);
    nActive(:, e) = sum(ok, 2);
end

% ---- outputs ------------------------------------------------------------
[x, ord] = sort(x);
B = struct();
B.x   = x(:);
B.z   = z;
B.iq  = iq(:, ord);
B.rf  = real(B.iq);
B.env = abs(B.iq);
B.nActive = nActive(:, ord);
B.params = struct('fNumber', o.FNumber, 'apod', string(o.Apod), ...
                  'dz', dz, 'c', c, 'fs', fs, 'fc', D.fc, 'fDemod', fd, ...
                  'nEvents', nEv, 'focus_m', D.tx(1).focus, ...
                  'dataFile', D.dataFile, 'preSetFile', D.preSetFile, ...
                  'frame', D.frame, 'cIsAssumed', D.cIsAssumed, ...
                  'tfSpread_ns', max(tfSpread)*1e9);

if o.Verbose
    fprintf('DAS focused\n');
    fprintf('  grid    : %d x %d  (dz %.1f um, line spacing %.1f um)\n', ...
        nz, numel(x), dz*1e6, median(diff(x))*1e6);
    fprintf('  extent  : x %.1f..%.1f mm,  z %.1f..%.1f mm\n', ...
        x(1)*1e3, x(end)*1e3, z(1)*1e3, z(end)*1e3);
    fprintf('  focus   : %.1f mm\n', D.tx(1).focus*1e3);
    fprintf('  f-number: %.2f (%s apodisation)\n', o.FNumber, o.Apod);
    fprintf('  aperture: %d..%d elements per pixel\n', ...
        min(nActive(:)), max(nActive(:)));
    fprintf('  focus check: max spread of arrival times %.1f ns', ...
        max(tfSpread)*1e9);
    if max(tfSpread) > 1/(4*D.fc)
        fprintf('   >> LARGE (> lambda/4). Delays, focus or c disagree.\n');
    else
        fprintf('   (ok, < lambda/4)\n');
    end
end
end
