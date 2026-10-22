% QUS_CUTOFF  Sweep the SWIFT detector threshold.
%
%   The weight map is built by flagging blocks whose first-pass delta-BSC
%   exceeds ratioCutOff, and down-weighting them in the second pass. The
%   inherited threshold of 10 dB comes from the lab's simulation script,
%   where contrasts are far larger than in a real phantom.
%
%   HOW THE RESULT IS JUDGED. Not by which threshold brings the ACS closest
%   to the certificate: choosing a parameter by its distance to the value
%   being validated invalidates the validation. The criterion here is
%   spatial. The detector is doing its job when the flagged blocks trace the
%   boundary of the inclusion, and it has been pushed too far when they
%   scatter over the whole map. That is read off the masks, not off the
%   final number, which is why every mask is saved and plotted.
%
%   The first pass does not depend on the threshold, so it is solved once
%   and only the weight map and the second pass repeat.
%
%   Run from the amejia folder, after qus_all has built bf_cache_PW.mat.
%
%   Plain ASCII on purpose.

clear; clc;

% ============================= SWEEP =====================================
CUTOFFS = [10 8 6 4];    % dB. First-pass delta-BSC reaches ~2.9
% =========================================================================

% ---- frozen parameters, identical to qus_all ----------------------------
P.fNumber = 2;    P.apod = 'hanning';
P.muB = 1e3;      P.muC = 1e1;
P.blockWL = [8 12];  P.overlap = 0.8;
P.freqL = 4e6;    P.freqH = 9e6;
P.acsRef = 0.5246;   P.gammaRef = 1.0019;
P.reject = 0.1;   P.extension = 3;
P.zBF = [3e-3 58e-3];   P.dz = 1515/(2*27.7778e6);   P.tol = 1e-3;
C_SOUND = 1542;

ANA = [-1.8 1.8 0.6 4.5];
TG  = struct('x', 0.6, 'z', [1.4 2.4]);

ROOT   = pwd;
NptodB = log10(exp(1))*20;
OUTDIR = fullfile(ROOT,'figs_cutoff');
if ~isfolder(OUTDIR), mkdir(OUTDIR); end

% T3 is the control: its certified backscatter contrast is -0.2 dB, so if
% the detector responds to contrast rather than to noise, this is the target
% where it should barely fire and where no arc should appear in the mask.
CASES = { 'T7', 'PLANE/RFdata_26_T7_-August-2026_15-53-50.mat', 0.95, +10.1
          'T8', 'PLANE/RFdata_26_T8_-August-2026_15-53-50.mat', 0.95,  -9.0
          'T3', 'PLANE/RFdata_26_T3_-August-2026_15-53-50.mat', 0.74,  -0.2 };

L0   = load(fullfile(ROOT,'bf_cache_PW.mat'), 'Bref');
Bref = L0.Bref;

R = struct([]);

for ic = 1:size(CASES,1)
    name = CASES{ic,1};
    cert = CASES{ic,3};
    fprintf('\n===== %s  (cert %.2f, delta-BSC %+.1f dB) =====\n', ...
            name, cert, CASES{ic,4});

    D  = load_vsx_rf(fullfile(ROOT, CASES{ic,2}), 'c', C_SOUND, 'Verbose', false);
    Bs = das_planewave(D, 'zLim', P.zBF, 'dz', P.dz, ...
                       'FNumber', P.fNumber, 'Apod', P.apod, 'Verbose', false);
    clear D

    sam1 = Bs.rf;
    xBm  = Bs.x(:)' * 100;
    zBm  = Bs.z(:)   * 100;
    fs   = Bs.params.fs;

    rfRef = zeros([size(sam1) numel(Bref)]);
    for j = 1:numel(Bref), rfRef(:,:,j) = Bref{j}.rf; end

    wl = C_SOUND / ((P.freqL + P.freqH)/2);
    bp.xInf = ANA(1); bp.xSup = ANA(2); bp.zInf = ANA(3); bp.zSup = ANA(4);
    bp.overlap = P.overlap;  bp.blocksize = P.blockWL * wl;

    [Sp, Sd, x_ACS, z_ACS, f] = getSpectrum(sam1,  xBm, zBm, fs, bp);
    [SpR,SdR,~,     ~,     ~] = getSpectrum(rfRef, xBm, zBm, fs, bp);

    dz_cm = zBm(2) - zBm(1);
    nz    = round(bp.blocksize(2)*100 / dz_cm);
    L     = (nz/2) * dz_cm;

    att_ref_map(1,1,:) = P.acsRef * f.^P.gammaRef / NptodB;
    comp = (log(SpR) - log(SdR)) - 4*L*att_ref_map;
    sld  =  log(Sp)  - log(Sd);

    rng = f > P.freqL/1e6 & f < P.freqH/1e6;
    b   = sld(:,:,rng) - comp(:,:,rng);
    ufr = f(rng);

    [m,n,p] = size(b);
    A1   = kron(4*L*ufr, speye(m*n));
    A2   = kron(ones(size(ufr)), speye(m*n));
    mask = ones(m*n*p,1);

    % ---- first pass: independent of the threshold, solved once ----------
    [~,Cn]  = optimAdmmTvTikhonov(A1,A2,b(:),P.muB,P.muC,m,n,P.tol,mask);
    bscMap1 = reshape(Cn*NptodB,m,n);
    fprintf('first-pass delta-BSC spans %.2f to %.2f dB\n', ...
            min(bscMap1(:)), max(bscMap1(:)));

    [Xa, Za] = meshgrid(x_ACS, z_ACS);
    sel = abs(Xa) < TG.x & Za > TG.z(1) & Za < TG.z(2);

    acs = nan(1,numel(CUTOFFS)); sdv = acs; pct = acs;
    masks = cell(1,numel(CUTOFFS)); maps = cell(1,numel(CUTOFFS));

    fprintf('%8s %8s %8s %9s\n','cutoff','ACS','sd','penal%');
    for iu = 1:numel(CUTOFFS)
        cut = CUTOFFS(iu);

        w    = (1-P.reject)*(abs(bscMap1) < cut) + P.reject;
        wExt = movmin(w, P.extension);
        pct(iu) = 100*mean(wExt(:) < 1);

        W = repmat(wExt,[1 1 p]);
        W = spdiags(W(:),0,m*n*p,m*n*p);
        [Bn2,~] = optimAdmmWeightedTvTikhonov(W*A1, W*A2, W*b(:), ...
                      P.muB, P.muC, m, n, P.tol, mask, wExt);
        acsMap = reshape(Bn2*NptodB,m,n);

        v = acsMap(sel); v = v(isfinite(v));
        acs(iu) = median(v); sdv(iu) = std(v);
        masks{iu} = wExt < 1;  maps{iu} = acsMap;

        fprintf('%8.1f %8.3f %8.3f %9.1f\n', cut, acs(iu), sdv(iu), pct(iu));
    end

        % ---- figure: the masks are the evidence, not the numbers ------------
    nShow = numel(CUTOFFS);
    fig = figure('Color','w','Position',[40 40 240*(nShow+1) 760]);
    t = tiledlayout(fig, 2, nShow+1, ...
                    'TileSpacing','tight', 'Padding','compact');

    nexttile;
    imagesc(x_ACS, z_ACS, bscMap1); axis image; colormap(gca,turbo);
    title('\DeltaBSC 1a pasada','FontSize',9);
    ylabel('Axial [cm]'); set(gca,'FontSize',8); colorbar;

    for iu = 1:nShow
        nexttile;
        imagesc(x_ACS, z_ACS, masks{iu}, [0 1]);
        axis image; colormap(gca, gray); set(gca,'FontSize',8,'YTickLabel',[]);
        title({sprintf('cut %.1f', CUTOFFS(iu)), ...
               sprintf('%.1f%%', pct(iu))}, 'FontSize',9);
    end

    nexttile; axis off;
    for iu = 1:nShow
        nexttile;
        imagesc(x_ACS, z_ACS, maps{iu}, [0 1.2]);
        axis image; colormap(gca, turbo); set(gca,'FontSize',8,'YTickLabel',[]);
        title(sprintf('ACS %.3f', acs(iu)),'FontSize',9);
        xlabel('Lateral [cm]','FontSize',8);
    end

    
    %exportgraphics(fig, fullfile(OUTDIR,['cutoff_' name '.png']), 'Resolution', 180);
    %close(fig);

    R(ic).name = name; R(ic).cut = CUTOFFS; R(ic).acs = acs;
    R(ic).sd = sdv; R(ic).pct = pct; R(ic).masks = masks;
    R(ic).x_ACS = x_ACS; R(ic).z_ACS = z_ACS; R(ic).bsc1 = bscMap1;

    clear Bs rfRef sld b A1 A2 W
end

save(fullfile(ROOT,'qus_cutoff.mat'), 'R', 'P', 'CUTOFFS', '-v7.3');

fprintf(['\n\nCOMO LEER LAS FIGURAS\n' ...
 '  Fila de arriba: los bloques que el detector castiga en cada umbral.\n' ...
 '  El umbral es adecuado cuando esos bloques dibujan el borde de la\n' ...
 '  inclusion. Es demasiado bajo cuando aparecen repartidos por todo el\n' ...
 '  mapa, porque entonces castiga ruido y no estructura.\n' ...
 '  Fila de abajo: el mapa de ACS resultante. Sirve para ver el efecto,\n' ...
 '  NO para elegir el umbral.\n']);
