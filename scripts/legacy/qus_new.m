% QUS_NEW  Process the DataTargets2809 block: 8 targets, 2 reference phantoms.
%
%   This block is self-contained: references and samples share the same
%   transducer, sampling rate and transmit scheme, so the reference method
%   is valid inside it. It cannot be mixed with the earlier acquisitions,
%   which used a different fc and fs.
%
%   TWO REFERENCES, NOT ONE SET. They are different phantoms, so their
%   certified ACS differs and they are never averaged. Every sample is
%   estimated twice, once against each. That is a free cross-check: if both
%   agree, the choice of reference introduces no bias; if they differ, the
%   difference is itself a measurement.
%
%   WITHOUT THE CERTIFICATES the script still works. What it reports is
%   delta-ACS, the difference against the reference, which needs no
%   certified value. Fill ACS_261 and ACS_544 when the certificates arrive
%   and it will also report absolute ACS.
%
%   Two detector thresholds are run in the same pass. The first SWIFT pass
%   does not depend on the threshold, so it is solved once per case.
%
%   Run from the amejia folder:  qus_new
%
%   Plain ASCII on purpose.

clear; clc;

ROOT = pwd;
BASE = fullfile(ROOT, 'DataTargets2809');

% ===================== TO FILL IN WHEN THE CERTIFICATES ARRIVE ===========
ACS_261 = NaN;     % certified ACS of reference phantom 261, dB/cm/MHz
ACS_544 = NaN;     % certified ACS of reference phantom 544
GAMMA   = 1.0;     % frequency exponent; 1.0 until the certificate says otherwise
% =========================================================================

% ---- frozen parameters --------------------------------------------------
% Band stays at 4-9 MHz: this is the same L11-4v as the earlier data, so the
% physical bandwidth of the element is unchanged even though the preSet
% reports a different reference frequency.
P.c        = 1540;          % as used by the lab for their own beamforming
P.fNumber  = 2;
P.apod     = 'hanning';
P.muB      = 1e3;
P.muC      = 1e1;
P.blockWL  = [8 12];        % [lateral axial], wavelengths
P.overlap  = 0.8;
P.freqL    = 4e6;           % restored to the current configuration (was 7e6)
P.freqH    = 9e6;
P.reject   = 0.1;
P.extension = 3;
P.zBF      = [3e-3 60e-3];
P.tol      = 1e-3;
ACS_261 = 0.54;
ACS_544 = 0.53;

CUTOFFS = [10 8 6];           % detector thresholds to run

ANA  = [-1.8 1.8 0.6 4.5];              % analysis region [xI xS zI zS], cm
TGT  = struct('x', 0.6, 'z', [1.4 2.4]);% target region, inside the inclusion
BGX  = 1.2;                             % background: |x| beyond this

OUTDIR = fullfile(ROOT,'figs_new');
if ~isfolder(OUTDIR), mkdir(OUTDIR); end
NptodB = log10(exp(1))*20;

% ---- file layout --------------------------------------------------------
% Each reference is reconstructed with ITS OWN speed of sound, from its
% certificate. Using 1540 for both put the 261 image on a depth axis 2 %
% off, since that phantom is 1509 m/s.
REFS = { '261', fullfile(BASE,'Reference','ChannelData','261'), ACS_261, 1509
         '544', fullfile(BASE,'Reference','ChannelData','544'), ACS_544, 1539 };

SAMPDIR = fullfile(BASE,'Samples','ChannelData');
SAMPLES = { 'T1','T1_F.mat'
            'T2','T2_F.mat'
            'T3','T3_F_3.mat'
            'T4','T4_F.mat'
            'T5','T5_F.mat'
            'T6','T6_F.mat'
            'T7','T7_F.mat'
            'T8','T8_F_3.mat' };

% Common axial grid. dz defaults to c/(2*fs); fixing it keeps every
% reconstruction on the same rows so sample and reference blocks sit at the
% same depths. The speed of sound still enters where it belongs, in the
% delays.
P.dz = P.c / (2 * 31.25e6);

bf = @(D) das_focused(D, 'zLim', P.zBF, 'dz', P.dz, ...
                      'FNumber', P.fNumber, 'Apod', P.apod, 'Verbose', false);

% ================= beamform the references once, cached ==================
for ir = 1:size(REFS,1)
    tag   = REFS{ir,1};
    cRef  = REFS{ir,4};                                                    % <- nueva
    cache = fullfile(ROOT, ['bf_new_ref' tag '_c' num2str(cRef) '.mat']);  % <- con sufijo
    if isfile(cache), fprintf('Cached references %s (c=%d)\n', tag, cRef); continue; end

    d = dir(fullfile(REFS{ir,2}, '*.mat'));
    d = d(~contains({d.name},'preSet') & ~contains({d.name},'preSave'));
    fprintf('Beamforming %d references, phantom %s (c=%d)...\n', numel(d), tag, cRef);

    Bref = cell(1, numel(d));
    for j = 1:numel(d)
        fprintf('  [%d/%d] %s\n', j, numel(d), d(j).name);
        Dj = load_vsx_rf(fullfile(d(j).folder, d(j).name), 'c', cRef, 'Verbose', false);  % <- cRef
        Bref{j} = bf(Dj);
        clear Dj
    end
    save(cache, 'Bref', '-v7.3');
end

% ============================= main loop =================================
R = struct([]);
row = 0;

for is = 1:size(SAMPLES,1)
    sName = SAMPLES{is,1};
    sFile = fullfile(SAMPDIR, SAMPLES{is,2});
    if ~isfile(sFile), warning('Missing %s -- skipped', SAMPLES{is,2}); continue; end

    fprintf('\n########## %s ##########\n', sName);
    Ds = load_vsx_rf(sFile, 'c', P.c, 'Verbose', false);
    Bs = bf(Ds); clear Ds

    sam1 = Bs.rf;
    xBm  = Bs.x(:)' * 100;
    zBm  = Bs.z(:)   * 100;
    fs   = Bs.params.fs;

    % ---- B-mode, saved once per sample ----------------------------------
    env = abs(hilbert(sam1)); bmode = 20*log10(env/max(env(:)));
    fg = figure('Visible','off','Position',[100 100 520 760],'Color','w');
    imagesc(xBm, zBm, bmode, [-60 0]); axis image; colormap gray; colorbar;
    xlabel('Lateral [cm]'); ylabel('Axial [cm]'); title(sName);
    hold on;
    rectangle('Position',[-TGT.x TGT.z(1) 2*TGT.x diff(TGT.z)], ...
              'EdgeColor',[0 1 0],'LineWidth',1.8);
    exportgraphics(fg, fullfile(OUTDIR,['bmode_' sName '.png']),'Resolution',180);
    close(fg);

    for ir = 1:size(REFS,1)
        tag    = REFS{ir,1};
        acsRef = REFS{ir,3};
        L0 = load(fullfile(ROOT, ['bf_new_ref' tag '_c' num2str(REFS{ir,4}) '.mat']), 'Bref');
        Bref = L0.Bref;

        if ~isequal(size(Bref{1}.rf), size(sam1))
            error('Grid mismatch vs ref %s: %dx%d vs %dx%d', tag, ...
                  size(sam1,1), size(sam1,2), size(Bref{1}.rf,1), size(Bref{1}.rf,2));
        end

        rfRef = zeros([size(sam1) numel(Bref)]);
        for j = 1:numel(Bref), rfRef(:,:,j) = Bref{j}.rf; end

        wl = P.c / ((P.freqL + P.freqH)/2);
        bp.xInf = ANA(1); bp.xSup = ANA(2); bp.zInf = ANA(3); bp.zSup = ANA(4);
        bp.overlap = P.overlap;  bp.blocksize = P.blockWL * wl;

        [Sp, Sd, x_ACS, z_ACS, f] = getSpectrum(sam1,  xBm, zBm, fs, bp);
        [SpR,SdR,~,     ~,     ~] = getSpectrum(rfRef, xBm, zBm, fs, bp);

        % axial half-length, from the block size in samples (the rounded
        % block step gave a 2.9 % scale error in the inherited code)
        dz_cm = zBm(2) - zBm(1);
        nz    = round(bp.blocksize(2)*100 / dz_cm);
        L     = (nz/2) * dz_cm;

        % With acsRef = NaN the reference is not compensated and the result
        % is delta-ACS. With a certified value it becomes absolute ACS.
        if isnan(acsRef)
            comp = log(SpR) - log(SdR);
            kind = 'dACS';
        else
            att_ref_map(1,1,:) = acsRef * f.^GAMMA / NptodB;
            comp = (log(SpR) - log(SdR)) - 4*L*att_ref_map;
            kind = 'ACS';
        end

        sld = log(Sp) - log(Sd);
        rngf = f > P.freqL/1e6 & f < P.freqH/1e6;
        b    = sld(:,:,rngf) - comp(:,:,rngf);
        ufr  = f(rngf);

        [m,n,p] = size(b);
        A1   = kron(4*L*ufr, speye(m*n));
        A2   = kron(ones(size(ufr)), speye(m*n));
        mask = ones(m*n*p,1);

        % first pass: does not depend on the threshold
        [~,Cn]  = optimAdmmTvTikhonov(A1,A2,b(:),P.muB,P.muC,m,n,P.tol,mask);
        bsc1 = reshape(Cn*NptodB,m,n);

        [Xa, Za] = meshgrid(x_ACS, z_ACS);
        selT = abs(Xa) < TGT.x & Za > TGT.z(1) & Za < TGT.z(2);
        selB = abs(Xa) > BGX;

        for iu = 1:numel(CUTOFFS)
            cut = CUTOFFS(iu);
            w    = (1-P.reject)*(abs(bsc1) < cut) + P.reject;
            wExt = movmin(w, P.extension);
            pct  = 100*mean(wExt(:) < 1);

            W = repmat(wExt,[1 1 p]);
            W = spdiags(W(:),0,m*n*p,m*n*p);
            [Bn2,Cn2] = optimAdmmWeightedTvTikhonov(W*A1, W*A2, W*b(:), ...
                            P.muB, P.muC, m, n, P.tol, mask, wExt);
            acsMap = reshape(Bn2*NptodB,m,n);
            bscMap = reshape(Cn2*NptodB,m,n);

            vT = acsMap(selT); vT = vT(isfinite(vT));
            vB = acsMap(selB); vB = vB(isfinite(vB));

            row = row + 1;
            R(row).sample = sName;   R(row).ref = tag;   R(row).cut = cut;
            R(row).kind   = kind;
            R(row).tg     = median(vT);  R(row).tgSd = std(vT);
            R(row).tgP10  = prctile(vT,10); R(row).tgP90 = prctile(vT,90);
            R(row).bg     = median(vB);  R(row).bgSd = std(vB);
            R(row).pct    = pct;
            R(row).acsMap = acsMap;  R(row).bscMap = bscMap;
            R(row).x_ACS  = x_ACS;   R(row).z_ACS  = z_ACS;

            fprintf('  ref %s  cut %4.1f  |  %s target %.3f +- %.3f  fondo %.3f +- %.3f  |  penal %.1f %%\n', ...
                    tag, cut, kind, median(vT), std(vT), median(vB), std(vB), pct);
        end

        clear rfRef sld b A1 A2 W
    end
    clear Bs
end

save(fullfile(ROOT,'qus_new_f49.mat'), 'R', 'P', 'CUTOFFS', '-v7.3');

% ============================== table ====================================
fid = fopen(fullfile(ROOT,'qus_new_f49_table.csv'),'w');
fprintf(fid,'sample,ref,cutoff,kind,target,target_sd,p10,p90,background,background_sd,penalised_pct\n');
fprintf('\n\n%-5s %5s %6s %6s %9s %8s %9s %8s %8s\n', ...
        'tgt','ref','cut','kind','target','sd','fondo','sd','penal%');
for i = 1:numel(R)
    fprintf('%-5s %5s %6.1f %6s %9.3f %8.3f %9.3f %8.3f %8.1f\n', ...
        R(i).sample, R(i).ref, R(i).cut, R(i).kind, ...
        R(i).tg, R(i).tgSd, R(i).bg, R(i).bgSd, R(i).pct);
    fprintf(fid,'%s,%s,%.1f,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f\n', ...
        R(i).sample, R(i).ref, R(i).cut, R(i).kind, ...
        R(i).tg, R(i).tgSd, R(i).tgP10, R(i).tgP90, R(i).bg, R(i).bgSd, R(i).pct);
end
fclose(fid);

% ---- the cross-check: do the two references agree? ----------------------
fprintf('\n\nDIFERENCIA ENTRE REFERENCIAS (261 menos 544, mismo target y umbral)\n');
fprintf('%-5s %6s %10s %10s %10s\n','tgt','cut','ref 261','ref 544','dif');
for is = 1:size(SAMPLES,1)
    for iu = 1:numel(CUTOFFS)
        a = R(strcmp({R.sample},SAMPLES{is,1}) & strcmp({R.ref},'261') & [R.cut]==CUTOFFS(iu));
        b2= R(strcmp({R.sample},SAMPLES{is,1}) & strcmp({R.ref},'544') & [R.cut]==CUTOFFS(iu));
        if isempty(a) || isempty(b2), continue; end
        fprintf('%-5s %6.1f %10.3f %10.3f %10.3f\n', ...
                SAMPLES{is,1}, CUTOFFS(iu), a.tg, b2.tg, a.tg - b2.tg);
    end
end

fprintf(['\nCOMO LEERLO\n' ...
 '  Sin certificados, la columna es delta-ACS: la diferencia contra el\n' ...
 '  phantom de referencia, no el valor absoluto. Al llenar ACS_261 y\n' ...
 '  ACS_544 pasa a ser ACS y la etiqueta cambia sola.\n' ...
 '  La ultima tabla es la validacion cruzada. Como 261 y 544 son phantoms\n' ...
 '  distintos, la diferencia entre sus delta-ACS deberia ser constante en\n' ...
 '  los ocho targets, e igual a la diferencia de ACS entre los dos\n' ...
 '  phantoms. Si varia por target, esa variacion es del metodo.\n']);
