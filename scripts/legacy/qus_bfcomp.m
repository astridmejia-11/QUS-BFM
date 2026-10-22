% QUS_BFCOMP  Two beamformers, the same raw echoes, the same estimator.
%
%   This is the direct form of the project's question. Every earlier test
%   varied settings inside one beamformer; here two independent
%   implementations reconstruct the same channel data, and each is divided
%   by references reconstructed by ITSELF. That pairing matters: the
%   reference method cancels the system response only when sample and
%   reference share a beam, so mixing beamformers across the ratio would
%   measure the mismatch instead of the tissue.
%
%   PREDICTION, WRITTEN BEFORE RUNNING. Their settings are f/2.5 with
%   rectangular apodisation, inside the range already swept in the
%   beamformer study (f/1 to f/3, Hanning and boxcar), where ACS varied less
%   than the spread of a single configuration. So the two pipelines should
%   agree within that spread. The one thing outside the sweep is their
%   dynamic aperture capped at 32 elements, which can change how resolution
%   varies with depth.
%
%   Stage 1 reports what is inside their files. If the fields cannot be
%   mapped automatically, it stops and prints the variable names so they can
%   be set by hand below.
%
%   Run from the amejia folder, after qus_new (it reuses the reference cache).
%
%   Plain ASCII on purpose.

clear; clc;

ROOT = pwd;
BASE = fullfile(ROOT, 'DataTargets2809');

% ---- manual overrides, only if stage 1 cannot map the fields ------------
FIELD_RF = '';    % name of the beamformed image variable, e.g. 'rf' or 'bf'
FIELD_X  = '';    % lateral axis variable
FIELD_Z  = '';    % axial axis variable

% ---- frozen estimator parameters, identical to qus_new ------------------
P.c = 1540;  P.fNumber = 2;  P.apod = 'hanning';
P.muB = 1e3; P.muC = 1e1;    P.blockWL = [8 12];  P.overlap = 0.8;
P.freqL = 4e6;  P.freqH = 9e6;
P.reject = 0.1; P.extension = 3;  P.ratioCut = 10;
P.zBF = [3e-3 60e-3];  P.dz = 1540/(2*31.25e6);  P.tol = 1e-3;

ANA = [-1.8 1.8 0.6 4.5];
TGT = struct('x', 0.6, 'z', [1.4 2.4]);
BGX = 1.2;

NptodB = log10(exp(1))*20;
OUTDIR = fullfile(ROOT,'figs_bfcomp');
if ~isfolder(OUTDIR), mkdir(OUTDIR); end

SAMPLES = { 'T1','T1_F'; 'T2','T2_F'; 'T3','T3_F_3'; 'T4','T4_F'
            'T5','T5_F'; 'T6','T6_F'; 'T7','T7_F'; 'T8','T8_F_3' };
REFTAG  = '261';                      % reference phantom used for both sides
REFNAMES = {'261_F','261_F_2','261_F_3'};

DIR_BF_S = fullfile(BASE,'Samples','BeamformedData');
DIR_BF_R = fullfile(BASE,'Reference','BeamformedData',REFTAG);
DIR_CH_S = fullfile(BASE,'Samples','ChannelData');

% ======================= STAGE 1: what is in their files =================
probe = fullfile(DIR_BF_S, [SAMPLES{1,2} '.mat']);
w = whos('-file', probe);
fprintf('Variables in %s:\n', SAMPLES{1,2});
for i = 1:numel(w)
    fprintf('  %-16s %-10s %s\n', w(i).name, w(i).class, mat2str(w(i).size));
end

% auto-map: the image is the largest 2D numeric array; the axes are the
% vectors whose lengths match its two dimensions
if isempty(FIELD_RF)
    is2d = arrayfun(@(s) numel(s.size)==2 && all(s.size>1), w);
    if ~any(is2d)
        error('No 2D array found. Set FIELD_RF by hand from the list above.');
    end
    cand = w(is2d);
    [~, k] = max(arrayfun(@(s) prod(s.size), cand));
    FIELD_RF = cand(k).name;
    sz = cand(k).size;
    isv = arrayfun(@(s) sum(s.size>1)==1, w);
    vec = w(isv);
    for i = 1:numel(vec)
        n = max(vec(i).size);
        if isempty(FIELD_Z) && n == sz(1), FIELD_Z = vec(i).name; end
        if isempty(FIELD_X) && n == sz(2), FIELD_X = vec(i).name; end
    end
end
fprintf('\nUsing: image = %s | x = %s | z = %s\n', FIELD_RF, ...
        ternstr(FIELD_X), ternstr(FIELD_Z));
if isempty(FIELD_X) || isempty(FIELD_Z)
    fprintf(['\nAxes not found in the file. They will be rebuilt from the\n' ...
             'preSet: z from startDepth with the native dz, x from the\n' ...
             'element positions.\n']);
end

% ======================= references, both sides ==========================
% theirs
RefT = cell(1,numel(REFNAMES));
for j = 1:numel(REFNAMES)
    S = load(fullfile(DIR_BF_R,[REFNAMES{j} '.mat']));
    RefT{j} = S.(FIELD_RF);
end
% mine, from the cache qus_new already built
Lc = load(fullfile(ROOT,['bf_new_ref' REFTAG '.mat']),'Bref');
RefM = Lc.Bref;

% axes for their grid
St = load(fullfile(DIR_BF_R,[REFNAMES{1} '.mat']));
if ~isempty(FIELD_X) && isfield(St,FIELD_X)
    xT = St.(FIELD_X)(:)'*100;  zT = St.(FIELD_Z)(:)*100;
    if max(abs(xT)) > 100, xT = xT/10; zT = zT/10; end   % guard: mm vs cm
else
    Dp = load_vsx_rf(fullfile(DIR_CH_S,[SAMPLES{1,2} '.mat']), ...
                     'c',P.c,'Verbose',false);
    nzT = size(RefT{1},1);
    dzT = P.c/(2*Dp.fs);
    zT  = (Dp.startDepth_m + (0:nzT-1)'*dzT)*100;
    xT  = Dp.elemX(:)'*100;
    clear Dp
end

xM = RefM{1}.x(:)'*100;  zM = RefM{1}.z(:)*100;
fsM = RefM{1}.params.fs;
fsT = fsM;   % same acquisition, so the sampling rate is the same

fprintf('\nTheir grid  %d x %d | z %.2f a %.2f cm\n', ...
        size(RefT{1},1), size(RefT{1},2), zT(1), zT(end));
fprintf('My grid     %d x %d | z %.2f a %.2f cm\n\n', ...
        size(RefM{1}.rf,1), size(RefM{1}.rf,2), zM(1), zM(end));

% ============================ compare ====================================
Rc = struct([]);

for is = 1:size(SAMPLES,1)
    nm = SAMPLES{is,1};
    fprintf('===== %s =====\n', nm);

    % --- their beamformed sample -----------------------------------------
    S = load(fullfile(DIR_BF_S,[SAMPLES{is,2} '.mat']));
    samT = S.(FIELD_RF);

    % --- my beamformed sample --------------------------------------------
    D = load_vsx_rf(fullfile(DIR_CH_S,[SAMPLES{is,2} '.mat']), ...
                    'c',P.c,'Verbose',false);
    Bm = das_focused(D,'zLim',P.zBF,'dz',P.dz,'FNumber',P.fNumber, ...
                     'Apod',P.apod,'Verbose',false);
    samM = Bm.rf;  clear D

    % --- estimate each side against ITS OWN references -------------------
    stT = estimateSide(samT, RefT,  xT, zT, fsT, P, ANA, TGT, BGX, NptodB, true);
    refMrf = cellfun(@(b) b.rf, RefM, 'UniformOutput', false);
    stM = estimateSide(samM, refMrf, xM, zM, fsM, P, ANA, TGT, BGX, NptodB, false);

    fprintf('  suyo   target %.3f +- %.3f | fondo %.3f +- %.3f | penal %.1f %%\n', ...
            stT.tg, stT.tgSd, stT.bg, stT.bgSd, stT.pct);
    fprintf('  mio    target %.3f +- %.3f | fondo %.3f +- %.3f | penal %.1f %%\n', ...
            stM.tg, stM.tgSd, stM.bg, stM.bgSd, stM.pct);
    fprintf('  dif    target %+.3f | fondo %+.3f\n', ...
            stM.tg-stT.tg, stM.bg-stT.bg);

    % --- B-modes side by side --------------------------------------------
    bT = 20*log10(abs(hilbert(samT))); bT = bT-max(bT(:));
    bM = 20*log10(abs(hilbert(samM))); bM = bM-max(bM(:));
    fg = figure('Visible','off','Position',[80 80 980 620],'Color','w');
    tl = tiledlayout(fg,1,2,'TileSpacing','compact','Padding','compact');
    nexttile; imagesc(xT,zT,bT,[-60 0]); axis image; colormap gray;
    title('laboratorio, f/2.5 rectangular'); xlabel('Lateral [cm]'); ylabel('Axial [cm]');
    nexttile; imagesc(xM,zM,bM,[-60 0]); axis image; colormap gray;
    title('propio, f/2 hanning'); xlabel('Lateral [cm]');
    title(tl, nm, 'FontWeight','bold');
    exportgraphics(fg, fullfile(OUTDIR,['bmode_' nm '.png']),'Resolution',180);
    close(fg);

    Rc(is).name = nm;  Rc(is).theirs = stT;  Rc(is).mine = stM;
end

save(fullfile(ROOT,'qus_bfcomp.mat'),'Rc','P','-v7.3');

% ============================== summary ==================================
dTg = arrayfun(@(r) r.mine.tg - r.theirs.tg, Rc);
dBg = arrayfun(@(r) r.mine.bg - r.theirs.bg, Rc);
sdI = arrayfun(@(r) max(r.mine.tgSd, r.theirs.tgSd), Rc);

fprintf('\n\n%-5s %9s %9s %9s %10s\n','tgt','suyo','mio','dif','sd indiv');
for i = 1:numel(Rc)
    fprintf('%-5s %9.3f %9.3f %+9.3f %10.3f\n', ...
            Rc(i).name, Rc(i).theirs.tg, Rc(i).mine.tg, dTg(i), sdI(i));
end
fprintf('\nDiferencia media en target %+.3f (desv %.3f)\n', mean(dTg), std(dTg));
fprintf('Diferencia media en fondo  %+.3f (desv %.3f)\n', mean(dBg), std(dBg));
fprintf('Desviacion individual tipica %.3f\n', mean(sdI));
fprintf(['\nLA PRUEBA\n' ...
 '  Si la diferencia media entre los dos beamformers queda por debajo de\n' ...
 '  la desviacion de una sola estimacion, las dos implementaciones son\n' ...
 '  indistinguibles y la conclusion del estudio se sostiene contra un\n' ...
 '  beamformer escrito por otras personas.\n']);


% ---------------------------------------------------------------- helpers
function s = ternstr(v)
    if isempty(v), s = '(no hallado)'; else, s = v; end
end

function st = estimateSide(sam, refCell, xBm, zBm, fs, P, ANA, TGT, BGX, NptodB, ~)
% One full estimation: spectra, first pass, weight map, second pass.
% Sample and references must come from the SAME beamformer.
    rfRef = zeros([size(sam) numel(refCell)]);
    for j = 1:numel(refCell), rfRef(:,:,j) = refCell{j}; end

    wl = P.c/((P.freqL+P.freqH)/2);
    bp.xInf=ANA(1); bp.xSup=ANA(2); bp.zInf=ANA(3); bp.zSup=ANA(4);
    bp.overlap=P.overlap; bp.blocksize=P.blockWL*wl;

    [Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam,  xBm, zBm, fs, bp);
    [SpR,SdR,~,~,~]       = getSpectrum(rfRef, xBm, zBm, fs, bp);

    dz_cm = zBm(2)-zBm(1);
    nz = round(bp.blocksize(2)*100/dz_cm);
    L  = (nz/2)*dz_cm;

    % no reference certificate available, so this is delta-ACS
    comp = log(SpR)-log(SdR);
    sld  = log(Sp)-log(Sd);
    rngf = f > P.freqL/1e6 & f < P.freqH/1e6;
    b = sld(:,:,rngf)-comp(:,:,rngf);
    ufr = f(rngf);

    [m,n,p] = size(b);
    A1 = kron(4*L*ufr, speye(m*n));
    A2 = kron(ones(size(ufr)), speye(m*n));
    mask = ones(m*n*p,1);

    [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),P.muB,P.muC,m,n,P.tol,mask);
    bsc1 = reshape(Cn*NptodB,m,n);
    w = (1-P.reject)*(abs(bsc1)<P.ratioCut)+P.reject;
    wE = movmin(w,P.extension);
    W = repmat(wE,[1 1 p]);
    W = spdiags(W(:),0,m*n*p,m*n*p);
    [Bn2,~] = optimAdmmWeightedTvTikhonov(W*A1,W*A2,W*b(:), ...
                  P.muB,P.muC,m,n,P.tol,mask,wE);
    acs = reshape(Bn2*NptodB,m,n);

    [Xa,Za] = meshgrid(x_ACS,z_ACS);
    vT = acs(abs(Xa)<TGT.x & Za>TGT.z(1) & Za<TGT.z(2));
    vB = acs(abs(Xa)>BGX);
    vT = vT(isfinite(vT)); vB = vB(isfinite(vB));

    st.tg = median(vT); st.tgSd = std(vT);
    st.bg = median(vB); st.bgSd = std(vB);
    st.pct = 100*mean(wE(:)<1);
    st.acs = acs; st.x_ACS = x_ACS; st.z_ACS = z_ACS;
end
