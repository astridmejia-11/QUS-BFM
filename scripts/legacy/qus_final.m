% QUS_FINAL  Main results, with the corrections the advisors asked for.
%
%   Changes with respect to qus_sweep:
%
%   1. Reports mean +/- STANDARD DEVIATION instead of IQR, to be comparable
%      with the NRMSE values in Merino's thesis and with the literature.
%   2. Region widened upward (zInf 1.0 -> 0.6 cm). The target spans roughly
%      1.0 to 3.2 cm, so its upper edge used to sit exactly on the boundary of
%      the analysis region, which is where the total-variation prior has the
%      fewest neighbours and therefore the least support.
%   3. Draws the analysis region and both measurement ROIs on the B-mode, so
%      the reported numbers can be checked against where they were taken.
%   4. Outlines the inclusion on the B-mode.
%   5. Outputs the delta-BSC map for every case, labelled with the method.
%   6. Adds a control run restricted to the background, to quantify how much
%      the presence of the target contaminates the background estimate through
%      the regularisation.
%
%   Run from the amejia folder:  qus_final
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT    = pwd;
SCHEMES = {'PW','FOC'};
TARGETS  = {'T1','T3','T5','T7','T8'};
ACS_CERT = [0.52 0.74 0.75 0.95 0.95];       % CIRS 049, serie 211270316
DBSC_CERT= [11.0 -0.2 -8.7 10.1 -9.0];       % dB respecto al fondo, a 3 MHz
ACS_BG   = 0.55;

C_SOUND = 1542; Z_BF = [3e-3 58e-3];
ACS_REF = 0.5246; GAMMA_REF = 1.0019;        % UW report 07/09/2025

X_INF = -1.8; X_SUP = 1.8;
Z_INF =  0.6; Z_SUP =  5.0;                  % AMPLIADA hacia arriba
X_TG  = 0.6;  Z_TG  = [1.4 2.4];
X_BG  = 1.2;
INCL_CENTER = [0.0 1.95]; INCL_RADIUS = 0.95; % delineado de la inclusion [cm]

FREQ_L = 4e6; FREQ_H = 9e6;
OVERLAP = 0.8; BLOCK_WL = [8 12];
FNUM = 2.0; APOD = 'hanning';
MU_RSLD = 10^3; MU_B = 10^3; MU_C = 10;
RATIO_CUT = 10; REJECT = 0.1; EXTENSION = 3;

DO_CONTROL = true;   % corrida restringida al fondo, sin el target en la region
% =========================================================================

dataOnly = @(d) d(~contains({d.name},'preSet') & ~contains({d.name},'preSave'));
Rall = struct([]); maps = {};

for iS = 1:numel(SCHEMES)
    scheme = SCHEMES{iS};
    switch scheme
        case 'PW'
            sampleDir=fullfile(ROOT,'PLANE'); refDir=fullfile(ROOT,'Referencia Plane');
            bf=@(D) das_planewave(D,'zLim',Z_BF,'FNumber',FNUM,'Apod',APOD,'Verbose',false);
            cacheFile=fullfile(ROOT,'bf_final_PW.mat');
        case 'FOC'
            sampleDir=fullfile(ROOT,'TARGETS'); refDir=fullfile(ROOT,'Ref Targets ATS');
            bf=@(D) das_focused(D,'zLim',Z_BF,'FNumber',FNUM,'Apod',APOD,'Verbose',false);
            cacheFile=fullfile(ROOT,'bf_final_FOC.mat');
    end

    if isfile(cacheFile)
        fprintf('[%s] referencias en cache\n', scheme);
        L0=load(cacheFile,'Bref'); Bref=L0.Bref;
    else
        rfl = dataOnly(dir(fullfile(refDir,'*.mat')));
        fprintf('[%s] beamformando %d referencias...\n', scheme, numel(rfl));
        Bref = cell(1,numel(rfl));
        for k=1:numel(rfl)
            Dk=load_vsx_rf(fullfile(rfl(k).folder,rfl(k).name),'c',C_SOUND,'Verbose',false);
            Bref{k}=bf(Dk); clear Dk
        end
        save(cacheFile,'Bref','-v7.3');
    end

    sf = dataOnly(dir(fullfile(sampleDir,'*.mat')));

    for iT = 1:numel(TARGETS)
        tg = TARGETS{iT};
        hit = find(contains({sf.name},['_' tg '_'])|contains({sf.name},['_' tg '-']),1);
        if isempty(hit), continue; end

        fprintf('\n--- %s / %s ---\n', scheme, tg);
        Ds=load_vsx_rf(fullfile(sf(hit).folder,sf(hit).name),'c',C_SOUND,'Verbose',false);
        Bs=bf(Ds); clear Ds

        O = runCase(Bs, Bref, C_SOUND, X_INF, X_SUP, Z_INF, Z_SUP, ...
            FREQ_L, FREQ_H, ACS_REF, GAMMA_REF, OVERLAP, BLOCK_WL, ...
            MU_RSLD, MU_B, MU_C, RATIO_CUT, REJECT, EXTENSION);

        [Xa,Za] = meshgrid(O.x_ACS, O.z_ACS);
        mBG = abs(Xa) > X_BG;
        mTG = abs(Xa) < X_TG & Za > Z_TG(1) & Za < Z_TG(2);
        ms  = @(v) [mean(v(isfinite(v))) std(v(isfinite(v)))];

        r.scheme=scheme; r.target=tg; r.cert=ACS_CERT(iT); r.dbscCert=DBSC_CERT(iT);
        v=ms(O.BSWIFT(mBG)); r.swBg=v(1); r.swBgSd=v(2);
        v=ms(O.BSWIFT(mTG)); r.swTg=v(1); r.swTgSd=v(2);
        v=ms(O.BRSLD(mBG));  r.rsBg=v(1); r.rsBgSd=v(2);
        v=ms(O.BRSLD(mTG));  r.rsTg=v(1); r.rsTgSd=v(2);
        r.dbsc = mean(O.CSWIFT(mTG),'omitnan') - mean(O.CSWIFT(mBG),'omitnan');
        r.ctrlBg = NaN;

        % ---- control: solo fondo en la region -----------------------------
        % Con regularizacion el resultado depende de que hay en la region, no
        % solo de donde se mide. Estimar sin el target adentro cuantifica
        % cuanto contamina su presencia al valor del fondo.
        if DO_CONTROL
            Oc = runCase(Bs, Bref, C_SOUND, X_BG, X_SUP, Z_INF, Z_SUP, ...
                FREQ_L, FREQ_H, ACS_REF, GAMMA_REF, OVERLAP, BLOCK_WL, ...
                MU_RSLD, MU_B, MU_C, RATIO_CUT, REJECT, EXTENSION);
            vv = Oc.BSWIFT(:); r.ctrlBg = mean(vv(isfinite(vv)));
        end
        Rall(end+1).scheme = r.scheme; %#ok<SAGROW>
        fn=fieldnames(r); for q=1:numel(fn), Rall(end).(fn{q})=r.(fn{q}); end

        maps{end+1} = struct('scheme',scheme,'target',tg,'O',O); %#ok<SAGROW>

        fprintf('  SWIFT fondo %.3f +/- %.3f | target %.3f +/- %.3f (cert %.2f)\n', ...
            r.swBg, r.swBgSd, r.swTg, r.swTgSd, r.cert);
        if DO_CONTROL
            fprintf('  control solo fondo: %.3f  (diferencia %+.3f)\n', ...
                r.ctrlBg, r.ctrlBg - r.swBg);
        end

        % ---- figura por caso ----------------------------------------------
        figure('Position',[40 40 1280 380]);
        subplot(1,3,1);
        imagesc(O.xBm,O.zBm,O.bMode,[-60 0]); axis image; colormap(gca,gray); colorbar;
        hold on;
        plot([X_INF X_SUP X_SUP X_INF X_INF],[Z_INF Z_INF Z_SUP Z_SUP Z_INF],...
             'y--','LineWidth',1.3);
        plot([-X_TG X_TG X_TG -X_TG -X_TG],[Z_TG(1) Z_TG(1) Z_TG(2) Z_TG(2) Z_TG(1)],...
             'c-','LineWidth',1.5);
        plot([-X_SUP -X_BG -X_BG -X_SUP -X_SUP],[Z_INF Z_INF Z_SUP Z_SUP Z_INF],...
             'g-','LineWidth',1.5);
        plot([X_BG X_SUP X_SUP X_BG X_BG],[Z_INF Z_INF Z_SUP Z_SUP Z_INF],...
             'g-','LineWidth',1.5);
        th=linspace(0,2*pi,200);
        plot(INCL_CENTER(1)+INCL_RADIUS*cos(th), INCL_CENTER(2)+INCL_RADIUS*sin(th), ...
             'r--','LineWidth',1.5);
        hold off;
        xlabel('Lateral [cm]'); ylabel('Axial [cm]');
        title(sprintf('%s %s | amarillo analisis, cyan target, verde fondo, rojo inclusion', ...
              scheme, tg), 'FontSize',8);

        subplot(1,3,2);
        imagesc(O.x_ACS,O.z_ACS,O.BSWIFT,[0 1.5]); axis image; colormap(gca,turbo); colorbar;
        xlabel('Lateral [cm]'); title(sprintf('SWIFT, ACS [dB/cm/MHz] - %s %s',scheme,tg));

        subplot(1,3,3);
        imagesc(O.x_ACS,O.z_ACS,O.CSWIFT,[-15 15]); axis image; colormap(gca,turbo); colorbar;
        xlabel('Lateral [cm]');
        title(sprintf('SWIFT, delta-BSC [dB] - cert %+.1f',DBSC_CERT(iT)));
    end
end

save(fullfile(ROOT,'qus_final_results.mat'),'Rall','maps','-v7.3');

% ---- tablas --------------------------------------------------------------
fprintf('\n\n========== TABLA PRINCIPAL (media +/- desv. estandar) ==========\n');
fprintf('Referencia: Wisconsin %.4f dB/cm/MHz (UW 07/09/2025)\n', ACS_REF);
fprintf('Fondo certificado CIRS: %.2f dB/cm/MHz\n\n', ACS_BG);
fprintf('%-4s %-4s %6s %16s %8s %16s %8s\n', ...
    'sch','tgt','cert','target SWIFT','Bias','fondo SWIFT','Bias');
for k=1:numel(Rall)
    R=Rall(k);
    fprintf('%-4s %-4s %6.2f %8.3f+/-%.3f %7.1f%% %8.3f+/-%.3f %7.1f%%\n', ...
        R.scheme,R.target,R.cert,R.swTg,R.swTgSd, ...
        100*(R.swTg-R.cert)/R.cert, R.swBg,R.swBgSd, ...
        100*(R.swBg-ACS_BG)/ACS_BG);
end

fprintf('\n========== RSLD contra SWIFT ==========\n');
fprintf('%-4s %-4s %16s %16s\n','sch','tgt','RSLD fondo','SWIFT fondo');
for k=1:numel(Rall)
    R=Rall(k);
    fprintf('%-4s %-4s %8.3f+/-%.3f %8.3f+/-%.3f\n', ...
        R.scheme,R.target,R.rsBg,R.rsBgSd,R.swBg,R.swBgSd);
end

if DO_CONTROL
    fprintf('\n========== CONTAMINACION POR REGULARIZACION ==========\n');
    fprintf('Fondo medido del mapa completo, contra fondo estimado solo.\n\n');
    fprintf('%-4s %-4s %12s %12s %10s\n','sch','tgt','mapa completo','solo fondo','dif');
    for k=1:numel(Rall)
        R=Rall(k);
        fprintf('%-4s %-4s %12.3f %12.3f %+10.3f\n', ...
            R.scheme,R.target,R.swBg,R.ctrlBg,R.ctrlBg-R.swBg);
    end
end

fprintf('\n========== CONTRASTE DE BSC ==========\n');
fprintf('%-4s %-4s %14s %14s\n','sch','tgt','medido [dB]','certificado');
for k=1:numel(Rall)
    R=Rall(k);
    fprintf('%-4s %-4s %14.2f %14.1f\n',R.scheme,R.target,R.dbsc,R.dbscCert);
end


% #########################################################################
function O = runCase(Bs, Bref, c0, xInf, xSup, zInf, zSup, freqL, freqH, ...
    acsRef, gammaRef, overlap, blockWl, muRsld, muB, muC, ratioCut, reject, extension)

sam1=Bs.rf; xBm=Bs.x(:)'*100; zBm=Bs.z(:)*100; fs=Bs.params.fs;
rfRef=zeros(size(sam1,1),size(sam1,2),numel(Bref));
for k=1:numel(Bref), rfRef(:,:,k)=Bref{k}.rf; end

wl=c0/((freqL+freqH)/2); NptodB=log10(exp(1))*20;
bp.xInf=xInf; bp.xSup=xSup; bp.zInf=zInf; bp.zSup=zSup;
bp.overlap=overlap; bp.blocksize=blockWl*wl;

evalc('[Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam1,xBm,zBm,fs,bp);');
evalc('[SpRef,SdRef,~,~,~]   = getSpectrum(rfRef,xBm,zBm,fs,bp);');

% Separacion entre ventana proximal y distal, tomada de la geometria de la
% ventana y no del paso entre bloques, que se redondea a muestras enteras y
% arrastra un 2.9% de error de escala directo a la atenuacion.
dz_cm=zBm(2)-zBm(1);
nz=2*round(bp.blocksize(2)/(dz_cm/100)/2);
L=(nz/2)*dz_cm;

att_ref_map(1,1,:)=acsRef*f.^gammaRef/NptodB;
sld=log(Sp)-log(Sd); sldRef=log(SpRef)-log(SdRef);
compensation=sldRef-4*L*att_ref_map;
range=f>freqL/1e6 & f<freqH/1e6;
b=sld(:,:,range)-compensation(:,:,range); ufr=f(range);

[m,n,p]=size(b);
A1=kron(4*L*ufr,speye(m*n)); A2=kron(ones(size(ufr)),speye(m*n));
mask=ones(m,n,p); tol=1e-3;

[Bn,~]=AlterOpti_ADMM(A1,A2,b(:),muRsld,muRsld,m,n,tol,mask(:));
O.BRSLD=reshape(Bn*NptodB,m,n);

[~,Cn]=optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
bscMap=reshape(Cn*NptodB,m,n);
w=(1-reject)*(abs(bscMap)<ratioCut)+reject;
wExt=movmin(w,extension);
W=repmat(wExt,[1 1 p]); W=spdiags(W(:),0,m*n*p,m*n*p);
[Bn2,Cn2]=optimAdmmWeightedTvTikhonov(W*A1,W*A2,W*b(:),muB,muC,m,n,tol,mask(:),wExt);

O.BSWIFT=reshape(Bn2*NptodB,m,n);
O.CSWIFT=reshape(Cn2*NptodB,m,n);
O.wExt=wExt;
O.x_ACS=x_ACS; O.z_ACS=z_ACS; O.xBm=xBm; O.zBm=zBm;
bM=20*log10(abs(hilbert(sam1))); O.bMode=bM-max(bM(:));
end
