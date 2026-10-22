% Script to test my simulations
clear,
clc,
close all;

dataDir = 'Q:\csoto\ACSSim\SOSPhantomSim\Variable';

% targetFiles = dir(fullfile(dataDir,"rf_prebf_sosI*Att*.mat"));
targetFiles = dir(fullfile(dataDir, '*.mat'));
targetFiles = targetFiles(~contains({targetFiles.name}, '.jpg'));


refDir = 'Q:\csoto\ACSSim\SOSPhantomSimRef\Ref1540';
% refFiles = dir(fullfile(refDir,'rf_prebf_sosI*Att*.mat'));
refFiles = dir(fullfile(refDir, '*.mat'));
refFiles = refFiles(~contains({refFiles.name}, '.jpg'));

resultsDir = 'Q:\csoto\ACSSim\ResultsACSVar\ResultsACS_Ref1540muB3.5';

mkdir(resultsDir)

%%
% idealAcs = 0.5;
c0 = 1540;
freqL = 1e6; freqH = 6e6; % [Original] freqL = 1e6; freqH = 5e6;
wl = c0/((freqL + freqH)/2);
alpha0Ref = 0.3; gammaRef = 1.00;
deadband = 0.05; % [cm]

blockParams.xInf = -1.8; % 0.8
blockParams.xSup = 1.8;
blockParams.zInf = 0.4;
blockParams.zSup = 5;
% blockParams.blocksize = [15 15]*wl;
% blockParams.blocksize = [15 20]*wl;
blockParams.blocksize = [8 12]*wl;
% blockParams.blocksize = [20 25]*wl;
blockParams.overlap = 0.8;

NptodB = log10(exp(1))*20;
% iAcq = 2;
acs_prom=0;
spectraResults = cell(length(targetFiles), 4);

%%
for iAcq = 1:length(targetFiles)
% for iAcq = 1
    %% Loading sample
    samName = targetFiles(iAcq).name(1:end);
    % blockParams.zInf = 0.5;
    % idealAcs = str2double(strrep(samName(19:21),'p','.'));
    % idealAcs= str2double(strrep(regexp(samName, '(?<=acs).{3}', 'match', 'once'), 'p', '.'));
    idealAcs = 0.5;

    fprintf("%i., %s\n",iAcq,targetFiles(iAcq).name);
    out = load(fullfile(dataDir,targetFiles(iAcq).name));
    % % out2 = load(fullfile(dataDir,targetFiles_noadm(iAcq).name));
    % xBm = out.lateral_positions_cm; % [cm]
    % zBm = out.axial_positions_cm'; % [cm]

    % Accomodations for CSG Simulations
    x = out.xr;
    z = out.zr;

    xBm = x'*1e2; % [cm]
    zBm = z'*1e2; % [cm]
    zValid = zBm>deadband;

    % sam1 = out.rf(zValid,:);
    % zBm = zBm(zValid);

    sam1 = out.rf;
    % zValid = zBm>deadband;
    % sam1 = sam1(zValid,:);
    % zBm = zBm(zValid);
    % sam=squeeze(sum(sam,2));
   
    % sam1=squeeze(sum(out.reconstructed_channel_data,2));
    % sam1 = sam1(zValid,:);
    % sam1(sam1 == 0) = noadm(sam1 == 0);
 
    env_adm=hilbert(sam1);
    env_adm=20*log10(abs(sam1));
    env_adm=env_adm-max(env_adm(:));
    % 
    % env_adm1=hilbert(sam1);
    % env_adm1=20*log10(abs(sam1(500:end,:)));
    % env_adm1=env_adm1-max(env_adm1(:));
    % 
    figure,
    imagesc(env_adm), colormap('gray'), clim([-55 0]), title('BMode')
    figure,
    plot(zBm, mean(env_adm, 2)), title('Envelope (dB)')
    % 
    % figure,
    % imagesc(env_adm1), colormap('gray'), clim([-60 0]), title('Bmode ADMIRE')

    fs = out.fs;

    % out = load(fullfile(mapsDir,samName(4:end)));
    % out.acsMap=idealAcs*ones([10390,64]);
    [xSim,zSim] = meshgrid(xBm/100,zBm/100);
    % attenuation_map = out.acsMap;
    % attenuation_map=attenuation_map(zSim);
      
    attenuation_map = idealAcs*ones(size(xSim));
    % attenuation_map(zSim<2e-2) = 1.3;
    % attenuation_map(zSim<1e-2) = 1;

    clear out

    % Plot region of interest B-mode image
    bMode = db(hilbert(sam1));
    bMode = bMode - max(bMode(:));

    % get Spectra
    [Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam1,xBm,zBm,fs,blockParams);

    % Plotting spectra
    spectrumSamzf = db(squeeze(mean(Sp/2+Sd/2, 2)));
    spectrumSamzf = spectrumSamzf - max(spectrumSamzf,[],2);
    figure,
    imagesc(f,z_ACS,spectrumSamzf, [-50 0]),
    ax = gca; ax.YDir = 'reverse';
    hold on
    xline(freqL/1e6, 'w--', 'LineWidth',2), xline(freqH/1e6, 'w--', 'LineWidth',2)
    hold off
    colorbar
    xlim([0 8])
    xlabel('f [MHz]')
    ylabel('z [cm]')
    title('Sample power spectrum by depth')

    %% Generating Diffraction compensation
    % Generating references
    clear att_ref_map rfRef

    for ff = 1:length(refFiles)
        out = load(fullfile(refDir,refFiles(ff).name));
        % out2=load(fullfile(refDir,files_noadm(ff).name));
        % rfRef1 = squeeze(sum(out.reconstructed_channel_data,2));
        % rfRef1=squeeze(sum(out2.beamformed_channel_data,2));
        % rfRef1(rfRef1 == 0) = noadm(rfRef1 == 0);
        % rfRef(:,:,ff)=rfRef1;
        rfRef(:,:,ff)=out.rf;
    end

    att_ref_map(1,1,:) = alpha0Ref*f.^gammaRef/NptodB;

    [SpRef,SdRef,~,~,~] = getSpectrum(rfRef,xBm,zBm,fs,blockParams);

    % Plotting spectra
    spectrumRefzf = db(squeeze(mean(SpRef/2+SdRef/2, 2)));
    spectrumRefzf = spectrumRefzf - max(spectrumRefzf,[],2);
    figure,
    imagesc(f,z_ACS,spectrumRefzf, [-50 0]),
    ax = gca; ax.YDir = 'reverse';
    hold on
    xline(freqL/1e6, 'w--', 'LineWidth',2), xline(freqH/1e6, 'w--', 'LineWidth',2)
    hold off
    colorbar
    xlim([0 8])
    xlabel('f [MHz]')
    ylabel('z [cm]')
    title('Reference power spectrum by depth')

    %% Setting up
    L = (z_ACS(2) - z_ACS(1))/(1 - blockParams.overlap)/2;   % (cm)
    sld = log(Sp) - log(Sd);
    sldRef = log(SpRef) - log(SdRef);
    compensation = sldRef - 4*L*att_ref_map;

    range = f>freqL/1e6 & f<freqH/1e6;
    b = sld(:,:,range) - compensation(:,:,range);
    ufr = f(range);

    [m,n,p] = size(b);
    A1 = kron( 4*L*ufr , speye(m*n) );
    A2 = kron( ones(size(ufr)) , speye(m*n) );
    mask = ones(m,n,p);
    tol = 1e-3;

    %% SLD by depth
    nLines = 3;
    depthPoints = floor(m/nLines);
    figure, hold on
    for ii = 1:nLines
        line = mean(b(ii*depthPoints-depthPoints+1:ii*depthPoints,:,:),[1 2]);
        zline = mean(z_ACS(ii*depthPoints-depthPoints+1:ii*depthPoints));
        plot(ufr,squeeze(line*NptodB/4/L), 'LineWidth',2)
        leg{ii} = "z = "+round(zline,2)+"cm";
    end
    plot(ufr,idealAcs*ufr, 'k--')
    hold off
    grid on
    leg{ii+1} = 'Ideal';
    legend(leg, 'Location','northwest')
    xlabel('Frequency [MHz]')
    ylabel('Attenuation [dB/cm]')


    %% RSLD-TV
    % muRsld = 10e2                      ;
    muRsld = 10e2                      ;

    tic
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muRsld,muRsld,m,n,tol,mask(:));
    toc
    BRSLD = (reshape(Bn*NptodB,m,n));
    CRSLD = (reshape(Cn*NptodB,m,n));

    % Creating masks and ideal map
    [Xq,Zq] = meshgrid(xBm,zBm);
    % attIdeal = interp2(xSim*100, zSim*100, attenuation_map, Xq, Zq);
    % attIdeal = interp2(xSim, zSim, attenuation_map, Xq, Zq);
    attIdeal=idealAcs;
    roi = ones(size(bMode));

    %%
    fit = [z_ACS,ones(m,1)]\mean(BRSLD,2);
    figure,
    plot(z_ACS,mean(BRSLD,2))
    grid on
    xlabel('Depth [cm]')
    ylabel('ACS [dB/cm/MHz]')
    yline(idealAcs,'k--')
    % title(sprintf("RSLD, ACS = %.2fz + %.2f",fit(1),fit(2)))
    title("RSLD")
    ylim([0 2])
    hold on,
    % plot(z_ACS,fit(1)*z_ACS + fit(2), 'r--')

    %% SWIFT
    % muBswift = 10^3; muCswift = 10^2;
    muBswift = 10^3.5; muCswift = 10^2;
    % muBswift = 10^4; muCswift = 10^1;
    % muBswift = 10^2; muCswift = 10^1;
    ratioCutOff = 10;
    reject = 0.1;
    extension = 3;

    % Plotting constants
    dynRange = [-60, 0];
    attRange = [0, 1.5];
    bsRange = [-10, 10];
    yLimits = [0, 5];

    % First iteration
    [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol,mask(:));
    bscMap = reshape(Cn*NptodB,m,n);
    
    % Weight map
    w = (1-reject)*(abs(bscMap)<ratioCutOff)+reject;
    % w = (1-reject)*(1 ./ (1 + (bscMap./ratioCutOff).^10))+reject;
    wExt = movmin(w,extension);
    
    % Weight matrices and new system
    W = repmat(wExt,[1 1 p]);
    W = spdiags(W(:),0,m*n*p,m*n*p);
    bw = W*b(:);        
    A1w = W*A1;
    A2w = W*A2;
    
    % Second iteration
    [Bn,cN] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
    BSWIFT = reshape(Bn*NptodB,m,n);
    CSWIFT = reshape(Cn*NptodB,m,n);

    fit = [z_ACS,ones(m,1)]\mean(BSWIFT,2);
    figure,
    plot(z_ACS,mean(BSWIFT,2))
    grid on
    xlabel('Depth [cm]')
    ylabel('ACS [dB/cm/MHz]')
    yline(idealAcs,'k--')
    % title(sprintf("SWIFT, ACS = %.2fz + %.2f",fit(1),fit(2)))
    title("SWIFT")
    ylim([0 2])
    hold on,
    yline(idealAcs, 'k--')
    hold off
    % plot(z_ACS,fit(1)*z_ACS + fit(2), 'r--')
%%
    figure('Units','centimeters', 'Position',[5 5 24 8]);
    tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

    t1 = nexttile;
    imagesc(xBm,zBm,bMode,dynRange)
    axis image
    c = colorbar(t1, 'westoutside');
    c.Label.String = 'dB';
    title('B-mode')
    ylabel('Axial [cm]')
    xlabel('Lateral [cm]')
    ylim(yLimits)

    % t2 = nexttile;
    % imagesc(xBm,zBm,attIdeal,attRange)
    % xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
    % axis image
    % title('Ideal')
    % ylim(yLimits)

    t3 = nexttile;
    [~,hB,hColor] = imOverlayInterp(bMode,BRSLD,dynRange,attRange,1,...
        x_ACS,z_ACS,roi,xBm,zBm);
    title('RSLD, ACS')
    xlabel('Lateral [cm]')
    ylim(yLimits)

    t4 = nexttile;
    [~,hB,hColor] = imOverlayInterp(bMode,BSWIFT,dynRange,attRange,1,...
        x_ACS,z_ACS,roi,xBm,zBm);
    title('SWIFT, ACS')
    xlabel('Lateral [cm]')
    ylim(yLimits)
    
    colormap(t1,gray)
    % colormap(t2,turbo)
    colormap(t3,turbo)
    colormap(t4,turbo)

    % acs_prom=[acs_prom mean(BSWIFT,'all')];
    acsResults{iAcq, 1} = samName;
    acsResults{iAcq, 2} = BSWIFT;
    acsResults{iAcq, 3} = x_ACS;
    acsResults{iAcq, 4} = z_ACS;
    acsResults{iAcq, 5} = xBm;
    acsResults{iAcq, 6} = zBm;
    acsResults{iAcq, 7} = bMode;
    acsResults{iAcq, 8} = dynRange;
    acsResults{iAcq, 9} = attRange;



    %%
    save_all_figures_to_directory(resultsDir,char(samName+"_fig"))

    close all
    save(fullfile(resultsDir,'acs_prom'),'acs_prom');
    
end
% save(fullfile(resultsDir,'acs_prom'),'acs_prom');
% figure, plot([2, 3, 4],acs_prom(2:end),'-o'), ylabel('ACS'), xlabel('Abdominal Wall (cm)'), hold on, plot([2, 3, 4],1.5*ones([3]),'--' ), ylim([0.3 2])
%  save_all_figures_to_directory(resultsDir,char(samName+"_prom"))