% ======================================================================
% ======================================================================
% phantoms_dt2809.m - Merino's estimator on the lab-beamformed DataTargets2809
% (Task 03, part A). Derived from merino/phantoms_repro.m. Changes against it,
% all required by the different data and nothing else:
%  a. paths: Samples/BeamformedData and Reference/BeamformedData/<ref>;
%  b. variable names: lab files store rf (2D), not RF;
%  c. region x_inf/x_sup = -1.85/+1.85 cm (lab grid is centred on 0);
%  d. c1x, c1z per target measured on the lab B-mode by circle fit
%     (merino/centres_dt2809.m -> results/merino_dt2809/centres.csv);
%  e. reference ACS 0.54 (261) or 0.53 (544), outer loop over both refs;
%  f. all eight targets; ground truth per target from the CIRS 049
%     certificate. Merino tuned mu only for T6, T7, T8 (his cases 1, 2, 3);
%     T1-T5 take his case-2 values (equal to case 3 for SWIFT). Not chosen by
%     closeness to the certificate: there is no published value for them;
%  g. metrics carry ref and file; figures saved per reference.
%  h. (Task 04) input set and SWIFT override from environment variables,
%     because setup clears the workspace:
%       QUS_BFSET  'lab' (default, lab BeamformedData) or 'own'
%                  (results/bf_own, das_focused with each phantom's c)
%       QUS_SWIFT  '[log10(muB) log10(muC) ratioCutOff]' applied to SWIFT
%                  for every target; empty = Merino's values
%       QUS_TAG    subfolder of the results folder for this run
%     ACS and first-pass BSC maps are saved to maps.mat.
% Everything else (band, windows, overlap, NFFT, weights, masks) unchanged.
%
% Inherited header of phantoms_repro.m:
% REPRODUCTION COPY of SWIFT_UFFC2024/phantoms.m (commit 2f37a01).
% Method unchanged. Changes are environment-only:
%  1. setup runs with pwd = Merino repo, so genpath(pwd) adds only Merino
%     code (running it from the project root would also add SWIFT_Simple/,
%     which holds same-name copies of spectra, AlterOpti_ADMM, etc.).
%  2. config.m is loaded AFTER setup, because setup calls "clear".
%  3. dataDir, refDir, resultsDir point to targets_edu via DATA_ROOT.
%  4. resultsDir is fixed to the project root, not the Merino repo.
%% PHANTOMSSS
cd(fullfile(fileparts(mfilename('fullpath')), '..', 'SWIFT_UFFC2024'));
setup,
cd(fullfile(fileparts(mfilename('fullpath')), '..'));
run(fullfile(pwd, 'config.m'));
% dataDir = ['C:\Users\sebas\Documents\Data\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
% refDir = ['C:\Users\sebas\Documents\Data\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];
% resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\UFFC2024results\phantoms';

bfSet = getenv('QUS_BFSET'); if isempty(bfSet), bfSet = 'lab'; end   % (h)
switch bfSet
    case 'lab'
        dataDir    = fullfile(DATA_ROOT, 'DataTargets2809', 'Samples', 'BeamformedData');
        refRoot    = fullfile(DATA_ROOT, 'DataTargets2809', 'Reference', 'BeamformedData');
        resultsDir = fullfile(pwd, 'results', 'merino_dt2809');
    case 'own'
        dataDir    = fullfile(DATA_ROOT, 'results', 'bf_own', 'Samples');
        refRoot    = fullfile(DATA_ROOT, 'results', 'bf_own', 'Reference');
        resultsDir = fullfile(pwd, 'results', 'merino_dt2809_own');
    otherwise
        error('Unknown QUS_BFSET %s', bfSet);
end
refList    = {'261', 0.54; '544', 0.53};               % (e)
centres    = readtable(fullfile(resultsDir, 'centres.csv'));   % (d)
swiftOverride = str2num(getenv('QUS_SWIFT')); %#ok<ST2NM>  (h)
runTag = getenv('QUS_TAG');
if ~isempty(runTag), resultsDir = fullfile(resultsDir, runTag); end
fprintf('Input set %s, SWIFT override %s, results %s\n', bfSet, ...
    mat2str(swiftOverride), resultsDir);

rawFiles = dir([dataDir,'\*.rf']);
targetFiles = dir([dataDir,'\*.mat']);
% (f) all eight targets: targetFiles = targetFiles(end-2:end);
if ~exist("resultsDir","dir"); mkdir(resultsDir); end
tableName = 'phantoms_dt2809.xlsx';

%% Constants
blocksize = 8;     % Block size in wavelengths
freq_L = 2.5e6; freq_H = 7.5e6;
freq_C = (freq_L + freq_H)/2;

overlap_pc      = 0.8;
ratio_zx        = 12/8;
x_inf = -1.85; x_sup = 1.85;   % (c)
z_inf = 0.2; z_sup = 3.5;
NptodB = log10(exp(1))*20;

% Weight parameters
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

groundTruthTargets = [0.52,0.55,0.74,0.81,0.75,0.97,0.95,0.95,0.55];   % (f) T1..T8, back

% Plotting constants
dynRange = [-50,0];
attRange = [0.2,1.2];

tol = 1e-3;

c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;
%% For looping each phantom

for iRef0 = 1:size(refList,1)                              % (e)
refTag = refList{iRef0,1};  acsRef = refList{iRef0,2};
refDir = fullfile(refRoot, refTag);
for iAcq = 1:numel(targetFiles)                              % (f)
iRow = (iRef0-1)*numel(targetFiles) + iAcq;
% (f) Merino's cases: T6 -> 1, T7 -> 2, T8 -> 3; T1..T5 -> 2
muCase = 2;
if iAcq >= 6, muCase = iAcq - 5; end
switch muCase
    % Optimal reg for BS 8x12, circular ROIs
    case 1
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^2.5;
        muBtvl1 = 10^3.5; muCtvl1 = 10^2;
        muBswift = 10^3.5; muCswift = 10^2;
    case 2
        muBtv = 10^3; muCtv = 10^3;
        muBswtv = 10^3; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
    case 3
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
end
if ~isempty(swiftOverride)                                   % (h)
    muBswift = 10^swiftOverride(1); muCswift = 10^swiftOverride(2);
    ratioCutOff = swiftOverride(3);
end

% (d) centre measured on the lab B-mode, replaces Merino's per-case values
c1x = centres.c1x(iAcq); c1z = centres.c1z(iAcq);

%%
fprintf("Ref %s, phantom no. %i, %s\n",refTag,iAcq,targetFiles(iAcq).name);
load(fullfile(dataDir,targetFiles(iAcq).name));
RF = rf;   % (b)

dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = RF(:,:,1);

%% Cropping and finding sample sizes

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
roi = ind_x.*ind_z';
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Wavelength size
c0 = 1540;
wl = c0/mean(freq_C);   % Wavelength (m)

% Lateral samples
wx = round(blocksize*wl*(1-overlap_pc)/dx);  % Between windows
nx = round(blocksize*wl/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blocksize*wl*(1-overlap_pc)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize*wl/dz /2 * ratio_zx); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);



%% Generating Diffraction compensation
if true %iAcq == 1
    % Generating references
    att_ref = acsRef*f/NptodB; % (e) From phantom especifications
    att_ref_map = zeros(m,n,p);
    for jj=1:n
        for ii=1:m
            att_ref_map(ii,jj,:) = att_ref;
        end
    end
    
    % Windows for spectrum
    windowing = hamming(nz/2);
    windowing = windowing*ones(1,nx);
    
    % For looping
    refFiles = dir([refDir,'\*.mat']);
    Nref = length(refFiles);
    
    % Memory allocation
    Sp_ref = zeros(m,n,p,Nref);
    Sd_ref = zeros(m,n,p,Nref);
    for iRef = 1:Nref
        out = load([refDir,'\',refFiles(iRef).name]);
        samRef = out.rf;   % (b)
        samRef = samRef(ind_z,ind_x); % Cropping
        % figure,imagesc(db(hilbert(samRef)))
        for jj=1:n
            for ii=1:m
                xw = x0(jj) ;   % x window
                zp = z0p(ii);
                zd = z0d(ii);
    
                sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
                sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
                [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
                [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
    
                Sp_ref(ii,jj,:,iRef) = (tempSp(rang));
                Sd_ref(ii,jj,:,iRef) = (tempSd(rang));
            end
        end
    end
    
    Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
    compensation = ( log(Sp) - log(Sd) ) - 4*L*att_ref_map;
    
    % Liberating memory to avoid killing my RAM
    clear Sp_ref Sd_ref
end

%% Spectrum
Sp = zeros(m,n,p);
Sd = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = sam1(zd:zd+nz/2-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

%% ROI selection
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
rInc = 0.95;
inc = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
back = ((Xq-c1x).^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;

x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
% [back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

% Setting up
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
mask = ones(m,n,p);

%% RSLD-TV

tic
[Bn,Cn,ite] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
exTime = toc;
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BR,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'TV';
r.ref = refTag; r.file = targetFiles(iAcq).name;   % (g)
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsTV(iRow) = r;

%% British Columbia Approach

% Weights
envelope = abs(hilbert(sam1));
SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);

        temp = [sub_block_p(:) sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% computation
tic
[Bn,Cn,ite] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv, ...
    m,n,tol,mask(:),wSNR);
exTime = toc;
BRBC = (reshape(Bn*NptodB,m,n));
CRBC = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BRBC,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'SWTV';
r.ref = refTag; r.file = targetFiles(iAcq).name;   % (g)
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsSWTV(iRow) = r;

%% SWIFT
% First iteration
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Weight map
w = (1-reject)*(abs(bscMap)<ratioCutOff)+reject;
wExt = movmin(w,extension);

% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
tic
[Bn,~,ite] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
exTime = toc;
BSWIFT = reshape(Bn*NptodB,m,n);

AttInterp = interp2(X,Z,BSWIFT,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'SWIFT';
r.ref = refTag; r.file = targetFiles(iAcq).name;   % (g)
r.sample = iAcq;
r.ite = ite;
r.exTime = exTime;
MetricsSWIFT(iRow) = r;
Maps(iRow) = struct('file', targetFiles(iAcq).name, 'ref', refTag, ...   % (h)
    'x_ACS', x_ACS, 'z_ACS', z_ACS, 'TV', BR, 'SWTV', BRBC, 'SWIFT', BSWIFT, ...
    'bscSWIFT', bscMap, 'w', w, 'c1x', c1x, 'c1z', c1z, 'x', x, 'z', z, ...
    'Bmode', Bmode, 'muBswift', muBswift, 'muCswift', muCswift, ...
    'ratioCutOff', ratioCutOff);

%%
figure('Units','centimeters', 'Position',[5 2 6 16]);
tl = tiledlayout(4,1, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
%subtitle(' ')
c = colorbar;
c.Label.String = 'dB';
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

% fontsize(gcf,8,'points')

t2 = nexttile;
imagesc(x_ACS,z_ACS,BR, attRange)
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('RSLD')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off

t3 = nexttile;
imagesc(x_ACS,z_ACS,BRBC, attRange)
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title('SWTV-ACE')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off

t5 = nexttile;
imagesc(x_ACS,z_ACS,BSWIFT, attRange)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('SWIFT')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';

% fontsize(gcf,8,'points')
hold on 
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
hold off
end
save_all_figures_to_directory(resultsDir,['phantom_ref' refTag '_'],'svg');   % (g)
close all
end

%%
%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results4 = struct2table(MetricsSWIFT);


T = [results1;results2;results4];
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);
save(fullfile(resultsDir,'maps.mat'),'Maps');   % (h)
%%
