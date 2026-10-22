function [Sp,Sd,x_ACS,z_ACS,f] = get_spectrum_opt(sam1,x,z,fs,blockParams)
%GET_SPECTRUM_OPT  Copy of lab/SWIFT_Simple/getSpectrum.m with two options exposed.
%
%   Identical to getSpectrum except:
%     blockParams.window    'tukey' (default, = getSpectrum: tukeywin 0.25)
%                           or 'hamming' (= Merino phantoms.m)
%     blockParams.nfftExtra 1 (default, = getSpectrum) or 2 (= Merino):
%                           NFFT = 2^(nextpow2(nz/2)+nfftExtra)
%     blockParams.verbose   false silences the block-size printout
%   With the defaults the output is bit-identical to getSpectrum (checked in
%   scripts/own_qus_edu.m).
if ~isfield(blockParams,'window'),    blockParams.window = 'tukey'; end
if ~isfield(blockParams,'nfftExtra'), blockParams.nfftExtra = 1;    end
if ~isfield(blockParams,'verbose'),   blockParams.verbose = true;   end

dx = (x(2)-x(1)) /100;
dz = (z(2)-z(1)) /100;

x_inf = blockParams.xInf;
x_sup = blockParams.xSup;
z_inf = blockParams.zInf;
z_sup = blockParams.zSup;
overlap_pc = blockParams.overlap;
blocksize = blockParams.blocksize;

%% Cropping and finding sample sizes
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x,:);

% Lateral samples
wx = round(blocksize(1)*(1-overlap_pc)/dx);  % Between windows
nx = round(blocksize(1)/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blocksize(2)*(1-overlap_pc)/dz); % Between windows
nz = 2*round(blocksize(2)/dz /2); % Window size
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+blockParams.nfftExtra);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
f  = band*1e-6; % [MHz]
p = length(f);

if blockParams.verbose
    fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
    fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
    fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);
end

%% Spectrum
switch lower(blockParams.window)
    case 'tukey',   windowing = tukeywin(nz/2,0.25);
    case 'hamming', windowing = hamming(nz/2);
    otherwise, error('get_spectrum_opt:window','Unknown window %s',blockParams.window);
end
windowing = windowing*ones(1,nx);

nSamples = size(sam1,3);
Sp_ref = zeros(m,n,p,nSamples);
Sd_ref = zeros(m,n,p,nSamples);
for iRef = 1:nSamples
    samRef = sam1(:,:,iRef);
    for jj=1:n
        for ii=1:m
            xw = x0(jj) ;   % x window
            zp = z0p(ii);
            zd = z0d(ii);

            sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
            sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
            [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
            [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);

            Sp_ref(ii,jj,:,iRef) = tempSp;
            Sd_ref(ii,jj,:,iRef) = tempSd;
        end
    end
end

Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
end
