% BF_ASTRID  Beamform the eight new CIRS acquisitions (L11-5v, focus 30 mm).
%
%   This block was acquired with a different transducer and a different
%   sampling rate than the rest of the dataset, and it comes with no
%   reference phantom. So this script does NOT estimate attenuation: the
%   reference method needs a sample and a reference that share the system
%   response, and there is none here.
%
%   What it does is reconstruct the eight files with the same beamformer
%   settings used throughout the study, so they can be compared, image by
%   image, against the externally beamformed version when it arrives.
%
%   Everything that depends on the transducer (t0, dz, element positions,
%   sampling rate) is read from each preSet, so no constant is hard-coded.
%
%   Run from the amejia folder:  bf_astrid
%
%   Outputs:
%       bf_astrid.mat        beamformed RF and axes for all eight
%       figs_astrid/*.png    B-mode per target
%
%   Plain ASCII on purpose.

clear; clc;

NEWDIR = fullfile(pwd, 'ASTRID');
OUTDIR = fullfile(pwd, 'figs_astrid');
OUTMAT = fullfile(pwd, 'bf_astrid.mat');

% ---- beamformer settings, same as the frozen set ------------------------
C_SOUND = 1542;          % CIRS 049 certificate
F_NUM   = 2;
APOD    = 'hanning';
Z_BF    = [3e-3 60e-3];  % focus is at 30 mm, so reconstruct deeper than before
DR      = 60;            % B-mode dynamic range, dB

if ~isfolder(OUTDIR), mkdir(OUTDIR); end

files = dir(fullfile(NEWDIR, '**', '*.mat'));
files = files(~[files.isdir]);

% keep only the data files
keep = false(numel(files),1);
for k = 1:numel(files)
    try
        w = whos('-file', fullfile(files(k).folder, files(k).name));
        keep(k) = any(strcmp({w.name},'RcvData'));
    catch
        keep(k) = false;
    end
end
files = files(keep);
fprintf('%d data files found\n\n', numel(files));

B = struct([]);

for k = 1:numel(files)
    f = fullfile(files(k).folder, files(k).name);
    [~, base] = fileparts(files(k).name);
    fprintf('[%d/%d] %s\n', k, numel(files), files(k).name);

    try
        D = load_vsx_rf(f, 'c', C_SOUND, 'Verbose', false);
    catch ME
        fprintf('   load failed: %s\n', ME.message);
        continue
    end

    % dz is left at its native value, c/(2*fs). It is NOT forced to match
    % the rest of the dataset: these files are never divided by a reference
    % from another transducer, so there is no grid to align with. When the
    % externally beamformed version arrives, resample whichever of the two
    % is needed onto the other's axes before comparing.
    try
        Bk = das_focused(D, 'zLim', Z_BF, 'FNumber', F_NUM, ...
                         'Apod', APOD, 'Verbose', false);
    catch ME
        fprintf('   beamforming failed: %s\n', ME.message);
        fprintf('   (128 focused transmits: check that das_focused handles\n');
        fprintf('    a line-by-line sweep with this TX layout)\n');
        continue
    end
    clear D

    dz_um = (Bk.z(2)-Bk.z(1))*1e6;
    fprintf('   %d x %d | dz %.2f um | x %.1f to %.1f mm | z %.1f to %.1f mm\n', ...
            size(Bk.rf,1), size(Bk.rf,2), dz_um, ...
            Bk.x(1)*1e3, Bk.x(end)*1e3, Bk.z(1)*1e3, Bk.z(end)*1e3);

    % ---- B-mode ---------------------------------------------------------
    env = abs(hilbert(Bk.rf));
    bm  = 20*log10(env / max(env(:)));

    fig = figure('Visible','off','Position',[100 100 520 760],'Color','w');
    imagesc(Bk.x*100, Bk.z*100, bm, [-DR 0]);
    axis image; colormap gray; colorbar;
    xlabel('Lateral [cm]'); ylabel('Axial [cm]');
    title(sprintf('%s   |   L11-5v, foco 30 mm, f/%.0f %s', ...
          strrep(base,'_','\_'), F_NUM, APOD));
    out = fullfile(OUTDIR, ['bmode_' base '.png']);
    exportgraphics(fig, out, 'Resolution', 200);
    close(fig);

    B(k).name  = base;
    B(k).rf    = Bk.rf;
    B(k).x     = Bk.x;
    B(k).z     = Bk.z;
    B(k).bmode = bm;
    B(k).params = Bk.params;
end

B = B(~arrayfun(@(s) isempty(s.name), B));
meta = struct('c', C_SOUND, 'fNumber', F_NUM, 'apod', APOD, 'zLim', Z_BF);
save(OUTMAT, 'B', 'meta', '-v7.3');

fprintf('\n%d reconstructed\n', numel(B));
fprintf('Saved: %s\n       %s\n', OUTMAT, OUTDIR);
fprintf(['\nNOTE\n' ...
 '  No attenuation is estimated here. The reference method divides the\n' ...
 '  sample spectrum by a reference acquired with the same system, and this\n' ...
 '  block has none. Dividing by the L11-4v Wisconsin would mostly measure\n' ...
 '  the difference between the two transducers.\n']);
