% FIG_BFCOMP  Two beamformers: the images differ, the ACS maps do not.
%
%   Top row: the B-modes. They look different, and they should: the two
%   reconstructions use different f-number, apodisation and aperture.
%   Bottom row: the ACS maps estimated from those same two reconstructions,
%   each divided by references reconstructed by itself. They agree.
%
%   That contrast is the whole conclusion of the beamformer study in one
%   figure: the beamformer determines how the image looks, not what the
%   method measures.
%
%   Needs qus_bfcomp.mat, written by qus_bfcomp.
%
%   Run from the amejia folder:  fig_bfcomp
%
%   Plain ASCII on purpose.

clear; clc;

TARGET  = 'T7';        % which target to show
ACS_REF = 0.54;        % reference 261 certificate, to plot absolute ACS
ACSLIM  = [0.4 1.6];   % colour scale for the ACS maps
DR      = 60;          % B-mode dynamic range, dB
OUTFILE = ['fig_bfcomp_' TARGET '.png'];

if ~isfile('qus_bfcomp.mat'), error('Run qus_bfcomp first.'); end
S = load('qus_bfcomp.mat','Rc');
k = find(strcmp({S.Rc.name}, TARGET), 1);
if isempty(k), error('%s not found in qus_bfcomp.mat', TARGET); end
r = S.Rc(k);

% The B-modes were saved as PNG by qus_bfcomp, not as arrays, so they are
% read back from disk rather than recomputed.
bmPng = fullfile('figs_bfcomp', ['bmode_' TARGET '.png']);
haveBm = isfile(bmPng);

fig = figure('Color','w','Position',[80 60 1000 820]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

% ---- top row: B-modes ---------------------------------------------------
if haveBm
    % the saved PNG already holds both B-modes side by side, so it spans
    % the top row as one tile
    ax = nexttile([1 2]);
    imshow(imread(bmPng),'Parent',ax);
    title(ax,'Modo B: laboratorio (izq.) y propio (der.)');
else
    nexttile; axis off; text(.5,.5,'falta figs\_bfcomp','HorizontalAlignment','center');
    nexttile; axis off;
end

% ---- bottom row: ACS maps ----------------------------------------------
ax1 = nexttile;
imagesc(r.theirs.x_ACS, r.theirs.z_ACS, r.theirs.acs + ACS_REF, ACSLIM);
axis image; colormap(ax1, turbo);
xlabel('Lateral [cm]'); ylabel('Axial [cm]');
title(sprintf('ACS laboratorio:  %.3f \\pm %.3f', ...
      r.theirs.tg + ACS_REF, r.theirs.tgSd));
c = colorbar; c.Label.String = 'dB/cm/MHz';

ax2 = nexttile;
imagesc(r.mine.x_ACS, r.mine.z_ACS, r.mine.acs + ACS_REF, ACSLIM);
axis image; colormap(ax2, turbo);
xlabel('Lateral [cm]');
title(sprintf('ACS propio:  %.3f \\pm %.3f', ...
      r.mine.tg + ACS_REF, r.mine.tgSd));
c = colorbar; c.Label.String = 'dB/cm/MHz';

title(tl, sprintf('%s   |   diferencia %+.3f dB/cm/MHz', ...
      TARGET, r.mine.tg - r.theirs.tg), 'FontWeight','bold');

exportgraphics(fig, OUTFILE, 'Resolution', 220);
fprintf('Written: %s\n', OUTFILE);
fprintf('  laboratorio %.3f +- %.3f\n', r.theirs.tg + ACS_REF, r.theirs.tgSd);
fprintf('  propio      %.3f +- %.3f\n', r.mine.tg + ACS_REF, r.mine.tgSd);
