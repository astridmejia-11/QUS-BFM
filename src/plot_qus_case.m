function fig = plot_qus_case(C, ttl, outFile)
%PLOT_QUS_CASE  B-mode with analysis and measurement regions, ACS map and
%   Delta-BSC map, all in physical coordinates [cm].
%
%   C fields: Bmode, x, z (pixel axes, cm, analysis region already cropped),
%             x_ACS, z_ACS, acs, bsc, centre [c1x c1z] cm (may be empty)
%   The analysis region is the extent of x, z. Measurement masks are
%   Merino's circles: inclusion r <= 0.85 cm, background r >= 1.05 cm.
fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 27 8.5], 'Color', 'w');
theme(fig, 'light');   % batch sessions default to the dark theme: grey text on white
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
th = linspace(0, 2*pi, 200);

ax = nexttile(tl);
imagesc(ax, C.x, C.z, C.Bmode, [-60 0]); axis(ax, 'image'); colormap(ax, gray);
hold(ax, 'on');
rectangle(ax, 'Position', [C.x(1) C.z(1) C.x(end)-C.x(1) C.z(end)-C.z(1)], ...
    'EdgeColor', [0.2 0.8 0.2], 'LineWidth', 1.2);
drawMasks(ax, C, th);
hold(ax, 'off');
title(ax, 'B-mode [dB]'); xlabel(ax, 'x [cm]'); ylabel(ax, 'z [cm]'); colorbar(ax);

ax = nexttile(tl);
imagesc(ax, C.x_ACS, C.z_ACS, C.acs, [0.2 1.4]); axis(ax, 'image'); colormap(ax, turbo);
xlim(ax, C.x([1 end])); ylim(ax, C.z([1 end]));
hold(ax, 'on'); drawMasks(ax, C, th); hold(ax, 'off');
title(ax, 'ACS [dB/cm/MHz]'); xlabel(ax, 'x [cm]'); colorbar(ax);

ax = nexttile(tl);
imagesc(ax, C.x_ACS, C.z_ACS, C.bsc, [-20 20]); axis(ax, 'image'); colormap(ax, parula);
xlim(ax, C.x([1 end])); ylim(ax, C.z([1 end]));
hold(ax, 'on'); drawMasks(ax, C, th); hold(ax, 'off');
title(ax, '\DeltaBSC [dB]'); xlabel(ax, 'x [cm]'); colorbar(ax);

title(tl, ttl, 'Interpreter', 'none', 'FontWeight', 'bold');
if nargin > 2 && ~isempty(outFile)
    od = fileparts(outFile); if ~isfolder(od), mkdir(od); end
    exportgraphics(fig, outFile, 'Resolution', 150);
    close(fig);
end
end

function drawMasks(ax, C, th)
if isfield(C, 'centre') && ~isempty(C.centre)
    plot(ax, C.centre(1) + 0.85*cos(th), C.centre(2) + 0.85*sin(th), 'w-', 'LineWidth', 1);
    plot(ax, C.centre(1) + 1.05*cos(th), C.centre(2) + 1.05*sin(th), 'w--', 'LineWidth', 1);
end
end
