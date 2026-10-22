% QUS_FIGS_ALL  One figure per case, straight out of qus_all_results.mat.
%
%   Three panels per case:
%       (1) B-mode, with the analysis region dashed and the measurement
%           region solid. This is what the advisors asked for: every number
%           in the table must be traceable to where it was measured.
%       (2) ACS map, same two regions drawn on it.
%       (3) BSC map, where the arcs at the target boundary show up.
%
%   Files are written to ./figs as fig_<case>.png, one per case.
%
%   Run from the amejia folder, after qus_all.
%
%   Plain ASCII on purpose.

clear; clc;

IN     = 'qus_all_results.mat';
OUTDIR = 'figs';
DPI    = 200;

if ~isfile(IN), error('Run qus_all first: %s not found.', IN); end
if ~isfolder(OUTDIR), mkdir(OUTDIR); end

S = load(IN, 'R');
R = S.R;

for i = 1:numel(R)
    if isempty(R(i).name), continue; end
    r  = R(i);
    st = r.stats;

    % ---- measurement region, as a box in physical coordinates -----------
    % Derived from the stored mask so it always matches what was measured,
    % whatever the region rule was ('all', 'bg', 'tg', 'box').
    [zi, xi] = find(r.selMask);
    mx = [min(r.x_ACS(xi)) max(r.x_ACS(xi))];
    mz = [min(r.z_ACS(zi)) max(r.z_ACS(zi))];

    % ---- colour scales --------------------------------------------------
    % Thyroid sits near 1.2 and the phantoms near 0.5. A single scale would
    % flatten the phantoms into one colour, so the range follows the case.
    if r.cert > 1, acsLim = [0 3.0]; else, acsLim = [0 1.2]; end
    q = prctile(r.bscMap(isfinite(r.bscMap)), [2 98]);
    bscLim = [-1 1] * max(2, max(abs(q)));

    % ---- figure ---------------------------------------------------------
    f = figure('Visible','off','Position',[100 100 1500 460],'Color','w');
    t = tiledlayout(f, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    % (1) B-mode
    ax = nexttile;
    imagesc(r.xBm, r.zBm, r.bMode, [-60 0]);
    colormap(ax, gray); axis image; hold on;
    drawRegions(r.ana, mx, mz);
    xlabel('Lateral [cm]'); ylabel('Axial [cm]');
    title('B-mode'); colorbar; hold off;

    % (2) ACS
    ax = nexttile;
    imagesc(r.x_ACS, r.z_ACS, r.acsMap, acsLim);
    colormap(ax, turbo); axis image; hold on;
    drawRegions(r.ana, mx, mz);
    xlabel('Lateral [cm]'); ylabel('Axial [cm]');
    title(sprintf('ACS = %.3f \\pm %.3f dB/cm/MHz', st.median, st.std));
    c = colorbar; c.Label.String = 'dB/cm/MHz'; hold off;

    % (3) BSC
    ax = nexttile;
    imagesc(r.x_ACS, r.z_ACS, r.bscMap, bscLim);
    colormap(ax, turbo); axis image; hold on;
    drawRegions(r.ana, mx, mz);
    xlabel('Lateral [cm]'); ylabel('Axial [cm]');
    title(sprintf('\\DeltaBSC   (detector fired on %.1f %% of blocks)', r.pctPenal));
    c = colorbar; c.Label.String = 'dB'; hold off;

    title(t, sprintf(['%s  |  %s  |  reference %.2f (%s)  |  bias %+.1f %%  |  ' ...
                      'p10-p90 [%.3f  %.3f]'], ...
            r.name, r.scheme, r.cert, r.certSrc, st.bias, st.p10, st.p90), ...
          'FontWeight','bold','Interpreter','none');

    % ---- save -----------------------------------------------------------
    safe = regexprep(r.name, '[^\w\-]', '_');
    out  = fullfile(OUTDIR, ['fig_' safe '.png']);
    exportgraphics(f, out, 'Resolution', DPI);
    close(f);
    fprintf('%-18s -> %s\n', r.name, out);
end

fprintf('\n%d figures written to ./%s\n', numel(R), OUTDIR);


% ---------------------------------------------------------------- helper
function drawRegions(ana, mx, mz)
% Dashed white: analysis region, where the estimator actually solves.
% Solid green: measurement region, where the reported number comes from.
% Both are drawn on every panel so the two are never confused.
    rectangle('Position', [ana(1) ana(3) ana(2)-ana(1) ana(4)-ana(3)], ...
              'EdgeColor','w', 'LineStyle','--', 'LineWidth', 1.2);
    rectangle('Position', [mx(1) mz(1) mx(2)-mx(1) mz(2)-mz(1)], ...
              'EdgeColor',[0 1 0], 'LineWidth', 1.8);
end
