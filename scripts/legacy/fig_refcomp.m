% FIG_REFCOMP  The same eight targets, measured against two reference phantoms.
%
%   Each target is one piece of material with one true attenuation, so the
%   two curves should lie on top of each other. They do not: they run
%   parallel, separated by a constant offset. The offset is the finding.
%
%   Reads qus_new_results.mat if it is there; otherwise the values below are
%   used, which are the ones already reported.
%
%   Run from the amejia folder:  fig_refcomp
%
%   Plain ASCII on purpose.

clear; clc;

CUT      = 10;        % threshold whose results are plotted
ACS_261  = 0.54;      % certificate, part 20960000, fibroglandular breast
ACS_544  = 0.53;      % certificate, part 20960001, liver
OUTFILE  = 'fig_refcomp.png';

tg = {'T1','T2','T3','T4','T5','T6','T7','T8'};

% ---- fallback values, already reported ----------------------------------
a261 = [0.859 0.651 0.782 1.251 0.827 1.000 1.133 1.021];
a544 = [0.720 0.472 0.644 1.083 0.671 0.874 0.965 0.827];
s261 = [0.049 0.081 0.082 0.091 0.104 0.074 0.075 0.068];
s544 = [0.060 0.060 0.083 0.079 0.105 0.060 0.074 0.083];
bg261 = 0.615;  bg544 = 0.506;

% ---- prefer the saved run, so the figure never drifts from the data -----
if isfile('qus_new_results.mat')
    S = load('qus_new_results.mat','R');
    R = S.R;
    for i = 1:numel(tg)
        k1 = find(strcmp({R.sample},tg{i}) & strcmp({R.ref},'261') & [R.cut]==CUT, 1);
        k2 = find(strcmp({R.sample},tg{i}) & strcmp({R.ref},'544') & [R.cut]==CUT, 1);
        if ~isempty(k1)
            off = 0; if strcmp(R(k1).kind,'dACS'), off = ACS_261; end
            a261(i) = R(k1).tg + off;  s261(i) = R(k1).tgSd;
        end
        if ~isempty(k2)
            off = 0; if strcmp(R(k2).kind,'dACS'), off = ACS_544; end
            a544(i) = R(k2).tg + off;  s544(i) = R(k2).tgSd;
        end
    end
    fprintf('Values taken from qus_new_results.mat (cut %g)\n', CUT);
else
    fprintf('qus_new_results.mat not found; using the reported values\n');
end

d = a261 - a544;
fprintf('Constant: mean %.3f, sd %.3f\n', mean(d), std(d));

% ---- figure -------------------------------------------------------------
x = 1:numel(tg);
C1 = [0.545 0.180 0.180];    % 261
C2 = [0.180 0.369 0.545];    % 544

fig = figure('Color','w','Position',[100 100 900 620]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

% --- top: the two series ---
ax1 = nexttile;
% error bars first so the markers sit on top
errorbar(x, a261, s261, 'o-', 'Color', C1, 'MarkerFaceColor', C1, ...
         'LineWidth', 1.6, 'MarkerSize', 7, 'CapSize', 8); hold on;
errorbar(x, a544, s544, 's-', 'Color', C2, 'MarkerFaceColor', C2, ...
         'LineWidth', 1.6, 'MarkerSize', 7, 'CapSize', 8);

% the background: one physical region, two answers
yline(bg261, ':', 'fondo vs 261', 'Color', C1, 'LineWidth', 1.2, ...
      'LabelHorizontalAlignment','left','FontSize',9);
yline(bg544, ':', 'fondo vs 544', 'Color', C2, 'LineWidth', 1.2, ...
      'LabelHorizontalAlignment','left','FontSize',9);

set(ax1,'XTick',x,'XTickLabel',tg,'XLim',[0.5 numel(tg)+0.5]);
grid on; box off;
ylabel('ACS [dB/cm/MHz]');
legend({'referencia 261  (cert. 0.54)','referencia 544  (cert. 0.53)'}, ...
       'Location','northwest','Box','off');
title('Los mismos targets, dos referencias');

% --- bottom: the offset ---
ax2 = nexttile;
bar(x, d, 0.55, 'FaceColor', [0.72 0.72 0.70], 'EdgeColor','none'); hold on;
yline(mean(d), '-', sprintf('media %.3f', mean(d)), ...
      'Color', [0.2 0.2 0.2], 'LineWidth', 1.6, 'FontSize', 10);
yline(ACS_261-ACS_544, '--', ...
      sprintf('diferencia certificada %.2f', ACS_261-ACS_544), ...
      'Color', C1, 'LineWidth', 1.4, 'FontSize', 10, ...
      'LabelHorizontalAlignment','right');

set(ax2,'XTick',x,'XTickLabel',tg,'XLim',[0.5 numel(tg)+0.5]);
grid on; box off;
ylabel('Diferencia [dB/cm/MHz]'); xlabel('Target');
title('Diferencia entre las dos estimaciones del mismo target');

exportgraphics(fig, OUTFILE, 'Resolution', 220);
fprintf('Written: %s\n', OUTFILE);
