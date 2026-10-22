% QUS_FIGS  Regenerates the three figures that were still missing, from the
%           .mat files already on disk. Nothing is recomputed.
%
%   1. The 5 x 2 grid of ACS maps          <- qus_final_results.mat
%   2. Bar chart, four materials           <- values typed in below
%   3. Beamformer sweep, both schemes      <- qus_bfsweep_PW.mat / _FOC.mat
%
%   For figure 3 you need the sweep saved under two names. If you overwrote
%   the file when you ran the focused case, re-run qus_bfsweep once per
%   scheme and save each result:
%       SCHEME='PW';  qus_bfsweep;  movefile('qus_bfsweep_results.mat','qus_bfsweep_PW.mat')
%       SCHEME='FOC'; qus_bfsweep;  movefile('qus_bfsweep_results.mat','qus_bfsweep_FOC.mat')
%
%   Run from the amejia folder:  qus_figs
%   Plain ASCII on purpose.

clear; clc;
ROOT = pwd;
ATT_RANGE = [0 1.5];

% =====================================================================
% 1. GRID OF ACS MAPS, 5 targets x 2 schemes
% =====================================================================
f1 = fullfile(ROOT,'qus_final_results.mat');
if isfile(f1)
    L = load(f1,'maps');
    maps = L.maps;
    schemes = unique(cellfun(@(m) string(m.scheme), maps), 'stable');
    targets = unique(cellfun(@(m) string(m.target), maps), 'stable');

    figure('Position',[30 30 240*numel(targets) 270*numel(schemes)]);
    tiledlayout(numel(schemes), numel(targets), ...
        'TileSpacing','compact','Padding','compact');
    for iS = 1:numel(schemes)
        for iT = 1:numel(targets)
            k = find(cellfun(@(m) string(m.scheme)==schemes(iS) & ...
                                  string(m.target)==targets(iT), maps), 1);
            nexttile
            if isempty(k)
                axis off; continue
            end
            O = maps{k}.O;
            imagesc(O.x_ACS, O.z_ACS, O.BSWIFT, ATT_RANGE);
            axis image; colormap turbo;
            title(sprintf('%s %s', schemes(iS), targets(iT)), 'FontSize',10);
            if iT == 1, ylabel('Axial [cm]'); end
            if iS == numel(schemes), xlabel('Lateral [cm]'); end
        end
    end
    cb = colorbar; cb.Label.String = 'ACS [dB/cm/MHz]';
    sgtitle('SWIFT, mapas de ACS', 'FontSize',12);

    % ---- pair T7 / T8, the controlled experiment ----------------------
    figure('Position',[40 40 900 560]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    for iS = 1:numel(schemes)
        for tg = ["T7","T8"]
            k = find(cellfun(@(m) string(m.scheme)==schemes(iS) & ...
                                  string(m.target)==tg, maps), 1);
            nexttile
            if isempty(k), axis off; continue; end
            imagesc(maps{k}.O.x_ACS, maps{k}.O.z_ACS, maps{k}.O.BSWIFT, ATT_RANGE);
            axis image; colormap turbo; colorbar;
            % Both have the same certified ACS: any difference between them
            % cannot come from attenuation.
            dbsc = 10.1*(tg=="T7") - 9.0*(tg=="T8");
            title(sprintf('%s %s   ACS cert 0.95, dBSC %+.1f dB', ...
                schemes(iS), tg, dbsc), 'FontSize',10);
            xlabel('Lateral [cm]'); ylabel('Axial [cm]');
        end
    end
    sgtitle('Par controlado: misma atenuacion, ecogenicidad opuesta','FontSize',12);
else
    fprintf('No se encontro %s, figura 1 y par T7/T8 omitidos.\n', f1);
end

% =====================================================================
% 2. BAR CHART, four materials
% =====================================================================
% Typed in rather than loaded: these come from three different scripts.
% Update if any value changes.
mat  = {'Wisconsin','ATS 539','Tiroides','CIRS fondo'};
meas = [0.530 0.539 1.194 0.456];
cert = [0.5246 0.50 1.34 0.55];
err  = [0     0.05  0.15  0    ];   % tolerancia declarada, cuando existe

figure('Position',[50 50 780 420]);
b = bar([meas; cert].', 'grouped');
b(1).FaceColor = [0.20 0.36 0.55];
b(2).FaceColor = [0.65 0.65 0.62];
hold on;
xt = (1:numel(mat)) + 0.145;        % centro de la segunda barra
errorbar(xt, cert, err, 'k.', 'LineWidth',1.1, 'CapSize',8);
hold off;
set(gca,'XTickLabel',mat);
ylabel('ACS [dB/cm/MHz]');
legend({'medido','referencia'},'Location','northwest');
grid on;
title('Medido contra referencia, cuatro materiales');
for i = 1:numel(mat)
    nb = 100*(meas(i)-cert(i))/cert(i);
    text(i, max(meas(i),cert(i))+0.06, sprintf('%+.1f%%',nb), ...
        'HorizontalAlignment','center','FontSize',9.5);
end
ylim([0 1.7]);

% =====================================================================
% 3. BEAMFORMER SWEEP, both schemes
% =====================================================================
fPW  = fullfile(ROOT,'qus_bfsweep_PW.mat');
fFOC = fullfile(ROOT,'qus_bfsweep_FOC.mat');
have = isfile(fPW) && isfile(fFOC);

if have
    A = load(fPW,'Rs');  B = load(fFOC,'Rs');
    sets = {A.Rs, B.Rs};
    names = {'Onda plana','Focalizado'};
    mk = {'o-','s--'};

    figure('Position',[40 40 1100 440]);

    subplot(1,2,1);
    for q = 1:2
        Rs = sets{q};
        apods = unique([Rs.apod],'stable');
        for a = 1:numel(apods)
            sel = [Rs.apod] == apods(a);
            sd = getSd(Rs(sel));
            errorbar([Rs(sel).fnum], [Rs(sel).acsBg], sd, mk{q}, ...
                'LineWidth',1.3,'MarkerSize',6, ...
                'DisplayName', sprintf('%s, %s', names{q}, apods(a)));
            hold on;
        end
    end
    yline(0.55,'k-','fondo certificado 0.55','LineWidth',1.2, ...
        'HandleVisibility','off');
    hold off; grid on; ylim([0 0.9]);
    xlabel('f-number de recepcion'); ylabel('ACS del fondo [dB/cm/MHz]');
    legend('Location','southeast','FontSize',8);
    title('La ACS no cambia');

    subplot(1,2,2);
    for q = 1:2
        Rs = sets{q};
        apods = unique([Rs.apod],'stable');
        for a = 1:numel(apods)
            sel = [Rs.apod] == apods(a);
            plot([Rs(sel).fnum], [Rs(sel).res_mm], mk{q}, ...
                'LineWidth',1.3,'MarkerSize',6, ...
                'DisplayName', sprintf('%s, %s', names{q}, apods(a)));
            hold on;
        end
    end
    hold off; grid on;
    xlabel('f-number de recepcion'); ylabel('Resolucion lateral [mm]');
    legend('Location','northwest','FontSize',8);
    title('La resolucion si cambia');
    sgtitle('El contraste entre estos dos paneles ES el resultado','FontSize',12);
else
    fprintf(['Faltan %s o %s.\n' ...
             'Corre qus_bfsweep una vez por esquema y renombra el .mat:\n' ...
             '  movefile(''qus_bfsweep_results.mat'',''qus_bfsweep_PW.mat'')\n'], ...
             fPW, fFOC);
end

fprintf('\nListo. Guarda cada figura con exportgraphics, por ejemplo:\n');
fprintf('  exportgraphics(figure(1), ''fig_grid.png'', ''Resolution'', 300)\n');


% =====================================================================
function sd = getSd(Rs)
%GETSD  Standard deviation field, whatever it ended up being called.
if isfield(Rs,'sdBg')
    sd = [Rs.sdBg];
elseif isfield(Rs,'iqrBg')
    sd = [Rs.iqrBg];      % corrida antigua: el campo conserva el nombre viejo
else
    sd = zeros(1,numel(Rs));
end
end
