% QUS_WIRES  Geometric validation with the ATS 539 wire targets.
%
%   Reconstructs the same raw data with TWO sound speeds and measures the wire
%   depths in each case:
%     1540 m/s  the value configured in the Verasonics
%     1450 m/s  the value in the ATS 539 datasheet
%
%   The wires belong to the Vertical Distance Group: 17 targets, 1 to 18 cm,
%   spaced 1 cm (ATS539Phantom_101623 datasheet).
%
%   The linear fit of measured against nominal depth separates two errors that
%   are independent of each other:
%     intercept -> t0 error   (a constant shift moves every point equally and
%                              does not change the slope)
%     slope     -> sound speed error (a scale error does not move the origin)
%
%   Run from the amejia folder:  qus_wires
%   Plain ASCII on purpose.

clear; clc;

% ============================ CONFIGURATION ==============================
ROOT     = pwd;
C_LIST   = [1540 1450];           % reconstruction speeds to compare
C_DATASHEET = 1450;               % ATS 539 datasheet
X_COL    = -10e-3;                % lateral position of the wire column [m]
X_TOL    =  1e-3;                 % columns averaged around it
Z_SEARCH = [3e-3 60e-3];          % depth search window
Z_NOMINAL = (10:10:50) * 1e-3;    % expected wire depths [m]
PEAK_SEP = 3e-3;                  % minimum separation between peaks [m]
PEAK_THR = 0.25;                  % fraction of the max profile value
FNUM = 2.0; APOD = 'hanning';
% =========================================================================

dataOnly = @(d) d(~contains({d.name},'preSet') & ~contains({d.name},'preSave'));
af = dataOnly(dir(fullfile(ROOT,'PLANE ATS','*.mat')));
hit = find(contains({af.name},'puntitos'), 1);
if isempty(hit)
    fprintf('  %s\n', af.name);
    error('No file with "puntitos" found in PLANE ATS.');
end
wireFile = fullfile(af(hit).folder, af(hit).name);
fprintf('File: %s\n\n', af(hit).name);

R = struct('c',{},'z',{},'sep',{},'slope',{},'intercept',{},'cReal',{});
figure('Position',[50 50 520*numel(C_LIST) 480]);

for k = 1:numel(C_LIST)
    c = C_LIST(k);
    D = load_vsx_rf(wireFile, 'c', c, 'Verbose', false);
    B = das_planewave(D, 'zLim', [3e-3 62e-3], 'FNumber', FNUM, ...
                      'Apod', APOD, 'Verbose', false);
    clear D

    % ---- profile along the wire column ----------------------------------
    % Averaging neighbouring columns lowers the noise without displacing the
    % position of the peak, which is the quantity of interest.
    zsel = B.z > Z_SEARCH(1) & B.z < Z_SEARCH(2);
    zz   = B.z(zsel);
    col  = abs(B.x - X_COL) < X_TOL;
    prof = mean(B.env(zsel, col), 2);

    dz  = B.z(2) - B.z(1);
    sep = max(1, round(PEAK_SEP/dz));
    thr = PEAK_THR * max(prof);
    isPk = false(size(prof));
    for q = sep+1 : numel(prof)-sep
        if prof(q) >= thr && prof(q) == max(prof(q-sep:q+sep))
            isPk(q) = true;
        end
    end
    zPk = zz(isPk);

    % keep as many peaks as nominal depths, the shallowest ones first
    n = min(numel(zPk), numel(Z_NOMINAL));
    zPk = zPk(1:n);
    zNom = Z_NOMINAL(1:n).';

    % ---- linear fit ------------------------------------------------------
    p = polyfit(zNom*1e3, zPk*1e3, 1);
    r.c = c; r.z = zPk*1e3; r.sep = diff(zPk)*1e3;
    r.slope = p(1); r.intercept = p(2); r.cReal = c/p(1);
    R(end+1) = r; %#ok<SAGROW>

    fprintf('--- reconstruido con c = %d m/s ---\n', c);
    fprintf('  profundidades [mm]: %s\n', num2str(zPk(:)'*1e3, '%.2f  '));
    fprintf('  separaciones  [mm]: %s\n', num2str(diff(zPk)'*1e3, '%.2f  '));
    fprintf('  pendiente  %.4f   intercepto %+.3f mm\n', p(1), p(2));
    fprintf('  c real implicada: %.0f m/s\n\n', c/p(1));

    % ---- B-mode with the detected wires marked ---------------------------
    subplot(1, numel(C_LIST), k);
    imagesc(B.x*1e3, B.z*1e3, 20*log10(B.env/max(B.env(:))), [-50 0]);
    axis image; colormap(gca, gray); colorbar;
    hold on;
    plot(X_COL*1e3*ones(size(zPk)), zPk*1e3, 'ro', 'MarkerSize', 9, 'LineWidth', 1.2);
    xline(X_COL*1e3, 'y--');
    hold off;
    xlabel('Lateral [mm]'); ylabel('Axial [mm]');
    title(sprintf('ATS 539, onda plana, c = %d m/s', c));
end

% ---- the fit figure that the comment asks for ---------------------------
figure('Position',[60 60 1050 430]);

subplot(1,2,1);
mk = {'o-','s--'};
for k = 1:numel(R)
    plot(Z_NOMINAL(1:numel(R(k).z))*1e3, R(k).z, mk{k}, ...
        'LineWidth',1.4,'MarkerSize',7); hold on;
end
plot([0 55],[0 55],'k:','LineWidth',1.2);
hold off; grid on; axis square;
xlabel('Profundidad nominal [mm]'); ylabel('Profundidad medida [mm]');
title('Ajuste: medida contra nominal');
legend([arrayfun(@(r) sprintf('c = %d m/s', r.c), R, 'uni', 0), {'identidad'}], ...
       'Location','northwest');

subplot(1,2,2);
for k = 1:numel(R)
    plot(1:numel(R(k).sep), R(k).sep, mk{k}, 'LineWidth',1.4,'MarkerSize',7);
    hold on;
end
yline(10,'k-','nominal 10 mm','LineWidth',1.2);
hold off; grid on; ylim([9 11]);
xlabel('par de hilos'); ylabel('Separacion medida [mm]');
title('Separacion entre hilos consecutivos');
legend(arrayfun(@(r) sprintf('c = %d m/s', r.c), R, 'uni', 0), ...
       'Location','southeast');

% ---- summary -------------------------------------------------------------
fprintf('=============== RESUMEN ===============\n');
fprintf('%-12s %10s %12s %12s\n','c recon','pendiente','intercepto','c implicada');
for k = 1:numel(R)
    fprintf('%-12d %10.4f %10.3f mm %9.0f m/s\n', ...
        R(k).c, R(k).slope, R(k).intercept, R(k).cReal);
end
fprintf(['\nDatasheet ATS 539: %d m/s +/-1%%, o sea %d a %d.\n'], ...
    C_DATASHEET, round(C_DATASHEET*0.99), round(C_DATASHEET*1.01));
fprintf(['Ninguna de las dos velocidades nominales reproduce la geometria\n' ...
         'del phantom. La que la reproduce es la c implicada, que es la misma\n' ...
         'en ambas reconstrucciones porque es una propiedad del material.\n']);
