% centres_dt2809.m - measure inclusion centre (c1x, c1z) of each target on
% the lab B-mode of DataTargets2809 (Task 03, part A). Output in cm, lab
% coordinates, which is the grid phantoms_dt2809 works on.
%
% Method: circle fit to edge points on the upper half and sides of the
% inclusion (fit_inclusion_circle). The sign (hyper/hypo) comes from the
% lateral contrast of measure_inclusion, which is reliable for all targets.
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(fullfile(DATA_ROOT, 'src'));
warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
% QUS_BFSET: 'lab' (default) or 'own' (results/bf_own), see phantoms_dt2809
bfSet = getenv('QUS_BFSET'); if isempty(bfSet), bfSet = 'lab'; end
if strcmp(bfSet, 'own')
    bfDir  = fullfile(DATA_ROOT, 'results', 'bf_own', 'Samples');
    outDir = fullfile(DATA_ROOT, 'results', 'merino_dt2809_own');
else
    bfDir  = fullfile(DATA_ROOT, 'DataTargets2809', 'Samples', 'BeamformedData');
    outDir = fullfile(DATA_ROOT, 'results', 'merino_dt2809');
end
if ~isfolder(outDir), mkdir(outDir); end

names = {'T1_F','T2_F','T3_F_3','T4_F','T5_F','T6_F','T7_F','T8_F_3'};
C = table('Size',[numel(names) 6], ...
    'VariableTypes',{'string','double','double','double','double','double'}, ...
    'VariableNames',{'file','c1x','c1z','r','rmsRes_mm','contrastLat'});
f = figure('Units','centimeters','Position',[1 1 32 16]);
for k = 1:numel(names)
    S = load(fullfile(bfDir, [names{k} '.mat']));
    env = abs(hilbert(S.rf));
    if strcmp(bfSet, 'own')
        % Own set: start from the lab centre mapped to own axes (same
        % acquisition; Task 02: own z = lab z - 1.50 mm, own x = lab x *
        % 0.300/0.29184). Sign from disc-vs-side contrast, which does not
        % depend on profile crossings. The centre still comes from the fit.
        Lc = readtable(fullfile(DATA_ROOT, 'results', 'merino_dt2809', 'centres.csv'));
        c0 = [Lc.c1x(k)/100*0.300/0.29184, Lc.c1z(k)/100 - 1.50e-3];
        M0.contrastLat = discContrast(env, S.x, S.z, c0);
        M = fit_inclusion_circle(env, S.x, S.z, c0, sign(M0.contrastLat));
        fprintf('%-8s start (%.3f, %.3f) cm -> fit (%.3f, %.3f) cm\n', names{k}, ...
            c0*100, M.xc*100, M.zc*100);
        % Rejection rule fixed beforehand: a fit that moves the centre more
        % than 2 mm from the start has locked onto something else (weak
        % edges of an isoechoic inclusion); keep the mapped lab centre.
        if hypot(M.xc - c0(1), M.zc - c0(2)) > 2e-3
            fprintf('         fit rejected (moved > 2 mm): using start\n');
            M.xc = c0(1); M.zc = c0(2); M.r = NaN; M.rmsRes = NaN;
            M.ex = NaN; M.ez = NaN;
        end
    else
        % initial lateral measure; retry with another depth guess if it fails
        for zg = [20e-3 19e-3 21e-3]
            M0 = measure_inclusion(env, S.x, S.z, [0 zg]);
            if isfinite(M0.contrastLat), break; end
        end
        sgn = sign(M0.contrastLat);
        M = fit_inclusion_circle(env, S.x, S.z, [M0.xc 20e-3], sgn);
    end
    C(k,:) = {names{k}, M.xc*100, M.zc*100, M.r*100, M.rmsRes*1e3, M0.contrastLat};
    subplot(2,4,k);
    imagesc(S.x*100, S.z*100, 20*log10(env/max(env(:))), [-60 0]); colormap gray; axis image
    ylim([0 4.5]); hold on
    plot(M.ex*100, M.ez*100, 'y.', 'MarkerSize', 6);
    rp = M.r; if ~isfinite(rp), rp = 0.95e-2; end
    rectangle('Position',[M.xc-rp, M.zc-rp, 2*rp, 2*rp]*100, 'Curvature',1, ...
        'EdgeColor','r','LineStyle','--');
    plot(M.xc*100, M.zc*100, 'r+'); hold off
    title(sprintf('%s  r=%.2f cm', names{k}, M.r*100), 'Interpreter','none');
end
exportgraphics(f, fullfile(outDir, 'centres.png'), 'Resolution', 120);
disp(C);
writetable(C, fullfile(outDir, 'centres.csv'));

% =========================================================================
function c = discContrast(env, x, z, c0)
% Median dB inside r < 6 mm minus median of the sides (12-16 mm from the
% centre, within +-6 mm in depth, away from shadow or enhancement below).
[X, Z] = meshgrid(x, z);
lg = 20*log10(env/max(env(:)) + eps);
r  = hypot(X - c0(1), Z - c0(2));
inside = r < 6e-3;
sides  = r >= 12e-3 & r <= 16e-3 & abs(Z - c0(2)) < 6e-3;
c = median(lg(inside)) - median(lg(sides));
end
