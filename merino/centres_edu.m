% centres_edu.m - Task 04, phase 2: inclusion centre (c1x, c1z) of the eight
% targets_edu targets, measured on the B-mode of frame 1 (the frame
% phantoms.m uses) by circle fit, as in centres_dt2809.m.
%
% Rule for low-contrast (isoechoic) targets, fixed before the fit was used:
% if |lateral contrast| < 3 dB, refit starting from the median centre of
% the targets with contrast >= 3 dB; if that fit moves more than 2 mm from
% its start it is rejected and the start is kept (for T6, Merino's published
% centre (1.8, 1.9) cm is used instead, since it exists).
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(fullfile(DATA_ROOT, 'src'));
warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
dataDir = fullfile(DATA_ROOT, 'targets_edu', 'data');
outDir  = fullfile(DATA_ROOT, 'results', 'calib_edu');
if ~isfolder(outDir), mkdir(outDir); end

names = arrayfun(@(k) sprintf('T%d', k), 1:8, 'UniformOutput', false);
C = table('Size',[8 6], ...
    'VariableTypes',{'string','double','double','double','double','double'}, ...
    'VariableNames',{'file','c1x','c1z','r','rmsRes_mm','contrastLat'});
f = figure('Units','centimeters','Position',[1 1 32 16]);
for k = 1:8
    S = load(fullfile(dataDir, [names{k} '.mat']), 'RF', 'x', 'z');
    rf = S.RF(:,:,1); x = S.x; z = S.z; clear S
    ENV{k} = abs(hilbert(rf)); X{k} = x; Z{k} = z; %#ok<SAGROW>
    env = abs(hilbert(rf));
    x0 = mean(x([1 end]));
    for zg = [19e-3 20e-3 18e-3]
        M0 = measure_inclusion(env, x, z, [x0 zg]);
        if isfinite(M0.contrastLat), break; end
    end
    M = fit_inclusion_circle(env, x, z, [M0.xc 19e-3], sign(M0.contrastLat));
    C(k,:) = {names{k}, M.xc*100, M.zc*100, M.r*100, M.rmsRes*1e3, M0.contrastLat};
    subplot(2,4,k);
    imagesc(x*100, z*100, 20*log10(env/max(env(:))), [-60 0]); colormap gray; axis image
    hold on
    plot(M.ex*100, M.ez*100, 'y.', 'MarkerSize', 6);
    rectangle('Position',[M.xc-M.r, M.zc-M.r, 2*M.r, 2*M.r]*100, 'Curvature',1, ...
        'EdgeColor','r','LineStyle','--');
    plot(M.xc*100, M.zc*100, 'r+'); hold off
    xlabel('x [cm]'); ylabel('z [cm]');
    title(sprintf('%s  r=%.2f cm', names{k}, M.r*100));
end
% ---- low-contrast rule ----
C.method = repmat("fit", 8, 1);
good = abs(C.contrastLat) >= 3;
cMed = [median(C.c1x(good)) median(C.c1z(good))]/100;
for k = find(~good).'
    M = fit_inclusion_circle(ENV{k}, X{k}, Z{k}, cMed, sign(C.contrastLat(k)));
    if hypot(M.xc - cMed(1), M.zc - cMed(2)) <= 2e-3
        C(k, {'c1x','c1z','r','rmsRes_mm'}) = {M.xc*100, M.zc*100, M.r*100, M.rmsRes*1e3};
        C.method(k) = "refit from median";
    elseif strcmp(names{k}, 'T6')
        C(k, {'c1x','c1z','r','rmsRes_mm'}) = {1.8, 1.9, NaN, NaN};
        C.method(k) = "Merino published";
    else
        C(k, {'c1x','c1z','r','rmsRes_mm'}) = {cMed(1)*100, cMed(2)*100, NaN, NaN};
        C.method(k) = "median of contrasted targets";
    end
    subplot(2,4,k); hold on
    rectangle('Position',[C.c1x(k)-0.95, C.c1z(k)-0.95, 1.9, 1.9], 'Curvature',1, ...
        'EdgeColor','c','LineStyle','-');
    plot(C.c1x(k), C.c1z(k), 'c+'); hold off
    title(sprintf('%s  (%s)', names{k}, C.method(k)));
end
exportgraphics(f, fullfile(outDir, 'centres.png'), 'Resolution', 120);
disp(C);
writetable(C, fullfile(outDir, 'centres.csv'));
