% old_pwfoc.m - Task 04, phase 3.1 and 3.4 (old data): PW versus FOC on the
% CIRS 049 (T1 T3 T5 T7 T8) with the SWIFT parameters frozen in phase 2,
% Merino's estimator (src/merino_qus), Wisconsin reference of each scheme
% (ACS 0.5246, exponent 1.0019), reconstructions of scripts/bf_own_old.m
% (samples 1542, Wisconsin 1535.7).
%
% Band rule (fixed beforehand, SUPUESTOS.md): if the frozen band 2.5-7.5 MHz
% is not inside the -20 dB band of the mean reference spectrum of a scheme,
% that scheme uses its -20 dB band.
% Threshold sweep (3.4): ratioCutOff 10, 8, 6 with muB, muC frozen.
% Centres: circle fit on each sample B-mode; low-contrast rule as in
% centres_edu (|contrast| < 3 dB -> refit from median of contrasted targets,
% rejected if it moves > 2 mm -> median).
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
addpath(genpath(fullfile(MERINO_DIR, 'AttenuationUtils')));
addpath(fullfile(DATA_ROOT, 'src'));
warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
bfDir  = fullfile(DATA_ROOT, 'results', 'bf_old');
outDir = fullfile(DATA_ROOT, 'results', 'old_pwfoc');
if ~isfolder(outDir), mkdir(outDir); end

F = load(fullfile(DATA_ROOT, 'results', 'calib_edu', 'calib_edu.mat'), 'best');
best = F.best;
fprintf('Frozen: muB 10^%.1f muC 10^%.1f cut %d\n', best.lmuB, best.lmuC, best.cut);
CUTS = unique([best.cut 10 8 6], 'stable');

tg   = {'T1','T3','T5','T7','T8'};
cert = [0.52 0.74 0.75 0.95 0.95];
schemes = {'PW','FOC'};
region = [-1.85 1.85 0.2 3.5];

res = struct([]);
for s = 1:2
    sd = fullfile(bfDir, schemes{s});
    rl = dir(fullfile(sd, 'ref', '*.mat'));
    for j = 1:numel(rl)
        Rj = load(fullfile(rl(j).folder, rl(j).name), 'rf', 'x', 'z', 'fs');
        if j == 1, ref = zeros([size(Rj.rf) numel(rl)]); x = Rj.x; z = Rj.z; fs = Rj.fs; end
        ref(:,:,j) = Rj.rf;
    end

    % ---- band rule ----
    [bandUse, b20] = bandRule(ref, z, fs, [2.5e6 7.5e6]);
    fprintf('\n%s: -20 dB band %.2f-%.2f MHz -> using %.2f-%.2f MHz\n', schemes{s}, ...
        b20/1e6, bandUse/1e6);

    % ---- samples and centres ----
    sl = dir(fullfile(sd, 'samples', '*.mat'));
    names = cell(1, numel(tg)); RF = cell(1, numel(tg)); cen = nan(numel(tg), 2); con = nan(1, numel(tg));
    for k = 1:numel(tg)
        hit = sl(contains({sl.name}, ['_' tg{k} '_']));
        S = load(fullfile(hit(1).folder, hit(1).name), 'rf');
        names{k} = hit(1).name; RF{k} = S.rf;
        env = abs(hilbert(S.rf));
        for zg = [20e-3 19e-3 21e-3 18e-3]
            M0 = measure_inclusion(env, x, z, [0 zg]);
            if isfinite(M0.contrastLat), break; end
        end
        con(k) = M0.contrastLat;
        if isfinite(con(k)) && abs(con(k)) >= 3
            M = fit_inclusion_circle(env, x, z, [M0.xc 20e-3], sign(con(k)));
            cen(k,:) = [M.xc M.zc]*100;
        end
    end
    good = all(isfinite(cen), 2);
    cMed = median(cen(good,:), 1);
    method = repmat("fit", numel(tg), 1);
    for k = find(~good).'
        env = abs(hilbert(RF{k}));
        sg = sign(con(k)); if ~isfinite(sg) || sg == 0, sg = -1; end
        M = fit_inclusion_circle(env, x, z, cMed/100, sg);
        if hypot(M.xc*100 - cMed(1), M.zc*100 - cMed(2)) <= 0.2
            cen(k,:) = [M.xc M.zc]*100; method(k) = "refit from median";
        else
            cen(k,:) = cMed; method(k) = "median";
        end
    end
    for k = 1:numel(tg)
        fprintf('  %s %-45s centre (%.3f, %.3f) cm  contrast %.1f dB  [%s]\n', tg{k}, ...
            names{k}, cen(k,:), con(k), method(k));
    end

    % ---- estimation, frozen mu, each threshold ----
    for k = 1:numel(tg)
        for cu = CUTS
            o = struct('band', bandUse, 'region', region, 'refAcs', 0.5246, 'refGamma', 1.0019, ...
                'lmuB', best.lmuB, 'lmuC', best.lmuC, 'cut', cu, 'centre', cen(k,:));
            R = merino_qus(RF{k}, ref, x, z, fs, o);
            res(end+1).scheme = string(schemes{s}); %#ok<SAGROW>
            res(end).target = string(tg{k}); res(end).cert = cert(k); res(end).cut = cu;
            res(end).frozen = (cu == best.cut);
            res(end).inc = R.incSWIFT; res(end).sdInc = R.sdIncSWIFT;
            res(end).back = R.backSWIFT; res(end).sdBack = R.sdBackSWIFT;
            res(end).incTV = R.incTV; res(end).backTV = R.backTV;
            res(end).pctPenal = R.pctPenal;
            res(end).bandL = bandUse(1); res(end).bandH = bandUse(2);
            res(end).c1x = cen(k,1); res(end).c1z = cen(k,2); res(end).centreMethod = method(k);
            res(end).map = R; res(end).file = string(names{k});
        end
        fprintf('  %s done\n', tg{k});
    end
end
T = struct2table(rmfield(res, 'map'));
writetable(T, fullfile(outDir, 'old_pwfoc_all.csv'));
save(fullfile(outDir, 'old_pwfoc.mat'), 'res', 'best', 'CUTS', '-v7.3');

% ---- PW vs FOC table at each threshold ----
for cu = CUTS
    fprintf('\nratioCutOff %d%s\n', cu, repmat(' (frozen)', 1, cu == best.cut));
    fprintf('%-4s %5s | %13s %13s %7s | %13s %13s %7s\n', 'tgt', 'cert', ...
        'PW inc', 'FOC inc', 'diff', 'PW back', 'FOC back', 'diff');
    for k = 1:numel(tg)
        a = res([res.scheme] == "PW"  & [res.target] == tg{k} & [res.cut] == cu);
        b = res([res.scheme] == "FOC" & [res.target] == tg{k} & [res.cut] == cu);
        fprintf('%-4s %5.2f | %6.3f+-%.3f %6.3f+-%.3f %+7.3f | %6.3f+-%.3f %6.3f+-%.3f %+7.3f\n', ...
            tg{k}, cert(k), a.inc, a.sdInc, b.inc, b.sdInc, a.inc-b.inc, ...
            a.back, a.sdBack, b.back, b.sdBack, a.back-b.back);
    end
end

% =========================================================================
function [bandUse, b20] = bandRule(ref, z, fs, bandFrozen)
% Mean power spectrum of the references over 1-3 cm depth, all columns,
% smoothed; -20 dB band around its peak.
sel = z >= 0.01 & z <= 0.03;
seg = ref(sel, :, :);
seg = reshape(seg, size(seg,1), []);
N = 2^nextpow2(size(seg,1));
P = mean(abs(fft(seg .* hamming(size(seg,1)), N)).^2, 2);
P = P(1:N/2); f = (0:N/2-1)'/N*fs;
P = movmean(P, 9);
PdB = 10*log10(P/max(P));
[~, ip] = max(PdB);
lo = find(PdB(1:ip) < -20, 1, 'last'); hi = ip - 1 + find(PdB(ip:end) < -20, 1, 'first');
b20 = [f(lo+1) f(hi-1)];
if bandFrozen(1) >= b20(1) && bandFrozen(2) <= b20(2)
    bandUse = bandFrozen;
else
    bandUse = b20;
end
end
