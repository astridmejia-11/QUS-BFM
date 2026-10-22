% phase3_dt_summary.m - Task 04, phase 3.2-3.4 on DataTargets2809 with the
% frozen SWIFT parameters (muB 10^4, muC 10^2): own beamformer per target
% and reference, lab versus own, and the threshold sweep (10, 8, 6).
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
R0 = fullfile(DATA_ROOT, 'results');
gt = [0.52 0.55 0.74 0.81 0.75 0.97 0.95 0.95];
get = @(set, cut) loadSwift(fullfile(R0, ['merino_dt2809' set], sprintf('frozen_cut%d', cut), ...
    'phantoms_dt2809.xlsx'));

% ---- 3.2: own beamformer, frozen threshold 8 ------------------------------
O = get('_own', 8);
fprintf('\n3.2 Own beamformer, SWIFT frozen (cut 8)\n');
fprintf('%-12s %5s | %7s %7s %7s | %7s %7s %7s\n', 'target', 'cert', 'inc261', 'inc544', ...
    'diff', 'bck261', 'bck544', 'diff');
for k = 1:8
    a = O(O.ref == "261", :); b = O(O.ref == "544", :);
    fprintf('%-12s %5.2f | %7.3f %7.3f %+7.3f | %7.3f %7.3f %+7.3f\n', a.file{k}, gt(k), ...
        a.meanInc(k), b.meanInc(k), a.meanInc(k)-b.meanInc(k), ...
        a.meanBack(k), b.meanBack(k), a.meanBack(k)-b.meanBack(k));
end
a = O(O.ref == "261", :); b = O(O.ref == "544", :);
fprintf('%-12s %5s | %7.3f %7.3f %+7.3f | %7.3f %7.3f %+7.3f  mean\n', 'all', '', ...
    mean(a.meanInc), mean(b.meanInc), mean(a.meanInc-b.meanInc), ...
    mean(a.meanBack), mean(b.meanBack), mean(a.meanBack-b.meanBack));
for rr = ["261" "544"]
    s = O(O.ref == rr, :);
    fprintf('ref %s: mean |rel bias| inc %.3f, back %.3f (back mean %.3f vs 0.55)\n', rr, ...
        mean(abs(s.meanInc - gt.')./gt.'), mean(abs(s.meanBack - 0.55)/0.55), mean(s.meanBack));
end

% ---- 3.3: lab versus own, same estimator, same frozen parameters ----------
L = get('', 8);
fprintf('\n3.3 Lab minus own (SWIFT frozen, cut 8)\n');
fprintf('%-12s | %8s %8s | %8s %8s\n', 'target', 'inc 261', 'inc 544', 'back 261', 'back 544');
for k = 1:8
    v = zeros(1,4); i = 0;
    for q = ["meanInc" "meanBack"]
        for rr = ["261" "544"]
            i = i + 1;
            l = L(L.ref == rr, :); o = O(O.ref == rr, :);
            v(i) = l.(q)(k) - o.(q)(k);
        end
    end
    nmL = L.file(L.ref == "261");
    fprintf('%-12s | %+8.3f %+8.3f | %+8.3f %+8.3f\n', nmL{k}, v);
end
for q = ["meanInc" "meanBack"]
    for rr = ["261" "544"]
        d = L.(q)(L.ref == rr) - O.(q)(O.ref == rr);
        sdI = mean(O.(strrep(q, 'mean', 'std'))(O.ref == rr));
        fprintf('  %s ref %s: lab - own mean %+.3f (std %.3f) | typical individual std %.3f\n', ...
            q, rr, mean(d), std(d), sdI);
    end
end

% ---- 3.4: threshold sweep -------------------------------------------------
fprintf('\n3.4 Threshold sweep, mean over 8 targets (SWIFT)\n');
fprintf('%-5s %-4s %5s | %8s %8s | %8s\n', 'set', 'ref', 'cut', 'inc', 'back', 'mean|relbias|');
for set = ["_own" ""]
    for rr = ["261" "544"]
        for cut = [10 8 6]
            T = get(char(set), cut);
            if isempty(T), continue; end
            s = T(T.ref == rr, :);
            e = mean([abs(s.meanInc - gt.')./gt.'; abs(s.meanBack - 0.55)/0.55]);
            nm = "own"; if set == "", nm = "lab"; end
            fprintf('%-5s %-4s %5d | %8.3f %8.3f | %8.3f\n', nm, rr, cut, mean(s.meanInc), mean(s.meanBack), e);
        end
    end
end

function T = loadSwift(f)
if ~isfile(f), T = []; return; end
T = readtable(f);
T.ref = string(T.ref); T.method = string(T.method);
T = T(T.method == "SWIFT", :);
end
