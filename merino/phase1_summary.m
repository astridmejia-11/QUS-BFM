% phase1_summary.m - Task 04, phase 1: per-target table and depth test.
% Reads results/<set>/phantoms_dt2809.xlsx and maps.mat for each input set
% (lab: lab beamformer, everything at 1540; own: das_focused, each phantom
% at its own c) and reports
%   1. per target: background and inclusion with 261 and 544, and the
%      261 - 544 constant, for TV and SWIFT;
%   2. depth test: background ACS by z band with each reference. Background
%      = ACS-grid points at >= 1.05 cm from the inclusion centre (Merino's
%      background radius). If the constant changes with z the cause is a
%      beam effect, not a certificate offset.
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
sets  = {'lab', 'merino_dt2809'; 'own', 'merino_dt2809_own'};
bands = [0.5 1.5; 1.5 2.5; 2.5 3.5];
out = struct();
for is = 1:size(sets,1)
    rd = fullfile(DATA_ROOT, 'results', sets{is,2});
    if ~isfile(fullfile(rd, 'maps.mat')), fprintf('\n[%s] no maps.mat, skipped\n', sets{is,1}); continue; end
    T = readtable(fullfile(rd, 'phantoms_dt2809.xlsx'));
    T.ref = string(T.ref); T.method = string(T.method);
    M = load(fullfile(rd, 'maps.mat')); M = M.Maps;
    fprintf('\n================ input set: %s ================\n', sets{is,1});
    for meth = ["TV" "SWIFT"]
        a = T(T.method == meth & T.ref == "261", :);
        b = T(T.method == meth & T.ref == "544", :);
        fprintf('\n%s %-9s %7s %7s %7s | %7s %7s %7s\n', meth, 'target', ...
            'bck261', 'bck544', 'diff', 'inc261', 'inc544', 'diff');
        for k = 1:height(a)
            fprintf('   %-12s %7.3f %7.3f %+7.3f | %7.3f %7.3f %+7.3f\n', a.file{k}, ...
                a.meanBack(k), b.meanBack(k), a.meanBack(k)-b.meanBack(k), ...
                a.meanInc(k), b.meanInc(k), a.meanInc(k)-b.meanInc(k));
        end
        dB = a.meanBack - b.meanBack; dI = a.meanInc - b.meanInc;
        fprintf('   %-12s %7.3f %7.3f %+7.3f | %7.3f %7.3f %+7.3f  mean\n', 'all', ...
            mean(a.meanBack), mean(b.meanBack), mean(dB), mean(a.meanInc), mean(b.meanInc), mean(dI));
        fprintf('   %-12s %7s %7s %7.3f | %7s %7s %7.3f  std\n', '', '', '', std(dB), '', '', std(dI));
        out.(sets{is,1}).(meth).dBack = dB; out.(sets{is,1}).(meth).dInc = dI;
        out.(sets{is,1}).(meth).back261 = a.meanBack; out.(sets{is,1}).(meth).back544 = b.meanBack;
    end

    % ---- depth test ----
    files = unique(string({M.file}), 'stable');
    for meth = ["TV" "SWIFT"]
        V = nan(numel(files), size(bands,1), 2);   % target x band x ref
        for k = 1:numel(files)
            for r = 1:2
                tag = ["261" "544"];
                mm = M(string({M.file}) == files(k) & string({M.ref}) == tag(r));
                [Xa, Za] = meshgrid(mm.x_ACS, mm.z_ACS);
                back = hypot(Xa - mm.c1x, Za - mm.c1z) >= 1.05;
                A = mm.(char(meth));
                for ib = 1:size(bands,1)
                    sel = back & Za >= bands(ib,1) & Za < bands(ib,2);
                    V(k, ib, r) = mean(A(sel), 'omitnan');
                end
            end
        end
        fprintf('\nDepth test, %s, background ACS by z band (mean over 8 targets +- std)\n', meth);
        fprintf('   %-10s %15s %15s %15s\n', 'z [cm]', 'ref 261', 'ref 544', '261 - 544');
        for ib = 1:size(bands,1)
            d = V(:,ib,1) - V(:,ib,2);
            fprintf('   %4.1f-%-5.1f %8.3f+-%.3f %8.3f+-%.3f %+8.3f+-%.3f\n', bands(ib,:), ...
                mean(V(:,ib,1)), std(V(:,ib,1)), mean(V(:,ib,2)), std(V(:,ib,2)), mean(d), std(d));
        end
        out.(sets{is,1}).(meth).depth = V;
    end
end
save(fullfile(DATA_ROOT, 'results', 'phase1_summary.mat'), 'out', 'bands');
