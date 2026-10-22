% summarise_dt2809.m - Task 03, part A: 261-minus-544 difference per target
run(fullfile(fileparts(mfilename('fullpath')), '..', 'config.m'));
T = readtable(fullfile(DATA_ROOT, 'results', 'merino_dt2809', 'phantoms_dt2809.xlsx'));
T.ref = string(T.ref); T.method = string(T.method);
for meth = ["TV" "SWTV" "SWIFT"]
    a = T(T.method == meth & T.ref == "261", :);
    b = T(T.method == meth & T.ref == "544", :);
    dI = a.meanInc - b.meanInc;  dB = a.meanBack - b.meanBack;
    fprintf('\n%s\n%-8s %8s %8s %8s | %8s %8s %8s\n', meth, 'target', ...
        'inc261', 'inc544', 'diff', 'bck261', 'bck544', 'diff');
    for k = 1:height(a)
        fprintf('%-8s %8.3f %8.3f %+8.3f | %8.3f %8.3f %+8.3f\n', a.file{k}, ...
            a.meanInc(k), b.meanInc(k), dI(k), a.meanBack(k), b.meanBack(k), dB(k));
    end
    fprintf('%-8s %8s %8s %+8.3f | %8s %8s %+8.3f   (mean)\n', 'all', '', '', mean(dI), '', '', mean(dB));
    fprintf('%-8s %8s %8s %8.3f | %8s %8s %8.3f   (std)\n', '', '', '', std(dI), '', '', std(dB));
end
