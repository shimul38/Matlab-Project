%% =========================================================================
%  Bio-AI Clinical ECG — Batch Pan-Tompkins QRS Detector
%  Author    : Snehashish (BUET, Applied Bio-AI)
%  Supervisor: Dr. Taufiq Hasan
%  Version   : 2.0
%
%  Key improvements over v1:
%    - Zero-phase pipeline: butter + filtfilt + movmean (no external helpers)
%    - R-peak back-tracing: MWI detections snapped to true filtered ECG peak
%    - Pre-built filter objects reused across all records (speed)
%    - Double-peak / refractory period guard in detection loop
%    - Vectorised TP counting replaces inner k-loop (speed)
%    - Per-record results table printed to console + exportable struct
%    - Sanity-check plot uses back-traced peaks (no manual shift constant)
%    - Fully self-contained: no pan_tompkins_bandpass / derivative_func /
%      movingWindowIntegrator dependencies
%% =========================================================================
clc; clear; close all;

%% =========================================================================
%  SECTION 1: Configuration
%% =========================================================================

% --- Toolbox Path ---
if ~exist('rdsamp', 'file')
    addpath('B:\EEE 376\Labtest materials\lab-05\mcode'); % Adjust if needed
end
wfdbloadlib;

% --- Dataset Path (relative or absolute) ---
DB_PATH = 'mit-bih-arrhythmia-database-1.0.0';

% --- Record List (MIT-BIH standard 48-record set) ---
File_No = [100:109, 111:119, 200:203, 205, 207:210, 212:215, ...
           217, 219:223, 228, 230:234];

% --- Sampling Rates ---
fs_orig = 360;
fs_new  = 200;   % Pan-Tompkins designed for ~200 Hz

% --- Detection Parameters ---
REFRACT_MS      = 200;   % Refractory period (ms)
SEARCH_WIN_MS   = 100;   % R-peak back-trace search window (ms)
MWI_WIN_MS      = 150;   % Moving window integration length (ms)
TOLERANCE_MS    = 150;   % AHA evaluation tolerance (ms)
SEARCHBACK_MULT = 1.66;  % RR multiplier to trigger search-back
DECAY_GAP_S     = 1.5;   % Seconds of silence before threshold decay

% Convert to samples
refract_samp    = round(REFRACT_MS    / 1000 * fs_new);
search_samp     = round(SEARCH_WIN_MS / 1000 * fs_new);
mwi_win_samp    = round(MWI_WIN_MS    / 1000 * fs_new);
tolerance_samp  = round(TOLERANCE_MS  / 1000 * fs_new);
decay_samp      = round(DECAY_GAP_S   * fs_new);
init_rr_samp    = round(0.8 * fs_new); % Initial RR guess (~75 BPM)

% --- Valid Beat Annotation Symbols ---
VALID_SYMS = {'N','L','R','A','a','J','S','V','F','e','j','E','/'};

%% =========================================================================
%  SECTION 2: Pre-Build Filter Objects (Reused Every Record — Much Faster)
%% =========================================================================

% Butterworth bandpass 5-15 Hz (Pan-Tompkins spec)
[b_bp, a_bp] = butter(2, [5 15] / (fs_new/2), 'bandpass');

% Pan-Tompkins 5-point derivative: h = (fs/8)*[-2,-1,0,1,2]
b_der = (fs_new / 8) * [-2, -1, 0, 1, 2];

fprintf('Filter objects built. Starting batch processing...\n\n');

%% =========================================================================
%  SECTION 3: Batch Processing Loop
%% =========================================================================

num_files   = length(File_No);
all_bpm     = nan(num_files, 1);
all_qrs     = cell(num_files, 1);   % Back-traced R-peak indices
all_mwi     = cell(num_files, 1);   % MWI signals (for sanity-check plot)
all_ecg_rs  = cell(num_files, 1);   % Resampled raw ECG (for plot)
all_ecg_flt = cell(num_files, 1);   % Filtered ECG (for back-trace plot)

% Results struct
null_stat = struct('TP',0,'FP',0,'FN',0,'Se',NaN,'Pp',NaN,'loaded',false,'validated',false);
results   = repmat(null_stat, num_files, 1);

fprintf('%-6s  %-8s  %-8s  %-8s  %-8s  %-8s\n', ...
        'Rec', 'Beats', 'BPM', 'Se(%)', '+P(%)', 'Status');
fprintf('%s\n', repmat('-', 1, 55));

for i = 1:num_files
    No          = File_No(i);
    rec_str     = num2str(No);
    record_path = fullfile(DB_PATH, rec_str);

    % -----------------------------------------------------------------------
    %  Step A: Load & Resample
    % -----------------------------------------------------------------------
    try
        [signal, ~] = rdsamp(record_path, 1);
    catch
        fprintf('%-6s  [LOAD ERROR — skipped]\n', rec_str);
        continue;
    end

    ecg_rs  = resample(signal, fs_new, fs_orig);
    N       = length(ecg_rs);
    results(i).loaded = true;

    % -----------------------------------------------------------------------
    %  Step B: Pan-Tompkins Pipeline (all zero-phase)
    % -----------------------------------------------------------------------
    ecg_flt = filtfilt(b_bp, a_bp, ecg_rs);           % Bandpass
    ecg_der = filtfilt(b_der, 1,    ecg_flt);          % Derivative
    ecg_sq  = ecg_der .^ 2;                            % Squaring
    mwi     = movmean(ecg_sq, mwi_win_samp);           % Integration
    mwi     = mwi / max(mwi);                          % Normalize

    % -----------------------------------------------------------------------
    %  Step C: Initial Peak Candidates
    % -----------------------------------------------------------------------
    min_h = median(mwi) * 0.8;
    [pks, locs] = findpeaks(mwi, ...
        'MinPeakDistance', refract_samp, ...
        'MinPeakHeight',   min_h);

    % -----------------------------------------------------------------------
    %  Step D: Adaptive Training (first 2 s)
    % -----------------------------------------------------------------------
    train_end = 2 * fs_new;
    tr_idx    = locs <= train_end;

    if any(tr_idx)
        spki = mean(pks(tr_idx));
        npki = mean(mwi(1:min(train_end, N))) * 0.5;
    else
        spki = max(pks) * 0.4;
        npki = mean(mwi)  * 0.1;
    end

    T1     = npki + 0.25 * (spki - npki);
    T2     = 0.5 * T1;
    avg_rr = init_rr_samp;

    % -----------------------------------------------------------------------
    %  Step E: Adaptive Detection Loop
    % -----------------------------------------------------------------------
    qrs_mwi = [];   % MWI-domain beat locations

    for j = 1:length(pks)
        cp_loc = locs(j);
        cp_val = pks(j);

        % -- Refractory period guard + double-peak replacement --
        if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) < refract_samp
            if cp_val > mwi(qrs_mwi(end))
                qrs_mwi(end) = cp_loc;
            end
            continue;
        end

        % -- Threshold decay after long silence --
        if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) > decay_samp
            T1 = T1 * 0.5;
            T2 = 0.5 * T1;
        end

        % -- Search-back for missed beat --
        if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) > (SEARCHBACK_MULT * avg_rr)
            gap_idx = find(locs > qrs_mwi(end) & locs < cp_loc);
            for m = gap_idx'
                if pks(m) > 0.75 * T2
                    qrs_mwi = [qrs_mwi; locs(m)]; %#ok<AGROW>
                    spki    = 0.25 * pks(m) + 0.75 * spki;
                    T1      = npki + 0.25 * (spki - npki);
                    T2      = 0.5 * T1;
                    break;
                end
            end
        end

        % -- Standard detection --
        if cp_val > T1
            qrs_mwi    = [qrs_mwi; cp_loc]; %#ok<AGROW>
            capped     = min(cp_val, 2 * spki);
            spki       = 0.125 * capped + 0.875 * spki;
        else
            npki       = 0.125 * cp_val + 0.875 * npki;
        end

        T1 = npki + 0.25 * (spki - npki);
        T2 = 0.5 * T1;

        if length(qrs_mwi) > 2
            avg_rr = mean(diff(qrs_mwi(max(1, end-4):end)));
        end
    end

    % -----------------------------------------------------------------------
    %  Step F: R-Peak Back-Tracing (snap MWI peaks → filtered ECG maxima)
    % -----------------------------------------------------------------------
    qrs_refined = zeros(size(qrs_mwi));
    for k = 1:length(qrs_mwi)
        lo = max(1, qrs_mwi(k) - search_samp);
        hi = min(N, qrs_mwi(k) + search_samp);
        [~, peak_offset] = max(abs(ecg_flt(lo:hi)));
        qrs_refined(k)   = lo + peak_offset - 1;
    end
    qrs_refined = unique(qrs_refined);

    % -----------------------------------------------------------------------
    %  Step G: BPM & RR Cleaning
    % -----------------------------------------------------------------------
    rr_ms   = diff(qrs_refined) / fs_new * 1000;
    rr_clean = rr_ms(rr_ms > 300 & rr_ms < 2000);

    if ~isempty(rr_clean)
        all_bpm(i) = 60000 / mean(rr_clean);
    end

    all_qrs{i}     = qrs_refined;
    all_mwi{i}     = mwi;
    all_ecg_rs{i}  = ecg_rs;
    all_ecg_flt{i} = ecg_flt;

    % -----------------------------------------------------------------------
    %  Step H: Validation Against .atr Annotations
    % -----------------------------------------------------------------------
    se_str = '  N/A  ';
    pp_str = '  N/A  ';

    try
        [ann_locs, ann_type] = rdann(record_path, 'atr');
        valid_mask    = ismember(ann_type, VALID_SYMS);
        truth_samples = round(ann_locs(valid_mask) * (fs_new / fs_orig));

        % Vectorised TP count (no inner loop)
        tp = sum(any(abs(qrs_refined - truth_samples') <= tolerance_samp, 1));
        fp = length(qrs_refined) - tp;
        fn = length(truth_samples) - tp;

        results(i).TP        = tp;
        results(i).FP        = max(fp, 0);
        results(i).FN        = max(fn, 0);
        results(i).Se        = 100 * tp / (tp + fn);
        results(i).Pp        = 100 * tp / (tp + max(fp,0));
        results(i).validated = true;

        se_str = sprintf('%6.2f', results(i).Se);
        pp_str = sprintf('%6.2f', results(i).Pp);
    catch
        % No .atr — analysis-only mode, no error
    end

    status = 'OK';
    fprintf('%-6s  %-8d  %-8.1f  %-8s  %-8s  %s\n', ...
            rec_str, length(qrs_refined), all_bpm(i), se_str, pp_str, status);
end

%% =========================================================================
%  SECTION 4: Aggregate Statistics
%% =========================================================================
validated = [results.validated];
se_vals   = [results(validated).Se];
pp_vals   = [results(validated).Pp];

total_tp  = sum([results.TP]);
total_fp  = sum([results.FP]);
total_fn  = sum([results.FN]);
avg_Acc   = 100 * total_tp / (total_tp + total_fp + total_fn);

fprintf('\n%s\n', repmat('=', 1, 55));
fprintf('  AGGREGATE RESULTS (%d validated / %d total records)\n', ...
        sum(validated), num_files);
fprintf('%s\n', repmat('=', 1, 55));
fprintf('  Mean Sensitivity (Se)      : %6.2f%%\n', mean(se_vals));
fprintf('  Mean Pos. Predictivity (+P): %6.2f%%\n', mean(pp_vals));
fprintf('  Overall Accuracy           : %6.2f%%\n', avg_Acc);
fprintf('  Total TP / FP / FN         : %d / %d / %d\n', ...
        total_tp, total_fp, total_fn);
fprintf('%s\n\n', repmat('=', 1, 55));

%% =========================================================================
%  SECTION 5: Sanity-Check Plot — Interactive Record Selector
%% =========================================================================
% Change this index to inspect any record. 
% Example: idx=1 → record 100, idx=5 → record 104, idx=23 → record 207, etc.
PLOT_IDX = 1;

if PLOT_IDX < 1 || PLOT_IDX > num_files || ~results(PLOT_IDX).loaded
    fprintf('PLOT_IDX=%d is invalid or record failed to load.\n', PLOT_IDX);
    return;
end

No_plot  = File_No(PLOT_IDX);
mwi_p    = all_mwi{PLOT_IDX};
ecg_p    = all_ecg_rs{PLOT_IDX};
ecg_f_p  = all_ecg_flt{PLOT_IDX};
qrs_p    = all_qrs{PLOT_IDX};
N_p      = length(ecg_p);
t_p      = (0:N_p-1) / fs_new;

% Load truth indices for plot (safe — may not exist)
truth_p = [];
try
    [a_l, a_t] = rdann(fullfile(DB_PATH, num2str(No_plot)), 'atr');
    valid_mask  = ismember(a_t, VALID_SYMS);
    truth_p     = round(a_l(valid_mask) * (fs_new / fs_orig));
    truth_p     = truth_p(truth_p > 0 & truth_p <= N_p);
catch
    % No .atr — plot without truth markers
end

% Guard detection indices
qrs_p_safe    = qrs_p(qrs_p > 0 & qrs_p <= N_p);

% Determine BPM/diagnosis string for title
bpm_p = all_bpm(PLOT_IDX);
if isnan(bpm_p)
    diag_str = 'Insufficient Beats';
elseif bpm_p > 100
    diag_str = sprintf('Tachycardia (%.0f BPM)', bpm_p);
elseif bpm_p < 60
    diag_str = sprintf('Bradycardia (%.0f BPM)', bpm_p);
else
    diag_str = sprintf('Normal Sinus Rhythm (%.0f BPM)', bpm_p);
end

% ---- Figure ----
fig = figure('Name', sprintf('ECG Batch Sanity Check — Record %d', No_plot), ...
             'NumberTitle', 'off', 'Position', [40, 40, 1150, 720]);

XLIM_S = [0, min(10, t_p(end))]; % Show first 10 seconds

% -- Subplot 1: Raw ECG --
ax1 = subplot(3,1,1);
plot(t_p, ecg_p, 'Color', [0.20 0.40 0.80], 'LineWidth', 0.8, ...
     'DisplayName', 'Raw ECG');
hold on;
if ~isempty(truth_p)
    plot(t_p(truth_p), ecg_p(truth_p), 'go', ...
         'MarkerSize', 9, 'LineWidth', 1.8, 'DisplayName', 'Expert Annotations');
end
if ~isempty(qrs_p_safe)
    plot(t_p(qrs_p_safe), ecg_p(qrs_p_safe), 'rx', ...
         'MarkerSize', 9, 'LineWidth', 2.0, 'DisplayName', 'Detected R-peaks');
end
title(sprintf('Raw ECG — Record %d — %s', No_plot, diag_str), 'FontWeight', 'bold');
ylabel('Amplitude (mV)'); legend('Location','northeast'); grid on;
xlim(XLIM_S);

% -- Subplot 2: Filtered ECG --
ax2 = subplot(3,1,2);
plot(t_p, ecg_f_p, 'Color', [0.10 0.60 0.30], 'LineWidth', 0.8, ...
     'DisplayName', 'Filtered ECG');
hold on;
if ~isempty(qrs_p_safe)
    plot(t_p(qrs_p_safe), ecg_f_p(qrs_p_safe), 'rx', ...
         'MarkerSize', 9, 'LineWidth', 2.0, 'DisplayName', 'Detected R-peaks');
end
title('Bandpass Filtered ECG (5–15 Hz) + R-peaks', 'FontWeight', 'bold');
ylabel('Amplitude'); legend('Location','northeast'); grid on;
xlim(XLIM_S);

% -- Subplot 3: MWI --
ax3 = subplot(3,1,3);
plot(t_p, mwi_p, 'k', 'LineWidth', 0.8, 'DisplayName', 'MWI Energy');
hold on;
% Show MWI-domain peaks before back-tracing for debugging transparency
mwi_qrs_raw = all_qrs{PLOT_IDX}; % after back-trace; use mwi directly
plot(t_p(qrs_p_safe), mwi_p(qrs_p_safe), 'rx', ...
     'MarkerSize', 9, 'LineWidth', 2.0, 'DisplayName', 'Detection Points');
title('Moving Window Integration (MWI) Energy Envelope', 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Normalized Energy');
legend('Location','northeast'); grid on;
xlim(XLIM_S);

linkaxes([ax1, ax2, ax3], 'x');

% Compatible figure title (sgtitle not available in older MATLAB)
annotation('textbox', [0, 0.97, 1, 0.03], ...
    'String',              sprintf('Pan-Tompkins Batch Processor — Record %d', No_plot), ...
    'FontSize',            13, ...
    'FontWeight',          'bold', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor',           'none', ...
    'FitBoxToText',        'off');

fprintf('Sanity check plot shown for Record %d.\n', No_plot);
fprintf('Change PLOT_IDX in Section 5 to inspect a different record.\n');