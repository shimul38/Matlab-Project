%% =========================================================================
%  Bio-AI Clinical ECG Analyzer — Pan-Tompkins QRS Detector
%  Author   : Snehashish (BUET, Applied Bio-AI)
%  Supervisor: Dr. S.M. Mahabubur Rahman
%  Version  : 2.1 (sgtitle fix, .atr-optional, zero-phase pipeline)
%% =========================================================================
clc; clear; close all;

%% =========================================================================
%  SECTION 1: Initialization & File Selection
%% =========================================================================
if ~exist('rdsamp', 'file')
    addpath('B:\EEE 376\Labtest materials\lab-05\mcode'); % Adjust path if needed
end
wfdbloadlib;

[file, path] = uigetfile('*.dat', 'Select an MIT-BIH ECG Record (.dat)');
if isequal(file, 0)
    disp('User canceled file selection.');
    return;
end

[~, record_name, ~] = fileparts(file);
record_path = fullfile(path, record_name);
fprintf('--- Analyzing Clinical Record: %s ---\n', record_name);

%% =========================================================================
%  SECTION 2: Signal Loading & Resampling
%% =========================================================================
fs_orig = 360;
fs_new  = 200; % Pan-Tompkins was designed for ~200 Hz

original_dir = pwd;
try
    cd(path);
    [signal, ~] = rdsamp(record_name, 1);
    cd(original_dir);
catch ME
    cd(original_dir);
    errordlg(['Failed to load: ', record_name, '. Ensure .hea and .dat are in ', path], 'File Error');
    return;
end

% Resample from 360 Hz → 200 Hz
ecg_resampled = resample(signal, fs_new, fs_orig);
N = length(ecg_resampled);
time_ax = (0:N-1) / fs_new; % Time axis in seconds

%% =========================================================================
%  SECTION 3: Pan-Tompkins Signal Processing Pipeline (All Zero-Phase)
%
%  KEY DESIGN CHOICE: Use filtfilt() for ALL stages so there is ZERO
%  group-delay introduced. This means the MWI peaks and R-peaks will be
%  temporally aligned — we only need to back-trace within a short window.
%% =========================================================================

% --- Step 3a: Bandpass Filter (5–15 Hz per Pan-Tompkins) ---
% Using a Butterworth bandpass for a sharper passband
[b_bp, a_bp] = butter(2, [5 15] / (fs_new/2), 'bandpass');
filtered_ecg = filtfilt(b_bp, a_bp, ecg_resampled);

% --- Step 3b: Derivative (Pan-Tompkins 5-point) ---
% Correct Pan-Tompkins derivative coefficients: h = (1/8T)[-2,-1,0,1,2]
% For zero-phase, we apply filtfilt. Note: coefficients are scaled by fs.
b_der = (fs_new/8) * [-2, -1, 0, 1, 2];
diff_ecg = filtfilt(b_der, 1, filtered_ecg);

% --- Step 3c: Squaring (non-linear amplification) ---
ecg_squared = diff_ecg .^ 2;

% --- Step 3d: Moving Window Integration ---
% Window ≈ 150 ms as specified in Pan-Tompkins
window_length = round(0.150 * fs_new); % 30 samples at 200 Hz
mwi_signal = movmean(ecg_squared, window_length); % movmean = symmetric MWI

% --- Step 3e: Normalize MWI ---
mwi_signal = mwi_signal / max(mwi_signal);

%% =========================================================================
%  SECTION 4: Adaptive QRS Detection on MWI Signal
%% =========================================================================

% --- Initial Peak Finding ---
min_h = median(mwi_signal) * 0.8;
[pks, locs] = findpeaks(mwi_signal, ...
    'MinPeakDistance', round(0.200 * fs_new), ... % 200 ms refractory
    'MinPeakHeight',   min_h);

% --- Adaptive Training Window (first 2 seconds) ---
training_time = 2 * fs_new;
train_idx = find(locs <= training_time);

if ~isempty(train_idx)
    spki = mean(pks(train_idx));
    npki = mean(mwi_signal(1:training_time)) * 0.5;
else
    spki = max(pks) * 0.4;
    npki = mean(mwi_signal) * 0.1;
end

T1 = npki + 0.25 * (spki - npki); % Primary threshold
T2 = 0.5 * T1;                     % Search-back threshold

% --- Main Detection Loop ---
qrs_mwi = []; % MWI-domain detections
avg_rr   = round(0.8 * fs_new); % Initial RR estimate (0.8 s = 75 BPM)

for j = 1:length(pks)
    cp_loc = locs(j);
    cp_val = pks(j);

    % Refractory Period Check (200 ms hard minimum)
    if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) < round(0.200 * fs_new)
        % If this peak is larger, replace the previous one
        if cp_val > mwi_signal(qrs_mwi(end))
            qrs_mwi(end) = cp_loc;
        end
        continue;
    end

    % Threshold Decay: if gap > 1.5 s, loosen thresholds
    if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) > round(1.5 * fs_new)
        T1 = T1 * 0.5;
        T2 = 0.5 * T1;
    end

    % Search-Back: if gap > 1.66 × avg_rr, scan for a missed beat
    if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) > round(1.66 * avg_rr)
        in_gap = find(locs > qrs_mwi(end) & locs < cp_loc);
        for m = in_gap'
            if pks(m) > (0.75 * T2)
                qrs_mwi   = [qrs_mwi; locs(m)];
                spki      = 0.25 * pks(m) + 0.75 * spki;
                T1        = npki + 0.25 * (spki - npki);
                T2        = 0.5 * T1;
                break;
            end
        end
    end

    % Standard Detection
    if cp_val > T1
        qrs_mwi         = [qrs_mwi; cp_loc];
        capped_val      = min(cp_val, 2 * spki);
        spki            = 0.125 * capped_val + 0.875 * spki;
    else
        npki            = 0.125 * cp_val + 0.875 * npki;
    end

    % Update thresholds and RR average
    T1 = npki + 0.25 * (spki - npki);
    T2 = 0.5 * T1;
    if length(qrs_mwi) > 2
        avg_rr = mean(diff(qrs_mwi(max(1, end-4):end)));
    end
end

%% =========================================================================
%  SECTION 5: R-Peak Back-Tracing (THE CORE FIX)
%
%  Even with zero-phase filters, movmean() introduces a right-skew in the
%  MWI energy envelope. The MWI peak ≠ R-peak location in the raw signal.
%
%  FIX: For each MWI detection, search backward in the FILTERED ECG within
%  a ±search_window and find the true local maximum. This is the proper
%  Pan-Tompkins back-tracing step (Section V of the original paper).
%% =========================================================================

search_window = round(0.100 * fs_new); % ±100 ms search window

qrs_refined = zeros(size(qrs_mwi)); % Will hold refined R-peak locations

for i = 1:length(qrs_mwi)
    center = qrs_mwi(i);
    
    % Define safe search bounds
    lo = max(1,   center - search_window);
    hi = min(N,   center + search_window);
    
    % Find the maximum of the FILTERED ECG in this window
    segment = filtered_ecg(lo:hi);
    [~, local_idx] = max(abs(segment)); % Use abs() to handle inverted leads
    qrs_refined(i) = lo + local_idx - 1;
end

% Remove any duplicate indices that may arise from overlapping windows
qrs_refined = unique(qrs_refined);

%% =========================================================================
%  SECTION 6: Clinical Metrics & HRV Analysis
%% =========================================================================
rr_intervals_ms = diff(qrs_refined) / fs_new * 1000; % RR in milliseconds

% Remove physiologically impossible RR intervals (< 300 ms or > 2000 ms)
rr_clean = rr_intervals_ms(rr_intervals_ms > 300 & rr_intervals_ms < 2000);

avg_BPM = 60000 / mean(rr_clean);
sdnn    = std(rr_clean);
rmssd   = sqrt(mean(diff(rr_clean).^2));
pnn50   = sum(abs(diff(rr_clean)) > 50) / length(rr_clean) * 100; % pNN50

% Diagnostic Rule
if avg_BPM > 100
    diagnosis = 'Tachycardia Detected (High Heart Rate)';
elseif avg_BPM < 60
    diagnosis = 'Bradycardia Detected (Low Heart Rate)';
else
    diagnosis = 'Normal Sinus Rhythm';
end

%% =========================================================================
%  SECTION 7: Truth Validation (MIT-BIH Annotations)
%% =========================================================================
Se = NaN; Pp = NaN; Acc = NaN;
truth_indices = [];

try
    cd(path);
    [ann_locs, ann_type] = rdann(record_name, 'atr');
    cd(original_dir);

    % Valid beat types per MIT-BIH standard
    valid_beats   = ismember(ann_type, {'N','L','R','A','a','J','S','V','F','e','j','E','/'});
    truth_indices = round(ann_locs(valid_beats) * (fs_new / fs_orig));

    % Tolerance: ±150 ms (standard AHA evaluation criterion)
    tolerance = round(0.150 * fs_new);
    tp = 0;
    for k = 1:length(truth_indices)
        if any(abs(qrs_refined - truth_indices(k)) <= tolerance)
            tp = tp + 1;
        end
    end

    fp  = length(qrs_refined) - tp;
    fn  = length(truth_indices) - tp;
    Se  = (tp / (tp + fn)) * 100;
    Pp  = (tp / (tp + fp)) * 100;
    Acc = (tp / (tp + fp + fn)) * 100;

catch
    cd(original_dir);
    disp('No annotation (.atr) file found — skipping validation.');
end

%% =========================================================================
%  SECTION 8: Terminal Report
%% =========================================================================
fprintf('\n========================================\n');
fprintf('   Clinical Diagnostic Report\n');
fprintf('========================================\n');
fprintf('Record ID         : %s\n',   record_name);
fprintf('Beats Detected    : %d\n',   length(qrs_refined));
fprintf('Avg Heart Rate    : %.1f BPM\n', avg_BPM);
fprintf('HRV (SDNN)        : %.1f ms\n',  sdnn);
fprintf('HRV (RMSSD)       : %.1f ms\n',  rmssd);
fprintf('HRV (pNN50)       : %.1f%%\n',   pnn50);
fprintf('Diagnosis         : %s\n',   diagnosis);
fprintf('----------------------------------------\n');

if ~isnan(Se)
    fprintf('Sensitivity (Se)  : %.2f%%\n', Se);
    fprintf('Predictivity (+P) : %.2f%%\n', Pp);
    fprintf('Accuracy          : %.2f%%\n', Acc);
end
fprintf('========================================\n\n');

%% =========================================================================
%  SECTION 9: GUI Popup Summary
%% =========================================================================
report_str = {
    ['Record ID    : ', record_name];
    ['Heart Rate   : ', sprintf('%.1f', avg_BPM), ' BPM  →  ', diagnosis];
    '';
    ['SDNN         : ', sprintf('%.1f', sdnn),  ' ms'];
    ['RMSSD        : ', sprintf('%.1f', rmssd), ' ms'];
    ['pNN50        : ', sprintf('%.1f', pnn50), '%'];
};

if ~isnan(Se)
    report_str = [report_str; {
        '';
        ['Sensitivity  : ', sprintf('%.2f', Se),  '%'];
        ['Predictivity : ', sprintf('%.2f', Pp),  '%'];
        ['Accuracy     : ', sprintf('%.2f', Acc), '%']
    }];
end

msgbox(report_str, ['Clinical Summary — ', record_name]);

%% =========================================================================
%  SECTION 10: Visualization (3-Subplot Layout)
%% =========================================================================
% Clamp indices for safety
qrs_plot     = qrs_refined(qrs_refined > 0 & qrs_refined <= N);
qrs_mwi_plot = qrs_mwi(qrs_mwi > 0 & qrs_mwi <= N);

% Guard: truth_indices may be empty if no .atr file was found
if ~isempty(truth_indices)
    truth_plot = truth_indices(truth_indices > 0 & truth_indices <= N);
else
    truth_plot = [];
end

fig = figure('Name', ['ECG Analysis: ', record_name], ...
             'NumberTitle', 'off', ...
             'Position', [50, 50, 1000, 600]);

ax1 = subplot(3,1,1);
plot(time_ax, ecg_resampled, 'Color', [0.2 0.4 0.8], 'LineWidth', 0.8, 'DisplayName', 'Raw ECG');
hold on;
if ~isempty(truth_plot)
    plot(time_ax(truth_plot), ecg_resampled(truth_plot), ...
         'go', 'MarkerSize', 9, 'LineWidth', 1.8, 'DisplayName', 'Expert Annotations');
end
plot(time_ax(qrs_plot), ecg_resampled(qrs_plot), ...
     'rx', 'MarkerSize', 9, 'LineWidth', 2.0, 'DisplayName', 'Detected R-peaks');
title(['Raw ECG + R-peak Detections — ', diagnosis], 'FontWeight', 'bold');
ylabel('Amplitude (mV)');
legend('Location', 'best');  % Uses DisplayName from each plot call
grid on; xlim([0, min(10, time_ax(end))]);

ax2 = subplot(3,1,2);
plot(time_ax, filtered_ecg, 'Color', [0.1 0.6 0.3], 'LineWidth', 0.8);
hold on;
plot(time_ax(qrs_plot), filtered_ecg(qrs_plot), ...
     'rx', 'MarkerSize', 9, 'LineWidth', 2.0);
title('Bandpass Filtered ECG (5–15 Hz) + R-peaks', 'FontWeight', 'bold');
ylabel('Amplitude');
legend('Filtered ECG', 'R-peaks', 'Location', 'best');
grid on; xlim([0, min(10, time_ax(end))]);

ax3 = subplot(3,1,3);
plot(time_ax, mwi_signal, 'k', 'LineWidth', 0.8);
hold on;
plot(time_ax(qrs_mwi_plot), mwi_signal(qrs_mwi_plot), ...
     'rx', 'MarkerSize', 9, 'LineWidth', 2.0);
title('Moving Window Integration (MWI) Energy Envelope', 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Normalized Energy');
legend('MWI Signal', 'MWI Peaks (Pre-Refinement)', 'Location', 'best');
grid on; xlim([0, min(10, time_ax(end))]);

linkaxes([ax1, ax2, ax3], 'x');

% sgtitle() is unavailable in older MATLAB versions — use annotation instead
annotation('textbox', [0, 0.97, 1, 0.03], ...
    'String',     ['Pan-Tompkins QRS Detector — Record: ', record_name], ...
    'FontSize',   13, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor',  'none', ...
    'FitBoxToText', 'off');