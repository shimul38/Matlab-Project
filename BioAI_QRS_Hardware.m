%% =========================================================================
%  Bio-AI Clinical ECG Analyzer — Pan-Tompkins QRS Detector
%  Author    : Snehashish (BUET, Applied Bio-AI)
%  Supervisor: Dr. S.M. Mahabubur Rahman
%  Version   : 3.0  — Hardware / Text-File Mode
%
%  Accepts plain-text ECG files (.txt, .csv, .dat) recorded from
%  any hardware ECG module.  No WFDB toolbox required.
%  Truth-annotation markers are removed (no .atr file needed).
%
%  Supported text formats (auto-detected):
%    • Single column of raw ADC / voltage values, one sample per line
%    • Two columns: [timestamp, value]  — only the value column is used
%    • Header lines starting with #, %, or letters are skipped
%    • Comma-separated or whitespace-separated
%% =========================================================================
clc; clear; close all;

%% =========================================================================
%  SECTION 1: Configuration — SET YOUR HARDWARE PARAMETERS HERE
%% =========================================================================

% --- Sampling rate of your ECG hardware module (Hz) ---
% Common values: 250, 360, 500, 1000
% If you are unsure, check your module's datasheet.
FS_HARDWARE = 360;      % <-- change this to match your hardware

% --- If your signal is in raw ADC counts, set a scale factor ---
% For voltage signals already in mV, use 1.0
% For Arduino 10-bit ADC at 5V reference: scale = 5000/1023  (gives mV)
% For Arduino 12-bit ADC at 3.3V:         scale = 3300/4095
% Set to 1.0 if your file already contains voltage in mV or arbitrary units
ADC_SCALE = 1.0;        % <-- adjust if needed

% --- Processing sample rate (Pan-Tompkins designed for ~200 Hz) ---
FS_NEW = 200;

%% =========================================================================
%  SECTION 2: File Selection & Smart Loading
%% =========================================================================
[file, fpath] = uigetfile( ...
    {'*.txt;*.csv;*.dat;*.TXT;*.CSV', 'ECG Text Files (*.txt, *.csv, *.dat)'; ...
     '*.*', 'All Files (*.*)'}, ...
    'Select your ECG text file');

if isequal(file, 0)
    disp('User cancelled file selection.');
    return;
end

full_path   = fullfile(fpath, file);
[~, record_name, ext] = fileparts(file);
fprintf('\n--- Loading Hardware ECG: %s ---\n', file);

% ── Smart reader: skip header lines, handle comma or whitespace delimiters
fid = fopen(full_path, 'r');
if fid == -1
    errordlg(['Cannot open file: ', full_path], 'File Error');
    return;
end

raw_lines  = textscan(fid, '%s', 'Delimiter', '\n', 'WhiteSpace', '');
fclose(fid);
raw_lines  = raw_lines{1};

% Filter out comment / header lines
data_lines = {};
for k = 1:length(raw_lines)
    ln = strtrim(raw_lines{k});
    if isempty(ln), continue; end
    first_char = ln(1);
    % Skip lines that start with non-numeric characters (headers/comments)
    if ismember(first_char, {'#','%','/','!'}) || isletter(first_char)
        fprintf('  [Skipped header line %d]: %s\n', k, ln);
        continue;
    end
    data_lines{end+1} = ln; %#ok<AGROW>
end

if isempty(data_lines)
    errordlg('No numeric data found in the file. Check the format.', 'Parse Error');
    return;
end

% Parse numeric values — handle comma or space/tab separation
num_cols = length(strsplit(strtrim(data_lines{1}), {' ','\t',','}));
raw_matrix = zeros(length(data_lines), num_cols);
parse_errors = 0;
for k = 1:length(data_lines)
    parts = strsplit(strtrim(data_lines{k}), {' ','\t',','}, 'CollapseDelimiters', true);
    vals  = str2double(parts);
    if any(isnan(vals)) || length(vals) ~= num_cols
        parse_errors = parse_errors + 1;
        continue;
    end
    raw_matrix(k,:) = vals;
end

if parse_errors > 0
    fprintf('  Warning: %d lines could not be parsed and were skipped.\n', parse_errors);
end

% Remove any all-zero rows from parse failures
raw_matrix = raw_matrix(any(raw_matrix ~= 0, 2), :);

% Extract the ECG column
if num_cols == 1
    % Single column — ECG samples directly
    ecg_raw = raw_matrix(:, 1);
    fprintf('  Detected: single-column format (%d samples)\n', length(ecg_raw));

elseif num_cols >= 2
    % Two+ columns — assume last column is ECG, first may be timestamp/index
    % If first column is monotonically increasing (timestamp/index), use last col
    col1 = raw_matrix(:,1);
    if all(diff(col1) > 0)
        ecg_raw = raw_matrix(:, end);
        fprintf('  Detected: multi-column format — using column %d as ECG (%d samples)\n', ...
                num_cols, length(ecg_raw));
    else
        % Not clearly a timestamp — use second column (common Arduino format)
        ecg_raw = raw_matrix(:, 2);
        fprintf('  Detected: multi-column format — using column 2 as ECG (%d samples)\n', length(ecg_raw));
    end
end

% Apply ADC scale factor
ecg_raw = double(ecg_raw) * ADC_SCALE;

% Basic sanity check
if length(ecg_raw) < FS_HARDWARE * 2
    warndlg( ...
        sprintf('Signal is only %.1f seconds long. Results may be unreliable.', ...
                length(ecg_raw)/FS_HARDWARE), ...
        'Short Recording Warning');
end

fprintf('  Hardware fs   : %d Hz\n', FS_HARDWARE);
fprintf('  Total samples : %d  (%.1f seconds)\n', length(ecg_raw), length(ecg_raw)/FS_HARDWARE);
fprintf('  Amplitude range: [%.2f, %.2f]\n', min(ecg_raw), max(ecg_raw));

%% =========================================================================
%  SECTION 3: Resampling to 200 Hz
%% =========================================================================
if FS_HARDWARE ~= FS_NEW
    ecg_resampled = resample(ecg_raw, FS_NEW, FS_HARDWARE);
    fprintf('  Resampled: %d Hz → %d Hz\n', FS_HARDWARE, FS_NEW);
else
    ecg_resampled = ecg_raw;
    fprintf('  No resampling needed (hardware fs = %d Hz)\n', FS_NEW);
end

% Ensure column vector
ecg_resampled = ecg_resampled(:);

N       = length(ecg_resampled);
time_ax = (0:N-1) / FS_NEW;   % time axis in seconds

%% =========================================================================
%  SECTION 4: Pan-Tompkins Pipeline (All Zero-Phase)
%
%  filtfilt() is used at every filtering stage → ZERO group delay.
%  This means detections will align with R-peaks in the raw signal
%  without any manual shift correction.
%% =========================================================================

% Step 4a: Bandpass filter 5–15 Hz (Pan-Tompkins specification)
[b_bp, a_bp] = butter(2, [5 15] / (FS_NEW/2), 'bandpass');
filtered_ecg  = filtfilt(b_bp, a_bp, ecg_resampled);

% Step 4b: 5-point derivative  h = (fs/8)·[-2,-1,0,1,2]
b_der    = (FS_NEW/8) * [-2, -1, 0, 1, 2];
diff_ecg = filtfilt(b_der, 1, filtered_ecg);

% Step 4c: Squaring
ecg_squared = diff_ecg .^ 2;

% Step 4d: Moving Window Integration (150 ms symmetric window)
win_len    = round(0.150 * FS_NEW);
mwi_signal = movmean(ecg_squared, win_len);

% Step 4e: Normalize
mwi_signal = mwi_signal / max(mwi_signal);

%% =========================================================================
%  SECTION 5: Adaptive QRS Detection on MWI Signal
%% =========================================================================

% Initial peak candidates
min_h = median(mwi_signal) * 0.8;
[pks, locs] = findpeaks(mwi_signal, ...
    'MinPeakDistance', round(0.200 * FS_NEW), ...
    'MinPeakHeight',   min_h);

% Adaptive training on first 2 seconds
training_end = 2 * FS_NEW;
train_idx    = locs <= training_end;

if any(train_idx)
    spki = mean(pks(train_idx));
    npki = mean(mwi_signal(1:min(training_end, N))) * 0.5;
else
    spki = max(pks) * 0.4;
    npki = mean(mwi_signal) * 0.1;
end

T1     = npki + 0.25 * (spki - npki);   % Primary threshold
T2     = 0.5 * T1;                       % Search-back threshold
avg_rr = round(0.8 * FS_NEW);            % Initial RR estimate (~75 BPM)

qrs_mwi = [];

for j = 1:length(pks)
    cp_loc = locs(j);
    cp_val = pks(j);

    % Refractory period guard (200 ms) with double-peak replacement
    if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) < round(0.200 * FS_NEW)
        if cp_val > mwi_signal(qrs_mwi(end))
            qrs_mwi(end) = cp_loc;
        end
        continue;
    end

    % Threshold decay after 1.5 s silence
    if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) > round(1.5 * FS_NEW)
        T1 = T1 * 0.5;
        T2 = 0.5 * T1;
    end

    % Search-back for missed beat (RR gap > 1.66× average)
    if ~isempty(qrs_mwi) && (cp_loc - qrs_mwi(end)) > round(1.66 * avg_rr)
        in_gap = find(locs > qrs_mwi(end) & locs < cp_loc);
        for m = in_gap'
            if pks(m) > 0.75 * T2
                qrs_mwi = [qrs_mwi; locs(m)]; %#ok<AGROW>
                spki    = 0.25 * pks(m) + 0.75 * spki;
                T1      = npki + 0.25 * (spki - npki);
                T2      = 0.5 * T1;
                break;
            end
        end
    end

    % Standard detection
    if cp_val > T1
        qrs_mwi = [qrs_mwi; cp_loc]; %#ok<AGROW>
        capped  = min(cp_val, 2 * spki);
        spki    = 0.125 * capped + 0.875 * spki;
    else
        npki    = 0.125 * cp_val + 0.875 * npki;
    end

    T1 = npki + 0.25 * (spki - npki);
    T2 = 0.5 * T1;
    if length(qrs_mwi) > 2
        avg_rr = mean(diff(qrs_mwi(max(1, end-4):end)));
    end
end

%% =========================================================================
%  SECTION 6: R-Peak Back-Tracing
%  Snaps each MWI detection to the true filtered-ECG peak within ±100 ms.
%  This eliminates the MWI envelope lag without any manual offset constant.
%% =========================================================================

search_window = round(0.100 * FS_NEW);
qrs_refined   = zeros(size(qrs_mwi));

for i = 1:length(qrs_mwi)
    lo = max(1,   qrs_mwi(i) - search_window);
    hi = min(N,   qrs_mwi(i) + search_window);
    [~, local_idx]  = max(abs(filtered_ecg(lo:hi)));
    qrs_refined(i)  = lo + local_idx - 1;
end

qrs_refined = unique(qrs_refined);

%% =========================================================================
%  SECTION 7: Clinical Metrics & HRV Analysis
%% =========================================================================
rr_ms    = diff(qrs_refined) / FS_NEW * 1000;
rr_clean = rr_ms(rr_ms > 300 & rr_ms < 2000);

if length(rr_clean) > 1
    avg_BPM   = 60000 / mean(rr_clean);
    sdnn_val  = std(rr_clean);
    rmssd_val = sqrt(mean(diff(rr_clean).^2));
    pnn50_val = 100 * sum(abs(diff(rr_clean)) > 50) / length(rr_clean);
else
    avg_BPM   = NaN;
    sdnn_val  = NaN;
    rmssd_val = NaN;
    pnn50_val = NaN;
    warning('Fewer than 2 valid RR intervals — HRV metrics cannot be computed.');
end

% Rhythm classification
if isnan(avg_BPM)
    diagnosis = 'Insufficient beats for classification';
elseif avg_BPM > 100
    diagnosis = 'Tachycardia Detected (High Heart Rate)';
elseif avg_BPM < 60
    diagnosis = 'Bradycardia Detected (Low Heart Rate)';
else
    diagnosis = 'Normal Sinus Rhythm';
end

%% =========================================================================
%  SECTION 8: Console Report
%% =========================================================================
fprintf('\n========================================\n');
fprintf('   Clinical Diagnostic Report\n');
fprintf('========================================\n');
fprintf('File              : %s\n',  file);
fprintf('Hardware fs       : %d Hz\n', FS_HARDWARE);
fprintf('Total Duration    : %.1f s\n', N / FS_NEW);
fprintf('Beats Detected    : %d\n',  length(qrs_refined));
if ~isnan(avg_BPM)
    fprintf('Avg Heart Rate    : %.1f BPM\n', avg_BPM);
    fprintf('HRV (SDNN)        : %.1f ms\n',  sdnn_val);
    fprintf('HRV (RMSSD)       : %.1f ms\n',  rmssd_val);
    fprintf('HRV (pNN50)       : %.1f%%\n',   pnn50_val);
end
fprintf('Diagnosis         : %s\n',  diagnosis);
fprintf('========================================\n\n');

%% =========================================================================
%  SECTION 9: GUI Popup Summary
%% =========================================================================
report_str = {
    ['File       : ', file];
    ['Heart Rate : ', sprintf('%.1f', avg_BPM), ' BPM'];
    ['Status     : ', diagnosis];
    '';
    ['SDNN       : ', sprintf('%.1f', sdnn_val),  ' ms'];
    ['RMSSD      : ', sprintf('%.1f', rmssd_val), ' ms'];
    ['pNN50      : ', sprintf('%.1f', pnn50_val), '%'];
    '';
    ['Beats      : ', num2str(length(qrs_refined))];
    ['Duration   : ', sprintf('%.1f', N/FS_NEW), ' s'];
};

msgbox(report_str, ['ECG Summary — ', record_name]);

%% =========================================================================
%  SECTION 10: Visualization (3-Subplot Layout)
%
%  Green truth-annotation markers are removed — this is hardware mode.
%  Only detected R-peaks (red ×) are shown on the signal plots.
%% =========================================================================
qrs_plot     = qrs_refined(qrs_refined > 0 & qrs_refined <= N);
qrs_mwi_plot = qrs_mwi(qrs_mwi > 0 & qrs_mwi <= N);

% Display up to 15 seconds (enough to see several beats clearly)
xlim_s = [0, min(15, time_ax(end))];

fig = figure('Name', ['Hardware ECG — ', file], ...
             'NumberTitle', 'off', ...
             'Position',    [40, 40, 1100, 680]);

% ── Subplot 1: Raw (resampled) ECG + detected R-peaks ───────────────────
ax1 = subplot(3,1,1);
plot(time_ax, ecg_resampled, 'Color', [0.15 0.35 0.75], 'LineWidth', 0.8, ...
     'DisplayName', 'Raw ECG (resampled)');
hold on;
if ~isempty(qrs_plot)
    plot(time_ax(qrs_plot), ecg_resampled(qrs_plot), 'rx', ...
         'MarkerSize', 9, 'LineWidth', 2.0, 'DisplayName', 'Detected R-peaks');
end

if ~isnan(avg_BPM)
    title_str = sprintf('Raw ECG  —  %s  (%.1f BPM)', diagnosis, avg_BPM);
else
    title_str = 'Raw ECG  —  Insufficient beats detected';
end
title(title_str, 'FontWeight', 'bold');
ylabel('Amplitude'); legend('Location', 'best');
grid on; xlim(xlim_s);

% ── Subplot 2: Bandpass-filtered ECG + R-peaks ──────────────────────────
ax2 = subplot(3,1,2);
plot(time_ax, filtered_ecg, 'Color', [0.08 0.55 0.28], 'LineWidth', 0.8, ...
     'DisplayName', 'Filtered ECG (5–15 Hz)');
hold on;
if ~isempty(qrs_plot)
    plot(time_ax(qrs_plot), filtered_ecg(qrs_plot), 'rx', ...
         'MarkerSize', 9, 'LineWidth', 2.0, 'DisplayName', 'R-peaks (back-traced)');
end
title('Bandpass Filtered ECG (5–15 Hz) — R-peak markers land on filtered peak', ...
      'FontWeight', 'bold');
ylabel('Amplitude'); legend('Location', 'best');
grid on; xlim(xlim_s);

% ── Subplot 3: MWI energy envelope + MWI detection points ───────────────
ax3 = subplot(3,1,3);
plot(time_ax, mwi_signal, 'k', 'LineWidth', 0.8, 'DisplayName', 'MWI Energy');
hold on;
if ~isempty(qrs_mwi_plot)
    plot(time_ax(qrs_mwi_plot), mwi_signal(qrs_mwi_plot), 'rx', ...
         'MarkerSize', 9, 'LineWidth', 2.0, 'DisplayName', 'MWI Detections (pre-refinement)');
end
title('Moving Window Integration (MWI) Energy Envelope', 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Normalized Energy');
legend('Location', 'best'); grid on; xlim(xlim_s);

linkaxes([ax1, ax2, ax3], 'x');

% Figure-level title (compatible with older MATLAB)
annotation('textbox', [0, 0.97, 1, 0.03], ...
    'String',              ['Pan-Tompkins QRS Detector  —  Hardware Recording: ', file], ...
    'FontSize',            12, ...
    'FontWeight',          'bold', ...
    'HorizontalAlignment', 'center', ...
    'EdgeColor',           'none', ...
    'FitBoxToText',        'off');

fprintf('Done. Plot window shows first %.0f seconds. Use axis zoom to inspect further.\n', xlim_s(2));
