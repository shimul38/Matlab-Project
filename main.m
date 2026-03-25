%% SECTION 1: Initialization & Data Loading
clc; clear; close all;

% Automatically add the toolbox to the path
if ~exist('rdsamp', 'file')
    addpath('B:\EEE 376\Labtest materials\lab-05\mcode');
end
wfdbloadlib; 

File_No = [100:109, 111:119, 200:234]; 
fs_orig = 360; 
fs_new = 200;  
num_files = length(File_No);
tf_ecg = cell(num_files, 1); 

%% SECTION 2: Pan-Tompkins Processing Loop
fprintf('--- Starting Signal Processing ---\n');
for i = 1:num_files
    No = File_No(i);
    record_path = ['mit-bih-arrhythmia-database-1.0.0\', num2str(No)]; 
    
    try
        [signal, ~] = rdsamp(record_path, 1);
        ecg_resampled = resample(signal, fs_new, fs_orig);
        filtered_ecg = pan_tompkins_bandpass(ecg_resampled, fs_new);
        
        b_der = [2 1 0 -1 -2] * (fs_new / 8); 
        diff_ecg = derivative_func(filtered_ecg, b_der, 1);
        
        ecg_squared = diff_ecg .^ 2; 
        ecg_squared = ecg_squared / max(abs(ecg_squared));
        
        window_length = 30; % Optimized window
        ecg_mwi = movingWindowIntegrator(ecg_squared, window_length);
        
        tf_ecg{i} = ecg_mwi;
        fprintf('Processed Record: %d\n', No);
    catch
        fprintf('Error loading Record: %d (Skipping)\n', No);
    end
end

%% SECTION 3: Adaptive Thresholding & Validation
fprintf('\n--- Starting Detection & Validation ---\n');
all_bpm = zeros(num_files, 1); 
all_qrs_indices = cell(num_files, 1); 
stats = struct('TP',0,'FP',0,'FN',0,'Se',0,'Pp',0);
results_log = repmat(stats, num_files, 1);

tolerance = 0.150 * fs_new; 
filter_delay_seconds = 0; 
delay_samples = 0;

for i = 1:num_files
    mwi_signal = tf_ecg{i};
    No = File_No(i);
    
    if isempty(mwi_signal), continue; end
    
    % --- Robust Normalization (Outlier Suppression) ---
    mwi_signal = mwi_signal / max(mwi_signal); 
    limit = prctile(mwi_signal, 95); 
    mwi_signal(mwi_signal > 1.5*limit) = 1.5*limit;
    mwi_signal = mwi_signal / max(mwi_signal);
    
    min_h = mean(mwi_signal) * 0.8; %if not working, change it back to 1
    [pks, locs] = findpeaks(mwi_signal, 'MinPeakDistance', 0.2*fs_new, 'MinPeakHeight', min_h);
    
    % --- SECTION 3: Adaptive Training (First 2 Seconds) ---
    training_time = 2 * fs_new; 
    train_idx = find(locs <= training_time);
    if ~isempty(train_idx)
        spki = mean(pks(train_idx)); 
        npki = mean(mwi_signal(1:training_time)) * 0.5; 
    else
        spki = max(pks) * 0.4; 
        npki = mean(mwi_signal) * 0.1;
    end
    
    T1 = npki + 0.25*(spki - npki);
    T2 = 0.5 * T1;
    
    qrs_indices = [];
    avg_rr = 0.8 * fs_new; 
    
    for j = 1:length(pks)
        current_peak_loc = locs(j);
        current_peak_val = pks(j);
        
        % --- ADDED: Threshold Decay (Fix for Record 203) ---
        % If no beat has been detected for 1.5 seconds, the threshold is likely
        % stuck on a large artifact or PVC. We drop the threshold to find smaller beats.
        if ~isempty(qrs_indices) && (current_peak_loc - qrs_indices(end) > 1.5 * fs_new)
            T1 = T1 * 0.5; 
            T2 = 0.5 * T1;
        end
        
        % 1. Search-Back Logic
        if ~isempty(qrs_indices)
            if (current_peak_loc - qrs_indices(end)) > (1.66 * avg_rr)
                potential_misses = find(locs > qrs_indices(end) & locs < current_peak_loc);
                for m = potential_misses'
                    if pks(m) > (0.75 * T2)
                        qrs_indices = [qrs_indices; locs(m)];
                        spki = 0.25*pks(m) + 0.75*spki; 
                        % Update thresholds immediately after finding a missed beat
                        T1 = npki + 0.25*(spki - npki);
                        T2 = 0.5 * T1;
                        break; 
                    end
                end
            end
        end
        
        % 2. Standard Detection
        % --- Inside the (current_peak_val > T1) block ---
        if current_peak_val > T1
            qrs_indices = [qrs_indices; current_peak_loc];
    
            % --- ADDED: Saturation Cap ---
            % Don't let a single peak be more than 2x the current spki
            % This stops massive PVCs from "skyrocketing" your threshold
            capped_val = min(current_peak_val, 2 * spki); 
            spki = 0.125*capped_val + 0.875*spki;
        else
            npki = 0.125*current_peak_val + 0.875*npki;
        end
        
        % Final Threshold Update
        T1 = npki + 0.25*(spki - npki);
        T2 = 0.5 * T1;
        
        if length(qrs_indices) > 2
            avg_rr = mean(diff(qrs_indices(max(1, end-4):end))); 
        end
    end
    
    all_qrs_indices{i} = qrs_indices; 
    
    % 3. BPM Calculation
    if length(qrs_indices) > 1
        all_bpm(i) = 60 / (mean(diff(qrs_indices)) / fs_new);
    else
        all_bpm(i) = NaN; 
    end
    
    % 4. Validation
    try
        [ann_locs, ann_type] = rdann(['mit-bih-arrhythmia-database-1.0.0\', num2str(No)], 'atr');
        valid_beats = ismember(ann_type, {'N','L','R','A','a','J','S','V','F','e','j','E','/'});
        truth_indices = round(ann_locs(valid_beats) * (fs_new / fs_orig)) + delay_samples;
        
        tp = 0;
        for k = 1:length(truth_indices)
            if any(abs(qrs_indices - truth_indices(k)) <= tolerance), tp = tp + 1; end
        end
        
        results_log(i).TP = tp;
        results_log(i).FP = length(qrs_indices) - tp;
        results_log(i).FN = length(truth_indices) - tp;
        results_log(i).Se = (tp / (tp + results_log(i).FN)) * 100;
        results_log(i).Pp = (tp / (tp + results_log(i).FP)) * 100;
    catch
        fprintf('Warning: Failed to validate Record %d\n', No);
    end
end

%% SECTION 4: Results & Export
avg_Se = mean([results_log.Se], 'omitnan');
avg_Pp = mean([results_log.Pp], 'omitnan');
avg_Acc = (sum([results_log.TP]) / (sum([results_log.TP]) + sum([results_log.FP]) + sum([results_log.FN]))) * 100;

fprintf('\n--- Final Pipeline Results ---\n');
fprintf('Average Sensitivity (Se): %.2f%%\n', avg_Se);
fprintf('Average Positive Predictivity (+P): %.2f%%\n', avg_Pp);
fprintf('Overall Detection Accuracy: %.2f%%\n', avg_Acc);

results_table = table(File_No', all_bpm, [results_log.Se]', [results_log.Pp]', ...
    'VariableNames', {'Record_ID', 'BPM', 'Sensitivity', 'Pos_Predictivity'});
writetable(results_table, 'ECG_Full_Analysis.xlsx');

%% SECTION 5: Visual Sanity Check
idx = 23; % Adjust this to check Record 203 (index 23 or 24 depending on File_No)
mwi_p = tf_ecg{idx} / max(tf_ecg{idx}); 
qrs_p = all_qrs_indices{idx};

[a_l, a_t] = rdann(['mit-bih-arrhythmia-database-1.0.0\', num2str(File_No(idx))], 'atr');
t_p = round(a_l(ismember(a_t, {'N','L','R','A','a','J','S','V','F','e','j','E','/'})) * (fs_new/fs_orig));

shift_constant = 11; 
t_p = t_p + shift_constant;

t_p_valid = t_p(t_p > 10 & t_p <= length(mwi_p) - 10);
t_p_heights = zeros(size(t_p_valid));
for k = 1:length(t_p_valid)
    t_p_heights(k) = max(mwi_p(t_p_valid(k)-10 : t_p_valid(k)+10));
end

figure; plot(mwi_p); hold on;
plot(t_p_valid, t_p_heights, 'go', 'MarkerSize', 12, 'LineWidth', 2); 
plot(qrs_p, mwi_p(qrs_p), 'rx', 'MarkerSize', 10);
legend('Signal', 'Truth', 'Detections'); 
title(['Validation: Record ', num2str(File_No(idx))]);
xlim([0 2000]); grid on;