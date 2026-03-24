% Automatically add the toolbox to the path
if ~exist('rdsamp', 'file')
    addpath('B:\EEE 376\Labtest materials\lab-05\mcode');
end
wfdbloadlib; 

File_No = [100:109, 111:119, 200:234]; 
fs_new = 200; 

% Use a Cell Array for safety if you want to keep everything in one variable
tf_ecg = cell(length(File_No), 1); 

for i = 1:length(File_No)
    No = File_No(i);
    
    % FIX 1: Use square brackets for path concatenation
    record_path = ['mit-bih-arrhythmia-database-1.0.0\', num2str(No)]; 
    
    try
        [signal, fs_orig] = rdsamp(record_path, 1);
        
        % Resampling
        ecg_resampled = resample(signal, fs_new, fs_orig);
        
        % Filtering
        filtered_ecg = pan_tompkins_bandpass(ecg_resampled, fs_new);
        
        % Derivative
        b_der = [2 1 0 -1 -2] * (fs_new / 8); 
        a_der = 1;
        % FIX 2: Variable name must match 'filtered_ecg'
        diff_ecg = derivative_func(filtered_ecg, b_der, a_der);
        
        % Squaring
        % FIX 3: Variable name must match 'diff_ecg'
        ecg_squared = diff_ecg .^ 2; 
        ecg_squared = ecg_squared / max(ecg_squared);
        
        % Moving Window integrator
        window_length = 20;
        ecg_mwi = movingWindowIntegrator(ecg_squared, window_length);
        
        % Store in Cell Array (prevents "Array Growth" speed issues)
        tf_ecg{i} = ecg_mwi;
        
        fprintf('Processed Record: %d\n', No);
    catch
        fprintf('Error loading Record: %d (Skipping)\n', No);
    end
end


% --- After the for loop ends ---

% --- Final Batch Analysis ---
num_files = length(tf_ecg);
all_bpm = zeros(num_files, 1); % Pre-allocate for speed


num_files = length(tf_ecg);
stats = struct('TP',0,'FP',0,'FN',0,'Se',0,'Pp',0);
results_log = repmat(stats, num_files, 1);

fs_old = 360; % MIT-BIH original sampling rate
tolerance = 0.150 * fs_new; % 150ms window in our 200Hz scale

for i = 1:num_files
    mwi_signal = tf_ecg{i};
    

    %---Check if data exists before calling findpeaks ---

    if isempty(mwi_signal)
        fprintf('Warning: Record %d is empty. Skipping...\n', File_No(i));
        all_bpm(i) = NaN; 
        continue; % Skip to the next iteration of the loop
    end

    % 1. Find potential peaks
    [pks, locs] = findpeaks(mwi_signal, 'MinPeakDistance', 0.2*fs_new);
    
    % 2. Apply Adaptive Thresholding
    spki = max(pks) * 0.5; 
    npki = mean(pks) * 0.1;
    threshold = npki + 0.25*(spki - npki);
    qrs_indices = [];
    
    for j = 1:length(pks)
        if pks(j) > threshold
            qrs_indices = [qrs_indices; locs(j)];
            spki = 0.125*pks(j) + 0.875*spki; % Update signal level
        else
            npki = 0.125*pks(j) + 0.875*npki; % Update noise level
        end
        threshold = npki + 0.25*(spki - npki);
    end
    
    % 3. Calculate BPM for this record
    if length(qrs_indices) > 1
        rr_intervals = diff(qrs_indices) / fs_new; % Seconds between beats
        all_bpm(i) = 60 / mean(rr_intervals);
    else
        all_bpm(i) = NaN; % Mark as 'Bad Data' if no beats found
    end

    % 2. Get Gold Standard (Physician Annotations)
    record_str = num2str(File_No(i));
    [ann_locs, ann_type] = rdann(record_str, 'atr');
    
    % Filter for QRS beats and Rescale to 200Hz
    valid_beats = ismember(ann_type, {'N','L','R','A','a','J','S','V','F','e','j','E','/'});
    truth_indices = round(ann_locs(valid_beats) * (fs_new / fs_old));
    
    % 3. Calculate TP, FP, FN
    tp_count = 0;
    for k = 1:length(truth_indices)
        if any(abs(qrs_indices - truth_indices(k)) <= tolerance)
            tp_count = tp_count + 1;
        end
    end
    
    fn_count = length(truth_indices) - tp_count;
    fp_count = length(qrs_indices) - tp_count;
    
    % 4. Store Results for this record
    results_log(i).TP = tp_count;
    results_log(i).FP = fp_count;
    results_log(i).FN = fn_count;
    results_log(i).Se = (tp_count / (tp_count + fn_count)) * 100;
    results_log(i).Pp = (tp_count / (tp_count + fp_count)) * 100;
end

% --- Final Average Calculation ---
avg_Se = mean([results_log.Se], 'omitnan');
avg_Pp = mean([results_log.Pp], 'omitnan');
avg_Acc = (sum([results_log.TP]) / (sum([results_log.TP]) + sum([results_log.FP]) + sum([results_log.FN]))) * 100;

fprintf('--- Final Pipeline Results ---\n');
fprintf('Average Sensitivity (Se): %.2f%%\n', avg_Se);
fprintf('Average Positive Predictivity (+P): %.2f%%\n', avg_Pp);
fprintf('Overall Detection Accuracy: %.2f%%\n', avg_Acc);

% Create a summary table
results_table = table(File_No', all_bpm, 'VariableNames', {'Record_ID', 'HeartRate_BPM'});
%disp(results_table);

% Export your results to the EEE 376 lab folder
writetable(results_table, 'B:\EEE 376\Labtest materials\lab-05\ECG_Analysis_Results.xlsx');
fprintf('Results exported successfully!\n');