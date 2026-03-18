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

% 1. Plot the FIRST record processed
figure('Name', 'First Record Verification');
first_data = tf_ecg{1}; 
fs_new = 200;
t1 = (0:length(first_data)-1) / fs_new;
zoom1 = (t1 >= 10) & (t1 <= 15); % 5-second window

subplot(2,1,1);
plot(t1(zoom1), first_data(zoom1), 'r', 'LineWidth', 1.2);
title(['Record ', num2str(File_No(1)), ': Integrated (MWI)']);
grid on; ylabel('Energy');

% 2. Plot the LAST record processed
figure('Name', 'Last Record Verification');
last_data = tf_ecg{end}; 
t2 = (0:length(last_data)-1) / fs_new;
zoom2 = (t2 >= 10) & (t2 <= 15);

subplot(2,1,2);
plot(t2(zoom2), last_data(zoom2), 'r', 'LineWidth', 1.2);
title(['Record ', num2str(File_No(end)), ': Integrated (MWI)']);
grid on; ylabel('Energy'); xlabel('Time (seconds)');