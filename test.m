record_path = 'mit-bih-arrhythmia-database-1.0.0\100';
[signal, fs_orig] = rdsamp(record_path, 1);

figure();
fs_new = 200;
% Resampling
ecg_resampled = resample(signal, fs_new, fs_orig);
dt = 1/fs_new;
t = 0:dt:2;
subplot(2,2, 1);
stem(t,ecg_resampled(1:length(t)));
xlabel("Time");
ylabel("Amplitude");
% Filtering
filtered_ecg = pan_tompkins_bandpass(ecg_resampled, fs_new);
subplot(2,2,2);
plot(t, filtered_ecg(1:length(t)))
xlabel("Time");
ylabel("Amplitude");
title('filtered ecg');
% Derivative
b_der = [2 1 0 -1 -2] * (fs_new / 8); 
a_der = 1;
% FIX 2: Variable name must match 'filtered_ecg'
diff_ecg = derivative_func(filtered_ecg, b_der, a_der);

%%Normalization

diff_ecg = diff_ecg / max(abs(diff_ecg));

subplot(2,2,3);
plot(t,diff_ecg(1:length(t)));
xlabel("Time");
ylabel("Amplitude");
title("derivative ecg");
% Squaring
% FIX 3: Variable name must match 'diff_ecg'
ecg_squared = diff_ecg .^ 2; 
        
% Moving Window integrator
window_length = 10;
ecg_mwi = movingWindowIntegrator(ecg_squared, window_length);
subplot(2,2,4);

plot(t, ecg_mwi(1:length(t)));
xlabel("Time");
ylabel("Amplitude");