% record_path should point to the folder 'mitdbdir'
record_path = 'mit-bih-arrhythmia-database-1.0.0\100'; 

% rdsamp reads the signal data (.dat) and header info (.hea)
% [signal, fs_orig, time] = rdsamp(record_name, channel_number, num_samples)
[signal, fs_orig] = rdsamp(record_path, 1); % Importing Channel 1 (usually Lead MLII)

t = (0:length(signal)-1)/fs_orig;
%%plot(t, signal);

%Resampling part

fs_new = 200; %resampling frequency

ecg_resampled = resample(signal, fs_new, fs_orig);

t_new = (0: length(ecg_resampled) - 1)/fs_new;
plot(t_new(1:1000), ecg_resampled(1:1000));


%%Filtering part
filtered_ecg = pan_tompkins_bandpass(ecg_resampled, fs_new);

figure(2);
plot(t_new(1:1000), filtered_ecg(1:1000));
title("Filted ECG Signal");
xlabel("Time(sec)");
ylabel("Voltage");



