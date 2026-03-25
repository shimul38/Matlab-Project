function y = pan_tompkins_bandpass(x, fs)
    % Butterworth Bandpass: 5Hz to 15Hz
    % A 3rd order filter becomes 6th order when using filtfilt
    [b, a] = butter(3, [5 15]/(fs/2), 'bandpass');
    
    % Convert to Second-Order Sections for numerical stability
    [sos, g] = tf2sos(b, a);
    
    % Apply Zero-Phase Filtering
    y = filtfilt(sos, g, x);
end