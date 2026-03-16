%%This is the derivative function for our ecg signal
function  ecg_deriv = derivative_func(ecg_signal, num, den)
% Apply the derivative to the bandpassed signal
    ecg_deriv = filter(num, den, ecg_signal);

% Normalize for visualization
    ecg_deriv = ecg_deriv / max(abs(ecg_signal));
end