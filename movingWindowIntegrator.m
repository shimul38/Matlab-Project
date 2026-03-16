function y = movingWindowIntegrator(ecg_signal, window_length)
    N = window_length;
    b_init = ones(1, N)/N; %% numerator coefficients 
    a_init = 1;         %denominator coefficients

    y = filter(b_init, a_init, ecg_signal)
end