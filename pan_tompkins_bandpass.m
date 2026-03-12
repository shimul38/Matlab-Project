function filtered_ecg = pan_tompkins_bandpass(x, fs)
%verifiying the sampling frequency
    if fs ~= 200
        warning("Pan-tompkins coefficients are designed for 200 Hz");
    end
    %%designing a low-pass filter 
    b_lp = [1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1]; %%coefficient of numerators: H(z)=1 - 2Z^-6 + Z^-12 
    a_lp = [1, -2, 1]; %% coefficients of denominators of transfer function H(z)

    y_lp = filter(b_lp, a_lp, x);

    %%designing a High-pass filter

    b_hp = [-1 zeros(1, 15), 32, zeros(1, 15), 1];
    a_hp = [1, 1];

    y_hp = filter(b_hp, a_hp, y_lp);

    %%Gain normalization
    filtered_ecg = y_hp / (36 * 32);
end