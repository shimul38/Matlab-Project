function y = derivative_func(x, b, a)
    % Use filtfilt to ensure the derivative doesn't shift the peaks
    % b = [2 1 0 -1 -2] * (fs/8) as defined in your main script
    y = filtfilt(b, a, x);
end