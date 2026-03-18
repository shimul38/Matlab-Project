function y = movingWindowIntegrator(x, N)
    % x is the squared signal
    % N is the window width (usually 30 for 200Hz)
    
    b = ones(1, N) / N; % The '1/N' is crucial to prevent the "clipping" you see
    a = 1;
    
    y = filter(b, a, x);
end