function Phi = STM(t)
    global A
    % These could be arguments but we are forced to only have t as an
    % argument...
    mu_star = 3.0034e-6;
    x_L1 = fzero(f_L1, 0.99);
    Phi = expm(A * t);
end