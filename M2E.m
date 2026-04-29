function E = M2E (M, e)
% on M2E implement newton at the end

f  = @(E) E - e.*sin(E) - M;
df = @(E) 1 - e.*cos(E); % df/dE

% Initial guess: E0 = M
E0 = M;

tol = 1e-8;
maxiter = 100;

% Call your generic Newton solver
E = functionshw2.newton(f, df, E0, tol, maxiter);
end


