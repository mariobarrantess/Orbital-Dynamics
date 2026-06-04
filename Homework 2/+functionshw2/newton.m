function xroot = newton(f,df,x,tol,maxiter)

for iter = 1:maxiter

    % Newton step
    x_new = x - f(x)/df(x);

    % absolute and relative error (like in slides)
    err_abs = abs(x_new - x);
    err_rel = err_abs/(1 + abs(x_new));

    % update current iterate
    x = x_new;

    % check convergence
    if (err_abs <= tol) && (err_rel <= tol)
        xroot = x;
        % fprintf('Convergence reached in %d iterations\n',iter);
        return
    end
end

% if we get here, no convergence
xroot = x;
warning('Maximum number of iterations (%d) reached',maxiter);
end