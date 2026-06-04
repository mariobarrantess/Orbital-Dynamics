function [a, e, Omega, i, omega, theta1, theta2, eT_real] = Lambert_solve(mu, r1, r2, Dt, k)

% INPUTS
% mu - Gravitational parameter of central body [km3/m2]
% r1 - Position vector at time t1 [km] 
% r2 - Position vector at time t2 [km]
% Dt - Time between position 1 and position 2 [s]
% k - Positive direction of the orbital plane normal

% OUTPUTS
% a - semi mayor axis [km]
% e - eccentricity [km]
% Omega - Right ascension of the ascending node = RAAN [rad]
% i - Inclination [rad]
% omega - Argument of the periapsis [rad]
% theta1 - True anomaly at t1 [rad]
% theta2 - True anomaly at t2 [rad]
% eT_real - transverse eccentricity component

% Position
R1 = norm(r1);
R2 = norm(r2);

% Chord vector
cc = r2 - r1;
c = norm(cc);

% Vector basis (unitary)
kk = k/norm(k);
ii = cc/c;
jj = cross(kk,ii);

% Eccentricity component
eF = (R1 - R2)/c;

aF = (R1+R2)/2; % semimajor axis of fundamental ellipse [km]
pF = aF*(1-eF^2); % semi-latum rectum of fundamental ellipse [km]

eTP = sqrt(1 - eF^2); % parabolic limit

% Search interval (elliptic only)
eT_lo = -eTP + 1e-6;
eT_hi =  eTP - 1e-6;

% Evaluate bounds
[a_lo, e_lo, ~, ~, ~, th1_lo, th2_lo] = functionshw2.Lambert_conic(r1, r2, eT_lo, k);
tof_lo = functionshw2.tof(mu, a_lo, e_lo, th1_lo, th2_lo);
err_lo = tof_lo - Dt;

[a_hi, e_hi, ~, ~, ~, th1_hi, th2_hi] = functionshw2.Lambert_conic(r1, r2, eT_hi, k);
tof_hi = functionshw2.tof(mu, a_hi, e_hi, th1_hi, th2_hi);
err_hi = tof_hi - Dt;

% If no sign change, no elliptic solution for this Dt
if ~isfinite(err_lo) || ~isfinite(err_hi) || err_lo*err_hi > 0
    a = NaN; e = NaN; Omega = NaN; i = NaN; omega = NaN; theta1 = NaN; theta2 = NaN; eT_real = NaN;
    return
end

% Bisection
tol = 1e-5;
max_iter = 200;

for iter = 1:max_iter
    eT_real = (eT_lo + eT_hi)/2;

    [a_mid, e_mid, ~, ~, ~, th1_mid, th2_mid] = functionshw2.Lambert_conic(r1, r2, eT_real, k);
    tof_mid = functionshw2.tof(mu, a_mid, e_mid, th1_mid, th2_mid);
    err_mid = tof_mid - Dt;

    if abs(err_mid) < tol
        break;
    end

    if err_lo * err_mid < 0
        eT_hi = eT_real;
        err_hi = err_mid;
    else
        eT_lo = eT_real;
        err_lo = err_mid;
    end
end

% Final solution
[a, e, Omega, i, omega, theta1, theta2] = functionshw2.Lambert_conic(r1, r2, eT_real, k);

% Diagnostics (optional)
Dtheta = mod(theta2 - theta1, 2*pi);
eTD = (pF*c)/(R1*R2*sin(Dtheta));
if eT_real < eTP && eT_real > -eTP
    % elliptic
elseif eT_real == -eTP || eT_real == eTP
    fprintf('\nOrbit is a parabola only allowing short arc solution\n');
elseif eT_real < -eTP
    fprintf('\nHyperbolic short arcs\n')
elseif eT_real > eTP && eT_real < eTD
    fprintf('\nHyperbolic long arcs\n')
elseif eT_real == eTD
    fprintf('\np=0 and F=F´--> Pair of lines that intersect at the origin\n')
elseif eT_real > eTD
    fprintf('\nInadmissible hyperbolic solutions\n')
end

end