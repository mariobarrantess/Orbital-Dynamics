function [a, e, Omega, i, omega, theta1, theta2, e_unit, kk] = Lambert_conic(rr1, rr2, eT, k)

kk = k/norm(k);

% Radial distances
r1 = norm(rr1);
r2 = norm(rr2);

% Chord vector
cc = rr2-rr1;
c = norm(cc);

% Parameters of the fundamental ellipse
eF = (r1-r2)/c;
aF = (r1+r2)/2;
pF = aF*(1-eF^2);

% Unit vectors in Lambert basis
ii = cc/c;
jj = cross(kk, ii);

% Eccentricity vector
eF_vect = eF*ii;
eT_vect = eT*jj;
ee = eF_vect + eT_vect;
e_unit = ee/norm(ee);
e = norm(ee);

% True anomaly increment
sin_Dtheta = dot(kk, cross(rr1, rr2)) / (r1*r2);
cos_Dtheta = dot(rr1, rr2) / (r1*r2);
Dtheta = atan2(sin_Dtheta, cos_Dtheta);
if Dtheta < 0
    Dtheta = Dtheta + 2*pi;
end

% Semi-latus rectum
p = pF - eT*r1*r2/c*sin(Dtheta);

% Semi-major axis
a = p/(1-e^2);

% Inclination
i = acos(kk(3));

% Line of nodes
N = cross([0;0;1], kk);
if norm(N) < 1e-10
    N = [1;0;0];
end
N = N/norm(N);

% RAAN
Omega = atan2(N(2), N(1));
if Omega < 0
    Omega = Omega + 2*pi;
end

% Argument of periapsis (omega) — robust quadrant
omega = atan2(dot(cross(N, e_unit), kk), dot(N, e_unit));
if omega < 0
    omega = omega + 2*pi;
end

% True anomalies — robust quadrant using atan2
q_hat = cross(kk, e_unit);  % 90 deg ahead in orbit plane

theta1 = atan2(dot(rr1, q_hat), dot(rr1, e_unit));
theta2 = atan2(dot(rr2, q_hat), dot(rr2, e_unit));

% Enforce sweep direction consistency with k
Dtheta_check = mod(theta2 - theta1, 2*pi);

% We detect which case we are in by checking sin_Dtheta sign
if sin_Dtheta >= 0
    % Short arc: sweep should be in (0, pi)
    % If Dtheta_check > pi, the atan2 wrapped wrong → subtract
    if Dtheta_check > pi
        theta2 = theta2 - 2*pi;  % pull it back so sweep < pi
        if theta2 < 0
            theta2 = theta2 + 2*pi;
        end
    end
else
    % Long arc: sweep should be in (pi, 2*pi)
    if Dtheta_check < pi
        theta2 = theta2 + 2*pi;  % push it forward so sweep > pi
    end
end

end