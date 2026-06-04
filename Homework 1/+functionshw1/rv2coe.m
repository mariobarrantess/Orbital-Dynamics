function [a, e, Omega, i, omega, theta_0] = rv2coe(r0, v0, mu_E)

% INPUTS
% r0 - position vector (as seen from ECI ref. frame) [km]
% v0 - velocity vector (as seen from ECI ref. frame) [km/s]
% mu - gravitational parameter of the Earth [km3/s2]
% 
% OUTPUTS
% a       - semi mayor axis [km]
% e       - eccentricity [km]
% Omega   - Right ascension of the ascending node = RAAN [rad]
% i       - Inclination [rad]
% omega   - Argument of the periapsis [rad]
% theta_0 - True anomaly (angle from pericenter to satellite)

r0_norm = norm(r0);
v0_norm = norm(v0);

% Semi-major axis (vis-viva)
energy = v0_norm^2/2 - mu_E/r0_norm;
a = -mu_E/(2*energy);

% Angular momentum
h_vect = cross(r0, v0);
h_norm = norm(h_vect);
h_hat = h_vect / h_norm;

% Eccentricity
e_vect = cross(v0, h_vect)/mu_E - r0/r0_norm;
e = norm(e_vect);

% Inclination 
% from third component of the angular momentum gives
i = acos( h_vect(3)/h_norm );

% Node vector
k = [0;0;1];
n_vect = cross(k, h_vect);
n_norm = norm(n_vect);

% RAAN
if n_norm > 1e-12
    n_hat = n_vect / n_norm;
    Omega = mod(atan2(n_hat(2), n_hat(1)), 2*pi);
else
    Omega = 0; % equatorial orbit
end

% Argument of periapsis and true anomaly
if e > 1e-12
    e_hat = e_vect / e;

    if n_norm > 1e-12
        g_hat = cross(h_hat, n_hat);
        omega = mod(atan2(dot(e_hat, g_hat), dot(e_hat, n_hat)), 2*pi);
    else
        omega = 0; % equatorial - undefined
    end
    % define vector q to make a base with h and e:
    q_hat = cross(h_hat, e_hat);
    theta_0 = mod(atan2(dot(r0, q_hat), dot(r0, e_hat)), 2*pi);

else
    % Circular orbit - omega undefined, using true anomaly
    omega = 0;
    if n_norm > 1e-12
        theta_0 = mod(atan2(dot(r0, cross(h_hat,n_hat)), dot(r0, n_hat)), 2*pi);
    else
        % Circular equatorial - only true anomaly
        theta_0 = mod(atan2(r0(2), r0(1)), 2*pi);
    end
end

end