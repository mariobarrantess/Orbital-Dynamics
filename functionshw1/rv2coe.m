function [a, e_norm, Omega, i, omega, theta_0] = rv2coe(r0,v0,mu_E)
% GOAL
% take vectors r0, v0 at a certain moment, as seen from Earth Centered Inertial reference frame
% and express the classical orbital parameters of that specific orbit 

% INPUTS
% r0 - position vector (as seen from ECI ref. frame) [km]
% v0 - velocity vector (as seen from ECI ref. frame) [km/s]
% mu_E - gravitational parameter of the Earth [km3/s2]
% 
% OUTPUTS
% a - semi mayor axis [km]
% e - eccentricity [km]
% Omega - Right ascension of the ascending node = RAAN [rad]
% i - Inclination [rad]
% omega - Argument of the periapsis [rad]
% theta_0 - True anomaly (angle from pericenter to satellite)

v0_norm = norm(v0);
r0_norm = norm(r0);

%from vis viva equation we find a:

a = (-mu_E/2)*(1/((v0_norm^2 /2)-(mu_E/r0_norm))); % semimajor axis [km] (from vis-viva equation)


h_vect = cross(r0,v0);
h_norm = norm(h_vect);
h_normalised = h_vect/h_norm;

e_vect = (cross(v0,h_vect)/mu_E - (r0/r0_norm));
e_norm = norm(e_vect);
e_normalised = e_vect./e_norm;
% e = sqrt(1-p/a);

p_vect = cross(h_normalised, e_normalised); %/ (h_norm * e_norm);
p = h_norm^2/mu_E;

i_vect = [1; 0; 0];
j_vect = [0; 1; 0];
k_vect = [0; 0; 1];

% the third component of the angular momentum gives the inclination of the orbit
i = acos(dot(h_vect,k_vect) / h_norm); % Inclination [rad]

n_vect = cross(k_vect,h_vect);
n_norm = norm(n_vect);
n_normalised = n_vect./n_norm;

Omega = mod(atan2(n_normalised(2), n_normalised(1)), 2*pi);

% define a vector g to make a base with h and n:
g_vect = cross(h_normalised,n_normalised);

% omega = acos(dot(n_vect, e_vect) / (n_norm * e));
% omega = mod(atan2(dot(e_normalised,n_normalised), dot(e_normalised, g_vect)), 2 * pi);
omega = mod(atan2(dot(e_normalised, g_vect), dot(e_normalised, n_normalised)), 2*pi);

% theta_0 = mod(atan2(dot(r0, e_normalised),dot(r0, p_vect)), 2*pi);
theta_0 = mod(atan2(dot(r0, p_vect), dot(r0, e_normalised)), 2*pi);

p_vector = e_vect/e_norm;
w_vector = h_vect/h_norm;
q_vector = cross(w_vector, p_vector);

end