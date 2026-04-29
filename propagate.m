% Propagate true anomaly from t0 to final time t0+Dt
function theta = propagate(mu, a, e, theta0, Dt)
% INPUTS
%   mu - gravitational parameter of central body
%   a - Semi-major axis of orbit [km]
%   e - Eccentricity of orbit [-]
%   theta - True anomaly of starting point [rad]
%   Dt - time increment [s] --> tf=t0+Dt
% OUTPUTS
%   theta - True anomaly of final point [rad]

n  = sqrt(mu/a^3);                      % mean motion [rad/s]
E0 = functionshw2.theta2E(theta0, e);   % initial eccentric anomaly
M0 = functionshw2.E2M(E0, e);           % initial mean anomaly
M  = M0 + n*Dt;                         % mean anomaly at t
M  = mod(M, 2*pi);                      % keep result in [0, 2pi]
E  = functionshw2.M2E(M, e);            % eccentric anomaly [rad]
theta = functionshw2.E2theta(E, e);     % convert back to true anomaly

end