function ap = J2_acceleration(r)
% Computes perturbation acceleration due to J2 oblateness effect
%
% INPUTS
%   r  - [km]     - position vector (3×1)
%
% OUTPUTS
%   ap - [km/s²]  - J2 perturbation acceleration (3×1)

muE = 398600;       % [km³/s²] Earth gravitational parameter
J2  = 0.00108263;   % [-]      J2 coefficient
RE  = 6371;         % [km]     Earth radius

k0  = [0; 0; 1];            % [-]  Z-axis unit vector
rn  = norm(r);              % [km] position magnitude
e_r  = r / rn;               % [-]  position unit vector
z   = dot(r, k0);           % [km] Z-component of position

ap = 1.5 * J2 * muE * RE^2 / rn^4 * ((5*z^2/rn^2 - 1)*e_r - 2*z/rn*k0);  % [km/s²]
end
