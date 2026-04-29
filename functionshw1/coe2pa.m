function [r_p, r_a] = coe2pa(a, e)

% INPUTS
% a - semimajor axis of orbit [km]
% e - eccentricity of orbit [-]

% OUTPUTS
% r_p - Radius of pericenter [km]
% r_a - Radius of apocenter [km]

r_p = a*(1-e);
r_a = a*(1+e);

end