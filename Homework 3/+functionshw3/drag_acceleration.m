function ap = drag_acceleration(zeta, v, bc)
% Computes perturbation acceleration due to atmospheric drag
%
% INPUTS
%  zeta - [km] altitude above Earth's surface
%  v    - [km/s] velocity vector (3×1)
%  bc   - [kg/m^2] ballistic coefficient
%
% OUTPUTS
%  ap   - [km/s^2] drag perturbation acceleration (3×1)

zeta0 = 100;    % [km]    reference altitude
rho0  = 1e-8;   % [kg/m^3] reference density at zeta0
H     = 50;     % [km]    scale height

if zeta < zeta0
    ap = nan(3, 1);
    return
end

rho = rho0*exp(-(zeta - zeta0) / H);           % [kg/m^3] Earth's atmosphere density
ap = -(rho * 1e3) / (2 * bc) * norm(v) * v;   % [km/s^2]

end