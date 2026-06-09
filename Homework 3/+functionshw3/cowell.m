function [r,v] = cowell(r0, v0, bc, t, flag_drag, flag_J2)
%
% INPUTS
%   - r0 [km] initial position vector
%   - v0 [km/s] initial velocity vector
%   - bc [log(kg/m^2)] ballistic coefficient
%   - t [s] time
%   - flag_drag (true or false) drag perturbations
%   - flag_J2 (true or false) J2 perturbations
%
% OUTPUTS
%   - r [km] position vector
%   - v [km/s] velocity vector
% 

X0 = [r0; v0]; % [km, km/s] initial position and velocity vector

options = odeset('relTol', 1E-6, 'AbsTol', 1E-9);

odefun = @(t,X) functionshw3.evolution(t, X, bc, flag_drag, flag_J2);

[~, X] = ode45(odefun, t, X0, options);

r = [X(:,1:3)]';            % [km] position vector
v = [X(:,4:6)]';            % [km/s] velocity vector
end