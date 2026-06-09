function [dX] = evolution(t, X, bc, flag_drag, flag_J2)
%Unzip

rr = X(1:3);            % [km] position vector
vv = X(4:6);            % [km/s] velocity vector
r  = norm(rr);          % [km] position magnitude

muE = 3.98600E5;  % muEarth [km3/s2]
RE = 6371; % Earth radius [km]

% Drag perturbation
if flag_drag
    f_drag = functionshw3.drag_acceleration(r - RE, vv, bc);  % [km/s^2]
else
    f_drag = 0;
end

% J2 perturbation
if flag_J2
    f_j2 = functionshw3.J2_acceleration(rr);   % [km/s^2]
else
    f_j2 = 0;  % [km/s^2]
end

% Compute RHS
dX(1:3,1) = vv; % [km/s]
dX(4:6,1) = -muE/r^3 * rr + f_j2 + f_drag; % [km/s^2]
end
