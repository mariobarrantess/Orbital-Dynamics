function [r0,v0,lambda,R_tot] = SEZ2ECI(r1,v1,rS,phi,t)
% INPUTS
%   r1      [km]     position in SEZ frame
%   v1      [km/s]   velocity in SEZ frame
%   rS      [km]     station distance from Earth center
%   phi     [rad]    observer latitude
%   t       [h]      time since Oxz plane crossing
%
% OUTPUTS
%   r0      [km]     position in ECI frame
%   v0      [km/s]   velocity in ECI frame

omega_earth = 72.921e-6; % [rad/s]
n_Earth = [0; 0; omega_earth];  % [rad/s] rotation vector

lambda = omega_earth * t;  % [rad] longitud S1 w/ respect t0 S0 

% Rotation about Z axis of ECI (longitude) %revise rotation matrix
R_z = [cos(lambda)  -sin(lambda)   0;
       sin(lambda)   cos(lambda)   0;
       0             0             1];

R_y = [sin(phi)   0    cos(phi); %
       0          1    0;
       -cos(phi)  0    sin(phi)];

% Total rotation matrix.
R_tot = R_z * R_y;

% Position transformation
rStaSpa_ECI = R_tot * r1;  % [km] spacecraft position relative to station in ECI
rSta0_ECI = rS * [cos(phi) * cos(lambda); cos(phi) * sin(lambda); sin(phi)];  % [km] station position in ECI

% Output position vector
r0 = rStaSpa_ECI + rSta0_ECI; % [km] spacecraft position in ECI frame

% Velocity components
v1_0 = R_tot * v1;  % [km/s] relative velocity in ECI

rStaSpa_rot = cross(n_Earth, rStaSpa_ECI);  % [km/s] velocity due to Earth rotation (relative position)
rSta0_rot = cross(n_Earth, rSta0_ECI);  % [km/s] velocity due to Earth rotation (station)

% Output velocity vector
v0 = v1_0 + rStaSpa_rot + rSta0_rot;  % [km/s] spacecraft velocity in ECI frame

end