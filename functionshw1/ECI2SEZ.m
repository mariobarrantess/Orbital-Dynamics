function [r1,v1,v1_0] = ECI2SEZ(r0,v0,rS,phi,t)
% INPUTS
%   r0      [km]     position in ECI frame
%   v0      [km/s]   velocity in ECI frame
%   rS      [km]     station distance from Earth center
%   phi     [rad]    observer latitude
%   t       [h]      time since Oxz plane crossing
%
% OUTPUTS
%   r1      [km]     position in SEZ frame
%   v1      [km/s]   velocity in SEZ frame
%   v1_0    [km/s]   relative velocity in ECI frame

omega_earth = 72.921e-6; % [rad/s]
n_Earth = [0; 0; omega_earth];  % [rad/s] rotation vector

lambda = omega_earth * t;  % [rad] longitud S1 w/ respect t0 S0 

R_z = [cos(lambda) -sin(lambda) 0;
        sin(lambda)  cos(lambda) 0;
        0            0           1];                              % Rotation about Z axis of ECI
R_y = [sin(phi) 0  cos(phi);
        0        1  0;
       -cos(phi) 0  sin(phi)];               % Rotation about Y axis of ECI

R_SEZ2ECI = R_z*R_y;
R_ECI2SEZ = R_SEZ2ECI'; % ECI2SEZ rotation matrix is the transpose of the SEZ2ECI 

rSta0_ECI = (rS)*[cos(phi)*cos(lambda) cos(phi)*sin(lambda) sin(phi)]';  % [km] station vector in ECI ref.frame
rStaSpa_ECI = r0 - rSta0_ECI;                         % [km] spacecraft position vector in ECI RF 
r1 = R_ECI2SEZ * rStaSpa_ECI;     % absolute position of spacecraft wrt to S1;                

% v_rel = v0 - (omega x r0)
v1_0 = v0 - cross(n_Earth, r0);        

v1 = R_ECI2SEZ * v1_0;  % absolute velocity of spacecraft wrt to S1;

end