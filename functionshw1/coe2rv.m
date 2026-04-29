function [r0, v0, p_vector, q_vector, w_vector, h_norm]  = coe2rv(a, e, Omega, i, omega, theta, mu_E)

% INPUTS
%  a = semi-major axis [km]
%  e = eccentricity
%  Omega = right ascension of the ascending node (RAAN) [rad]
%  i = inclination [rad] 
%  omega = argument of perigee [rad]
%  theta = true anomaly [rad]
%  mu_E = gravitational parameter of Earth [km^3/s^2]
%
% OUTPUTS
%  r0 = position vector wrt S_0 [km]
%  v0 = velocity vector wrt S_0 [km/s]

    p = a*(1-e^2); % [km]
    
    h_norm = sqrt(p*mu_E); % [km^2/s]

    r_norm = p/(1+e*cos(theta)); % [km]
    r_vect = r_norm.*[cos(theta); sin(theta); 0];       % [km] position in p directions

    v_norm = mu_E/h_norm;                          
    v_vect = v_norm.*[-sin(theta); e+cos(theta); 0];    % [km/s] velocity in p directions

    % Rotation matrix of argument of perigee - 3rd axis rotation
    % R_omega = [cos(omega) sin(omega) 0; 
    %           -sin(omega) cos(omega) 0; 
    %            0          0          1];
    % 
    % % Rotation matrix of inclination - 1st axis rotation
    % R_i = [1  0      0; 
    %        0  cos(i) sin(i); 
    %        0 -sin(i) cos(i)];
    % 
    % % Rotation matrix of RAAN - 3rd axis rotation
    % R_Omega = [cos(Omega) sin(Omega) 0; 
    %           -sin(Omega) cos(Omega) 0; 
    %            0          0          1];
    % 
    % %Rotation matrix from orbit plane to S0
    % R_P20 = R_omega*R_i*R_Omega;
    % 
    % % r0 = (r_vect*R_P20)'; % position vector in S0(ECI) [km]
    % % v0 = (v_vect*R_P20)'; % velocity vector in S0(ECI) [km/s]
    % r0 = R_P20'*r_vect; % position vector in S0(ECI) [km]
    % v0 = R_P20'*v_vect; % velocity vector in S0(ECI) [km/s]
    
    % Rotating from ECI to Orbital plane

    % Rotation matrix of argument of perigee - 3rd axis rotation
    R_omega = [cos(omega) -sin(omega) 0; 
               sin(omega)  cos(omega) 0; 
               0           0          1];

    % Rotation matrix of inclination - 1st axis rotation
    R_i = [1  0       0; 
           0  cos(i) -sin(i); 
           0  sin(i)  cos(i)];

    % Rotation matrix of RAAN - 3rd axis rotation
    R_Omega = [cos(Omega) -sin(Omega) 0; 
               sin(Omega)  cos(Omega) 0; 
               0           0          1];

    % Rotation matrix from orbit plane to S0
    % R_ECI2ORB = R_omega*R_i*R_Omega;
    R_ORB2ECI = R_Omega*R_i*R_omega;

    % r0 = (r_vect*R_P20)'; % position vector in S0(ECI) [km]
    % v0 = (v_vect*R_P20)'; % velocity vector in S0(ECI) [km/s]
    r0 = R_ORB2ECI*r_vect; % position vector in S0(ECI) [km]
    v0 = R_ORB2ECI*v_vect; % velocity vector in S0(ECI) [km/s]

    % h_vect = cross(r0,v0);
    % h_norm = norm(h_vect);
    % e_vect = cross(v0,h_vect)/mu_E - (r0/norm(r0));
    % 
    % p_vector = e_vect/(norm(e_vect));
    % w_vector = h_vect/norm(h_vect);
    % q_vector = cross(w_vector, p_vector);

    p_vector = R_ORB2ECI(:, 1);  % p̂ — along periapsis
    q_vector = R_ORB2ECI(:, 2);  % q̂ — 90° ahead in orbit
    w_vector = R_ORB2ECI(:, 3);  % ŵ — angular momentum direction

end