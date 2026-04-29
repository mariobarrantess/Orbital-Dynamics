function [t_f] = deltatheta2t(theta_0, delta_theta, mu_E, a, e, t_i)

theta_f = theta_0 + delta_theta;

n = sqrt(mu_E/a^3);

%initial state (theta_0=perigee)
E_0 = functionshw1.theta2E(theta_0, e);
M_0 = functionshw1.E2M(E_0, e);

%final state (theta_0+pi/2)
E_f = functionshw1.theta2E(theta_f, e);
M_f = functionshw1.E2M(E_f, e);

t_f = t_i + ((M_f - M_0) /(n*3600)); % [s]->[h]

% t_f = t_i + t_travel(i)./3600; 

end