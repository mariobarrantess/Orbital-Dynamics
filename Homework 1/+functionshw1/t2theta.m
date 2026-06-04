function [theta_f] = t2theta(theta_0, mu_E, a, e, t2)

E_0 = functionshw1.theta2E (theta_0 , e); 
M_0 = functionshw1.E2M (E_0, e);

n = sqrt(mu_E/a^3);

Dt = t2*3600; % [s]

DM = Dt*n;

M_f= M_0 + DM;
M_f = mod(M_f, 2*pi);

E_f = functionshw1.M2E (M_f, e);

theta_f = functionshw1.E2theta (E_f, e);

theta_f = mod(theta_f, 2*pi);

end