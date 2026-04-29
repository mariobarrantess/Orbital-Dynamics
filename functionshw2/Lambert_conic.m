% % function [a, e, Omega, i, omega, theta1, theta2] = Lambert_conic(r1, r2, eT, k, tm)
% 
% % ii.- [a,e,Omega,i,omega,theta1,theta2] = Lambert conic(r1,r2,eT,k) 
% % to compute the COE of the transfer trajectory (a [km] and angles [rad]) between position
% % vectors r1 and r2 (in km), as a function of the transverse eccentricity component (eT).
% % Vector k defines the positive direction of the orbital plane normal and therefore selects
% % between short and long arc solutions.
% 
% % INPUTS
% % r1 - Position vector at time t1 [km] 
% % r2 - Position vector at time t2 [km]
% % eT - Transverse eccentricity component
% % k - Positive direction of the orbital plane normal (used to select between short & long arc solutions)
% 
% % OUTPUTS
% % a - semi mayor axis [km]
% % e - eccentricity [km]
% % Omega - Right ascension of the ascending node = RAAN [rad]
% % i - Inclination [rad]
% % omega - Argument of the periapsis [rad]
% % theta1 - True anomaly at t1 [rad]
% % theta2 - True anomaly at t2 [rad]
% 
% %%% version 1 - not 100% correct %%%
% % % %%%%%%%%%%%%%%%%% FROM MAIN %%%%%%%%%%%%%%%%%%%%%
% % % r1 = [15945.34;5000;3000]; %km
% % % r2 = [12214.83899;2000;10249.64731]; %km
% % % eT = 0.7; %
% % % k = cross(r2,r1); % Normal vector to orbit plane
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % 
% % % tm = +1; % short arc solution - deltatheta<pi
% % % tm = -1; % long arc solution - deltatheta>pi
% % % kk = (tm*k)/norm(tm*k);
% % 
% % kk = k/norm(k);     % normalized k vector
% % 
% % R1 = norm(r1);  % radial distance at t1 [km]
% % R2 = norm(r2);  % radial distance at t2 [km]
% % 
% % c_vect = r2 - r1; % chord vector between r1 and r2 [km]
% % c = norm(c_vect); % chord length between r1 and r2 [km]
% % % ic = c_vect./c; % normalised chord vector
% % 
% % ii = c_vect/c;
% % jj = cross(kk,ii);
% % 
% % eF = (R1 - R2)/norm(c_vect); % Eccentricity component along the chord
% % eF_vect = eF*ii; % eF vector
% % e_vect = (eF*ii + eT*jj);% /norm(eF*ii + eT*jj);
% % e = sqrt(eF^2 + eT^2); % e = norm(e_vect);
% % e_unit = e_vect/norm(e_vect);
% % 
% % a = (R1+R2)/2; % Semimajor axis of fundamental ellipse [km]
% % p = a/(1-norm(e_vect)^2); % Semilatus rectum of fundamental ellipse [km]
% % 
% % 
% % % COE
% % i = acos(kk(3));                       % inclination [rad]
% % Omega = atan2(kk(2), kk(1));           % RAAN [rad]
% % omega = atan2(eF_vect(3), eF_vect(1)); % argument of periapsis [rad]
% % 
% % % only valid if orbit is equatorial (i=0)
% % % theta1 = atan2(r1(2), r1(1));          % true anomaly at t1 [rad]
% % % theta2 = atan2(r2(2), r2(1));          % true anomaly at t2 [rad]
% % 
% % theta1 = acos(max(-1, min(1, dot(e_unit, r1/R1))));
% % if dot(kk, cross(e_unit, r1)) < 0 
% %     theta1 = mod(-theta1, 2*pi); %if r1 is in the lower half of orbit, reflect it into [pi, 2pi]
% % end
% % 
% % theta2 = acos(max(-1, min(1, dot(e_unit, r2/R2))));
% % if dot(kk, cross(e_unit, r2)) < 0
% %     theta2 = mod(-theta2, 2*pi); %if r1 is in the lower half of orbit, reflect it into [pi, 2pi]
% % end
% % 
% % Dtheta = theta2 - theta1;    % [rad] difference in true anomaly
% % if tm < 0
% %   Dtheta = 2*pi - Dtheta;
% % end
% % 
% % fprintf('\ne_unit = [%.6f, %.6f, %.6f]\n', e_unit(1), e_unit(2), e_unit(3));
% % fprintf('kk     = [%.6f, %.6f, %.6f]\n', kk(1), kk(2), kk(3));
% % fprintf('dot(e_unit, kk) = %.6f  (should be ~0)\n', dot(e_unit, kk));
% % end

% function [a, e, Omega, i, omega, theta1, theta2] = Lambert_conic(rr1, rr2, eT, k, tm)
function [a, e, Omega, i, omega, theta1, theta2, e_unit, kk] = Lambert_conic(rr1, rr2, eT, k)

kk = k/norm(k);

% Radial distances
r1 = norm (rr1);    % [km] normalized position vector 1
r2 = norm (rr2);    % [km] normalized position vector 2

% Chord vector
cc = rr2-rr1;       % [km] chord vector
c = norm(cc);       % [km] normalized chord vector

% Parameters of the fundamental ellipse
eF = (r1-r2)/c;     % [-] projection of the eccentricity vector along the chord vector
aF = (r1+r2)/2;     % [km] semi-major axis of the fundamental ellipse
pF = aF*(1-eF^2);   % [km] semi-latus rectum of the fundamental ellipse

% Compute the rest of unit vectors used in Lambert's problem
ii = cc/c;
jj = cross(kk, ii);

% Eccentricity
eF_vect = eF*ii;
eT_vect = eT*jj;
ee = eF_vect + eT_vect; % Eccentricity vector
e_unit = ee/norm(ee); % Unitary eccentricity vector
e = norm(ee);       % Eccentricity value

% Compute increment of true anomaly
tm = dot(cross(rr1,rr2)/norm(cross(rr1,rr2)), kk);

Dtheta = acos(dot(rr1,rr2)/(r1*r2)) + 2*pi*tm;   % [rad] difference in true anomaly
% from class notes
%   m = 0 - short-arc solution - m > 0 waiting orbits in direction of k
%   m = -1 - long-arc solution - m <-1 waiting orbits in direction of -k

if tm<0 % long-arc solution
    Dtheta = 2*pi - Dtheta; % [rad] difference in true anomaly
end

% Semi-latus rectum
p = pF - eT*r1*r2/c*sin(Dtheta);  % [km] semi-latus rectum of the generic conic

% Semi-major axis
a = p/(1-e^2); % [km] 

% COE
i = acos(kk(3));  % inclination [rad]

% only valid if orbit is equatorial (i=0)
% theta1 = atan2(r1(2), r1(1));    % true anomaly at t1 [rad]
% theta2 = atan2(r2(2), r2(1));    % true anomaly at t2 [rad]

theta1 = acos(max(-1, min(1, dot(e_unit, rr1/r1))));
if dot(kk, cross(e_unit, rr1)) < 0 
    theta1 = mod(-theta1, 2*pi); %if r1 is in the lower half of orbit, reflect it into [pi, 2pi]
end

theta2 = acos(max(-1, min(1, dot(e_unit, rr2/r2))));
if dot(kk, cross(e_unit, rr2)) < 0
    theta2 = mod(-theta2, 2*pi); %if r1 is in the lower half of orbit, reflect it into [pi, 2pi]
end

% revisar si esto etá bien
N = cross([0;0;1], kk); % line of nodes
if norm(N) < 1e-10
    N = [1;0;0];
end

N = N/norm(N); % normalised vector

Omega = atan2(N(2), N(1));  % RAAN [rad]

omega = acos(max(-1, min(1, dot(N, e_unit)))); % true anomaly [rad]
if e_unit(3) < 0
    omega = 2*pi - omega;
end

% Dtheta_2 = mod(theta2 - theta1, 2*pi)   % [rad] difference in true anomaly
% 
% if tm<0 % long-arc solution
%     Dtheta_2 = 2*pi - Dtheta_2 % [rad] difference in true anomaly
% end

end