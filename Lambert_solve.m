function [a, e, Omega, i, omega, theta1, theta2, eT_real] = Lambert_solve(mu, r1, r2, Dt, k)

% INPUTS
% mu - Gravitational parameter of central body [km3/m2]
% r1 - Position vector at time t1 [km] 
% r2 - Position vector at time t2 [km]
% eT - Transverse eccentricity component
% Dt - Time between position 1 and position 2 [s]
% k - Positive direction of the orbital plane normal (used to select between short & long arc solutions)

% OUTPUTS
% a - semi mayor axis [km]
% e - eccentricity [km]
% Omega - Right ascension of the ascending node = RAAN [rad]
% i - Inclination [rad]
% omega - Argument of the periapsis [rad]
% theta1 - True anomaly at t1 [rad]
% theta2 - True anomaly at t2 [rad]

% Solve Lambert's problem using the previous two functions. Internally, the function should use the bisection method. 
% Restrict the search to elliptic solutions - 0<e<1
% The units of input and output arguments must be the same as in previous functions.

% Position
R1 = norm(r1);
R2 = norm(r2);

% Chord vector
cc = r2 - r1;
c = norm(cc);

% Vector basis (unitary)
kk = k/norm(k);
ii = cc/c;
jj = cross(kk,ii);

% Eccentricity component
eF = (R1 - R2)/c;

aF = (R1+R2)/2; % semimajor axis of fundamental ellipse [km]

pF = aF*(1-eF^2); % semi-latum rectum of fundamental ellipse [km]

eTP = sqrt(1 - eF^2); % parabolic limit - from notes p.47
eT_lo = -eTP + 1e-6;
eT_hi =  eTP - 1e-6;

% bisection method to find eT that makes tof(eT)=Dt --> tof(eT)-Dt=0
tol = 1e-8;
max_iter = 200;

for iter = 1:max_iter
    eT_real = (eT_lo + eT_hi) / 2;

    % Evaluate tof at midpoint
    [a_mid, e_mid, ~, ~, ~, theta1_mid, theta2_mid] = functionshw2.Lambert_conic(r1, r2, eT_real, k);
    tof_real = functionshw2.tof(mu, a_mid, e_mid, theta1_mid, theta2_mid);
    
    err = tof_real - Dt; 
    if abs(err) < tol %if error<tol we have finished
        break;
    end

    % to determine if tof is monotonic
    [a_lo, e_lo, ~, ~, ~, th1_lo, th2_lo] = functionshw2.Lambert_conic(r1, r2, eT_lo, k);
    tof_lo = functionshw2.tof(mu, a_lo, e_lo, th1_lo, th2_lo);

    if (tof_lo - Dt) * err < 0
        eT_hi = eT_real;
    else
        eT_lo = eT_real;
    end
end

if iter == max_iter && abs(err) > tol
    warning('Lambert_solve: bisection did not converge. err=%.2e, Dt=%.2f s', err, Dt);
end

% Eccentricity vector
ee = eF*ii + eT_real*jj;
e = norm(ee); %sqrt(eF^2 + eT^2);
% e_unit = ee/e;

% Final solution with converged eT
[a, e, Omega, i, omega, theta1, theta2] = functionshw2.Lambert_conic(r1, r2, eT_real, k);

% % alpha1 = acos(dot(r1,ii));
% alpha1 = atan2(dot(cross(ii, r1/R1), kk), dot(ii, r1/R1));  
% % alpha2 = acos(dot(r2,ii));
% alpha2 = atan2(dot(cross(ii, r2/R2), kk), dot(ii, r2/R2));  
% 
% delta = atan2(eT, eF); % from simplifying e in delta=atan2(eT/e, eF/e);

% theta1 = delta + alpha1;
% theta2 = delta + alpha2;
Dtheta = theta2 - theta1;

% % Semi-latum rectum of ellipse w/ given eT
% p = pF - eT * (R1*R2)/c * sin(Dtheta);

% % Semi-major axis of ellipse w/ given eT
% a = p/(1 - e^2);

% from page 47 in class notes
% eTP = sqrt(1 - eF^2);
eTD = (pF*c)/(R1*R2*sin(Dtheta));
if eT_real < eTP && eT_real > -eTP
    % fprintf('\nOrbit is elliptical. Correct\n');
elseif eT_real == -eTP || eT_real == eTP
    fprintf('\nOrbit is a parabola only allowing short arc solution\n');
elseif eT_real < -eTP
    fprintf('\nHyperbolic short arcs\n')
elseif eT_real > eTP && eT_real < eTD
    fprintf('\nHyperbolic long arcs\n')
elseif eT_real == eTD
    fprintf('\np=0 and F=F´--> Pair of lines that intersect at the origin\n')
elseif eT_real > eTD
    fprintf('\nInadmissible hyperbolic solutions\n')
end

end