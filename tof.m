function Dt = tof(mu,a,e,theta1,theta2)

theta = [theta1, theta2];

n = sqrt(mu/a^3);

E = zeros(size(theta));
M = zeros(size(theta));
t = zeros(size(theta));

for i = 1:length(theta)
    E(i) = functionshw2.theta2E(theta(i), e);
    M(i) = functionshw2.E2M(E(i), e);
    t(i) = M(i)/n;  % time from periapsis [s]
end

Dt = t(end) - t(1);

if Dt < 0 % If Dt<0 the arc crosses periapsis → add one full period
    T = 2*pi/n;   % orbital period [s]
    Dt = Dt + T;
end

end