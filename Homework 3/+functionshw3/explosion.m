function [v,bc] = explosion(mu_v,sigma_v,mu_bc,sigma_bc)

% INPUTS
%   mu_v     - [km/s]           Mean velocity
%   sigma_v  - [km/s]           Velocity standard deviation
%   mu_bc    - [log(kg/m^2)]    Logarithmic mean of the ballistic coefficient
%   sigma_bc - [-]              Ballistic coefficient standard deviation
%
% OUTPUTS
%   v        - [km/s]           velocity array
%   bc       - [log(kg/m^2)]    Logarithmic ballistic coefficient

N = 30;                              % number of pieces of debris

rng('default')
v  = mu_v + sigma_v * randn(3, N);   % [km/s]        (3×N)
bc = lognrnd(mu_bc, sigma_bc, 1, N); % [kg/m^2]      (1×N)
end