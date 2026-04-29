function [theta_burn, Dv] = hyperbola(mu, a, e, Omega, i, omega, vh)
% INPUTS
% mu    - Gravitational parameter of central body [km^3/s^2]
% a     - Semi-major axis of parking orbit [km]
% e     - Eccentricity of parking orbit [-]
% Omega - RAAN of parking orbit [rad]
% i     - Inclination of parking orbit [rad]
% omega - Argument of periapsis of parking orbit [rad]
% vh    - Velocity vector at infinity [km/s]
% OUTPUTS
% theta_burn - True anomaly of burn point [rad]
% Dv         - Injection impulse vector [km/s]

theta_vec = linspace(0, 2*pi, 360);

% Semi-latus rectum of parking orbit [km] - use p_park to avoid conflict
p_park = a*(1-e^2);

v_inf = norm(vh);
a_hyp = -mu/v_inf^2;

% Rotation matrix from orbit to ECI
R_ORB2ECI = [cos(Omega)*cos(omega)-sin(Omega)*sin(omega)*cos(i),  -cos(Omega)*sin(omega)-sin(Omega)*cos(omega)*cos(i),  sin(Omega)*sin(i);
             sin(Omega)*cos(omega)+cos(Omega)*sin(omega)*cos(i),  -sin(Omega)*sin(omega)+cos(Omega)*cos(omega)*cos(i), -cos(Omega)*sin(i);
             sin(omega)*sin(i),                                    cos(omega)*sin(i),                                   cos(i)];

% Perifocal basis vectors in ECI
p_hat = R_ORB2ECI(:,1);   % eccentricity direction
q_hat = R_ORB2ECI(:,2);   % 90 deg ahead in orbit plane
w_hat = R_ORB2ECI(:,3);   % orbital plane normal

% Pre-allocate
dv_mag = inf(1, 360);
Dv_all = zeros(3, 360);

for ii = 1:length(theta_vec)
    th = theta_vec(ii);

    % Position and velocity on parking orbit
    r_mag = p_park / (1 + e*cos(th));
    r_vec = r_mag * (cos(th)*p_hat + sin(th)*q_hat);
    r_hat = r_vec / r_mag;

    v_park_mag = sqrt(mu*(2/r_mag - 1/a));
    vr_hat = (-sin(th)*p_hat + (e+cos(th))*q_hat) / sqrt(1+e^2+2*e*cos(th));
    v_park = v_park_mag * vr_hat;

    % Two possible hyperbola plane normals
    cross_r_vh = cross(r_vec, vh);
    if norm(cross_r_vh) < 1e-10
        fprintf('theta = %.2f deg: r_vec and vh are parallel, skipping.\n', rad2deg(th));
        continue;
    end
    normals = [cross_r_vh/norm(cross_r_vh), -cross_r_vh/norm(cross_r_vh)];

    dv_best = inf;
    Dv_best = zeros(3,1);

    for jj = 1:2
        kk = normals(:,jj);

        v_hyp_mag = sqrt(v_inf^2 + 2*mu/r_mag);

        vh_hat = vh/v_inf;
        r_along_vh = dot(r_vec, vh_hat);
        r_perp_vh = r_vec - r_along_vh*vh_hat;
        r_perp_in_plane = r_perp_vh - dot(r_perp_vh, kk)*kk;
        B_eff = norm(r_perp_in_plane);

        h_hyp = B_eff * v_inf;
        if h_hyp < 1e-6
            continue;
        end

        p_hyp_val = h_hyp^2 / mu;
        x = p_hyp_val/r_mag - 1;
        v_t_hyp = h_hyp/r_mag;

        v_r_sq = v_hyp_mag^2 - v_t_hyp^2;
        if v_r_sq < 0
            continue;
        end

        v_r_hyp = sqrt(v_r_sq);
        y = v_r_hyp*h_hyp/mu;
        e_hyp_val = sqrt(x^2 + y^2);
        if e_hyp_val <= 1
            continue;
        end

        v_t_dir = cross(kk, r_hat);
        v_t_dir = v_t_dir/norm(v_t_dir);

        for sign_vr = [+1, -1]
            v_hyp_vec = v_t_hyp*v_t_dir + sign_vr*v_r_hyp*r_hat;
            if dot(v_hyp_vec, vh_hat) < 0
                continue;
            end
            dv_vec = v_hyp_vec - v_park;
            dv_val = norm(dv_vec);
            if dv_val < dv_best
                dv_best = dv_val;
                Dv_best = dv_vec;
            end
        end
    end

    dv_mag(ii)   = dv_best;
    Dv_all(:,ii) = Dv_best;
end

[~, best_idx] = min(dv_mag);
theta_burn = theta_vec(best_idx);
Dv = Dv_all(:, best_idx);

end