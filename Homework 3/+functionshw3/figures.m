function [] = figures(zeta_a, zeta_p, tau, a_M, e_M, Omega_M, case_names, c)

% INPUTS
%   fig_num [-]   figure counter
%   tau     [s]   (M+1 x N) orbital period
%   r_p     [km]  (M+1 x N) pericenter radius
%   r_a     [km]  (M+1 x N) apocenter radius
%   a_M     [km]  (M+1 x N) semi-major axis
%   e_M     [-]   (M+1 x N) eccentricity
%   Omega_M [rad] (M+1 x N) RAAN
%
% OUTPUT
%   fig_num [-]   updated figure counter
%
RE   = 6371;             % [km]
M    = size(tau, 1) - 1; % number of time steps
% Define colormap for the weeks (M+1 colors for t=0 and M weeks)
    cmap = parula(M+1);

    % Colorbar tick settings (one label per month = every 4 weeks)
    cb_ticks  = 0:4:M;
    cb_labels = arrayfun(@(w) sprintf('%d',w), cb_ticks, 'UniformOutput', false);

    
    % (i) Gabbard plot: perigee and apogee altitude vs orbital period
   
    figure; hold on;
    for k = 1:(M+1)
        % Plot the cloud at week k (using HandleVisibility 'off' so they don't flood the legend)
        plot(tau(k,:)/60, zeta_p(k,:), '*', 'Color', cmap(k,:), 'HandleVisibility', 'off');
        plot(tau(k,:)/60, zeta_a(k,:), '^', 'Color', cmap(k,:), 'HandleVisibility', 'off');
    end

    % Create dummy plots just for the legend
    plot(NaN, NaN, '*k', 'DisplayName', 'Perigee \zeta_p');
    plot(NaN, NaN, '^k', 'DisplayName', 'Apogee \zeta_a');
    legend('Location', 'best');

    title(['(i) Gabbard Plot | Case: ', case_names{c}]);
    xlabel('Orbital period \tau [min]');
    ylabel('Altitude \zeta [km]');
    colormap(gca, cmap); clim([0 M]);
    cb = colorbar; cb.Label.String = 'Time [weeks]';
    cb.Ticks = cb_ticks; cb.TickLabels = cb_labels;
    grid on; hold off;

    
    % (ii) Eccentricity (e) vs Semi-major axis (a)
    
    figure; hold on;
    for k = 1:(M+1)
        plot(a_M(k,:), e_M(k,:), '*', 'Color', cmap(k,:));
    end
    title(['(ii) Eccentricity vs Semi-major axis | Case: ', case_names{c}]);
    xlabel('Semi-major axis a [km]');
    ylabel('Eccentricity e [-]');
    colormap(gca, cmap); clim([0 M]);
    cb = colorbar; cb.Label.String = 'Time [weeks]';
    cb.Ticks = cb_ticks; cb.TickLabels = cb_labels;
    grid on; hold off;

    
    % (iii) RAAN (Omega) vs Semi-major axis (a)
    
    figure; hold on;
    for k = 1:(M+1)
        plot(a_M(k,:), rad2deg(Omega_M(k,:)), '*', 'Color', cmap(k,:));
    end
    title(['(iii) RAAN \Omega vs Semi-major axis | Case: ', case_names{c}]);
    xlabel('Semi-major axis a [km]');
    ylabel('RAAN \Omega [deg]');
    colormap(gca, cmap); clim([0 M]);
    cb = colorbar; cb.Label.String = 'Time [weeks]';
    cb.Ticks = cb_ticks; cb.TickLabels = cb_labels;
    grid on; hold off;

    drawnow; % This is important(!) for the .mlx workspace.

end

