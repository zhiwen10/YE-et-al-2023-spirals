% Kuramoto model starting from Euclidean grid points
% Maps Euclidean (x,y) points to polar (theta, r) coordinates
clear; close all;

%% Parameters
% Network geometry in Euclidean space
N_euclidean = 80;  % Grid points per dimension in Euclidean space
r_max = 1;         % Maximum radius (only use points within this radius)

% Kuramoto parameters
K = 1.2;           % Coupling strength
r_coupling = 0.5;  % Coupling radius in POLAR coordinate space
p_coupling = 0.25; % probability of forming a link if within radius
omega_mean = 5;    % Mean natural frequency
omega_std = 0.5;   % Std dev of natural frequencies

% Time parameters
dt = 0.01;         % Time step
T_total = 30;     % Total simulation time
t_plot_interval = 25;  % Save/plot every N steps
t_analyze_start = 400;  % Start analyzing after transients

%% Setup Euclidean grid and filter points within radius
x_euclidean = linspace(-r_max, r_max, N_euclidean);
y_euclidean = linspace(-r_max, r_max, N_euclidean);
[X_euclidean, Y_euclidean] = meshgrid(x_euclidean, y_euclidean);

% Find points within the circle
R_euclidean = sqrt(X_euclidean.^2 + Y_euclidean.^2);
mask = R_euclidean <= r_max;

% Extract valid points
x_points = X_euclidean(mask);
y_points = Y_euclidean(mask);
n_points = length(x_points);

fprintf('Total points within radius: %d\n', n_points);

% Convert to polar coordinates
[theta_points, r_points] = cart2pol(x_points, y_points);
% theta is in [-pi, pi], convert to [0, 2*pi]
theta_points = mod(theta_points, 2*pi);

%% Visualize the point distribution
figure('Position', [100 100 1200 500]);

subplot(1, 3, 1);
scatter(x_points, y_points, 5, 'filled');
axis equal;
xlabel('x');
ylabel('y');
title('Points in Euclidean Space');
grid on;
viscircles([0 0], r_max, 'Color', 'r', 'LineStyle', '--');

subplot(1, 3, [2,3]);
scatter(theta_points, r_points, 5, 'filled');
xlabel('\theta (radians)');
ylabel('r');
title('Points in Polar Coordinates');
grid on;
xlim([0 2*pi]);
ylim([0 r_max]);

%% Initialize natural frequencies for each point
omega_points = omega_mean + omega_std * randn(n_points, 1);

%% Compute coupling matrix based on polar coordinates
fprintf('Computing coupling matrix based on polar distances...\n');

% Polar distance parameters
theta_scale = 1;  % scaling for angular distance
r_scale = 1;      % scaling for radial distance
Lx = 2*pi;        % Periodic in theta

coupling = cell(n_points, 1);

for i = 1:n_points
    neighbors_idx = [];
    weights = [];
    
    for j = 1:n_points
        if i ~= j
            % Distance calculation in polar coordinates (periodic in theta, regular in r)
            dtheta_periodic = theta_scale * min(abs(theta_points(i) - theta_points(j)), ...
                                                Lx - abs(theta_points(i) - theta_points(j)));
            dr_regular = r_scale * abs(r_points(i) - r_points(j));
            dist = sqrt(dtheta_periodic^2 + dr_regular^2);
            
            if dist <= r_coupling && rand(1) < p_coupling
                neighbors_idx = [neighbors_idx; j];
                weights = [weights; 1.0];
            end
        end
    end
    
    % Normalize weights
    if ~isempty(weights)
        coupling{i}.neighbors = neighbors_idx;
        coupling{i}.weights = weights / sum(weights);
    else
        coupling{i}.neighbors = [];
        coupling{i}.weights = [];
    end
end

fprintf('Coupling matrix computed (polar metric).\n');
%% Visualize coupling structure
figure('Position', [100 100 800 800]);
pPlot = 0.0003;  % Plot a fraction of connections
hold on;

for i = 1:n_points
    if ~isempty(coupling{i}.neighbors)
        neighbors_idx = coupling{i}.neighbors;
        thisX = x_points(i);
        thisY = y_points(i);
        
        for qq = 1:numel(neighbors_idx)
            if rand(1) < pPlot
                thatX = x_points(neighbors_idx(qq));
                thatY = y_points(neighbors_idx(qq));
                plot([thisX, thatX], [thisY, thatY], 'k'); 
            end
        end
    end
end

axis equal;
xlabel('x');
ylabel('y');
title('Coupling Structure (Cartesian coordinates)');
% viscircles([0 0], r_max, 'Color', 'r', 'LineStyle', '--');
% grid on;
%% Run Kuramoto model simulation
[phase_history, time_points, grad_theta_all, grad_r_all] = ...
    run_kuramoto_points(x_points, y_points, theta_points, r_points, ...
                        coupling, K, omega_points, dt, T_total, ...
                        t_analyze_start, t_plot_interval, ...
                        true, X_euclidean, Y_euclidean, mask);
%% Plot final phase distribution
fprintf('Creating final visualizations...\n');

figure('Position', [100 100 1200 500]);

% Euclidean space
subplot(1, 2, 1);
final_phase = phase_history(:, end);
scatter(x_points, y_points, 50, final_phase, 'filled');
colorbar;
colormap(hsv);
caxis([0 2*pi]);
axis equal;
xlabel('x');
ylabel('y');
title(sprintf('Final Phase (Euclidean) at t = %.1f', time_points(end)));
viscircles([0 0], r_max, 'Color', 'k', 'LineStyle', '--');

% Polar space
subplot(1, 2, 2);
scatter(theta_points, r_points, 50, final_phase, 'filled');
colorbar;
colormap(hsv);
caxis([0 2*pi]);
xlabel('\theta (radians)');
ylabel('r');
title('Final Phase (Polar)');
xlim([0 2*pi]);
ylim([0 r_max]);

fprintf('Done!\n');
