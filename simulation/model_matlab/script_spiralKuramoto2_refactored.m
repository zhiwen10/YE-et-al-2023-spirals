% Kuramoto model on a cylinder - REFACTORED VERSION
% Uses modular functions for coupling computation and model iteration
clear; close all;

%% Parameters
% Network geometry
Lx = 2*pi;  % Length/circumference (periodic)
Ly = 1;     % Width/height (free boundaries)
Nx = 200;   % Grid points in x (long axis)
Ny = 32;    % Grid points in y (short axis)

% Kuramoto parameters
K = 1.2;           % Coupling strength
r_coupling = 0.5;  % Coupling radius
p_coupling = 0.25; % probability of forming a link if within radius
omega_mean = 5;    % Mean natural frequency
omega_std = 0.5;   % Std dev of natural frequencies
theta_scale = 1;   % scaling for theta distance
r_scale = 1;       % scaling for radius distance

% Time parameters
dt = 0.01;         % Time step
T_total = 30;     % Total simulation time
t_plot_interval = 25;  % Save every N steps
t_analyze_start = 400;  % Start analyzing after transients

%% Setup grid
x = linspace(0, Lx, Nx); % really this is theta
y = linspace(0, Ly, Ny); % this is radius
[X, Y] = meshgrid(x, y);

% 'real space' grid, with nearest neighbor indices
Nr = 100;
xr = linspace(-1, 1, Nr);
yr = linspace(-1, 1, Nr);
[Xr, Yr] = meshgrid(xr, yr); 
[xCart, yCart] = pol2cart(X, Y);

nearestIdx = NaN(size(Xr)); 
for p = 1:numel(xr)
    for q = 1:numel(yr)
        if xr(p)^2 + yr(q)^2 < 1 % within circle
            dists = ((Xr(p,q) - xCart).^2 + (Yr(p,q) - yCart).^2).^(0.5);
            [~, idx] = min(dists(:));
            nearestIdx(p,q) = idx;
        end
    end
end

%% Initialize natural frequencies
omega = omega_mean + omega_std * randn(Ny, Nx);

%% Compute coupling matrix using function
coupling = compute_coupling_matrix(X, Y, r_coupling, p_coupling, theta_scale, r_scale, Lx);

%% Visualize coupling structure
figure('Position', [100 100 600 600]);
pPlot = 0.0003;
hold on; 
for i = 5:Ny
    for j = 1:Nx
        if ~isempty(coupling{i,j}.i)
            neighbors_i = coupling{i,j}.i;
            neighbors_j = coupling{i,j}.j;
            [thisX, thisY] = pol2cart(x(j), y(i));
            
            for qq = 1:numel(neighbors_i)
                if rand(1) < pPlot
                    [thatX, thatY] = pol2cart(x(neighbors_j(qq)), y(neighbors_i(qq)));
                    plot([thisX, thatX], [thisY, thatY], 'k'); 
                end
            end
        end
    end
end
title('Coupling Structure');
axis equal;
drawnow;

%% Run Kuramoto model simulation using function (with real-time animation)
[phase_history, time_points, grad_x_all, grad_y_all] = ...
    run_kuramoto_model(X, Y, coupling, K, omega, dt, T_total, t_analyze_start, ...
                       t_plot_interval, true, nearestIdx, xr, yr);

%% Plot results
fprintf('Creating visualizations...\n');

% Create figure for phase evolution
figure('Position', [100 100 1200 800]);

% Determine which time points to plot
n_plots = min(9, size(phase_history, 3));
plot_indices = round(linspace(1, size(phase_history, 3), n_plots));

for idx = 1:n_plots
    subplot(3, 3, idx);
    
    phase_snapshot = phase_history(:, :, plot_indices(idx));
    imagesc(x, y, phase_snapshot);
    colorbar;
    colormap(hsv);
    caxis([0 2*pi]);
    xlabel('theta');
    ylabel('radius');
    title(sprintf('t = %.1f', time_points(plot_indices(idx))));
    axis xy;
end
sgtitle('Phase Field Evolution');

% Create figure for final state in both coordinate systems
figure('Position', [100 100 1200 500]);

% Final phase in cylindrical coordinates
subplot(1, 2, 1);
final_phase = phase_history(:, :, end);
imagesc(x, y, final_phase);
colorbar;
colormap(hsv);
caxis([0 2*pi]);
xlabel('theta');
ylabel('radius');
title(sprintf('Final Phase Field (t = %.1f)', time_points(end)));
axis xy;

% Final phase in real space
subplot(1, 2, 2);
realTheta = zeros(size(Xr)); 
idxs = nearestIdx(~isnan(nearestIdx)); 
realTheta(~isnan(nearestIdx)) = final_phase(idxs);
imagesc(xr, yr, realTheta);
colorbar;
colormap(hsv);
caxis([0 2*pi]);
xlabel('x');
ylabel('y');
title('Final Phase in Real Space');
axis equal;

fprintf('Done!\n');
