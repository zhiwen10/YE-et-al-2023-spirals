function [phase_history, time_points, grad_theta_all, grad_r_all] = ...
    run_kuramoto_points(x_points, y_points, theta_points, r_points, ...
                        coupling, K, omega, dt, T_total, ...
                        t_analyze_start, save_interval, ...
                        show_animation, X_grid, Y_grid, mask)
% RUN_KURAMOTO_POINTS Simulate Kuramoto oscillators on point cloud
%
% Inputs:
%   x_points, y_points - Euclidean coordinates of oscillators
%   theta_points, r_points - Polar coordinates of oscillators
%   coupling - coupling structure (cell array)
%   K - coupling strength
%   omega - natural frequencies (n_points x 1 vector)
%   dt - time step
%   T_total - total simulation time
%   t_analyze_start - time to start collecting statistics
%   save_interval - save phase every N steps
%   show_animation - boolean to show real-time animation
%   X_grid, Y_grid - meshgrid for Euclidean visualization
%   mask - logical mask for points within radius
%
% Outputs:
%   phase_history - 2D array (n_points x n_saves) of phases
%   time_points - vector of times
%   grad_theta_all, grad_r_all - collected gradients

n_points = length(x_points);

% Initialize phases
theta = 2*pi * rand(n_points, 1);  % Random initial phases

%% Convert coupling structure to sparse matrix format
fprintf('Converting to sparse matrix format...\n');

% Count total connections
total_connections = 0;
for i = 1:n_points
    if ~isempty(coupling{i}.neighbors)
        total_connections = total_connections + length(coupling{i}.neighbors);
    end
end

% Pre-allocate arrays for sparse matrix
source_idx = zeros(total_connections, 1);
target_idx = zeros(total_connections, 1);
coupling_weights = zeros(total_connections, 1);

% Fill arrays
connection_counter = 1;
for i = 1:n_points
    if ~isempty(coupling{i}.neighbors)
        n_neighbors = length(coupling{i}.neighbors);
        
        source_idx(connection_counter:connection_counter+n_neighbors-1) = coupling{i}.neighbors;
        target_idx(connection_counter:connection_counter+n_neighbors-1) = i;
        coupling_weights(connection_counter:connection_counter+n_neighbors-1) = K * coupling{i}.weights;
        
        connection_counter = connection_counter + n_neighbors;
    end
end

% Create sparse adjacency matrix
W_sparse = sparse(target_idx, source_idx, coupling_weights, n_points, n_points);

fprintf('Sparse matrix created: %d total connections\n', total_connections);
fprintf('Average connections per oscillator: %.1f\n', total_connections/n_points);

%% Setup visualization if requested
if show_animation
    fig_anim = figure('Position', [100 100 1200 500]);
    
    % Euclidean coordinates subplot
    ax1 = subplot(1, 3, 1);
    % Create image representation
    phase_grid1 = NaN(size(X_grid));
    phase_grid1(mask) = theta;
    im1 = imagesc(X_grid(1,:), Y_grid(:,1), phase_grid1);
    colorbar;
    colormap(hsv);
    caxis([0 2*pi]);
    xlabel('x');
    ylabel('y');
    title_handle1 = title('Phase field at t = 0 (Euclidean)');
    axis xy equal;
    axis off;
    
    % Polar coordinates subplot
    ax2 = subplot(1, 3, [2,3]);
    scatter_handle = scatter(theta_points, r_points, 50, theta, 'filled');
    colorbar;
    colormap(hsv);
    caxis([0 2*pi]);
    xlabel('\theta (radians)');
    ylabel('r');
    title_handle2 = title('Phase at t = 0 (Polar)');
    xlim([0 2*pi]);
    ylim([0 max(r_points)*1.1]);
end

%% Time integration
num_steps = round(T_total / dt);
n_saves = floor(num_steps / save_interval);

% Pre-allocate storage for phase history
phase_history = zeros(n_points, n_saves);
time_points = zeros(n_saves, 1);
save_counter = 0;

% Storage for gradients (placeholder for now)
grad_theta_all = [];
grad_r_all = [];

fprintf('Starting simulation...\n');
for step = 1:num_steps
    t = step * dt;
    
    % Natural frequency term
    dtheta_omega = omega * dt;
    
    % Compute coupling term using sparse matrix operations
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    
    W_sin = W_sparse * sin_theta;
    W_cos = W_sparse * cos_theta;
    
    coupling_term = W_sin .* cos_theta - W_cos .* sin_theta;
    dtheta_coupling = coupling_term * dt;
    
    % Update phases
    theta = theta + dtheta_coupling + dtheta_omega;
    theta = mod(theta, 2*pi);  % Keep in [0, 2pi]
    
    % Save phase at specified intervals
    if mod(step, save_interval) == 0
        save_counter = save_counter + 1;
        phase_history(:, save_counter) = theta;
        time_points(save_counter) = t;
        
        % Update visualization if animation is enabled
        if show_animation
            % Update Euclidean coordinates plot
            phase_grid1 = NaN(size(X_grid));
            phase_grid1(mask) = theta;
            set(im1, 'CData', phase_grid1);
            set(title_handle1, 'String', sprintf('Phase at t = %.1f (Euclidean)', t));
            
            % Update polar coordinates plot
            set(scatter_handle, 'CData', theta);
            set(title_handle2, 'String', sprintf('Phase at t = %.1f (Polar)', t));
            
            drawnow;
        end
    end
    
    % Progress update
    if mod(step, 1000) == 0
        fprintf('Step %d/%d (t = %.1f)\n', step, num_steps, t);
    end
end

fprintf('Simulation complete!\n');

end
