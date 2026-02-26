function [phase_history, time_points, grad_x_all, grad_y_all] = ...
    run_kuramoto_model(X, Y, coupling, K, omega, dt, T_total, t_analyze_start, save_interval, show_animation, nearestIdx, xr, yr)
% RUN_KURAMOTO_MODEL Simulate Kuramoto oscillators and return phase history
%
% Inputs:
%   X - meshgrid of x (theta) coordinates
%   Y - meshgrid of y (radius) coordinates  
%   coupling - coupling structure from compute_coupling_matrix
%   K - coupling strength
%   omega - natural frequencies (Ny x Nx matrix)
%   dt - time step
%   T_total - total simulation time
%   t_analyze_start - time to start collecting gradient statistics
%   save_interval - save phase every N steps (default: 25)
%   show_animation - boolean to show real-time animation (default: false)
%   nearestIdx - nearest neighbor indices for polar to cartesian mapping (required if show_animation=true)
%   xr - x coordinates for real space grid (required if show_animation=true)
%   yr - y coordinates for real space grid (required if show_animation=true)
%
% Outputs:
%   phase_history - 3D array (Ny x Nx x n_saves) of phases at saved time points
%   time_points - vector of times corresponding to saved phases
%   grad_x_all - collected x-gradients after t_analyze_start
%   grad_y_all - collected y-gradients after t_analyze_start

if nargin < 9
    save_interval = 25;
end

if nargin < 10
    show_animation = false;
end

[Ny, Nx] = size(X);
Lx = max(X(:)) - min(X(:));
dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);

% Initialize phases
theta = 2*pi * rand(Ny, Nx);  % Random initial phases

%% Convert coupling structure to sparse matrix format
fprintf('Converting to sparse matrix format...\n');

% Count total connections
total_connections = 0;
for i = 1:Ny
    for j = 1:Nx
        if ~isempty(coupling{i,j}.i)
            total_connections = total_connections + length(coupling{i,j}.i);
        end
    end
end

% Pre-allocate arrays for sparse matrix
source_idx = zeros(total_connections, 1);
target_idx = zeros(total_connections, 1);
coupling_weights = zeros(total_connections, 1);

% Fill arrays
connection_counter = 1;
for i = 1:Ny
    for j = 1:Nx
        if ~isempty(coupling{i,j}.i)
            linear_self = sub2ind([Ny, Nx], i, j);
            n_neighbors = length(coupling{i,j}.i);
            linear_neighbors = sub2ind([Ny, Nx], coupling{i,j}.i, coupling{i,j}.j);
            
            source_idx(connection_counter:connection_counter+n_neighbors-1) = linear_neighbors;
            target_idx(connection_counter:connection_counter+n_neighbors-1) = linear_self;
            coupling_weights(connection_counter:connection_counter+n_neighbors-1) = K * coupling{i,j}.w;
            
            connection_counter = connection_counter + n_neighbors;
        end
    end
end

% Create sparse adjacency matrix
N_total = Ny * Nx;
W_sparse = sparse(target_idx, source_idx, coupling_weights, N_total, N_total);

fprintf('Sparse matrix created: %d total connections\n', total_connections);
fprintf('Average connections per oscillator: %.1f\n', total_connections/N_total);
%% Setup visualization if requested
if show_animation
    if nargin < 11 || isempty(nearestIdx) || isempty(xr) || isempty(yr)
        error('nearestIdx, xr, and yr are required when show_animation is true');
    end
    
    x_coords = X(1, :);
    y_coords = Y(:, 1);
    
    fig_anim = figure('Position', [100 100 900 300]);
    
    % Polar coordinates subplot
    ax1 = subplot(1,3, [1,2]);
    im1 = imagesc(x_coords, y_coords, zeros(Ny, Nx));
    colorbar;
    colormap(hsv);
    caxis([0 2*pi]);
    xlabel('theta');
    ylabel('radius');
    title_handle1 = title('Phase field at t = 0');
    axis xy;
    
    % Euclidean coordinates subplot
    ax2 = subplot(1,3, 3);
    Xr_plot = xr;
    Yr_plot = yr;
    im2 = imagesc(Xr_plot, Yr_plot, zeros(length(yr), length(xr)));
    colorbar;
    colormap(hsv);
    caxis([0 2*pi]);
    xlabel('x');
    ylabel('y');
    title('Phase in Real Space');
    axis equal;
    axis off;
end
%% Time integration
num_steps = round(T_total / dt);
n_saves = floor(num_steps / save_interval);

% Pre-allocate storage for phase history
phase_history = zeros(Ny, Nx, n_saves);
time_points = zeros(n_saves, 1);
save_counter = 0;

% Storage for phase gradients
grad_x_all = [];
grad_y_all = [];

fprintf('Starting simulation...\n');
for step = 1:num_steps
    t = step * dt;
    
    % Natural frequency term
    dtheta_omega = omega * dt;
    
    % Compute coupling term using sparse matrix operations
    theta_vec = theta(:);
    sin_theta = sin(theta_vec);
    cos_theta = cos(theta_vec);
    
    W_sin = W_sparse * sin_theta;
    W_cos = W_sparse * cos_theta;
    
    coupling_term_vec = W_sin .* cos_theta - W_cos .* sin_theta;
    dtheta_coupling = reshape(coupling_term_vec * dt, [Ny, Nx]);
    
    % Update phases
    theta = theta + dtheta_coupling + dtheta_omega;
    theta = mod(theta, 2*pi);  % Keep in [0, 2pi]
    
    % Save phase at specified intervals
    if mod(step, save_interval) == 0
        save_counter = save_counter + 1;
        phase_history(:, :, save_counter) = theta;
        time_points(save_counter) = t;
        
        % Update visualization if animation is enabled
        if show_animation
            % Update polar coordinates plot
            set(im1, 'CData', theta);
            set(title_handle1, 'String', sprintf('Phase field at t = %.1f', t));
            
            % Update euclidean coordinates plot
            realTheta = zeros(size(nearestIdx)); 
            idxs = nearestIdx(~isnan(nearestIdx)); 
            realTheta(~isnan(nearestIdx)) = theta(idxs);
            set(im2, 'CData', realTheta);
            
            drawnow;
        end
    end
    
    % Collect phase gradients after transients
    if t >= t_analyze_start
        [grad_x, grad_y] = compute_phase_gradient(theta, dx, dy, Lx);
        grad_x_all = [grad_x_all; grad_x(:)];
        grad_y_all = [grad_y_all; grad_y(:)];
    end
    
    % Progress update
    if mod(step, 1000) == 0
        fprintf('Step %d/%d (t = %.1f)\n', step, num_steps, t);
    end
end

fprintf('Simulation complete!\n');

end
