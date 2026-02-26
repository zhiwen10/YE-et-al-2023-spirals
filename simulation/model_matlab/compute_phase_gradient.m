function [grad_x, grad_y] = compute_phase_gradient(theta, dx, dy, Lx)
% COMPUTE_PHASE_GRADIENT Compute spatial gradients of phase field
%
% Inputs:
%   theta - phase field (Ny x Nx)
%   dx - grid spacing in x
%   dy - grid spacing in y
%   Lx - length in x (for periodic boundary)
%
% Outputs:
%   grad_x - gradient in x direction
%   grad_y - gradient in y direction

[Ny, Nx] = size(theta);

% Compute gradients with periodic boundary in x, regular in y
grad_x = zeros(Ny, Nx);
grad_y = zeros(Ny, Nx);

% X-gradient (periodic boundary)
for i = 1:Ny
    for j = 1:Nx
        j_next = mod(j, Nx) + 1;  % Periodic wrap
        j_prev = mod(j - 2, Nx) + 1;
        
        % Central difference with angle wrapping
        dtheta = angle_diff(theta(i, j_next), theta(i, j_prev));
        grad_x(i, j) = dtheta / (2 * dx);
    end
end

% Y-gradient (free boundary)
for i = 2:Ny-1
    for j = 1:Nx
        % Central difference
        dtheta = angle_diff(theta(i+1, j), theta(i-1, j));
        grad_y(i, j) = dtheta / (2 * dy);
    end
end

% One-sided differences at boundaries
for j = 1:Nx
    dtheta = angle_diff(theta(2, j), theta(1, j));
    grad_y(1, j) = dtheta / dy;
    
    dtheta = angle_diff(theta(Ny, j), theta(Ny-1, j));
    grad_y(Ny, j) = dtheta / dy;
end

end

function diff = angle_diff(theta1, theta2)
% Compute theta1 - theta2, accounting for circular nature
diff = theta1 - theta2;
diff = mod(diff + pi, 2*pi) - pi;  % Map to [-pi, pi]
end
