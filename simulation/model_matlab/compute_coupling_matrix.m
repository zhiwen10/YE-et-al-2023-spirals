function coupling = compute_coupling_matrix(X, Y, r_coupling, p_coupling, theta_scale, r_scale, Lx)
% COMPUTE_COUPLING_MATRIX Compute coupling structure for Kuramoto oscillators
%
% Inputs:
%   X - meshgrid of x (theta) coordinates
%   Y - meshgrid of y (radius) coordinates
%   r_coupling - coupling radius
%   p_coupling - probability of forming a link if within radius
%   theta_scale - scaling factor for theta distance
%   r_scale - scaling factor for radius distance
%   Lx - length/circumference (for periodic boundary in x)
%
% Output:
%   coupling - cell array where coupling{i,j} contains:
%              .i - neighbor row indices
%              .j - neighbor column indices
%              .w - normalized weights

[Ny, Nx] = size(X);

fprintf('Computing coupling matrix...\n');
coupling = cell(Ny, Nx);  % Store neighbors and weights for each oscillator

for i = 1:Ny
    for j = 1:Nx
        neighbors_i = [];
        neighbors_j = [];
        weights = [];
        
        for ii = 1:Ny
            for jj = 1:Nx
                % Distance calculation (periodic in x, regular in y)
                dx_periodic = theta_scale * min(abs(X(i,j) - X(ii,jj)), Lx - abs(X(i,j) - X(ii,jj)));
                dy_regular = r_scale * abs(Y(i,j) - Y(ii,jj));
                dist = sqrt(dx_periodic^2 + dy_regular^2);
                
                if dist > 0 && dist <= r_coupling && rand(1) < p_coupling
                    neighbors_i = [neighbors_i; ii];
                    neighbors_j = [neighbors_j; jj];
                    weights = [weights; 1.0];  % Uniform coupling within radius
                end
            end
        end
        
        % Normalize weights
        if ~isempty(weights)
            coupling{i,j}.i = neighbors_i;
            coupling{i,j}.j = neighbors_j;
            coupling{i,j}.w = weights / sum(weights);
        else
            coupling{i,j}.i = [];
            coupling{i,j}.j = [];
            coupling{i,j}.w = [];
        end
    end
end

fprintf('Coupling matrix computed.\n');

end
