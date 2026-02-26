"""
kuramoto_models.py

Module containing functions for Kuramoto oscillator simulations with
Cartesian and Polar coupling schemes.
"""

import numpy as np
from scipy.sparse import csr_matrix
import time


def setup_euclidean_grid(N_euclidean, r_max=1.0):
    """
    Create an evenly spaced Euclidean grid within a circular region.
    
    Parameters:
    -----------
    N_euclidean : int
        Number of grid points per dimension
    r_max : float
        Maximum radius for circular boundary
        
    Returns:
    --------
    x_points : ndarray
        x coordinates of points within circle
    y_points : ndarray
        y coordinates of points within circle
    theta_points : ndarray
        Angular coordinates (0 to 2π)
    r_points : ndarray
        Radial coordinates (0 to r_max)
    mask : ndarray
        Boolean mask for grid points within circle
    X_grid : ndarray
        Meshgrid X coordinates
    Y_grid : ndarray
        Meshgrid Y coordinates
    """
    x_euclidean = np.linspace(-r_max, r_max, N_euclidean)
    y_euclidean = np.linspace(-r_max, r_max, N_euclidean)
    X_grid, Y_grid = np.meshgrid(x_euclidean, y_euclidean)
    
    # Find points within circle
    R_euclidean = np.sqrt(X_grid**2 + Y_grid**2)
    mask = R_euclidean <= r_max
    
    # Extract valid points
    x_points = X_grid[mask]
    y_points = Y_grid[mask]
    
    # Convert to polar coordinates
    theta_points = np.arctan2(y_points, x_points)
    r_points = np.sqrt(x_points**2 + y_points**2)
    theta_points = np.mod(theta_points, 2*np.pi)  # Convert to [0, 2π]
    
    print(f'Grid created: {len(x_points)} points within radius {r_max}')
    
    return x_points, y_points, theta_points, r_points, mask, X_grid, Y_grid


def compute_coupling_cartesian(x_points, y_points, r_coupling, p_coupling):
    """
    Compute coupling structure based on Cartesian/Euclidean distances.
    
    Parameters:
    -----------
    x_points, y_points : ndarray
        Cartesian coordinates of oscillators
    r_coupling : float
        Coupling radius
    p_coupling : float
        Probability of forming connection within radius
        
    Returns:
    --------
    coupling : list of dict
        Each element contains 'neighbors' and 'weights' arrays
    """
    n_points = len(x_points)
    coupling = []
    
    print('Computing coupling matrix (Cartesian distances)...')
    
    for i in range(n_points):
        neighbors_idx = []
        weights = []
        
        for j in range(n_points):
            if i != j:
                # Euclidean distance
                dx = x_points[i] - x_points[j]
                dy = y_points[i] - y_points[j]
                dist = np.sqrt(dx**2 + dy**2)
                
                if dist <= r_coupling and np.random.rand() < p_coupling:
                    neighbors_idx.append(j)
                    weights.append(1.0)
        
        # Normalize weights
        if len(weights) > 0:
            weights = np.array(weights) / np.sum(weights)
            coupling.append({
                'neighbors': np.array(neighbors_idx, dtype=int),
                'weights': weights
            })
        else:
            coupling.append({
                'neighbors': np.array([], dtype=int),
                'weights': np.array([])
            })
        
        if (i + 1) % 500 == 0:
            print(f'  Processed {i+1}/{n_points} points')
    
    n_connections = sum([len(c['neighbors']) for c in coupling])
    print(f'Coupling computed: {n_connections} total connections')
    print(f'Average connections per oscillator: {n_connections/n_points:.1f}')
    
    return coupling


def compute_coupling_polar(theta_points, r_points, r_coupling, p_coupling, 
                          theta_scale=1.0, r_scale=1.0, Lx=2*np.pi):
    """
    Compute coupling structure based on polar distances with periodic boundaries.
    
    Parameters:
    -----------
    theta_points, r_points : ndarray
        Polar coordinates of oscillators
    r_coupling : float
        Coupling radius in polar metric
    p_coupling : float
        Probability of forming connection within radius
    theta_scale : float
        Scaling factor for angular distance
    r_scale : float
        Scaling factor for radial distance
    Lx : float
        Period length for angular coordinate (default 2π)
        
    Returns:
    --------
    coupling : list of dict
        Each element contains 'neighbors' and 'weights' arrays
    """
    
    n_points = len(theta_points)
    coupling = []
    
    print('Computing coupling matrix (Polar distances)...')
    
    for i in range(n_points):
        neighbors_idx = []
        weights = []
        
        for j in range(n_points):
            if i != j:
                # Distance in polar coordinates (periodic in theta)
                dtheta = abs(theta_points[i] - theta_points[j])
                dtheta_periodic = theta_scale * min(dtheta, Lx - dtheta)
                dr_regular = r_scale * abs(r_points[i] - r_points[j])
                dist = np.sqrt(dtheta_periodic**2 + dr_regular**2)
                
                if dist <= r_coupling and np.random.rand() < p_coupling:
                    neighbors_idx.append(j)
                    weights.append(1.0)
        
        # Normalize weights
        if len(weights) > 0:
            weights = np.array(weights) / np.sum(weights)
            coupling.append({
                'neighbors': np.array(neighbors_idx, dtype=int),
                'weights': weights
            })
        else:
            coupling.append({
                'neighbors': np.array([], dtype=int),
                'weights': np.array([])
            })
        
        if (i + 1) % 500 == 0:
            print(f'  Processed {i+1}/{n_points} points')
    
    n_connections = sum([len(c['neighbors']) for c in coupling])
    print(f'Coupling computed: {n_connections} total connections')
    print(f'Average connections per oscillator: {n_connections/n_points:.1f}')
    
    return coupling


def run_kuramoto_simulation(coupling, omega, dt, T_total, save_interval, 
                           K=1.2, noise_scale=0.0, seed=None, verbose=True):
    """
    Run Kuramoto model simulation with given coupling structure.
    
    The dynamics follow:
    dθ_i/dt = ω_i + K * Σ_j w_ij * sin(θ_j - θ_i) + u_i * I(t)
    
    where u_i * I(t) is a noise term with:
    - u_i: fixed random coefficient for each oscillator (heterogeneous noise sensitivity)
    - I(t): time-varying common noise (affects all oscillators but weighted by u_i)
    
    Parameters:
    -----------
    coupling : list of dict
        Coupling structure from compute_coupling_* functions
    omega : ndarray
        Natural frequencies for each oscillator
    dt : float
        Time step
    T_total : float
        Total simulation time
    save_interval : int
        Save phase every N steps
    K : float
        Coupling strength
    noise_scale : float
        Standard deviation of noise term I(t). If 0, no noise is added.
    seed : int, optional
        Random seed for initial phases and noise
    verbose : bool
        Print progress updates
        
    Returns:
    --------
    phase_history : ndarray
        Array of shape (n_points, n_saves) with phase values
    time_points : ndarray
        Array of time points corresponding to saved phases
    """
    if seed is not None:
        np.random.seed(seed)
    
    n_points = len(coupling)
    num_steps = int(T_total / dt)
    
    # Initialize phases
    theta = 2 * np.pi * np.random.rand(n_points)
    
    # Initialize noise terms (if noise_scale > 0)
    if noise_scale > 0:
        # u: heterogeneous noise sensitivity for each oscillator (fixed over time)
        u = np.random.normal(0, 1, n_points)
        # I: common noise fluctuation over time
        I = np.random.normal(0, noise_scale, num_steps)
        if verbose:
            print(f'Noise enabled: noise_scale={noise_scale:.3f}')
    
    # Convert coupling to sparse matrix
    if verbose:
        print('Converting to sparse matrix format...')
    
    total_connections = sum([len(c['neighbors']) for c in coupling])
    
    # Pre-allocate arrays
    source_idx = np.zeros(total_connections, dtype=int)
    target_idx = np.zeros(total_connections, dtype=int)
    coupling_weights = np.zeros(total_connections)
    
    # Fill arrays
    counter = 0
    for i in range(n_points):
        if len(coupling[i]['neighbors']) > 0:
            n_neighbors = len(coupling[i]['neighbors'])
            source_idx[counter:counter+n_neighbors] = coupling[i]['neighbors']
            target_idx[counter:counter+n_neighbors] = i
            coupling_weights[counter:counter+n_neighbors] = K * coupling[i]['weights']
            counter += n_neighbors
    
    # Create sparse matrix
    W_sparse = csr_matrix((coupling_weights, (target_idx, source_idx)), 
                          shape=(n_points, n_points))
    
    if verbose:
        print(f'Sparse matrix created: {total_connections} connections')
        print(f'Average connections: {total_connections/n_points:.1f}')
    
    # Time integration
    n_saves = num_steps // save_interval
    
    phase_history = np.zeros((n_points, n_saves))
    time_points = np.zeros(n_saves)
    save_counter = 0
    
    if verbose:
        print('Starting simulation...')
    start_time = time.time()
    
    for step in range(num_steps):
        t = step * dt
        
        # Natural frequency term
        dtheta_omega = omega * dt
        
        # Coupling term using sparse matrix
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        
        W_sin = W_sparse.dot(sin_theta)
        W_cos = W_sparse.dot(cos_theta)
        
        coupling_term = W_sin * cos_theta - W_cos * sin_theta
        dtheta_coupling = coupling_term * dt
        
        # Noise term (if enabled)
        if noise_scale > 0:
            dtheta_noise = u * I[step] * dt
        else:
            dtheta_noise = 0
        
        # Update phases
        theta = theta + dtheta_coupling + dtheta_omega + dtheta_noise
        theta = np.mod(theta, 2*np.pi)
        
        # Save phase
        if step % save_interval == 0:
            phase_history[:, save_counter] = theta
            time_points[save_counter] = t
            save_counter += 1
        
        # Progress update
        if verbose and (step + 1) % 1000 == 0:
            elapsed = time.time() - start_time
            print(f'Step {step+1}/{num_steps} (t={t:.1f}), elapsed: {elapsed:.1f}s')
    
    if verbose:
        print('Simulation complete!')
    
    return phase_history, time_points


def compute_order_parameter(phase_history):
    """
    Compute Kuramoto order parameter over time.
    
    Parameters:
    -----------
    phase_history : ndarray
        Phase evolution array (n_points, n_times)
        
    Returns:
    --------
    order_parameter : ndarray
        Order parameter R at each time point
    """
    n_times = phase_history.shape[1]
    order_parameter = np.zeros(n_times)
    
    for i in range(n_times):
        z = np.mean(np.exp(1j * phase_history[:, i]))
        order_parameter[i] = np.abs(z)
    
    return order_parameter


def plot_connectivity_map(x_points, y_points, coupling, r_max, title, pPlot=0.0003):
    """
    Plot connectivity structure in Cartesian space.
    
    Parameters:
    -----------
    x_points, y_points : ndarray
        Cartesian coordinates
    coupling : list of dict
        Coupling structure
    r_max : float
        Domain radius
    title : str
        Plot title
    pPlot : float
        Fraction of connections to display
        
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(8, 8))
    n_points = len(x_points)
    
    # Plot connections
    for i in range(n_points):
        if len(coupling[i]['neighbors']) > 0:
            neighbors_idx = coupling[i]['neighbors']
            thisX = x_points[i]
            thisY = y_points[i]
            
            for neighbor_idx in neighbors_idx:
                if np.random.rand() < pPlot:
                    thatX = x_points[neighbor_idx]
                    thatY = y_points[neighbor_idx]
                    ax.plot([thisX, thatX], [thisY, thatY], 'k-', linewidth=0.5, alpha=0.4)
    
    # Add circle boundary
    circle = plt.Circle((0, 0), r_max, color='r', fill=False, linestyle='--', linewidth=2)
    ax.add_patch(circle)
    
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    return fig, ax


def plot_phase_comparison(x_points, y_points, phase_history_1, phase_history_2,
                         time_points, r_max, title1, title2, n_frames=6):
    """
    Plot side-by-side comparison of phase evolution for two models.
    
    Parameters:
    -----------
    x_points, y_points : ndarray
        Cartesian coordinates
    phase_history_1, phase_history_2 : ndarray
        Phase evolution arrays for two models
    time_points : ndarray
        Time array
    r_max : float
        Domain radius
    title1, title2 : str
        Titles for left and right columns
    n_frames : int
        Number of time frames to display
        
    Returns:
    --------
    fig, axes : matplotlib figure and axes objects
    """
    import matplotlib.pyplot as plt
    
    plot_indices = np.linspace(0, len(time_points)-1, n_frames, dtype=int)
    
    fig, axes = plt.subplots(n_frames, 2, figsize=(12, 4*n_frames))
    
    for row, plot_idx in enumerate(plot_indices):
        t = time_points[plot_idx]
        
        # Model 1 (left column)
        phase_1 = phase_history_1[:, plot_idx]
        scatter1 = axes[row, 0].scatter(x_points, y_points, c=phase_1, 
                                        s=40, cmap='hsv', vmin=0, vmax=2*np.pi)
        circle1 = plt.Circle((0, 0), r_max, color='k', fill=False, linestyle='--', linewidth=1.5)
        axes[row, 0].add_patch(circle1)
        axes[row, 0].set_xlabel('x', fontsize=11)
        axes[row, 0].set_ylabel('y', fontsize=11)
        axes[row, 0].set_title(f'{title1} - t = {t:.1f}', fontsize=12, fontweight='bold')
        axes[row, 0].set_aspect('equal')
        axes[row, 0].grid(True, alpha=0.2)
        plt.colorbar(scatter1, ax=axes[row, 0], label='Phase')
        
        # Model 2 (right column)
        phase_2 = phase_history_2[:, plot_idx]
        scatter2 = axes[row, 1].scatter(x_points, y_points, c=phase_2, 
                                        s=40, cmap='hsv', vmin=0, vmax=2*np.pi)
        circle2 = plt.Circle((0, 0), r_max, color='k', fill=False, linestyle='--', linewidth=1.5)
        axes[row, 1].add_patch(circle2)
        axes[row, 1].set_xlabel('x', fontsize=11)
        axes[row, 1].set_ylabel('y', fontsize=11)
        axes[row, 1].set_title(f'{title2} - t = {t:.1f}', fontsize=12, fontweight='bold')
        axes[row, 1].set_aspect('equal')
        axes[row, 1].grid(True, alpha=0.2)
        plt.colorbar(scatter2, ax=axes[row, 1], label='Phase')
    
    plt.tight_layout()
    return fig, axes


def compute_spiral_index(x_points, y_points, phases, center=(0, 0), direction=1):
    """
    Compute how well the phase pattern matches a perfect spiral.
    
    A perfect spiral has phase = atan2(y - cy, x - cx) + offset.
    The spiral index measures correlation between observed phases and expected spiral phases.
    
    Parameters:
    -----------
    x_points, y_points : ndarray
        Cartesian coordinates of oscillators
    phases : ndarray
        Observed phases at each point
    center : tuple
        (cx, cy) center of the spiral
    direction : int
        1 for counter-clockwise, -1 for clockwise
        
    Returns:
    --------
    spiral_index : float
        Value between 0 and 1, where 1 indicates perfect spiral
    """
    cx, cy = center
    
    # Compute ideal spiral phase pattern
    ideal_spiral = np.arctan2(y_points - cy, x_points - cx)
    ideal_spiral = np.mod(ideal_spiral, 2*np.pi)
    
    if direction == -1:
        ideal_spiral = 2*np.pi - ideal_spiral
    
    # Use circular correlation to measure match
    # Convert phases to complex numbers on unit circle
    observed_complex = np.exp(1j * phases)
    ideal_complex = np.exp(1j * ideal_spiral)
    
    # Correlation: |mean(observed * conj(ideal))|
    correlation = np.abs(np.mean(observed_complex * np.conj(ideal_complex)))
    
    return correlation


def compute_spiral_index_timeseries(x_points, y_points, phase_history, center=(0, 0)):
    """
    Compute spiral index over time for entire phase history.
    
    Parameters:
    -----------
    x_points, y_points : ndarray
        Cartesian coordinates
    phase_history : ndarray
        Phase evolution array (n_points, n_times)
    center : tuple
        Center of spiral
        
    Returns:
    --------
    spiral_indices : ndarray
        Spiral index at each time point
    """
    n_times = phase_history.shape[1]
    spiral_indices = np.zeros(n_times)
    
    for i in range(n_times):
        # Try both directions and take the maximum
        idx_ccw = compute_spiral_index(x_points, y_points, phase_history[:, i], center, direction=1)
        idx_cw = compute_spiral_index(x_points, y_points, phase_history[:, i], center, direction=-1)
        spiral_indices[i] = max(idx_ccw, idx_cw)
    
    return spiral_indices

