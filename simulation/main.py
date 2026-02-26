import numpy as np
import scipy.spatial as ss
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.spatial import Voronoi
import colorcet
import tqdm
import scipy.interpolate as scint


def generate_network_connectivity(Nacross=35, connectivity_type='circular_bias', sigma=1000, ell=0.05, Wscale=1, seed=0, 
                                    yhole=.5, xhole=.5):
    """
    Generate network connectivity matrix with additional features.
    
    Args:
    Nacross (int): Number of nodes across one dimension (total nodes = Nacross^2)
    connectivity_type (str): Type of connectivity ('circular_bias', 'uniform', or 'random')
    sigma (float): Parameter for circular bias
    ell (float): Length scale for exponential decay
    Wscale (float): Scaling factor for weights
    seed (int): Random seed for reproducibility
    xhole (float): x-coordinate of the hole center
    yhole (float): y-coordinate of the hole center
    rhole (float): Radius of the hole
    
    Returns:
    tuple: (W, r) - Weight matrix and node positions
    """
    np.random.seed(seed)
    Ncells = Nacross**2
    # setup spatial connectivity
    r = np.random.normal(0, 0.01, (Ncells, 2))
    xs = np.linspace(0, 1, Nacross)
    ys = np.linspace(0, 1, Nacross) 
    Xs, Ys = np.meshgrid(xs, ys)
    r[:, 0] += Xs.ravel()
    r[:, 1] += Ys.ravel()    
    delpoints = ss.distance.cdist(r, r)

    # crete coupling weights J
    pinhib = 0.001
    scale = 1
    inhib = np.random.choice([0,1], Ncells,pinhib)
    inhib = np.zeros(inhib.shape)
    J = scale*(-1)**inhib[:,None]*np.exp(-delpoints/ell)*(delpoints<(5*ell))
    np.fill_diagonal(J, 0)
    xs = r[:, 0]-xhole # center to 0
    ys = r[:, 1]-yhole # center to 0
    rx = (xs[:, None] + xs[None, :]) / 2
    dx = xs[:, None] - xs[None, :]
    ry = (ys[:, None] + ys[None, :]) / 2
    dy = ys[:, None] - ys[None, :]
    anglebw = np.arctan2(dy, dx) - np.arctan2(ry, rx)
    if connectivity_type == 'circular_bias':
        J = J * (np.cos(anglebw)**2 + sigma*np.sin(anglebw)**2)
    elif connectivity_type == 'uniform':
        J = J * (np.cos(anglebw)**2 + 1*np.sin(anglebw)**2)
        
    W = .001* Wscale * (J / np.mean(J))
    
    return W, r
  


def visualize_network(W, r, title="Network Visualization"):
    """
    Create a graph visualization scatterplot for W and r.
    
    Args:
    W (np.array): Weight matrix
    r (np.array): Node positions
    title (str): Title for the plot
    
    Returns:
    matplotlib.figure.Figure: The created figure
    """
    fig, ax = plt.subplots(figsize=(10, 10))
    cmap = cm.bone_r(np.linspace(0,1,100))#cmap = colorcet.cm.CET_C6(np.linspace(0,1,100))

    # Plot connections
    for i in range(len(r)):
        for j in range(i+1, len(r)):
            if np.abs(W[i, j]) > 0.1:  # Only plot strong connections
                ax.plot([r[i,0], r[j,0]], [r[i,1], r[j,1]], color=cmap[int(99* np.min([W[i,j], 0.3])/0.3), :], alpha = 0.5, lw = 10*np.min([W[i,j], 0.4]))

    # Plot nodes
    ax.scatter(r[:, 0], r[:, 1], c='red', s=100)
    
    ax.set_aspect('equal')
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    
    return fig

def run_network_simulation(W, r, Ncells, num_steps=2000, dt=0.01, psi=None, disorder=0,
                            meanfreq=10, noise=0, seed=1, simulation_type='phase_oscillator',
                            sample_freq = 25, Ninside = 100):
    """
    Run network simulation evolving psi.
    
    Args:
    W (np.array): Weight matrix
    r (np.array): Node positions
    Ncells (int): Number of cells
    numsteps (int): Number of simulation steps
    dt (float): Time step
    psi (np.array): Initial phase values (optional)
    disorder (float): Disorder parameter for frequency distribution
    meanfreq (float): Mean frequency for oscillators
    noise (float): Noise level
    seed (int): Random seed for reproducibility
    sample_freq (int): frequency at which to interpolate grid samples
    Ninside (int): nxn dimension for interpolated phase grid
    
    Returns:
    np.array: psi values over time
    np.array: omega values (frequencies)
    """
    np.random.seed(seed)
    
    if psi is None:
        psi = np.random.random((Ncells, 1)) * 2 * np.pi
    
    omega = np.random.normal(meanfreq, disorder*meanfreq, (Ncells, 1)) # no disorder, all cell use omega at 10
    u = np.random.normal(0, 1, (Ncells, 1))
    I = np.random.normal(0, noise, (num_steps, 1))
    
    ##################
    # settings to interpolate nXn grid of phase  
    ##################
    numsamples = int(num_steps/sample_freq)+1
    xin = np.linspace(0,1, Ninside)
    yin = np.linspace(0,1, int(Ninside))
    Xpt, Ypt = np.meshgrid(xin,yin)
    psitimegrid = np.zeros([Xpt.shape[0], Xpt.shape[1], numsamples])
    cts = 0
    ##################
    
    psi_all = np.zeros((Ncells, num_steps))
    
    for t in tqdm.tqdm(range(num_steps)):
        if simulation_type == 'phase_oscillator':
            psi = psi + dt*(omega + np.sum(W*np.sin((psi - psi[:,None])[:,:,0]), axis=1)[:,None] + u*I[t])
            psi = np.mod(psi, 2*np.pi)
        elif simulation_type == 'linear':
            psi += dt * (omega + np.dot(W, psi) + u * I[t])
            
        psi_all[:, t] = psi.flatten()
        if t%sample_freq==0:
            interp = scint.NearestNDInterpolator((r[:,0], r[:,1]), psi)
            psiinterp = interp(Xpt,Ypt)
            psitimegrid[:,:,cts] = psiinterp[:,:,0]
            cts = cts+1
    
    return psi_all, omega, psitimegrid




def visualize_network_simulation(psi_history, r, plot_every_n_frames=1, save_path=None, x_lim=[0,1], y_lim=[0,1]):
    """
    Visualize network simulation using Voronoi diagrams.
    
    Args:
    psi_history (np.array): History of psi values, shape (Ncells, num_timesteps)
    r (np.array): Node positions, shape (Ncells, 2)
    plot_every_n_frames (int): Plot every n-th frame to reduce computation time
    save_path (str): Path to save figures. If None, figures are displayed but not saved.
    
    Returns:
    list: List of figure objects
    """
    vor = Voronoi(r)
    regions = vor.regions
    vertices = vor.vertices
    cmap = colorcet.cm.CET_C6(np.linspace(0,1,100))
    figures = []
    
    for t in range(0, psi_history.shape[1], plot_every_n_frames):
        fig = plt.figure(figsize=[10,7.5])
        ax = fig.add_subplot(111)
        ax.set_facecolor('k')
        
        for i in range(len(r)):
            region = regions[vor.point_region[i]]
            if region and not any([v == -1 for v in region]):
                polygon = vertices[region]
#                 polygon = np.clip(polygon, 0,1)
                normedphase = int(100 * (psi_history[i, t] / (2*np.pi)))
                ax.fill(*zip(*polygon), alpha=0.4, color=cmap[normedphase,:])
        
        # case where did flipping stuff and want to draw blue line
        if np.min(x_lim) < 0:
            plt.plot([0,0], [0,1.5], '--k')
            plt.plot([-1,1], [1,1], '--k')
            
        plt.ylim(*y_lim)
        plt.xlim(*x_lim)
        
        if save_path:
            plt.savefig(f"{save_path}/frame_{t:05d}.png")
            plt.close(fig)
        else:
            plt.show()
        
        figures.append(fig)
   
    
    return figures


def visualize_network_simulation2(psi_history, r, plot_every_n_frames=1, save_path=None, x_lim=[0,1], y_lim=[0,1]):
    """
    Visualize network simulation using Voronoi diagrams with subplots in a row.

    Args:
    psi_history (np.array): History of psi values, shape (Ncells, num_timesteps)
    r (np.array): Node positions, shape (Ncells, 2)
    plot_every_n_frames (int): Plot every n-th frame to reduce computation time
    save_path (str): Path to save figures. If None, figures are displayed but not saved.
    
    Returns:
    list: List of figure objects

    """
    cmap = colorcet.cm.CET_C6(np.linspace(0,1,100))
    vor = Voronoi(r)
    regions = vor.regions
    vertices = vor.vertices

    fig = plt.figure(figsize=(10,7.5), layout='constrained')
    count1 = 1
    for t in range(0, psi_history.shape[1], plot_every_n_frames):   
        ax = plt.subplot(1, 10, count1)
        ax.set_facecolor('k')
        for i in range(len(r)):
            region = regions[vor.point_region[i]]
            if -1 not in region:
                polygon = [vertices[p] for p in region]
                polygon = [[_el if _el < 1 else 1 for _el in _ar] for _ar in polygon]
                polygon = [[_el if _el > 0 else 0 for _el in _ar] for _ar in polygon]
                normedphase = int(100*(psi_history[i,t]/(2*np.pi)))
                ax.fill(*zip(*polygon), alpha=0.4, color=cmap[normedphase,:])

        plt.ylim(0, 1)
        plt.xlim(0, 1)
        plt.yticks([])
        plt.xticks([])
        ax.set_aspect('equal')
        count1 += 1
      
    plt.show()
    # return fig

def mirror_sim(seed, Xpt, Ypt, mirror = True):

    # generate the random configuration (could expand to hyperuniform or foam like if we like)
    Ncells = 1000
    np.random.seed(seed)
    r = np.random.random([Ncells,2])
    randomshift = 0.001
    r = np.vstack([r,r+np.random.normal(0,randomshift,r.shape),r+np.random.normal(0,randomshift,r.shape),r+np.random.normal(0,randomshift,r.shape)])

    # filp along the x so that we don't have connectivity there.
    # now we can shift 50% of the cells across the mirror in x

      #let's calculate the synpatic matrix
    if mirror:
        delpoints = ss.distance.cdist(r,r)
        ell = 0.05
        pinhib = 0.001
        scale = 1
        inhib = np.random.choice([0,1], Ncells,pinhib)
        inhib = np.zeros(inhib.shape)
        W = scale*np.exp(-delpoints/ell)*(delpoints<(5*ell))
        for i in range(Ncells):
            W[i,i] = 0

        r[Ncells:2*Ncells, 0] = -1*r[Ncells:2*Ncells, 0]
        r[2*Ncells:3*Ncells, 1] = 2-r[2*Ncells:3*Ncells, 1]
        r[3*Ncells:4*Ncells, 1] = 2-r[3*Ncells:4*Ncells, 1]
        r[3*Ncells:4*Ncells, 0] = -r[3*Ncells:4*Ncells, 0]

    else:
        r[Ncells:2*Ncells, 0] = -1*r[Ncells:2*Ncells, 0]
        r[2*Ncells:3*Ncells, 1] = 2-r[2*Ncells:3*Ncells, 1]
        r[3*Ncells:4*Ncells, 1] = 2-r[3*Ncells:4*Ncells, 1]
        r[3*Ncells:4*Ncells, 0] = -r[3*Ncells:4*Ncells, 0]

        delpoints = ss.distance.cdist(r,r)
        ell = 0.05
        pinhib = 0.001
        scale = 1
        inhib = np.random.choice([0,1], Ncells,pinhib)
        inhib = np.zeros(inhib.shape)

        W = scale*np.exp(-delpoints/ell)*(delpoints<(5*ell))
        for i in range(Ncells):
            W[i,i] = 0

    scale2 = 0.003/np.mean(W)
    W = scale2*W

    # initialize a quick simulation
    dt = 0.01
    numsteps = 1000
    # sample = 25 # 5 is about 20Hz sampling
    sample  = 10
    numsamples = int(numsteps/sample)+1
    N = Ncells
    disorder = 0
    meanfreq = 10 #2pi/30 = 0.2 sec or 5Hz
    noise = 0
    psi = np.random.random([4*N,1])*2*np.pi
    omega = np.random.normal(meanfreq,disorder*meanfreq, [4*N,1])
    u = np.random.normal(0,1,[4*N, 1]) # unstructured input for now
    I = np.random.normal(0,0,[numsteps,1])  #noise for now
    # cmap = cm.twilight(np.linspace(0,1,100))
    cmap = colorcet.cm.CET_C6(np.linspace(0,1,100))
    vor = Voronoi(r)
    regions = vor.regions
    vertices = vor.vertices
    psiall = []
    psitimegrid = np.zeros([Xpt.shape[0], Xpt.shape[1], numsamples])
    cts = 0
    for t in range(numsteps):
        psi = psi + dt*(omega + np.sum(W*np.sin(psi - psi[:,None])[:,:,0], axis=1)[:,None] + u*I[t]+ np.random.normal(0, noise, psi.shape))
        psi = np.mod(psi, np.pi*2)
        if len(psiall)>0:
            psiall = np.hstack([psiall,psi])
        else:
            psiall = psi
        if t%sample==0:
            # print('.', end='')
            interp = scint.NearestNDInterpolator((r[:,0], r[:,1]), psi)
            psiinterp = interp(Xpt,Ypt)
            psitimegrid[:,:,cts] = psiinterp[:,:,0]
            cts = cts+1

    return psitimegrid


def calc_sync_spir(psitimegrid, n_across=100):
    """
    calculate synchronization and spiral metrics from phase grid time series
    
    Args:
    psitimegrid (np.array): grid of phases over time
    n_across (int): number of grid points along each dimension
    
    Returns:
    tuple: (sync, spir, anglein, phaseall) - synchronization and spiral metrics, angle field, final phase field
    """
    x = np.linspace(-1,1,n_across)
    y = np.linspace(-1,1,n_across)
    X,Y = np.meshgrid(x,y)
    anglein = np.arctan2(Y,X)
    
    sync, spir = [], []
    for t in range(0, psitimegrid.shape[-1]):
        phaseall = psitimegrid[:,:,t] 
        # phaseall = phaseall - np.pi
        sync.append(np.abs(np.sum(np.exp(1j*phaseall)))/(100**2))
        ccw_spir = np.abs(np.sum(np.exp(1j*(phaseall - anglein)))/(100**2))
        cw_spir = np.abs(np.sum(np.exp(1j*(phaseall + anglein)))/(100**2))
        spir.append(max(ccw_spir, cw_spir))
    return sync, spir, anglein, phaseall
