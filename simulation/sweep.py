import matplotlib.pyplot as plt
from main import generate_network_connectivity, run_network_simulation, visualize_network, visualize_network_simulation
import numpy as np
from scipy.spatial import Voronoi
import tqdm
import scipy.interpolate as scint
import multiprocessing
import time
import pickle
from enum import Enum
from dataclasses import dataclass
from typing import List, Dict, Any, Optional
from pathlib import Path
import os

class SweepType(Enum):
    """
    enumeration of different types of parameter sweeps
    
    Args:
    none
    
    Returns:
    none
    """
    NOISE_WEIGHT = "noise_weight"
    SIGMA_WEIGHT = "sigma_weight"
    SEED = "seed"

@dataclass
class SweepConfig:
    """
    configuration class for parameter sweeps
    
    Args:
    sweep_type (SweepType): type of parameter sweep to perform
    grid_size (int): number of points along each dimension in parameter grid
    n_across (int): number of nodes across one dimension of network
    num_steps (int): number of simulation timesteps
    noise_range (tuple): range of noise values to sweep
    output_dir (str): directory to save results
    connectivity_types (List[str]): list of connectivity types to test
    seeds (List[int]): list of random seeds for network initialization
    wscales (List[float]): list of weight scaling factors
    sigmas (List[float]): list of sigma values for circular bias
    seed_range (tuple): range of seeds to search
    default_wscale (float): default weight scale for seed search
    default_sigma (float): default sigma for seed search
    
    Returns:
    none
    """
    sweep_type: SweepType
    grid_size: int = 10
    n_across: int = 32
    num_steps: int = 3000
    noise_range: tuple = (0, 30)
    output_dir: str = 'outputs/combined_sweeps'
    connectivity_types: List[str] = None
    seeds: Optional[List[int]] = None
    wscales: Optional[List[float]] = None
    sigmas: Optional[List[float]] = None
    seed_range: Optional[tuple] = None
    default_wscale: float = 1.0
    default_sigma: float = 50.0

    def __post_init__(self):
        """
        validate configuration parameters and create output directory
        
        Args:
        none
        
        Returns:
        none
        """
        if self.connectivity_types is None:
            self.connectivity_types = ['circular_bias', 'uniform']
            
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        
        if self.sweep_type in [SweepType.NOISE_WEIGHT, SweepType.SIGMA_WEIGHT]:
            if not self.seeds or not self.wscales or len(self.seeds) != len(self.wscales):
                raise ValueError("NOISE_WEIGHT and SIGMA_WEIGHT sweeps require matching seeds and wscales lists")
            
            if self.sweep_type == SweepType.SIGMA_WEIGHT and (not self.sigmas or len(self.sigmas) != len(self.seeds)):
                raise ValueError("SIGMA_WEIGHT sweep requires sigmas list of same length as seeds")
                
        elif self.sweep_type == SweepType.SEED:
            if not self.seed_range or len(self.seed_range) != 2:
                raise ValueError("SEED sweep requires seed_range tuple (start, end)")

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
        phaseall = phaseall - np.pi
        sync.append(np.abs(np.sum(np.exp(1j*phaseall)))/(100**2))
        ccw_spir = np.abs(np.sum(np.exp(1j*(phaseall - anglein)))/(100**2))
        cw_spir = np.abs(np.sum(np.exp(1j*(phaseall + anglein)))/(100**2))
        spir.append(max(ccw_spir, cw_spir))
    return sync, spir, anglein, phaseall

def run_spiral_worker(params: Dict[str, Any]) -> tuple:
    """
    run single spiral simulation with given parameters
    
    Args:
    params (dict): dictionary of simulation parameters including:
        w_scale (float): weight scaling factor
        seed (int): random seed
        n_across (int): number of nodes across
        mode (str): connectivity type
        sigma (float): circular bias parameter
        noise (float): noise strength
        num_steps (int): number of timesteps
    
    Returns:
    tuple: (seed, spiral_metric) for seed search or spiral_metric for other sweeps
    """
    w_scale = params['w_scale']
    seed = params['seed']
    n_across = params['n_across']
    mode = params['mode']
    sigma = params.get('sigma', 50)
    noise = params.get('noise', 0)
    
    W, r = generate_network_connectivity(
        Nacross=n_across, 
        connectivity_type=mode, 
        Wscale=w_scale, 
        sigma=sigma, 
    )
    
    psi_all, omega, psitimegrid = run_network_simulation(
        W, r, n_across*n_across, 
        num_steps=params.get('num_steps', 3000), 
        seed=seed, 
        noise=noise
    )
    
    psitimegrid = psitimegrid[:,:,:-1]
    sync, spir, _, _ = calc_sync_spir(psitimegrid)
    
    if params.get('sweep_type') == SweepType.SEED:
        return params['seed'], np.mean(spir[-10:])
    else:
        return np.mean(spir[-10:])

def generate_parameter_sets(config: SweepConfig) -> Dict[str, List[Dict[str, Any]]]:
    """
    generate parameter combinations for sweep based on configuration
    
    Args:
    config (SweepConfig): sweep configuration object
    
    Returns:
    dict: parameter sets organized by connectivity type
    """
    parameter_sets = {}
    
    for mode in config.connectivity_types:
        spiral_paramsets = []
        if config.sweep_type == SweepType.NOISE_WEIGHT:
            noise_sweep = np.linspace(*config.noise_range, config.grid_size)
            for seed_idx, seed in enumerate(config.seeds):
                for i, noise in enumerate(noise_sweep):
                    w_scale_sweep = np.linspace(0, config.wscales[seed_idx], config.grid_size)
                    for j, w_scale in enumerate(w_scale_sweep):
                        params = {
                            'noise': noise,
                            'w_scale': w_scale,
                            'seed': seed,
                            'n_across': config.n_across,
                            'i': i,
                            'j': j,
                            'mode': mode,
                            'num_steps': config.num_steps
                        }
                        spiral_paramsets.append(params)
        elif config.sweep_type == SweepType.SIGMA_WEIGHT:
            for seed_idx, seed in enumerate(config.seeds):
                w_scale_sweep = np.linspace(0, config.wscales[seed_idx], config.grid_size)
                for i, w_scale in enumerate(w_scale_sweep):
                    sigma_sweep = np.linspace(0, config.sigmas[seed_idx], config.grid_size)
                    for j, sigma in enumerate(sigma_sweep):
                        params = {
                            'sigma': sigma,
                            'w_scale': w_scale,
                            'seed': seed,
                            'n_across': config.n_across,
                            'i': i,
                            'j': j,
                            'mode': mode,
                            'num_steps': config.num_steps
                        }
                        spiral_paramsets.append(params)
        elif config.sweep_type == SweepType.SEED:
            start_seed, end_seed = config.seed_range
            for seed in range(start_seed, end_seed):
                params = {
                    'seed': seed,
                    'w_scale': config.default_wscale,
                    'sigma': config.default_sigma,
                    'n_across': config.n_across,
                    'mode': mode,
                    'num_steps': config.num_steps,
                    'sweep_type': SweepType.SEED
                }
                spiral_paramsets.append(params)
                
        parameter_sets[mode] = spiral_paramsets
    
    return parameter_sets

def run_sweep(config: SweepConfig):
    """
    run parameter sweep with given configuration
    
    Args:
    config (SweepConfig): sweep configuration object
    
    Returns:
    none - saves results to files in output directory
    """
    parameter_sets = generate_parameter_sets(config)
    # iterate over circular / iso, for loop len should be small
    for mode, paramsets in parameter_sets.items():
        print(f"Running sweep for {mode} connectivity")
        
        start = time.time()
        with multiprocessing.Pool(20) as pool:
            results = pool.map(run_spiral_worker, paramsets)
        end = time.time()
        print(f'Sweep took: {end - start:.2f} seconds')
        
        os.makedirs(config.output_dir, exist_ok=True)
        base_filename = f"{config.sweep_type.value}_{mode}"

        np.save(f'{config.output_dir}/results_{base_filename}.npy', np.array(results))
        with open(f'{config.output_dir}/paramsets_{base_filename}.pkl', 'wb') as f:
            pickle.dump(paramsets, f)

def sweep_report(seed_search_config, criterion='sigma'):
    circ_spir = np.load(os.path.join(seed_search_config.output_dir, 'results_seed_circular_bias.npy' ))[:,1]
    iso_spir = np.load(os.path.join(seed_search_config.output_dir, 'results_seed_uniform.npy' ))[:,1]
    with open(os.path.join(seed_search_config.output_dir,'paramsets_seed_circular_bias.pkl'),'rb') as f:
        circ_params = pickle.load(f)
    with open(os.path.join(seed_search_config.output_dir,'paramsets_seed_uniform.pkl'),'rb') as f:
        iso_params = pickle.load(f)

    # SELECTION CRITERIA FOR NOISE SWEEP
    if criterion == 'noise':
        valid_seed_inds = np.where( (iso_spir > .65) & (circ_spir > .65))[0].astype(int)

    # SELECTION CRITERIA FOR SIGMA 
    if criterion == 'sigma':
        valid_seed_inds = np.where( (iso_spir > .4) & (circ_spir > .6))[0].astype(int)
    valid_seeds = [param['seed'] for param in np.array(circ_params)[valid_seed_inds]]
    print(f"List of valid seeds for {seed_search_config.sweep_type}: {valid_seeds}")


if __name__ == '__main__':
    multiprocessing.freeze_support()
    
    seed_search_config = SweepConfig(
        sweep_type=SweepType.SEED,
        seed_range=(0, 4),
        default_wscale=1.0,
        default_sigma=50.0,
        n_across=32,
        num_steps=3000,
        output_dir='outputs/seed_search',
    )
    
    noise_weight_config = SweepConfig(
        sweep_type=SweepType.NOISE_WEIGHT,
        seeds=[22, 31, 89, 100, 148, 264],
        wscales=[7, 10, 10, 10, 5, 8],
        grid_size=10,
        n_across=32,
        output_dir='outputs/noise_weight_sweeps'
    )
    
    sigma_weight_config = SweepConfig(
        sweep_type=SweepType.SIGMA_WEIGHT,
        seeds=[251, 252, 313, 364, 377, 390],
        wscales=[10, 10, 10, 10, 10, 10],
        sigmas=[10, 10, 10, 10, 10, 10],
        grid_size=10,
        n_across=32,
        output_dir='outputs/sigma_weight_sweeps'
    )
    
    run_sweep(seed_search_config)
    sweep_report(seed_search_config, criterion='sigma')
    run_sweep(noise_weight_config)
    run_sweep(sigma_weight_config)