githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\cbrewer2'));
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals-test')));       % paper repository
%% Figure 2c,e
save_folder = fullfile(data_folder,'axons');
T1 = getAxonBiasTable(data_folder,save_folder);                            % extract soma and axon info for the 435 sensory neurons 
%% Figure 2d,f
getSpiralsPhaseMap(T,data_folder,save_folder);                             % save spirals phase maps for all sessions