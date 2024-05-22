githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\cbrewer2'));
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals-test')));       % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
T1 = T(1,:);
%% Figure 3c
save_folder = fullfile(data_folder, 'spirals_mirror\regression_ap');
getReducedRankRegressionAP(T,data_folder,save_folder);                     % save regression coeffs, kernels and R2 from AP regression
save_folder = fullfile(data_folder, 'spirals_mirror\regression_hemi');
getReducedRankRegressionHEMI(T,data_folder,save_folder);                   % save regression coeffs, kernels and R2 from hemi regression
%% Figure 3h-k
save_folder = fullfile(data_folder, 'spirals_mirror\matching_index');
getExampleKernelAP(T, data_folder,save_folder);                            % save kernel maps for example 8 pixels, AP
getExampleKernelHEMI(T, data_folder,save_folder);                          % save kernel maps for example 8 pixels, hemi
getAxonMapAP(data_folder,save_folder);                                     % get axon projection maps in the anterior cortex (MO)
getAxonMapHEMI(data_folder,save_folder);                                   % get axon projection maps in the left hemipshere
getAxonMapInjection(data_folder,save_folder);                              % get viral injection maps in the sensory cortex
%% Extended Data Fig.11
save_folder = fullfile(data_folder, 'spirals_mirror\regression_kernels');
getMapsSession(T,data_folder,save_folder);                                 % get 8 example kernels for all 15 sessions 