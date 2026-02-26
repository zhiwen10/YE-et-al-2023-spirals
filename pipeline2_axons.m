githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% data folders 
data_folder = 'E:\spiral_data_share\data'; 
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%% Figure 2c,e
save_folder = fullfile(data_folder,'axons');
getAxonBiasTable(data_folder,save_folder);                            % extract soma and axon info for the 435 sensory neurons 
%% Figure 2d,f
save_folder = fullfile(data_folder,'axons','spirals_70pixels_mean_flow_left');
getSpiralsPhaseMap(T,data_folder,save_folder);                             % save spirals phase maps for all sessions
getSpiralsPhaseMapLeft(T,data_folder,save_folder);
%% axons revision
save_folder = fullfile(data_folder,'revision','axons');
getAxonBiasTableMO(data_folder,save_folder);                               % get table for all cells in left hemisphere 
getAxonBiasTableMO2(data_folder,save_folder);                              % get table for all MO cells
getMOroi(data_folder,save_folder);
getAxonBiasTableSSp2(data_folder,save_folder);
save_folder = fullfile(data_folder,'axons','spirals_100pixels_mean_flow');
getSpiralsPhaseMap2(T,data_folder,save_folder);                            % save spirals phase maps for all sessions