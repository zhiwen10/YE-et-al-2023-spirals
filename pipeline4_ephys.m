githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder, 'tables',...
    'spirals_ephys_sessions_new2.csv'));
%% Figure 4g
save_folder = fullfile(data_folder,'ephys','dv_prediction');
getdVPrediction(T,data_folder,save_folder);                                % predict dV from spiking data, cross validated
save_folder = fullfile(data_folder,'ephys','dv_permute');
getdVPredictionPermute(T,data_folder,save_folder);                         % predict dV from spiking data, shuffled
save_folder = fullfile(data_folder,'ephys','roi');
getEphysROI(T,data_folder,save_folder);                                    
save_folder = fullfile(data_folder,'ephys','spirals_raw');
getSpiralsRaw(T,data_folder,save_folder);                                  % detect widefield spirals within the brain ROI
save_folder = fullfile(data_folder,'ephys','spirals_raw_fftn');
getEphysSpiralsGrouping(T,data_folder,save_folder);                        % group spials based on spatiotemporal structure
save_folder = fullfile(data_folder,'ephys','spirals_predict');
getSpiralsPrediction(T,data_folder,save_folder);                           % detect spirals in predicted data, no need to group
save_folder = fullfile(data_folder,'ephys','spirals_compare');
getSpiralComparePredict(T,data_folder,save_folder);                        % assess spiral pairs in raw and predicted data
getSpiralComparePermute(T,data_folder,save_folder);                        % assess spiral pairs in raw and permuted data
%% Figure 4h,i
save_folder = fullfile(data_folder,'ephys','flow_var');
getPhaseFlowMatchingIndex(T,data_folder,save_folder);                      % get matching index for phase and flow
%% Extended Data Fig.13d
save_folder = fullfile(data_folder,'ephys','var_ordered');
getVarOrderedByNeuron(T,data_folder,save_folder);                          % sort variance explained by neuron contribution
