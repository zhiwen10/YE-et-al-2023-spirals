githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures';
%%
save_folder = fullfile(data_folder,'whisker','whisker_mean_maps');         
getWhiskerMeanMaps(data_folder,save_folder);                               % get whisker mean maps across 5 mice 
getSpiralsWhiskerMeanMaps(data_folder,save_folder);                        % detect spirals from mean maps
%%
save_folder = fullfile(data_folder,'whisker','spirals_peri_stim');  
getSpiralsPrePost(data_folder,save_folder);                                % concatenate spirals pre and post stimulus
getSpiralsPeriStim(data_folder,save_folder);                               % concatente sprials over time across trials
%%
save_folder = fullfile(data_folder,'whisker','single_trials'); 
getWhiskerSingleTrials(data_folder,save_folder);                           % single trial maps