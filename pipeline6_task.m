%% pipeline task
githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals_20250218')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
%% get psychometric curve and trial outcome
save_folder = fullfile(data_folder,'task','psychometric_curve');
getPsychometricCurve(data_folder,save_folder);                             % get psychometric curves for all subjects
save_folder = fullfile(data_folder,'task','task_outcome');
getTaskTrialOutcome(data_folder, save_folder);                             % get task trial outcome for each trial
%% get mean maps for different trial types
save_folder = fullfile(data_folder,'task','task_mean_maps');
getMeanMapSession(data_folder, save_folder);                               % get mean map for all contrast and trial types across sessions
getMeanMapsAll(data_folder,save_folder);                                   % average mean maps at [-2,2]s around stim onset for 3 trial types 
%% get sprials for different trial types
save_folder = fullfile(data_folder,'task','spirals');
getCorrectSpiralDensity(data_folder,save_folder);
%% get phase and amplitude around stim onset in VISp
freq1 = [0.05,2];
getTaskOnsetPhase(data_folder,save_folder,freq1);                          % get phase and amplitude around onset time
freq2 = [2, 8];
getTaskOnsetPhase(data_folder,save_folder,freq2);                          % get phase and amplitude around onset time
%% spirals during each trial
save_folder = fullfile(data_folder,'task','spirals');
getTaskSpirals(data_folder,save_folder);                                   % task spirals at[-2,2]s around onset time in each trial
getPassiveSpirals(data_folder,save_folder);                                % passive spirals at[-2,2]s around onset time in each trial
getSprialCountByRadiusTime(data_folder,save_folder);                       % sort task spirals by radius and time around stim onset
%%
getPassiveSpiralPrePost(data_folder,save_folder);                          % concatenate all passive spirals before and after stim onset
getCorrectSpiralPrePost(data_folder,save_folder);                          % concatenate all task spirals before and after stim onset