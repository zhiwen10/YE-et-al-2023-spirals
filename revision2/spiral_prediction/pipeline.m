githubdir = '/home/dome/Documents/git';                                    % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%% prediction use right SSp
save_folder = fullfile(data_folder,'revision2','V_prediction'); 
getdVpredictionWF(T,data_folder,save_folder);
%% prediction use all of right hemisphere
save_folder = fullfile(data_folder,'revision2','V_prediction2'); 
getdVpredictionWF2(T,data_folder,save_folder);
%% prediction use right MO
save_folder = fullfile(data_folder,'revision2','V_prediction_MO'); 
getdVpredictionWF_MO(T,data_folder,save_folder);
%%
save_folder = fullfile(data_folder,'revision2','V_permute'); 
getdVpredictionWFPermute(T,data_folder,save_folder);
%%
save_folder = fullfile(data_folder,'revision2','spirals_predict');         % run time estimate: 1h for each session
getSpiralDetectionPredict(T,data_folder, save_folder);
%% spatiotemporal clustering of spirals
save_folder = fullfile(data_folder,'revision2','spirals_grouping');
getSpiralGrouping2(T,data_folder,save_folder);                             % group spirals by spatotemporal proximity  
%% compare spirals
save_folder = fullfile(data_folder,'revision2','spirals_compare');
getSpiralComparePredictWF(T,data_folder,save_folder);
getSpiralComparePredictWFPermute(T,data_folder,save_folder);
%% compare maps
compareSpiralMaps(T,data_folder,save_folder);
%% prediction use both MO
save_folder = fullfile(data_folder,'revision2','V_prediction_MO_both'); 
getdVpredictionWF_MO_both(T,data_folder,save_folder);
save_folder = fullfile(data_folder,'revision2','V_permute_MO_both'); 
getdVpredictionWFPermuteMOboth(T,data_folder,save_folder);
%% spatiotemporal clustering of spirals
save_folder = fullfile(data_folder,'revision2','MO_predict','spirals_grouping');
getSpiralGrouping2(T,data_folder,save_folder);                             % group spirals by spatotemporal proximity
%% compare spirals
save_folder = fullfile(data_folder,'revision2','MO_predict','spirals_compare');
getSpiralComparePredictWF_MOboth(T,data_folder,save_folder);
getSpiralComparePredictWFPermute_MOboth(T,data_folder,save_folder);
%%
spirals_match_rate;