githubdir = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\git';          % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
%% Figure 5.
save_folder = fullfile(figure_folder, 'Fig5');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
h5b = plotWhiskerMeanMap(data_folder,save_folder);                         % plot mean maps across 5 mice
h5c = plotWhiskerMeanTrace(data_folder,save_folder);                       % lot mean traces around SSp center
[h5d,h5f] = plotDensityPrePost(data_folder,save_folder);                   % plot spiral density map pre/post stim
h5e = plotSpiralRatePeriStimulus(data_folder,save_folder);                 % plot spiral density over time peri stim
%%
save_folder = fullfile(figure_folder, 'FigS15');
hs15a = plotWhiskerSingleTrials(data_folder,save_folder);                  % plot example trials with spirals
