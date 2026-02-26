githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load widefield session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
T = readtable(fullfile(data_folder, 'tables',...
    'spirals_ephys_sessions_new2.csv'));
%% Figure 4.
save_folder = fullfile(figure_folder,'Fig4');
mkdir(save_folder);
[h4bc,h4de] = plotSpiralPredictionExample(T,data_folder,save_folder);      % plot example spiral prediction
h4f = plotProbeLocation(T,data_folder,save_folder);                        % plot all probe locations in atlas space
h4g = plotSpiralsMatchingRate(T,data_folder,save_folder);                  % plot spiral matching ratio vs shuffle condition
[h4h,h4i] = plotMatchingIndex(T,data_folder,save_folder);                  % plot phase and wave matching index vs shuffle condition 
close all;
%% Extented Data Fig.12
save_folder = fullfile(figure_folder,'FigS12');
mkdir(save_folder);
[hs12a,hs12b] = plotSpiralPredictionExample1(T,data_folder,save_folder);   % plot example spiral prediction from striatal spiking data
[hs12c,hs12d] = plotSpiralPredictionExample2(T,data_folder,save_folder);   % plot example spiral prediction from midbrain spiking data
[hs12e,hs12f] = plotWavePredictionExample1(T,data_folder,save_folder);     % plot example wave prediction from striatal spiking data
[hs12g,hs12h] = plotWavePredictionExample2(T,data_folder,save_folder);     % plot example wave prediction from midbrain spiking data
close all;
%% Extented Data Fig.13
save_folder = fullfile(figure_folder,'FigS13');
mkdir(save_folder);
hs13a = plotVarMap(T,data_folder,save_folder);                             % plot variance explained maps from prediction for all sessions
hs13be = plotVarSummary(T,data_folder,save_folder);                        % plot variance explained summary and relationship with neuron N
hs13f = plotWaveMatchingSession(T,data_folder,save_folder);                % plot wave matching index with 2-8Hz amp for all sessions
hs13g = sort_sprials_by_arousal(T,data_folder,save_folder);
hs13g = plotSpiralsMatchingRate_Arousal(T,data_folder,save_folder);
close all;