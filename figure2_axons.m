githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load widefield session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session info
%% Figure 2.
save_folder = fullfile(figure_folder,'Fig2');
h2abc = plotAxonOrientation(data_folder,save_folder);                      % plot axon bias of all senory neurons
h2e = plotBiasHitogram(data_folder,save_folder);                           % plot bias histogram and permutation
h2df = plotAxonFlowMatch2(T,data_folder,save_folder);                       % plot mean optical flow and axon-flow matching index
close all;
%% Extended Data Fig.10
save_folder = fullfile(figure_folder,'FigS10');
[hs10ad] = plotAxonOrientationMO5(data_folder,save_folder);
hs10e = plotAxonFlowMatch_large2_s(T,data_folder,save_folder);
hs10f = plotAxonFlowMatch_large_s(T,data_folder,save_folder);
