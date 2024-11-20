githubdir = '/Users/zhiwenye/Documents/git';                               % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load widefield session table
data_folder = '/Users/zhiwenye/Documents/data';                            % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = '/Users/zhiwenye/Documents/data/figures';                  % specify a folder to save figures
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session info
%% Figure 2.
save_folder = fullfile(figure_folder,'Fig2');
mkdir(save_folder);
h2abc = plotAxonOrientation(data_folder,save_folder);                      % plot axon bias of all senory neurons
h2e = plotBiasHitogram(data_folder,save_folder);                           % plot bias histogram and permutation
h2df = plotAxonFlowMatch2(T,data_folder,save_folder);                      % plot mean optical flow and axon-flow matching index
close all;
%% Extended Data Fig.10
save_folder = fullfile(figure_folder,'FigS10');
mkdir(save_folder);
hs10ad = plotAxonOrientationMO5(data_folder,save_folder);                  % plot axon bias of all neurons in left hemisphere
hs10e = plotAxonFlowMatch_large2_s(T,data_folder,save_folder);             % plot mean optical flow and axon-flow matching index for PPC
hs10f = plotAxonFlowMatch_large_s(T,data_folder,save_folder);              % plot mean optical flow and axon-flow matching index for PPC2
close all;
