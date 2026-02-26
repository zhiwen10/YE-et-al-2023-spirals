githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
% load session table
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session info
%% Figure 3.
save_folder = fullfile(figure_folder,'Fig3');
mkdir(save_folder);
h3b = plotCortexDivision(data_folder,save_folder);                         % plot division (AP, hemi)
h3c = plotVarianceExplained(T,data_folder,save_folder);                    % plot variance explained from regression
h3df = plotExampleKernelHEMI(data_folder,save_folder);                     % plot example hemi regression kernel
h3eg = plotExampleKernelAP(data_folder,save_folder);                       % plot example AP regression kernel
h3h = plotKernelMapsHEMI(data_folder,save_folder);                         % plot hemi kernel maps
h3j = plotKernelMapsAP(data_folder,save_folder);                           % plot AP kernel maps
h3i = plotMatchingIndexHEMI(data_folder,save_folder);                      % plot hemi matching index
h3k = plotMatchingIndexAP(data_folder,save_folder);                        % plot AP matching index
close all;
%% Extended Data Fig.11
save_folder = fullfile(figure_folder,'FigS11');
mkdir(save_folder);
[hs11a, hs11b] = plotMapsSession(T,data_folder,save_folder);               % plot kernel maps across sessions
close all;