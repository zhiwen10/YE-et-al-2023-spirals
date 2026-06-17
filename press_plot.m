githubdir = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
% addpath(genpath(fullfile(githubdir, 'matplotlib'))); 
%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session table
%% Figure 1.
save_folder = fullfile(figure_folder, 'press');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
makePressVideoExample(data_folder, save_folder);
close all;