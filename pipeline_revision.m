githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(githubdir, 'spikes'))); 
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\cbrewer2'));
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
%%
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
T1 = T(1,:);
%% Figure revision.
save_folder = fullfile(data_folder,'revision','planar_wave');
h5b = getBorderPlanarWave(data_folder,save_folder);
%% axons
save_folder = fullfile(data_folder,'revision','axons');
getAxonBiasTableMO2(data_folder,save_folder);
%% power spectrum
save_folder = fullfile(data_folder, 'spirals\spirals_power_spectrum2');    % power spectrum map for each pixel in each session  
getTaperPowerMap2(T,data_folder,save_folder);                              % calculate tapered power spectrum
%%
save_folder = fullfile(data_folder, 'spirals\spectrum');                   % power spectrum map for each pixel in each session 
getExamplePixelTrace_005_8Hz(T,data_folder,save_folder);
%%
getPowerBandRatio(T,data_folder,save_folder);
