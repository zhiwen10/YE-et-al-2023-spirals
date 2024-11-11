githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(githubdir, 'spikes'))); 
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\cbrewer2'));
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
%% Figure revision plane wave
save_folder = fullfile(data_folder,'revision','plane_wave');
getPlaneWaveRoi(data_folder,save_folder);
getBrainMask(data_folder,save_folder);
getPlaneWaveMirror(data_folder,save_folder);
getPlaneWaveMirrorAmp(data_folder,save_folder);
getBorderPlaneWave2(data_folder,save_folder);
%% Figure revision new g8 sessions
save_folder = fullfile(data_folder,'revision','spirals_new\rfmap_8x');
getCCFregistration_8x(data_folder,save_folder);
%%
save_folder = fullfile(data_folder,'revision','spirals_new');
% getSpiralsDensityLine_newsession(data_folder,save_folder);
getPlaneWaveMirror_newsession(data_folder,save_folder);
%% Figure revision axons
save_folder = fullfile(data_folder,'revision','axons');
getAxonBiasTableMO(data_folder,save_folder);
getAxonBiasTableMO2(data_folder,save_folder);
getMOBorderLineroi(data_folder,save_folder);
getMOroi(data_folder,save_folder);
getMOLineroi(data_folder,save_folder);
getSSproi(data_folder,save_folder);
%%