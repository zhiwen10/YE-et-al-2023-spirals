%% figure revision
githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(githubdir, 'spikes'))); 
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\cbrewer2'));
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%% Figure revision plane wave
save_folder = fullfile(figure_folder,'Revision','plane_wave');
hr1a = plotExamplePlaneWaveSeries2(data_folder,save_folder);
hr1b = plotWaveIndexAmp(data_folder,save_folder);
hr1c = plotWaveRatio2(data_folder,save_folder);
hr1d = plotExamplePlaneWave(data_folder,save_folder);
hr1e = plotSymmetry3(data_folder,save_folder);
hr1h = plotBorderPlaneWave(data_folder,save_folder);
%%
save_folder = fullfile(figure_folder,'Revision','axons');
[hr2,p] = plotAxonOrientationMO3(data_folder,save_folder);
%% power spectrum
save_folder = fullfile(figure_folder,'Revision','power_spectrum');
hr3a = plotPowerSpectrum3(T,data_folder,save_folder);                       % plot mean power map across 15 sessions
hr3b = plotExamplePowerSpectrum(T,data_folder,save_folder);
hr3c = plotPowerRatio(data_folder,save_folder);
hr3d = plotAlphaSpectrum(T,data_folder,save_folder);                       % plot alpha epoch bumps and ratio
%%
[hr3e,hr3f] = plotPowerSpectrumParams(T,data_folder,save_folder);