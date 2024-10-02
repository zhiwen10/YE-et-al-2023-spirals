githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures'; 
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
T1 = T(1,:);
%%
save_folder = fullfile(data_folder, 'spirals\spectrum');
getPowerSpectrumPoints(T,data_folder,save_folder);
%%
save_folder = fullfile(figure_folder, 'Fig_spectrum');
plotPowerSpectrumLog(T,data_folder,save_folder);
%%
save_folder = fullfile(data_folder, 'spirals\spectrum');
freq = [2,8];
getScatterDataVsFftn(T,data_folder,save_folder,freq);
%%
save_folder = fullfile(figure_folder, 'Fig_spectrum');
plotSpiralvsSpectrum(data_folder,save_folder);
%% save downsampled spectrum across sessions
save_folder = fullfile(data_folder, 'spirals\spectrum\sessions');
getPowerSpectrumSession(T,data_folder,save_folder);
%%
save_folder = fullfile(figure_folder, 'Fig_spectrum\sessions');
plotFooofParams(T,data_folder,save_folder);
%% detect alpha epochs and get spiral density for alpha and non-alpha epochs
freq = [2 8];
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_fft_bump');
getAlphaSprialMap(T,freq,data_folder,save_folder);
%% get mean fftn across sessions by alpha epochs and nonalpha epochs
freq = [2 8];
save_folder = fullfile(data_folder, 'spirals\spectrum\alpha_fft');
getAlphafft(T,freq,data_folder,save_folder);
%%
save_folder = fullfile(data_folder, 'spirals\spectrum\example_traces_05_8Hz');
getExamplePixelTrace_05_8Hz(T,data_folder,save_folder);
%%
save_folder = fullfile(data_folder, 'spirals\spectrum\example_traces_01_8Hz');
getExamplePixelTrace_01_8Hz(T,data_folder,save_folder);