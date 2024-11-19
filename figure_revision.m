githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%% spiral in different frequency band
save_folder = fullfile(figure_folder,'Revision','frequency_band');
[hr1a,hr1b] = plotScatterDataVsFftn_Freq(T,data_folder,save_folder);
hr1c = plotSpiralDensity_freq(data_folder,save_folder);
%%  plane wave
save_folder = fullfile(figure_folder,'Revision','plane_wave');
hr2a = plotExamplePlaneWaveSeries2(data_folder,save_folder);
hr2b = plotWaveIndexAmp(data_folder,save_folder);
hr2c = plotWaveRatio2(data_folder,save_folder);
hr2d = plotExamplePlaneWave(data_folder,save_folder);
hr2e = plotSymmetry4(data_folder,save_folder);
hr2f = plotBorderPlaneWave(data_folder,save_folder);
%% axons
save_folder = fullfile(figure_folder,'Revision','axons');
hr3 = plotAxonOrientationMO4(data_folder,save_folder);
hr3e = plotAxonFlowMatch_large2(T,data_folder,save_folder);
hr3f = plotAxonFlowMatch_large(T,data_folder,save_folder);
%% power spectrum1
save_folder = fullfile(figure_folder,'Revision','power_spectrum');
hr4a = plotPowerSpectrum3(T,data_folder,save_folder);                      % plot mean power map across 15 sessions
hr4b = plotExamplePowerSpectrum(T,data_folder,save_folder);
%% power spectrum2
hr5a = plotPowerRatio2(data_folder,save_folder);
hr5beh = plotPowerRatioRegression(T,data_folder,save_folder);
hr5fg = plotExampleAlphaEpoch(T,data_folder,save_folder);