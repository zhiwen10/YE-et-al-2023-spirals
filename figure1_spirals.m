githubdir = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
% addpath(genpath(fullfile(githubdir, 'matplotlib'))); 
%% load session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));      % load session table
%% Figure 1.
save_folder = fullfile(figure_folder, 'Fig1');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
h1a = plotSpiralTimeSeries3d(data_folder,save_folder);                    % plot example spiral time series & frame with optical flow
h1c = plotSpiralTimeSeries(data_folder,save_folder);                      % plot example spiral frame
h1b = plotSpiralSequence3(data_folder,save_folder);                        % plot example spiral sequence
h1d = plotSpiralDuration(T,data_folder,save_folder);                       % plot spiral duration ratio vs scrambled distribution
h1e = plotSpiralDensityAllSessions(data_folder,save_folder);               % plot spiral desnity map (combine all sessions)
h1fg = plotSpiralSpeedSummary(T,data_folder,save_folder);                  % plot spiral speeds for all spiral radii
close all;
%% Extended Data Fig.1
save_folder = fullfile(figure_folder, 'FigS1');
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
hs1ab = plotExampleOscillation(data_folder,save_folder);                   % plot example horizontal view and time series 
hs1c = plotPowerSpectrum4(T,data_folder,save_folder);                      % plot mean power map across 15 sessions
hs1d = plotExamplePowerSpectrum2(T,data_folder,save_folder);
hs1e = plotPowerRatio3(data_folder,save_folder);
hs1f = plotPowerRatioRegression3(T,data_folder,save_folder);
%%
[hs1g,hs1h] = plotExampleSpiral2b(T,data_folder,save_folder);              % plot LK_0003 example spirals
[hs1g2] = plotExampleSpiralSpectrum2(T,data_folder,save_folder);           % plot LK_0003 power spectrum
[hs1i,hs1j] = plotExampleSpiral3b(T,data_folder,save_folder);              % plot ZYE_0067 example spirals
[hs1i2] = plotExampleSpiralSpectrum3(T,data_folder,save_folder);           % plot ZYE_0067 power spectrum
% close all;
%% Extended Data Fig.2
save_folder = fullfile(figure_folder, 'FigS2');
mkdir(save_folder);
hs2 = plotSpiralDetectionPipeline(T,data_folder,save_folder);              % plot spiral detection pipeline illustration
% close all;
%% Extended Data Fig.3
save_folder = fullfile(figure_folder, 'FigS3');
mkdir(save_folder);
freq = [2 8];
hs3ab = plotExampleDataVsFft(T,data_folder,save_folder);                   % plot example epoch of data and 3d-fft
hs3c = plotMapDataVsFftn(data_folder,save_folder,freq);                    % plot density map (combine all sessions) for data and 3d-fft
[hs3d,hs3e] = plotScatterDataVsFftn(T,data_folder,save_folder,freq);       % plot peak desnity across sessions for data and 3d-fft
% close all;
%% Extended Data Fig.4
save_folder = fullfile(figure_folder, 'FigS4');
mkdir(save_folder);
[hs4bc,hs4de] = plotLFPspirals(data_folder, save_folder);                  % plot example spirals in cortical LFP
% close all;
%% Extended Data Fig.5
save_folder = fullfile(figure_folder, 'FigS5');
mkdir(save_folder);
session_rows = 1:6;                                                        % plot 6/15 sessions in one figure, to avoid crowding
hs5a1 = plotSpiralsBySession1(T,session_rows,data_folder,save_folder);     % plot visual rf mapping by session
hs5a2 = plotSpiralsBySession2(T,session_rows,data_folder,save_folder);     % plot spiral distribution by session
% close all;
%% Extended Data Fig.6
save_folder = fullfile(figure_folder, 'FigS6');
mkdir(save_folder);
hs6ab = plotSpiralDensityByRadius(data_folder,save_folder);                % plot spiral density maps across different radius
hs6cd = plotSpiralDensityByDuration(data_folder,save_folder);              % plot spiral density maps across different durations
close all;
%% Extended Data Fig.7
save_folder = fullfile(figure_folder, 'FigS7');
mkdir(save_folder);
hs7a = plotExampleSpiralTrajectory(T,data_folder,save_folder);             % plot example grouped sprial sequences
hs7b = plotSpiralDirectionRatio2(T,data_folder,save_folder);               % plot CCW spirals ratio across sessions
hs7cd = plotSpiralDensitySessionsMeanSEM(T,data_folder,save_folder);       % plot spiral density map mean and SEM across sessions
hs7ef = plotSpiralsSymmetryRatio(T, data_folder, save_folder);             % plot spiral symmetry ratio across sessions
close all;
%% Extended Data Fig.8
save_folder = fullfile(figure_folder, 'FigS8');
mkdir(save_folder);
[hs8e, hs8ad,hs8f] = plotSpiralSyncIndex(T, data_folder,save_folder);      % plot example time series and example frame index
hs8gh = plotMotionEnergyAmpX(T,data_folder,save_folder);                   % plot 2-8Hz amp vs Motion energy relationship
hs8i = plotAmpIndex(data_folder,save_folder);                              % plot 2-8Hz amp vs index across sessions
hs8j = plotMotionEnergyIndex(data_folder,save_folder);                     % plot motion energy vs index across sessions
hs8k = plotExamplePlaneWaveSeries2(data_folder,save_folder);               % plot example plane wave
hs8l = plotWaveIndexAmp(data_folder,save_folder);                          % plot plane wave index vs 2-8 Hz amp
hs8m = plotWaveRatio2(data_folder,save_folder);                            % plot plane wave vs sprial wave ratio
hs8n = plotExamplePlaneWave(data_folder,save_folder);                      % plot example plane wave symmetry 
hs8o = plotSymmetry4(data_folder,save_folder);                             % plot plane wave angle distribution across frames
hs8p = plotBorderPlaneWave(data_folder,save_folder);                       % plot plane wave angle distribution on the border
close all;
%% Extended Data Fig.9
save_folder = fullfile(figure_folder, 'FigS9');
mkdir(save_folder);
h9af = plotSpiralSpeedExample(data_folder,save_folder);                    % plot example spiral speed calculation illustration
radius = 50;
h9gj = plotSpeedForRadius(radius,data_folder,save_folder);                 % plot all spiral speeds at 50-pixel radius
radius = 100;
h9kn = plotSpeedForRadius(radius,data_folder,save_folder);                 % plot all spiral speeds at 100-pixel radius
close all;