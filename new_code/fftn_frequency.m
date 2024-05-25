githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\cbrewer2'));
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals-test')));       % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
T1 = T(1,:);
%% pipeline
freq = [4 8]; label1 = 'control'; label2 = 'fftn';
save_folder = fullfile(data_folder, 'spirals\spirals_fftn');
% getSpiralDetectionFftnRaw(T,freq,data_folder, save_folder);                % spiral detection in raw data
% getSpiralDetectionFftnPermute(T,freq,data_folder,save_folder);             % spiral detection in fftn data
getFFTNSpiralsMap(T,freq,label1,data_folder,save_folder);                  % spiral maps in raw data
getFFTNSpiralsStats(T,freq,label1,data_folder,save_folder);                % spiral stats in raw data
getFFTNSpiralsMap(T,freq,label2,data_folder,save_folder);                  % spiral maps in fftn data
getFFTNSpiralsStats(T,freq,label2,data_folder,save_folder);                % sprial stats in fftn data
% Extended Data Fig.3
save_folder = fullfile(figure_folder, 'Fig_fftn');
hs3c = plotMapDataVsFftn(data_folder,save_folder,freq);                    % plot density map (combine all sessions) for data and 3d-fft
[hs3d,hs3e] = plotScatterDataVsFftn(T,data_folder,save_folder,freq);       % plot peak desnity across sessions for data and 3d-fft
close all;
%%
freq = [0.5 2];
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
save_folder = fullfile(data_folder, ...
    'spirals/spirals_freq/spirals_scrambled',freq_folder);
getSpiralGroupingScrambled_freq(T,freq,data_folder,save_folder);              
%% spirals grouping in frame scrambled data, 10x
save_folder = fullfile(data_folder, 'spirals\spirals_duration');
getSpiralDurationRatio(T,data_folder,save_folder);                         % ratio calculation