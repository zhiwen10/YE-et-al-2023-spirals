githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures'; 
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
T1 = T(1,:);
%% pipeline
freq = [0.2 0.5]; label1 = 'control'; label2 = 'fftn';
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
freqs = [0.1 0.2;2 4; 3 6; 4 8];
for ifreq = 1
    freq = freqs(ifreq,:);
    freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
    save_folder = fullfile(data_folder, ...
        'spirals\spirals_freq\spirals_scrambled',freq_folder);
    getSpiralGroupingScrambled_freq(T,freq,data_folder,save_folder);      
end
%% spirals grouping in frame scrambled data, 10x
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_duration');
freqs = [0.1 0.2;2 4; 3 6; 4 8];
for ifreq = 1
    freq = freqs(ifreq,:);
    getSpiralDurationRatio_freq(T,freq,data_folder,save_folder);
end
%% plot duration of all frequencies
save_folder = fullfile(figure_folder, 'Fig_fftn');
freqs = [0.1 0.2; 0.5 2; 2 8];
h1d = figure('Renderer', 'painters', 'Position', [100 100 800 300]);
for ifreq = 1:3
    ax(ifreq) = subplot(1,3,ifreq);
    freq = freqs(ifreq,:);
    plotSpiralDuration_freq(ax(ifreq),T,freq,data_folder);
end
print(h1d, fullfile(save_folder,'Spiral_duration_scramble_all.pdf'), ...
    '-dpdf', '-bestfit', '-painters');
%% get spiral density maps
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_density');
freqs = [0.2 0.5; 0.5 2; 2 8];
for ifreq = 1
    freq = freqs(ifreq,:);
    getSpiralDensityMap_freq(T,freq,data_folder,save_folder);
end
%% get spiral density lines
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_density_line');
freqs = [0.1 0.2;0.5 2; 2 8];
for ifreq = 3
    freq = freqs(ifreq,:);
    getSpiralsDensityLine_freq(T,freq,data_folder,save_folder);
end
%% plot spiral density maps
save_folder = fullfile(figure_folder, 'Fig_fftn');
freqs = [0.1 0.2;0.5 2;2 8];
for ifreq = 3
    freq = freqs(ifreq,:);
    h1e = plotSpiralDensityAllSessions_freq(freq,data_folder,save_folder);
end
%% get spiral density line
save_folder = fullfile(figure_folder, 'Fig_fftn');
freqs = [0.1 0.2;0.5 2;2 8];
for ifreq = 1:3
    freq = freqs(ifreq,:);
    h1e = plotSpiralDensitySessionsMeanSEM_freq(freq,data_folder,save_folder);
end
%%
save_folder = fullfile(data_folder, 'spirals\spirals_fftn');
label1 = 'control'; label2 = 'fftn';
freqs = [0.2 0.5; 0.5 2; 2 4; 3 6; 4 8];
for ifreq = 1
    freq = freqs(ifreq,:);
    getFFTNSpiralsMap(T,freq,label1,data_folder,save_folder);                  % spiral maps in raw data
    getFFTNSpiralsStats(T,freq,label1,data_folder,save_folder);                % spiral stats in raw data
    getFFTNSpiralsMap(T,freq,label2,data_folder,save_folder);                  % spiral maps in fftn data
    getFFTNSpiralsStats(T,freq,label2,data_folder,save_folder);                % sprial stats in fftn data
end
%% spirals raw vs fftn
save_folder = fullfile(figure_folder, 'Fig_fftn');
freqs = [0.1 0.2;0.5 2; 2 8];
for i = 1
    freq = freqs(i,:);
    save_folder = fullfile(figure_folder, 'Fig_fftn');
    % hs3c = plotMapDataVsFftn(data_folder,save_folder,freq);                    % plot density map (combine all sessions) for data and 3d-fft
    [hs3d,hs3e] = plotScatterDataVsFftn(T,data_folder,save_folder,freq);       % plot peak desnity across sessions for data and 3d-fft
    close all;
end
%%
for i = 1:15
    mouseID = T.MouseID{i};
    save_folder = fullfile(figure_folder, 'Fig_fftn',mouseID);
    if ~exist(save_folder,'dir')        
        mkdir(save_folder);
    end
    freq = [4 8];
    plotSpiralSequence_freq(T,mouseID,freq,data_folder,save_folder);
end
%%
for i = 7
    mouseID = T.MouseID{i};
    save_folder = fullfile(figure_folder, 'Fig_fftn',mouseID);
    if ~exist(save_folder,'dir')        
        mkdir(save_folder);
    end
    freq = [0.5 2];
    plotSpiralSequence_freq(T,mouseID,freq,data_folder,save_folder);
end
%%
for i = 7
    mouseID = T.MouseID{i};
    save_folder = fullfile(figure_folder, 'Fig_fftn/final',mouseID);
    if ~exist(save_folder,'dir')        
        mkdir(save_folder);
    end
    freq = [2 8];
    plotSpiralSequence_freq4(T,mouseID,freq,data_folder,save_folder);
    freq = [0.5 2];
    plotSpiralSequence_freq4(T,mouseID,freq,data_folder,save_folder);
end