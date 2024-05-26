githubdir = 'dome/Documents/git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'npy-matlab')));                       % https://github.com/kwikteam/npy-matlab
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = '/dome/Downloads/data';     
figure_folder = '/dome/Downloads/figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
T1 = T(1,:);
%% pipeline
freq = [4 8]; label1 = 'control'; label2 = 'fftn';
save_folder = fullfile(data_folder, 'spirals/spirals_fftn');
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
freqs = [0.5 2;2 4; 3 6; 4 8];
for ifreq = 1:4
    freq = freqs(ifreq,:);
    freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
    save_folder = fullfile(data_folder, ...
        'spirals/spirals_freq/spirals_scrambled',freq_folder);
    getSpiralGroupingScrambled_freq(T,freq,data_folder,save_folder);      
end
%%
freqs = [0.5 2;2 4; 3 6; 4 8];
save_folder = fullfile(data_folder,'spirals/spirals_freq','spirals_fftn');
for ifreq = 1:4
    freq = freqs(ifreq,:);
    freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
    save_folder1 = fullfile(save_folder, freq_folder);
    mkdir(save_folder1);
    for kk = 1:size(T,1)
        mn = T.MouseID{kk};
        tda = T.date(kk);
        en = T.folder(kk);
        td = datestr(tda,'yyyy-mm-dd');
        tdb = datestr(tda,'yyyymmdd');
        %% load raw detected spirals
        fname = [mn '_' tdb '_' num2str(en)];
        load(fullfile(data_folder,'spirals/spirals_freq/raw',...
            freq_folder,[fname '_spirals_all.mat']));
        filteredSpirals = pwAll(pwAll(:,3)>=40,:);                             % only use sprials with radius >40 pixels, based on 3d-fft
        %% temporal grouping
        filteredSpirals =unique(filteredSpirals, 'rows');                      % get rid of duplication, in case any
        filteredSpirals = sortrows(filteredSpirals,5);                         % sort based on frame number
        [archiveCell,test_stats] = getGroupingAlgorithm(filteredSpirals);       % main grouping algorithm
        save(fullfile(save_folder1,[fname '_spirals_group_fftn.mat']),...       % save all archived Cells  
            'archiveCell');
    end
end
%% spirals grouping in frame scrambled data, 10x
save_folder = fullfile(data_folder, 'spirals/spirals_duration');
getSpiralDurationRatio(T,data_folder,save_folder);                         % ratio calculation