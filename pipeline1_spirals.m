githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';     
figure_folder = 'E:\spiral_data_share\figures';
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%% spiral detection algorithm 
save_folder = fullfile(data_folder,'spirals','spirals_raw');                 % run time estimate: 1h for each session
getSpiralDetection(T1,data_folder, save_folder);
%% spatiotemporal clustering of spirals
save_folder = fullfile(data_folder, 'spirals','spirals_grouping');
getSpiralGrouping(T1,data_folder,save_folder);                             % group spirals by spatotemporal proximity    
%% Figure 1d
save_folder = fullfile(data_folder, 'spirals','spirals_scrambled');
getSpiralGroupingScrambled(T,data_folder,save_folder);                     % spirals grouping in frame scrambled data, 10x
save_folder = fullfile(data_folder, 'spirals','spirals_duration');
getSpiralDurationRatio(T,data_folder,save_folder);                         % ratio calculation
%% Figure 1e, Extended Data Fig.7c,d
save_folder = fullfile(data_folder, 'spirals','spirals_density');
getSpiralDensityMap(T,data_folder,save_folder);                            % calculate spiral density               
getSpiralsDensityLine(T,data_folder,save_folder);                          % calculate sprial density lines
%% Extended Data Fig.1c,d,e,f, revision
save_folder = fullfile(data_folder, 'spirals','spirals_power_spectrum2');    % power spectrum map for each pixel in each session  
getTaperPowerMap3(T,data_folder,save_folder);                              % calculate tapered power spectrum
getExamplePixelTrace_005_8Hz(T,data_folder,save_folder);
getPowerBandRatio(T,data_folder,save_folder);
setAlphaThreshold;
%% Extended Data Fig.3 
freq = [2,8]; label1 = 'control'; label2 = 'fftn';
save_folder = fullfile(data_folder, 'spirals','spirals_fftn');
getSpiralDetectionFftnRaw(T,freq,data_folder, save_folder);                % spiral detection in raw data
getSpiralDetectionFftnPermute(T,freq,data_folder,save_folder);             % spiral detection in fftn data
getFFTNSpiralsMap(T,freq,label1,data_folder,save_folder);                  % spiral maps in raw data
getFFTNSpiralsStats(T,freq,label1,data_folder,save_folder);                % spiral stats in raw data
getFFTNSpiralsMap(T,freq,label2,data_folder,save_folder);                  % spiral maps in fftn data
getFFTNSpiralsStats(T,freq,label2,data_folder,save_folder);                % sprial stats in fftn data
%% Extended Data Fig.3 revision: different frequency bands
save_folder = fullfile(data_folder,'spirals','spirals_freq','spirals_fftn');
getSpiralsGroup_freq(T, data_folder,save_folder);
save_folder = fullfile(data_folder, 'spirals','spirals_freq\spirals_scrambled');
getSpiralsScrambled_freq(T,data_folder,save_folder);
save_folder = fullfile(data_folder, 'spirals','spirals_freq\spirals_duration');
getDurationFreq(T,data_folder,save_folder);
save_folder = fullfile(data_folder, 'spirals','spirals_freq\spirals_density_line');
getSpiralDensityLine_all(T,data_folder,save_folder);
save_folder = fullfile(data_folder, 'spirals','spirals_fftn');
getSpiralsFFTnFreq(T, data_folder, save_folder);
%% Extended Data Fig.6
save_folder = fullfile(data_folder, 'spirals','spirals_density_duration');
getSpiralDensityByDuration(T,data_folder,save_folder);                     % spiral maps across durations
save_folder = fullfile(data_folder, 'spirals','spirals_density_radius');     
getSpiralDensityByRadius(T,data_folder,save_folder);                       % sprials maps across radius
%% Extended Data Fig.8
save_folder = fullfile(data_folder, 'spirals','spirals_index');
getMotionEnergyIndex(data_folder,save_folder);                             % bin spirality index based on motion energy
getAmpIndex(data_folder,save_folder);                                      % bin spirality index based on 2-8Hz amplitude
%% Extended Data Fig.9
save_folder = fullfile(data_folder, 'spirals','spirals_speed');
getSpiralSpeed(T,data_folder,save_folder);                                 % calculate spiral speed for each session
getSpiralSpeedConcat(T,data_folder,save_folder);                           % sort spiral speed based on radius for each session
radius = 50;
getSpeedForRadius(T,radius,data_folder,save_folder);                       % get speed for all sprials with radius of 50 pixels
radius = 100;
getSpeedForRadius(T,radius,data_folder,save_folder);                       % get speed for all sprials with radius of 100 pixels





