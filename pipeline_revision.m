githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';   
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%% sprials in different frequency band
save_folder = fullfile(data_folder,'spirals\spirals_freq','spirals_fftn');
getSpiralsGroup_freq(T, data_folder,save_folder);
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_scrambled');
getSpiralsScrambled_freq(T,data_folder,save_folder);
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_duration');
getDurationFreq(T,data_folder,save_folder);
save_folder = fullfile(data_folder, 'spirals\spirals_freq\spirals_density_line');
getSpiralDensityLine_all(T,data_folder,save_folder);
save_folder = fullfile(data_folder, 'spirals\spirals_fftn');
getSpiralsFFTnFreq(T, data_folder, save_folder);
%% Figure plane wave
save_folder = fullfile(data_folder,'revision','plane_wave');
getBorderPlaneWave3(data_folder,save_folder);
getPlaneWaveMirror_newsession(data_folder,save_folder);
getPlaneWaveMirror2(data_folder,save_folder);
getSpiralsDensityLine_newsession(data_folder,save_folder);
save_folder = fullfile(data_folder,'revision','plane_wave\rfmap_8x');
getCCFregistration_8x(data_folder,save_folder);
%% axons
save_folder = fullfile(data_folder,'revision','axons');
getAxonBiasTableMO(data_folder,save_folder);                               % get table for all cells in left hemisphere 
getAxonBiasTableMO2(data_folder,save_folder);                              % get table for all MO cells
getMOroi(data_folder,save_folder);
getAxonBiasTableSSp2(data_folder,save_folder);
save_folder = fullfile(data_folder,'axons','spirals_100pixels_mean_flow');
getSpiralsPhaseMap2(T,data_folder,save_folder);                            % save spirals phase maps for all sessions
%% power spectrum
save_folder = fullfile(data_folder, 'spirals\spirals_power_spectrum2');    % power spectrum map for each pixel in each session  
getTaperPowerMap3(T,data_folder,save_folder);                              % calculate tapered power spectrum
%%
save_folder = fullfile(data_folder, 'revision','spectrum');                   % power spectrum map for each pixel in each session 
getExamplePixelTrace_005_8Hz(T,data_folder,save_folder);
getPowerBandRatio(T,data_folder,save_folder);
setAlphaThreshold;

