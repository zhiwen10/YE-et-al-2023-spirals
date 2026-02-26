githubdir = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\git';          % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load widefield session table
data_folder = 'E:\spiral_data_share\data';                                 % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                            % specify a folder to save figures
%% Figure 5.
save_folder = fullfile(figure_folder,'Fig5');
mkdir(save_folder);
h5b = plotCorrectMapsFlow(data_folder,save_folder);                        % plot mean widefield map in correct trials
h5b = getMeanSpiralsDetection(data_folder,save_folder);
h5c = plotCorrectMeanTrace2(data_folder,save_folder);                       % plot mean traces in correct trials
h5d = plotCorrectSpiralDensity(data_folder,save_folder);                   % plot spiral density maps before and after onset
h5e = plotSpiralRateAll(data_folder,save_folder);                          % plot spiral rate change around onset time
h5f = plotMeanTraceExampleSlow(data_folder,save_folder);                   % plot mean traces between 0.05-2Hz
h5i = plotTraceExample2_8Hz(data_folder,save_folder);                      % plot mean traces between 2-8Hz
h5gh = plotHitRateSlow(data_folder,save_folder);                           % Psychometric curve with pre-stim 0.05-2Hz power
h5jk = plotHitRate2_8Hz(data_folder,save_folder);                          % Psychometric curve with pre-stim 2-8Hz phase
close all;
%% Extented Data Fig.14
save_folder = fullfile(figure_folder,'FigS14');
mkdir(save_folder);
hs14a = plotPsychometricCurve(data_folder,save_folder);                    % plot psychometric curve
hs14b = plotTaskSessionEpochExample(data_folder,save_folder);              % plot example session with low-arousal states and miss
hs14cd = plotMeanTraceAcrossContrasts(data_folder,save_folder);            % plot mean responses across contrasts 
hs14efg = plotExamplePhase2_8Hz(data_folder,save_folder);                  % plot example 2-8Hz around onset
close all;
%% Extented Data Fig.15
save_folder = fullfile(figure_folder,'FigS15');
mkdir(save_folder);
hs15bdf = plotMeanMapsAll(data_folder,save_folder);                        % plot mean wf maps in different trial types 
hs15ace = plotMeanTraceAll(data_folder,save_folder);                       % plot traces in different trial types
hs15gh = plotSpiralTrialExample(data_folder, save_folder);                  % plot single trial maps
hs15i = plotPassiveSpiralPrePost(data_folder,save_folder);                 % plot spiral desnity maps in passive viewing
hs15j = plotCorrectSpiralPrePost(data_folder,save_folder);                 % plot spiral density maps in active task
hs15k = plotSpiralCountOverTime(data_folder,save_folder);                  % plot sprial density maps across radius in task
close all;