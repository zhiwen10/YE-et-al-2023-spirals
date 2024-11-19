githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% load session table
data_folder = 'E:\spiral_data_share\data';   
figure_folder = 'E:\spiral_data_share\figures';
%% Figure 5.
save_folder = fullfile(figure_folder,'Fig5');
h5b = plotCorrectMapsFlow(data_folder,save_folder);
h5c = plotCorrectMeanTrace(data_folder,save_folder);
h5d = plotCorrectSpiralDensity(data_folder,save_folder);
h5e = plotSpiralRateAll(data_folder,save_folder);
h5f = plotMeanTraceExampleSlow(data_folder,save_folder);
h5i = plotTraceExample2_8Hz(data_folder,save_folder);
h5gh = plotHitRateSlow(data_folder,save_folder);
h5jk = plotHitRate2_8Hz(data_folder,save_folder);
%% Extented Data Fig.14
save_folder = fullfile(figure_folder,'FigS14');
hs14a = plotPsychometricCurve(data_folder,save_folder);
hs14b = plotTaskSessionEpochExample(data_folder,save_folder);
hs14cd = plotMeanTraceAcrossContrasts(data_folder,save_folder);
hs14efg = plotExamplePhase2_8Hz(data_folder,save_folder);
%% Extented Data Fig.15
save_folder = fullfile(figure_folder,'FigS15');
hs15bdf = plotMeanMapsAll(data_folder,save_folder);
hs15ace = plotMeanTraceAll(data_folder,save_folder);
hs15gh = plotSpiralTrialExample(data_folder, save_folder);
hs15i = plotPassiveSpiralPrePost(data_folder,save_folder); 
hs15j = plotCorrectSpiralPrePost(data_folder,save_folder);
hs15k = plotSpiralCountOverTime(data_folder,save_folder);