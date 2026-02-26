
%% load widefield session table
data_folder = 'E:\spiral_data_share\data';                            % https://doi.org/10.6084/m9.figshare.27850707
figure_folder = 'E:\spiral_data_share\figures';                  % specify a folder to save figures
%%
save_folder1 = fullfile(data_folder,'task','spirals');
getCorrectSpiralPrePost2(data_folder,save_folder1);
%%
save_folder2 = fullfile(figure_folder,'Fig5');
h1 = plotCorrectSpiralPrePost2(data_folder,save_folder2);
%%
h3 = plotSpiralHemisphere(data_folder,save_folder);
%%
plotSpiralDirection(data_folder,save_folder);
%%
plotSpiralDirection2(data_folder,save_folder);
%% large spirals
h5e = plotSpiralRateAll7(data_folder,save_folder);
%% small spirals
h5e = plotSpiralRateAll8(data_folder,save_folder);
%% all spirals comparison
h5e = plotSpiralRateAll9(data_folder,save_folder);