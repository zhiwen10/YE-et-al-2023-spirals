function hr1c = plotSpiralDensity_freq(data_folder,save_folder)
%% get spiral density line
freqs = [0.05 0.5;0.5 2;2 8];
for ifreq = 1:3
    freq = freqs(ifreq,:);
    hr1c = plotSpiralDensitySessionsMeanSEM_freq(freq,data_folder,save_folder);
    close all;
end