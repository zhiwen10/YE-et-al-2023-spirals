function hs3h = plotSpiralDensity_freq2(data_folder,save_folder)
%% get spiral density line
freqs = [0.1 0.2;0.5 2];
for ifreq = 1:2
    freq = freqs(ifreq,:);
    hs3h = plotSpiralDensitySessionsMeanSEM_freq2(freq,data_folder,save_folder);
    close all;
end