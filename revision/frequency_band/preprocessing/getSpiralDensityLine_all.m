function getSpiralDensityLine_all(T,data_folder,save_folder)
%% get spiral density lines
freq_all = [0.05 0.5;0.5 2; 2 8];
for ifreq = 1
    freq = freq_all(ifreq,:);
    getSpiralsDensityLine_freq(T,freq,data_folder,save_folder);
end