function getSpiralDensityLine_all(T,data_folder,save_folder)
%% get spiral density lines
freq_all = [0.1 0.2;0.5 2; 2 8];
for ifreq = 1:3
    freq = freq_all(ifreq,:);
    getSpiralsDensityLine_freq(T,freq,data_folder,save_folder);
end