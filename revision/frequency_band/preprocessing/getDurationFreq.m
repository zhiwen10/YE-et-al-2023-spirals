function getDurationFreq(T,data_folder,save_folder)
%% spirals grouping in frame scrambled data, 10x
freq_all = [0.1,0.2;0.5,2;2, 8];
for ifreq = 1
    freq = freq_all(ifreq,:);
    getSpiralDurationRatio_freq(T,freq,data_folder,save_folder);
end