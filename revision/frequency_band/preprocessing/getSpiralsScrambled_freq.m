function getSpiralsScrambled_freq(T,data_folder,save_folder)
%%
freq_all = [0.1 0.2;0.5,2; 2, 8];
for ifreq = 1:3
    freq = freq_all(ifreq,:);
    freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
    save_folder1 = fullfile(save_folder,freq_folder);
    getSpiralGroupingScrambled_freq(T,freq,data_folder,save_folder1);      
end