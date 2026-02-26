function getSpiralsFFTnFreq(T, data_folder, save_folder)
freq_all = [0.05 0.5; 0.5 2; 2 8];
label1 = 'control'; label2 = 'fftn';
for i = 1
    freq = freq_all(i,:);
    getSpiralDetectionFftnRaw(T,freq,data_folder, save_folder);            % spiral detection in raw data
    getSpiralDetectionFftnPermute(T,freq,data_folder,save_folder);         % spiral detection in fftn data
    getFFTNSpiralsMap(T,freq,label1,data_folder,save_folder);              % spiral maps in raw data
    getFFTNSpiralsStats(T,freq,label1,data_folder,save_folder);            % spiral stats in raw data
    getFFTNSpiralsMap(T,freq,label2,data_folder,save_folder);              % spiral maps in fftn data
    getFFTNSpiralsStats(T,freq,label2,data_folder,save_folder);            % sprial stats in fftn data
end