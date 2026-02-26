function [hr1a,hr1b] = plotScatterDataVsFftn_Freq(T,data_folder,save_folder)
%% spirals raw vs fftn
freqs = [0.1 0.2;0.5 2];
for i = 1:2
    freq = freqs(i,:);
    [hr1a,hr1b] = plotScatterDataVsFftn1(T,data_folder,save_folder,freq);   % plot peak desnity across sessions for data and 3d-fft
    close all;
end