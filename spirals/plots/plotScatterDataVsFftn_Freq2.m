function [hs3f,hs3g] = plotScatterDataVsFftn_Freq2(T,data_folder,save_folder)
%% spirals raw vs fftn
freqs = [0.1 0.2;0.5 2];
for i = 1:2
    freq = freqs(i,:);
    [hs3f,hs3g] = plotScatterDataVsFftn2(T,data_folder,save_folder,freq);   % plot peak desnity across sessions for data and 3d-fft
    close all;
end