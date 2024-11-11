function hr1d = plotSpiralDuration_freq_all(T,data_folder, save_folder)
%% plot duration of all frequencies
freqs = [0.1 0.2; 0.5 2; 2 8];
hr1d = figure('Renderer', 'painters', 'Position', [100 100 800 300]);
for ifreq = 1:3
    ax(ifreq) = subplot(1,3,ifreq);
    freq = freqs(ifreq,:);
    plotSpiralDuration_freq(ax(ifreq),T,freq,data_folder);
end
print(hr1d, fullfile(save_folder,'FigR1d_Spiral_duration_scramble_all.pdf'), ...
    '-dpdf', '-bestfit', '-painters');