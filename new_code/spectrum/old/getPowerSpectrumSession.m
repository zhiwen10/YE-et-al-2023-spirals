function getPowerSpectrumSession(T,data_folder,save_folder)
%%
for kk = 1:15
    clear powerV
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals\spirals_power_spectrum2',...
        [fname '_fftSpectrum2.mat']));
    %%
    powerV = psdxMeanTransformed(1:8:end,1:8:end,:);
    %%
    writeNPY(powerV,fullfile(save_folder,[fname '_spectrum.npy']));
    writeNPY(freq,fullfile(save_folder,'frequency.npy'));
end