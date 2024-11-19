function getConcatPowerSpectrum(T, data_folder,save_folder)
%% load tapered power sprectrum for all sessions
psdxAllNorm = zeros(1320,1140,29,15);
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals\spirals_power_spectrum2',...
        [fname '_fftSpectrum2.mat']));
    psdxAllNorm(:,:,:,kk) = psdxMeanTransformed;
end
save(fullfile(save_folder,'fftSpectrumAllNew.mat'),'psdxAllNorm','freq');