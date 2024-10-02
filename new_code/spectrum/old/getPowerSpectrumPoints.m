function getPowerSpectrumPoints(T,data_folder,save_folder)
%%
pixels(1,:) = [845,835]; % VISp
pixels(2,:) = [775,650]; % RSP
pixels(3,:) = [590,750]; % SSp-ul
pixels(4,:) = [520,850]; % SSp-ll
pixels(5,:) = [480,950]; % SSp-m
pixels(6,:) = [550,950]; % SSp-n
pixels(7,:) = [675,905]; % SSp-bfd
area_names = {'VISp','RSP','SSp-ul','SSp-ll','SSp-m','SSp-n','SSp-bfd'};
%%
powerV = nan(29,7,15);
% params
for kk = 1:15
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals\spirals_power_spectrum2',...
        [fname '_fftSpectrum2.mat']));
    pixels_copy = pixels;
    if kk == 3
        pixels_copy(5,:) = [560,920];
        pixels_copy(6,:) = [655,890];
    end
    for i = 1:7
        ax(i) = subplot(2,4,1+i);
        powerV(:,i,kk) = squeeze(psdxMeanTransformed(pixels_copy(i,1),pixels_copy(i,2),:))';
    end
end
writeNPY(powerV,fullfile(save_folder,'powerSpectrum.npy'));
writeNPY(freq,fullfile(save_folder,'frequency.npy'));
