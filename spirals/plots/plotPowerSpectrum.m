function hs1c = plotPowerSpectrum(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%% load tapered power sprectrum for all sessions
psdxAllNorm = zeros(1320,1140,29,15);
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' td '_' num2str(en)];
    load(fullfile(data_folder,'spirals\spirals_power_spectrum',...
        [fname '_fftSpectrum.mat']));
    maxPower = max(psdxMeanTransformed(:));
    psdxAllNorm(:,:,:,kk) = psdxMeanTransformed/maxPower;
end
%%
[~,indx1] = min(abs(freq-2));
[~,indx2] = min(abs(freq-8));
%%
psdxAllNormMean = squeeze(mean(psdxAllNorm,4));
imax = max(psdxAllNormMean(:));
imin = min(psdxAllNormMean(:));
%% sum power between 2-8Hz
scale = 1;
BW = logical(projectedAtlas1);
powerSum = sum(psdxAllNormMean(:,:,indx1:indx2),3);
hs1c = figure; 
im1 = imagesc(powerSum);
axis image; axis off;
colormap(hot)
colorbar;
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,scale,'w');
axis off; axis image;
%%
print(hs1c, fullfile(save_folder,'Figs1c_2-8hz_power_map_15mice'),...
    '-dpdf', '-bestfit', '-painters');