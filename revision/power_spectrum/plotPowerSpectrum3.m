function hr4a = plotPowerSpectrum3(T,data_folder,save_folder)
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
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals\spirals_power_spectrum2',...
        [fname '_fftSpectrum2.mat']));
    psdxAllNorm(:,:,:,kk) = psdxMeanTransformed;
end
%%
BW = logical(projectedAtlas1);
powerRatio1 = zeros(size(BW,1),size(BW,2),size(T,1));
for kk = 1:size(T,1)
    psdxAllNormMean = squeeze(psdxAllNorm(:,:,:,kk));
    psdxAllNormMean2 = reshape(psdxAllNormMean,...
        size(psdxAllNormMean,1)*size(psdxAllNormMean,2),size(psdxAllNormMean,3));
    powerAlpha = bandpower(psdxAllNormMean2',freq,[2,7.9],'psd');
    powerTotal = bandpower(psdxAllNormMean2',freq,[0.3,7.9],'psd');
    powerRatio = powerAlpha./powerTotal;
    powerRatio(powerRatio(:)>1) = nan;
    powerRatio(powerRatio(:)<0) = nan;
    powerRatio1(:,:,kk) = reshape(powerRatio,size(BW));
end
%%
powerRatio2 = mean(powerRatio1,3,'omitnan');
%%
[~,indx1] = min(abs(freq-2));
[~,indx2] = min(abs(freq-8));
psdxAllNormMean = squeeze(mean(psdxAllNorm,4));
% sum power between 2-8Hz
scale = 1;
powerAlpha = sum(psdxAllNormMean(:,:,indx1:indx2),3);
hr4a = figure('Renderer', 'painters', 'Position', [100 100 900 300]);
subplot(1,2,1);
im1 = imagesc(powerAlpha);
axis image; axis off;
colormap(hot)
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,scale,'w');
axis off; axis image;
cb1 = colorbar;
%
subplot(1,2,2);
im1 = imagesc(powerRatio2);
axis image; axis off;
colormap(hot)
colorbar;
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,scale,'w');
axis off; axis image;
caxis([0,0.4]);
cb2 = colorbar;
cb2.Ticks = 0:0.1:0.4;
cb2.TickLabels = {'0','0.1','0.2','0.3','0.4'};
%%
print(hr4a, fullfile(save_folder,'FigR4a_2-8hz_power_map_15mice'),...
    '-dpdf', '-bestfit', '-painters');