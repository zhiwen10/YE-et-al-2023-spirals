function hr4a = plotPowerSpectrum3(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
%%
pixel(1,:) = [845,835]; % VISp
pixel(2,:) = [775,650]; % RSP
% pixel(3,:) = [590,750]; % SSp-ul
pixel(3,:) = [520,850]; % SSp-ll
pixel(4,:) = [290,700]; % MOs
scale = 8;
pixel = round(pixel/scale);
%% load tapered power sprectrum for all sessions
load(fullfile(data_folder,'spirals','spirals_power_spectrum2','fftSpectrumAllNew3.mat'));
BW = logical(projectedAtlas1(1:scale:end,1:scale:end));
powerRatio1 = zeros(size(BW,1),size(BW,2),size(T,1));
for kk = 1:size(T,1)
    psdxAllNormMean = squeeze(psdxAllNorm(:,:,:,kk));
    psdxAllNormMean2 = reshape(psdxAllNormMean,...
        size(psdxAllNormMean,1)*size(psdxAllNormMean,2),size(psdxAllNormMean,3));
    powerAlpha = bandpower(psdxAllNormMean2',freq,[2,7.9],'psd');
    powerTotal = bandpower(psdxAllNormMean2',freq,[0.05,7.9],'psd');
    powerRatio = powerAlpha./powerTotal;
    powerRatio(powerRatio(:)>1) = nan;
    powerRatio(powerRatio(:)<0) = nan;
    powerRatio1(:,:,kk) = reshape(powerRatio,size(BW));
end
powerRatio2 = mean(powerRatio1,3,'omitnan');
%%
[~,indx1] = min(abs(freq-2));
[~,indx2] = min(abs(freq-8));
psdxAllNormMean = squeeze(mean(psdxAllNorm,4));
% sum power between 2-8Hz

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
for i = 1:4
    hold on;
    scatter(pixel(i,2),pixel(i,1),12,'k','filled');
end
%
subplot(1,2,2);
im1 = imagesc(powerRatio2);
axis image; axis off;
colormap(hot)
colorbar;
set(im1, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,scale,'w');
for i = 1:4
    hold on;
    scatter(pixel(i,2),pixel(i,1),12,'k','filled');
end
axis off; axis image;
caxis([0,0.2]);
cb2 = colorbar;
cb2.Ticks = 0:0.1:0.2;
cb2.TickLabels = {'0','0.1','0.2','0.3','0.4'};
%%
print(hr4a, fullfile(save_folder,'FigR4a_2-8hz_power_map_15mice'),...
    '-dpdf', '-bestfit', '-painters');