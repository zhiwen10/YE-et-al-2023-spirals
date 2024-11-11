%% planar wave
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
%%
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
freq = [2,8];
%% use 600s to 635s as epoch examples, 600x35 = 21000 samples
epoch = 21000:22250;
%%
kk = 15;
mn = T.MouseID{kk};
tda = T.date(kk);
en = T.folder(kk);    
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals\svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root); 
dV = [zeros(size(V,1),1) diff(V,[],2)];
% registration
load(fullfile(data_folder,'spirals\rf_tform_8x',[fname '_tform_8x.mat']));   % load atlas transformation matrix tform;
% get flow field from unregistered frame frist, then transform to
% registered, to avoid interpolation problem with circular phase 
U1 = U(1:params.downscale:end,1:params.downscale:end,:);
mimg1 = mimg(1:params.downscale:end,1:params.downscale:end);
U1 = U1./mimg1;
[~,traceAmp1,tracePhase1] = spiralPhaseMap_freq(U1(:,:,1:50),dV(1:50,:),t,params,freq,rate);
% clear U U1 V dV
tracePhase1 = permute(tracePhase1,[3,1,2]);
useGPU = 1;
 tracePhase2 = tracePhase1(21000:22250,:,:);
[vxRaw,vyRaw] = HS_flowfield(tracePhase2,useGPU);
% now let's transform flow field to registered atlas
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
vxRaw1 = permute(vxRaw,[2,3,1]);
vyRaw1 = permute(vyRaw,[2,3,1]);
vxRawt = imwarp(vxRaw1,tform,'OutputView',imref2d(size(BW1)));
vyRawt = imwarp(vyRaw1,tform,'OutputView',imref2d(size(BW1)));
%% get phase map from artlas transformed U space
Ut = imwarp(U1(:,:,1:50),tform,'OutputView',imref2d(size(BW1)));
[~,traceAmp1t,tracePhase1t] = spiralPhaseMap_freq(Ut,dV(1:50,21000:22250),t,params,freq,rate);
tracePhase1t = permute(tracePhase1t,[3,1,2]);
%% save boarder rectangle roi
% h1 = figure('Renderer', 'painters', 'Position', [50 50 800 350]);
% frame = 699;
% ax1 = subplot(1,1,1);
% im1 = imagesc(squeeze(tracePhase1t(frame,:,:)));
% colormap(ax1,colorcet('C06'));
% axis image; axis off;
% set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
% hold on;
% overlayOutlines(coords,params.downscale);
% hold on;
% plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
% roi = drawpolygon(ax1);
% bwroi = createMask(roi);
% save('SSP_MOp_boarder_roi.mat','roi','bwroi');
%%
load('SSP_MOp_boarder_roi.mat');
lineColor = 'k'; lineColor1 = 'w';
hemi = [];
scale3 = 5/8;
h1 = figure('Renderer', 'painters', 'Position', [50 50 800 350]);
frame = 699;
ax1 = subplot(1,3,1);
im1 = imagesc(squeeze(tracePhase1t(frame,:,:)));
colormap(ax1,colorcet('C06'));
axis image; axis off;
set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,params.downscale);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
title('Frame 1');

ax2 = subplot(1,3,2);
im1 = imagesc(squeeze(tracePhase1t(frame+1,:,:)));
colormap(ax2,colorcet('C06'));
axis image; axis off;
set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,params.downscale);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
title('Frame 2');

ax3 = subplot(1,3,3);
im1 = imagesc(squeeze(tracePhase1t(frame,:,:)));
colormap(ax3,colorcet('C06'));
hold on;
vxRaw1a = squeeze(vxRawt(:,:,frame));
vyRaw1a = squeeze(vyRawt(:,:,frame));
vxRaw2a = nan(size(vxRaw1a));
vyRaw2a = nan(size(vyRaw1a));
skip = 3;
zoom_scale = 2;
vxRaw2a(1:skip:end,1:skip:end) = vxRaw1a(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2a(1:skip:end,1:skip:end) = vyRaw1a(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2a(not(BW1)) = nan;
vyRaw2a(not(BW1)) = nan;
imH1Raw4 = quiver(vxRaw2a,vyRaw2a,'k','lineWidth',1,'autoScale','off');
hold on;
axis image; axis off;
set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,params.downscale);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
% mean flow angle in rectangle roi
vxRaw3 = vxRaw1a(bwroi);
vyRaw3 = vyRaw1a(bwroi);
vxRaw3_mean = mean(vxRaw3);
vyRaw3_mean = mean(vyRaw3);
% mean rectangle center
[row,col] = find(bwroi);
row_mean = mean(row);
col_mean = mean(col);
positions = roi.Position;
positions(end+1,:) = positions(1,:);
plot(positions(:,1),positions(:,2),'r');
hold on;
scatter(col_mean,row_mean,10,'w','filled');
hold on;
quiver(col_mean,row_mean,vxRaw3_mean*10,vyRaw3_mean*10,'color','w','LineWidth',1,'AutoScale','off');
title('Flow field');
%%
% mean flow angle in rectangle roi
vxRaw3_mean = zeros(size(vxRawt,3),1);
vyRaw3_mean = zeros(size(vyRawt,3),1);
for frame = 1:size(vxRawt,3)
    vxRaw1a = squeeze(vxRawt(:,:,frame));
    vyRaw1a = squeeze(vyRawt(:,:,frame));
    vxRaw3 = vxRaw1a(bwroi);
    vyRaw3 = vyRaw1a(bwroi);
    vxRaw3_mean(frame,1) = mean(vxRaw3);
    vyRaw3_mean(frame,1) = mean(vyRaw3);
end
%%
figure;
posx = zeros(size(vxRaw3_mean,1),1)+col_mean;
posy = zeros(size(vyRaw3_mean,1),1)+row_mean;
quiver(posx,posx,vxRaw3_mean,vyRaw3_mean,'k');
%%
figure;
ax1 = subplot(1,1,1);
frame = 699;
im1 = imagesc(squeeze(tracePhase1t(frame,:,:)));
colormap(ax1,colorcet('C06'));
set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
overlayOutlines(coords,params.downscale);
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
% mean rectangle center
[row,col] = find(bwroi);
row_mean = mean(row);
col_mean = mean(col);
positions = roi.Position;
positions(end+1,:) = positions(1,:);
hold on;
plot(positions(:,1),positions(:,2),'r');
hold on;
scatter(col_mean,row_mean,10,'w','filled');
hold on;
posx = zeros(size(vxRaw3_mean,1),1)+col_mean;
posy = zeros(size(vyRaw3_mean,1),1)+row_mean;
quiver(posx,posy,vxRaw3_mean*4,vyRaw3_mean*4,'k');
axis image; axis off;
title('Flow field');
%%
vxy_mean = complex(vxRaw3_mean,vyRaw3_mean);
angle_mean = angle(vxy_mean);
%%
figure;
histogram(angle_mean);