function hr1 = plotPlanarWave(data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
params.downscale = 8;
BW1 = BW(1:params.downscale:end,1:params.downscale:end);
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
load(fullfile(data_folder,'revision','plane_wave','SSP_MOp_border_roi.mat'));
% load(fullfile(data_folder,'revision','planar_wave','RSP_roi.mat'));
load(fullfile(data_folder,'revision','plane_wave','example_flow_field.mat'));
lineColor = 'k';
hemi = [];
scale3 = 5/8;
hr1 = figure('Renderer', 'painters', 'Position', [50 50 800 350]);
frame = 699;
ax1 = subplot(1,3,1);
im1 = imagesc(squeeze(tracePhase1t(frame,:,:)));
colormap(ax1,colorcet('C06'));
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
title('Example flow field');
%
vxyRaw3_mean = complex(vxRaw3_mean,vyRaw3_mean);
vxyRaw3_mean_angle = angle(vxyRaw3_mean);
%
% mean flow angle in rectangle roi
vxRaw3_mean1 = zeros(size(vxRawt,3),1);
vyRaw3_mean1 = zeros(size(vyRawt,3),1);
for frame = 1:size(vxRawt,3)
    vxRaw1a = squeeze(vxRawt(:,:,frame));
    vyRaw1a = squeeze(vyRawt(:,:,frame));
    vxRaw3 = vxRaw1a(bwroi);
    vyRaw3 = vyRaw1a(bwroi);
    vxRaw3_mean1(frame,1) = mean(vxRaw3);
    vyRaw3_mean1(frame,1) = mean(vyRaw3);
end
ax2 = subplot(1,3,2);
frame = 699;
im1 = imagesc(squeeze(tracePhase1t(frame,:,:)));
colormap(ax2,colorcet('C06'));
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
posx = zeros(size(vxRaw3_mean1,1),1)+col_mean;
posy = zeros(size(vyRaw3_mean1,1),1)+row_mean;
quiver(posx,posy,vxRaw3_mean1*5,vyRaw3_mean1*5,'k');
hold on;
scatter(col_mean,row_mean,10,'w','filled');
hold on;
plot(positions(:,1),positions(:,2),'r');
axis image; axis off;
title('Example mean flow');
%
ax3 = subplot(1,3,3);
load(fullfile(data_folder,'revision','plane_wave','border_angle_mean.mat'));
edges = [-pi:pi/6:pi];
for i = 1:15
    [N(i,:),edges] =histcounts(angle_mean(i,:),edges);
end
N1 = N/size(angle_mean,2);
N_mean = mean(N1,1);
N_sem = std(N1,1)./sqrt(15);
ax3 = subplot(1,3,3);
N_mean(end+1) = N_mean(1);
N_sem(end+1) = N_sem(1);
bar(edges,N_mean,1,'faceColor',[0.5,0.5,0.5]);

% [xb,yb] = stairs(edges,N_mean,'k');
% [xb,yb] = stairs(edges(1:end-1),N_mean,'k');
% im1 = area(xb,yb,0,'faceColor','w','edgeColor','k','showBaseline',1);
hold on;
er = errorbar(edges,N_mean,N_sem);
% er = errorbar(edges+pi/12,N_mean,N_sem);
er.Color = 'k';
er.LineStyle = 'None';
% xlim([-pi,pi]);
xticks([-pi:pi/2:pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
title('Mean flow across mice');
%%
print(hr1, fullfile(save_folder,'FigR1_border_planar_wave.pdf'),...
    '-dpdf', '-bestfit', '-painters');
