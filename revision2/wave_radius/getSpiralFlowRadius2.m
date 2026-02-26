githubdir = 'C:\Users\Steinmetz lab\Documents\git';                        % folder where repositories are hosted 
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals')));            % paper repository
%% data folders 
data_folder = 'E:\spiral_data_share\data'; 
T = readtable(fullfile(data_folder,'tables','spiralSessions3.xlsx'));
%%
save_folder = fullfile(data_folder,'revision2','spiral_flow_radius');
radius_all = [40:10:100];
for k = 1:7
    radius = radius_all(k);
    getSpiralsPhaseMapRadius(T,data_folder,save_folder, radius);
end
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
colora2= colorcet('C06','N',180);
BW = logical(projectedAtlas1);
root1 = '/997/';
ctx = '/997/8/567/688/';
scale = 8;
BW1 = BW(1:scale:end,1:scale:end);
% mask and Kernel regression map for right and left
BW1_left = BW1;
BW1_left(:,72:end) = 0;
th2 = 1:5:360; 
pixSize = 3.45/1000/0.6*3;  % mm / pix for spiral radius
%% load spiral phase maps for all sessions and cancatenate
radius_all = [40:10:100];
vxRaw_session_all = zeros(165,143,15,7);
vyRaw_session_all = zeros(165,143,15,7);
k = 5;
radius = radius_all(k);
spiral_phase_all = [];
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(save_folder,['radius_' num2str(radius) 'mm'],...
        [fname '_mean_flow.mat']),'spiral_phase_all_norm');
    spiral_phase_all = cat(4,spiral_phase_all,spiral_phase_all_norm);
end
% average the phase maps and get mean spiral optical flow field
spiral_phase_mean = circ_mean(spiral_phase_all, [], 4);
spiral_phase_mean1 = permute(spiral_phase_mean,[3 1 2]);   
spiral_phase_all1 = squeeze(spiral_phase_all(:,:,1,:));
spiral_phase_var = circ_var(spiral_phase_all1, [],[], 3);
useGPU = 0;
[vxRaw,vyRaw] = HS_flowfield(spiral_phase_mean1,useGPU);
vxRaw = squeeze(vxRaw); vyRaw = squeeze(vyRaw);
%%
points2(1,:) = [78,80]; %out
points2(2,:) = [73,101]; %in
points2(3,:) = [107,93]; %out
points2(4,:) = [81,107];  %in 
%% load spiral phase maps for all sessions and cancatenate
radius_all = [40:10:100];
k = 5;
radius = radius_all(k);
spiral_phase_all = [];
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(save_folder,['radius_' num2str(radius) 'mm'],...
        [fname '_mean_flow.mat']),'spiral_phase_all_norm');
    %%
    spiral_phase_mean2 = squeeze(circ_mean(spiral_phase_all_norm, [], 4));
    spiral_phase_mean2 = permute(spiral_phase_mean2,[3 1 2]); 
    [vxRawa,vyRawa] = HS_flowfield(spiral_phase_mean2,useGPU);
    vxRawAll(:,:,kk) = squeeze(vxRawa); vyRawAll(:,:,kk) = squeeze(vyRawa);
end
%%            
hemi = [];
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
px1 = 110;
py1 = 72;
skip = 5;
zoom_scale = 2;
h1b = figure('Renderer', 'painters', 'Position', [10 10 1000 400]);
radius_all = [40:10:100];
radius = radius_all(5);
    
ax1 = subplottight(2,3,1);
framea = squeeze(spiral_phase_mean1(1,:,:));
im_phase = imagesc(framea);
colormap(ax1,colorcet('C06'));
axis image; 
axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
r = radius*pixSize/(0.01*8);
cx2 = r*cosd(th2)+px1;
cy2 = r*sind(th2)+py1;
hold on;
color1 = 'w';                                                               % clockwise, then color black
plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);               % draw the circle at max radius
hold on;
scatter(px1,py1,8,color1,'filled');
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
cb1 = colorbar;
pos1 = cb1.Position;
cb1.Position = [pos1(1)+0.05,pos1(2)+0.1,pos1(3),pos1(4)/2];

ax2 = subplottight(2,3,2);
framea = squeeze(spiral_phase_mean1(2,:,:));
im_phase = imagesc(framea);
colormap(ax2,colorcet('C06'));
axis image; 
axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
r = radius*pixSize/(0.01*8);
cx2 = r*cosd(th2)+px1;
cy2 = r*sind(th2)+py1;
hold on;
color1 = 'w';                                                               % clockwise, then color black
plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);               % draw the circle at max radius
hold on;
scatter(px1,py1,8,color1,'filled');
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
cb2 = colorbar;
pos1 = cb2.Position;
cb2.Position = [pos1(1)+0.05,pos1(2)+0.1,pos1(3),pos1(4)/2];
%
ax3 = subplottight(2,3,3);
% framea = squeeze(spiral_phase_mean1(1,:,:));

blank1 = ones(size(framea));
vxRaw2 = nan(size(vxRaw));
vyRaw2 = nan(size(vyRaw));
vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
im_phase = imagesc(blank1);
% colormap(ax3,colorcet('C06'));
colormap(ax3,ones(10,3));
axis image; 
axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
hold on;
imH1Raw4 = quiver(vxRaw2,vyRaw2,'k','lineWidth',1,'autoScale','off');
r = radius*pixSize/(0.01*8);
cx2 = r*cosd(th2)+px1;
cy2 = r*sind(th2)+py1;
hold on;
color1 = 'w';                                                               % clockwise, then color black
plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);               % draw the circle at max radius
hold on;
scatter(px1,py1,8,color1,'filled');
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
for i = 1:4
    hold on;
    scatter(points2(i,2),points2(i,1),24,'MarkerEdgeColor',color3(order1(i),:),'MarkerFaceColor','None');
end
cb3 = colorbar;
pos1 = cb3.Position;
cb3.Position = [pos1(1)+0.05,pos1(2)+0.1,pos1(3),pos1(4)/2];

ax4 = subplottight(2,3,4);
framea = squeeze(spiral_phase_var);
im_phase = imagesc(framea);
color2 = cbrewer2('seq','OrRd',100);
colormap(ax4,color2);
axis image; 
axis off;
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);

r = radius*pixSize/(0.01*8);
cx2 = r*cosd(th2)+px1;
cy2 = r*sind(th2)+py1;
hold on;
color1 = 'w';                                                               % clockwise, then color black
plot([cx2 cx2(1)],[cy2 cy2(1)],'color',color1,'LineWidth',2);               % draw the circle at max radius
hold on;
scatter(px1,py1,8,color1,'filled');
set(im_phase, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
cb4 = colorbar;
pos1 = cb4.Position;
cb4.Position = [pos1(1)+0.05,pos1(2)+0.1,pos1(3),pos1(4)/2];
%
ax5 = subplottight(2,3,5);
color3 = cbrewer2('qual','Paired',10);
order1 = [1,2,5,6];
color1a = {'r','b','g','m'};
mn = 'ZYE_0052';
td = '2021-12-18';
tdb = datestr(td,'yyyymmdd');
en = 2;
subfolder = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'spirals','svd',subfolder);
[U,V,t,mimg] = loadUVt1(session_root);                                     % load U,V, t
fname = [mn '_' tdb '_' num2str(en)];
load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));               % load atlas transformation matrix tform;
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));

mimgt = mimgt(1:8:end,1:8:end);
im1 = imagesc(mimgt);
colormap(ax5,gray);
hold on;
plot([cx2 cx2(1)],[cy2 cy2(1)],'color','w','LineWidth',2);               % draw the circle at max radius
hold on;
scatter(px1,py1,8,'w','filled');
set(im1, 'AlphaData', BW1, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
axis off; axis image;
for i = 1:4
    hold on;
    scatter(points2(i,2),points2(i,1),24,'MarkerEdgeColor',color3(order1(i),:),'MarkerFaceColor','None');
end
cb5 = colorbar;
pos1 = cb5.Position;
cb5.Position = [pos1(1)+0.05,pos1(2)+0.1,pos1(3),pos1(4)/2];
%
ax6 = subplottight(2,3,6);
for i = 1:4
    vector2(i,1) = vxRaw(points2(i,1),points2(i,2));
    vector2(i,2) = vyRaw(points2(i,1),points2(i,2));
end
vector2(:,2) = -vector2(:,2);
c1 = compass(squeeze(vector2(4,1)),squeeze(vector2(4,2)),color1a{4});
c1.LineWidth = 2; c1.Color = color3(order1(4),:);
hold on;
c2 = compass(squeeze(vector2(3,1)),squeeze(vector2(3,2)),color1a{4});
c2.LineWidth = 2; c2.Color = color3(order1(3),:);
hold on;
c3 = compass(squeeze(vector2(2,1)),squeeze(vector2(2,2)),color1a{4});
c3.LineWidth = 2; c3.Color = color3(order1(2),:);
hold on;
c4 = compass(squeeze(vector2(1,1)),squeeze(vector2(1,2)),color1a{4});
c4.LineWidth = 2; c4.Color = color3(order1(1),:);

clear vector3
for i = 1:4
    vector3(i,1,:) = vxRawAll(points2(i,1),points2(i,2),:);
    vector3(i,2,:) = vyRawAll(points2(i,1),points2(i,2),:);
end
vector3(:,2,:) = -vector3(:,2,:);
color1 = {'r','b','g','m'};
for k = 1:15
    c1a = compass(squeeze(vector3(4,1,k)),squeeze(vector3(4,2,k)),color1a{4});
    c1a.LineWidth = 0.5; c1a.Color = color3(order1(4),:); 
    hold on;
    c2a = compass(squeeze(vector3(3,1,k)),squeeze(vector3(3,2,k)),color1a{3});
    c2a.LineWidth = 0.5; c2a.Color = color3(order1(3),:);
    hold on;
    c3a = compass(squeeze(vector3(2,1,k)),squeeze(vector3(2,2,k)),color1a{2});
    c3a.LineWidth = 0.5; c3a.Color = color3(order1(2),:);
    hold on;
    c4a = compass(squeeze(vector3(1,1,k)),squeeze(vector3(1,2,k)),color1a{1});
    c4a.LineWidth = 0.5; c4a.Color = color3(order1(1),:);
end
%%
print(h1b, 'mean_rotating_wave2.pdf','-dpdf', '-bestfit', '-painters');