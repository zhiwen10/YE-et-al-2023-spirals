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
for k = 1:7
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
        %%
%         spiral_phase_all_norm1 = squeeze(spiral_phase_all_norm(:,:,1,:));
%         spiral_phase_session_mean = squeeze(circ_mean(spiral_phase_all_norm1,[],3));
%         spiral_phase_session_var = squeeze(circ_var(spiral_phase_all_norm1, [], [], 3));
%         spiral_phase_session_mean1 = permute(spiral_phase_session_mean,[3 1 2]);   
%         useGPU = 0;
%         [vxRaw_session,vyRaw_session] = HS_flowfield(spiral_phase_session_mean1,useGPU);
%         vxRaw_session_all(:,:,kk,k) = squeeze(vxRaw_session); 
%         vyRaw_session_all(:,:,kk,k) = squeeze(vyRaw_session);
%         phase_session_mean_all(:,:,kk,k) = squeeze(spiral_phase_session_mean(:,:,1));
%         phase_session_var_all(:,:,kk,k) = squeeze(spiral_phase_session_mean(:,:,1));
    end
    % average the phase maps and get mean spiral optical flow field
    spiral_phase_mean = circ_mean(spiral_phase_all, [], 4);
    spiral_phase_mean1 = permute(spiral_phase_mean,[3 1 2]);   
    spiral_phase_all1 = squeeze(spiral_phase_all(:,:,1,:));
    spiral_phase_var = circ_var(spiral_phase_all1, [],[], 3);
    useGPU = 0;
    [vxRaw,vyRaw] = HS_flowfield(spiral_phase_mean1,useGPU);
    vxRaw = squeeze(vxRaw); vyRaw = squeeze(vyRaw);
    vxRaw_all(:,:,k) = vxRaw;
    vyRaw_all(:,:,k) = vyRaw;
    spiral_phase_mean_all(:,:,k) = squeeze(spiral_phase_mean1(1,:,:));
    spiral_phase_var_all(:,:,k) = squeeze(spiral_phase_var);  
end
%%
hemi = [];
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
px1 = 110;
py1 = 72;
skip = 5;
zoom_scale = 2;
a = 60; b = 135; % height
c = 65; d = 140; % width
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
h1b = figure('Renderer', 'painters', 'Position', [10 10 1000 400]);
radius_all = [40:10:100];
for k = 1:7
    %%
    radius = radius_all(k);
    ax1 = subplottight(2,7,k);
    framea = squeeze(spiral_phase_mean_all(:,:,k));
    
    vxRaw = vxRaw_all(:,:,k);
    vyRaw = vyRaw_all(:,:,k);
    vxRaw2 = nan(size(vxRaw));
    vyRaw2 = nan(size(vyRaw));
    vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
    vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
    vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
    im_phase = imagesc(framea);
    colormap(ax1,colorcet('C06'));
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
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','w','lineWidth',2);
    %%
    ax2(k) = subplottight(2,7,7+k);
    ax2(k).Position(1) = ax2(k).Position(1)+0.005;
    ax2(k).Position(2) = ax2(k).Position(2);
    ax2(k).Position(3) = ax2(k).Position(3)-0.015;
    ax2(k).Position(4) = ax2(k).Position(4)-0.03;
    
    im_phase = imagesc(framea);
    colormap(ax2(k),colorcet('C06'));
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
    ylim([a,b]);
    xlim([c,d]);
end
%%
print(h1b, 'mean_rotating_wave.pdf','-dpdf', '-bestfit', '-painters');
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
for k = 1:7
    %%
    radius = radius_all(k);
    %%
    ax1 = subplottight(3,7,k);
    framea = squeeze(spiral_phase_mean_all(:,:,k));
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
    %%
    ax1 = subplottight(3,7,7+k);
    framea = squeeze(spiral_phase_var_all(:,:,k));
    im_phase = imagesc(framea);
    % colormap(ax1,colorcet('C06'));
    color2 = cbrewer2('seq','OrRd',100);
    % color2 = flipud(color2);
    colormap(ax1,color2);
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
    %%
    ax1 = subplottight(3,7,14+k);
    framea = squeeze(spiral_phase_mean_all(:,:,k));
    
    vxRaw = vxRaw_all(:,:,k);
    vyRaw = vyRaw_all(:,:,k);
    vxRaw2 = nan(size(vxRaw));
    vyRaw2 = nan(size(vyRaw));
    vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
    vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
    vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
    im_phase = imagesc(framea);
    colormap(ax1,colorcet('C06'));
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
end
%%
points2(1,:) = [72,105];
% points2(1,:) = [85,115]; % SSp
points2(2,:) = [71,85];
points2(3,:) = [81,103];
points2(4,:) = [108,88];
%%
h1= figure('Renderer', 'painters', 'Position', [100 100 800 200]);
color1 = {'r','b','g','m'};
x = zeros(15,1); y = zeros(15,1);
for k = 1:7
    subplot(1,7,k);
    clear vector2
    for i = 1:4
        vector2(1,:) = squeeze(vxRaw_session_all(points2(i,1),points2(i,2),:,k));
        vector2(2,:) = squeeze(vyRaw_session_all(points2(i,1),points2(i,2),:,k));
        vector2(2,:) = -vector2(2,:);
        vector3 = complex(squeeze(vector2(1,:)),squeeze(vector2(2,:)));
        compass(vector3,color1{i});
        hold on;
    end 
end
%%
color1 = cbrewer2('seq','greys',10);
figure;
for k  = 1:4
    subplot(1,4,k)
    vx1 = squeeze(vxRaw_all(points2(k,1),points2(k,2),:));
    vy1 = squeeze(vxRaw_all(points2(k,1),points2(k,2),:));
    for i = 1:7
        quiver(0,0,vx1(i),vy1(i),'color',color1(2+i,:));
        hold on;
    end
end
%% load spiral phase maps for all sessions and cancatenate
hemi = [];
scale3 = 5/8;
lineColor = 'k'; lineColor1 = 'w';
px1 = 103;
py1 = 90;
skip = 5;
zoom_scale = 2;
radius = 100;
spiral_phase_all = [];
for kk = 1:size(T,1)
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    %% load raw detected spirals
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(save_folder,['radius_' num2str(radius) 'mm-2'],...
        [fname '_mean_flow.mat']),'spiral_phase_all_norm');
    spiral_phase_all = cat(4,spiral_phase_all,spiral_phase_all_norm);
end
% average the phase maps and get mean spiral optical flow field
spiral_phase_mean = circ_mean(spiral_phase_all, [], 4);
spiral_phase_mean1 = permute(spiral_phase_mean,[3 1 2]);   
useGPU = 0;
[vxRaw,vyRaw] = HS_flowfield(spiral_phase_mean1,useGPU);
vxRaw = squeeze(vxRaw); vyRaw = squeeze(vyRaw);
%%
h1b = figure('Renderer', 'painters', 'Position', [10 10 400 400]);
ax1 = subplottight(1,1,1);
framea = squeeze(spiral_phase_mean1(1,:,:));
vxRaw2 = nan(size(vxRaw));
vyRaw2 = nan(size(vyRaw));
vxRaw2(1:skip:end,1:skip:end) = vxRaw(1:skip:end,1:skip:end)*zoom_scale;
vyRaw2(1:skip:end,1:skip:end) = vyRaw(1:skip:end,1:skip:end)*zoom_scale;  
vxRaw2(~BW1) = nan; vyRaw2(~BW1) = nan;
im_phase = imagesc(framea);
colormap(ax1,colorcet('C06'));
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