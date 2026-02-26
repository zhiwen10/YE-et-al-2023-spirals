githubdir2 = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir2, 'YE-et-al-2023-spirals')));            % paper repository
%% load atlas brain horizontal projection and outline
data_folder = 'E:\spiral_data_share\data'; 
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
%%
mn = 'ZYE_0092';
tda = '2025-06-21';
en = 1;
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
data_folder = 'E:\task2\task_svd';
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,fname);
[U1,V1,t1,mimg1] = loadUVt2(session_root); 
folder_rf = 'E:\task2';
load(fullfile(folder_rf,'rfmap',mn,[fname '.mat']));
U1 = U1(:,:,1:50)./mimg1;
U_raw = imwarp(U1,tform,'OutputView',imref2d(size(projectedAtlas1)));
U_raw = U_raw(1:8:end,1:8:end,:);
dV_raw = [zeros(size(V1,1),1) diff(V1,[],2)];
dV_raw = dV_raw(1:50,:);
mimgt1 = imwarp(mimg1,tform,'OutputView',imref2d(size(projectedAtlas1)));
%%
load('ZYE_0090_20260621_1_match_nomatch_spirals_all.mat');
load('ZYE_0090_20260621_1_craniotomy.mat');
%% predicted trace (Unew1, V_preidct)
reg_folder = 'E:\task2\rfmap\ZYE_0092';
load('E:\muscimol_prediction\prediction_svd\ZYE_0092_prediction_MO_bilateral2.mat');
%%
centers1 = [];
edges1 = [];
figure;
ax = subplot(1,1,1);
imagesc(mimg1);
colormap('gray');
axis image;
for k = 1:2
    roi = drawpoint(ax);
    centers1(k,:) = roi.Position;
    roi = drawpoint(ax);
    edges1(k,:) = roi.Position;
end
for k = 1:2
    radius1(k,1) = hypot(edges1(k,1)-centers1(k,1),edges1(k,2)-centers1(k,2));
end    
centers1 = round(centers1);
radius1 = round(radius1);
%% load spontaneous mimg for registration of predicted trace
fname2 = 'ZYE_0092_20250619_2';
session_root2 = fullfile(data_folder,fname2);
load(fullfile(session_root2, 'meanImage.mat')); %mimg
mimg2 = mimg;
load(fullfile(reg_folder,[fname2 '.mat']));
Unew2 = Unew1./mimg2;
U_predict = imwarp(Unew1,tform,'OutputView',imref2d(size(projectedAtlas1)));
U_predict = U_predict(1:8:end,1:8:end,:);
dV_predict = [zeros(size(V_predict,1),1) diff(V_predict,[],2)];
%%
Unew2_ds = Unew1(1:8:end,1:8:end,:);
params.downscale = 8;
params.lowpass = 0;
params.gsmooth = 0;
rate = 1;
freq = [2,8];
[trace2d1,~,tracePhase1] = spiralPhaseMap_freq(Unew2_ds,dV_predict,t1,params,freq,rate);
%%
for k = 1:2
    radius(k,1) = hypot(edges(k,1)-centers(k,1),edges(k,2)-centers(k,2));
end    
centers = round(centers);
radius = round(radius);
%%
hemi = [];
lineColor = 'k';
th2 = 1:5:360; 
scale1 = 5;
BW = logical(projectedAtlas1);
th2 = 1:5:360; 
pixSize = 3.45/1000/0.6*3;  % mm / pix
pixArea = pixSize^2; %mm^2
hist_bin = 40;
figure;
ax1 = subplot(1,2,1);
im  = imagesc(zeros(size(BW)));
hold on;
[unique_spirals,scolor,low_color_bound,high_color_bound] = density_color_plot(spiral_match_all,hist_bin);
frame_all = numel(t1);
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_all*35; % spirals/(mm^2*s)
obj2 = scatter(unique_spirals(:,1),unique_spirals(:,2),3,unique_spirals_unit,'filled');
set(gca,'Ydir','reverse')
hold on;
for kk = 1:2
    px1 = centers(kk,1);
    py1 = centers(kk,2);
    r = radius(kk,1);
    cx2 = round(r*cosd(th2)+px1);
    cy2 = round(r*sind(th2)+py1);
    hold on;
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
end
set(im, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale1,lineColor);
axis image;

ax1 = subplot(1,2,2);
im  = imagesc(zeros(size(BW)));
hold on;
[unique_spirals,scolor,low_color_bound,high_color_bound] = density_color_plot(spiral_nomatch_all2,hist_bin);
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_all*35; % spirals/(mm^2*s)
obj2 = scatter(unique_spirals(:,1),unique_spirals(:,2),3,unique_spirals_unit,'filled');
set(gca,'Ydir','reverse')
hold on;
for kk = 1:2
    px1 = centers(kk,1);
    py1 = centers(kk,2);
    r = radius(kk,1);
    cx2 = round(r*cosd(th2)+px1);
    cy2 = round(r*sind(th2)+py1);
    hold on;
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
end
set(im, 'AlphaData', BW, 'AlphaDataMapping', 'scaled');
hold on;
plotOutline(maskPath([1:11]),st,atlas1,hemi,scale1,lineColor);
axis image;
%%
index = (spiral_nomatch_all2(:,1)>250 & spiral_nomatch_all2(:,1)<400 &...
    spiral_nomatch_all2(:,2)>600 & spiral_nomatch_all2(:,2)<750 & spiral_nomatch_all2(:,4)==1);
spiral_select = spiral_nomatch_all2(index,:);
frame_select = spiral_select(:,5);
tracePhase_select = tracePhase1(:,:,frame_select);
%% predicted 
point1 = [32,14];
% point1 = [32,25];
phase1 = squeeze(tracePhase_select(point1(1),point1(2),:));
% indx = find(phase1>-0.5 & phase1<0.5);
indx = find(phase1>3*pi/4 | phase1<-3*pi/4);
frame_select2 = frame_select(indx);
tracePhase_select2 = tracePhase1(:,:,frame_select2);
tracePhase_mean = circ_mean(tracePhase_select2, [], 3);
%% raw
U1_ds = U1(1:8:end,1:8:end,:);
[trace2d_raw1,~,tracePhase_raw1] = spiralPhaseMap_freq(U1_ds,dV_raw,t1,params,freq,rate);
%%
tracePhase_select_raw2 = tracePhase_raw1(:,:,frame_select2);
tracePhase_raw_mean = circ_mean(tracePhase_select_raw2, [], 3);
%%
figure;
ax1 = subplot(1,2,1);
im = imagesc(tracePhase_mean);
colormap(ax1,colorcet('C06'));
axis image;
hold on;
for kk = 1:2
    px1 = centers(kk,1)/8;
    py1 = centers(kk,2)/8;
    r = radius(kk,1)/8;
    cx2 = r*cosd(th2)+px1;
    cy2 = r*sind(th2)+py1;
    hold on;
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
end
hold on;
scatter(point1(2),point1(1));
title('Predicted mean');
ax2 = subplot(1,2,2);
im = imagesc(tracePhase_raw_mean);
colormap(ax2,colorcet('C06'));
axis image;
hold on;
for kk = 1:2
    px1 = centers(kk,1)/8;
    py1 = centers(kk,2)/8;
    r = radius(kk,1)/8;
    cx2 = r*cosd(th2)+px1;
    cy2 = r*sind(th2)+py1;
    hold on;
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
end
hold on;
scatter(point1(2),point1(1));
title('Raw mean');
%%
figure;
for i = 1:8
    ax1 = subplot(2,8,i);
    im = imagesc(squeeze(tracePhase_select2(:,:,i)));
    colormap(ax1,colorcet('C06'));
    axis image; axis off;
    hold on;
    for kk = 1:2
        px1 = centers(kk,1)/8;
        py1 = centers(kk,2)/8;
        r = radius(kk,1)/8;
        cx2 = r*cosd(th2)+px1;
        cy2 = r*sind(th2)+py1;
        hold on;
        plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
    end
    ax2 = subplot(2,8,i+8);
    im = imagesc(squeeze(tracePhase_select_raw2(:,:,i)));
    colormap(ax2,colorcet('C06'));
    axis image; axis off;
    hold on;
    for kk = 1:2
        px1 = centers(kk,1)/8;
        py1 = centers(kk,2)/8;
        r = radius(kk,1)/8;
        cx2 = r*cosd(th2)+px1;
        cy2 = r*sind(th2)+py1;
        hold on;
        plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
    end
end