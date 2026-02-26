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
% area_index:{MO_L_index,MO_R_index,SSp_L_index,SSp_R_index};
area_index = getAreaIndex(data_folder);
BW = logical(projectedAtlas1);
%%
data_folder = 'E:\task2';
list_folder = 'C:\Users\Steinmetz lab\Documents\git\YE-et-al-2023-spirals\revision2\cutting';
T1 = readtable(fullfile(list_folder,'cutting_session_list.xlsx'));
T2 = T1(ismember(T1.area,'SSp') & ismember(T1.hemisphere,'bilateral'),:);
sessions = unique(T2.session_id);
%%
folder1 = 'E:\task2';
spiral_folder = 'E:\task2\spirals\spirals_grouping';
spiral_folder_predicted = 'E:\manipulation_prediction\cutting_prediction\spirals\spirals_grouping';
svd_predict_folder = 'E:\manipulation_prediction\cutting_prediction\prediction_svd';
%%
i = 3;
%% load raw image with cutting
Ti = T2(T2.session_id == sessions(i) & strcmp(T2.label,'cutting'),:);
mn = Ti.MouseID{1};
tda = Ti.date(1);
en = Ti.folder(1);
tdb = datestr(tda,'yyyymmdd');
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(folder1,'task_svd',fname);
[U1,V1,t1,mimg1] = loadUVt2(session_root); 
dV1 = [zeros(size(V1,1),1) diff(V1,[],2)];
dV1 = dV1(1:50,:);
load(fullfile(session_root, 'meanImage.mat')); %mimg
load(fullfile(folder1,'rfmap',mn,[fname '.mat']));
mimgt1 = imwarp(mimg1,tform,'OutputView',imref2d(size(projectedAtlas1)));
%%
U1t = imwarp(U1,tform,'OutputView',imref2d(size(projectedAtlas1)));
U1t = U1t(1:8:end,1:8:end,1:50);
%% get craniotomy ROI
[centers,radius] = getCraniotomyROI(mimgt1);
%% load raw spirals
clear archiveCell spirals_raw
load(fullfile(spiral_folder,[fname '_spirals_group_fftn.mat']));
spirals_raw = transformGroupedSprials(archiveCell,tform);
spirals_grouped_raw = [];
for k = 1:4
    clear brain_index_temp liaL locbL
    brain_index_temp = area_index{k}; 
    [liaL,locbL] = ismember(spirals_raw(:,1:2),brain_index_temp,'rows');
    spirals_grouped_raw{k} = spirals_raw(liaL,:);
end
%% load spirals from prediction
clear archiveCell tform spirals_predict
load(fullfile(spiral_folder_predicted,[fname '_spirals_group_fftn_MO_cutting.mat']));
% load first tform in merged_data (spontaneous)
Ti = T2(T2.session_id == sessions(i) & strcmp(T2.label,'spontaneous'),:);
mn = Ti.MouseID{1}; tda = Ti.date(1); en = Ti.folder(1);
tdb = datestr(tda,'yyyymmdd');
fname1 = [mn '_' tdb '_' num2str(en)];
session_root0 = fullfile(folder1,'task_svd',fname1);
[U0,V0,t0,mimg0] = loadUVt2(session_root0); 
load(fullfile(folder1,'rfmap',mn,[fname1 '.mat']));
%%
figure;
subplot(1,2,1);
imagesc(mimg0);
axis image; axis off;
subplot(1,2,2);
imagesc(mimg1);
axis image; axis off;
%% outline binary image
clear cpstruct
mimg1 = mimg1/max(mimg1(:));
mimg0 = mimg0/max(mimg0(:));
h = cpselect(mimg0(1:8:end,1:8:end),mimg1(1:8:end,1:8:end));
%%
[movingPoints,fixedPoints] = cpstruct2pairs(cpstruct); 
tform_bwn = fitgeotrans(movingPoints,fixedPoints,'affine');
sizeTemplate = size(mimg1);
mimg0_tf = imwarp(mimg0,tform_bwn,'OutputView',imref2d(sizeTemplate));
h1 = figure('Renderer', 'painters', 'Position', [100 100 800 600]);
ax1 = subplot(1,1,1); 
C = imfuse(mimg1,mimg0_tf,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
imshow(C);
%%
spirals_predict = transformGroupedSprials(archiveCell,tform);
spirals_grouped_predict = [];
for k = 1:4
    clear brain_index_temp liaL locbL
    brain_index_temp = area_index{k}; 
    [liaL,locbL] = ismember(spirals_predict(:,1:2),brain_index_temp,'rows');
    spirals_grouped_predict{k} = spirals_predict(liaL,:);
end
%% find spirals within 100 pixels distance from each raw spiral
spiral_match = [];
spiral_nomatch = [];
for k = 1:4
    clear spirals_raw1 spirals_predict1
    spirals_raw1 = spirals_grouped_raw{k};
    spirals_predict1 = spirals_grouped_predict{k};
    spiral_match{k} = getMatchedSprials(spirals_raw1,spirals_predict1);
    spiral_nomatch{k} = getNoMatchedSprials(spirals_raw1,spirals_predict1);
end
%%
N_match = cellfun(@(x) size(x,1), spiral_match);
N_nomatch = cellfun(@(x) size(x,1), spiral_nomatch);
N_raw = cellfun(@(x) size(x,1), spirals_grouped_raw);
N_predict = cellfun(@(x) size(x,1), spirals_grouped_predict);
N_array = [N_raw;N_predict;N_match;N_nomatch];
%%
Ta = array2table(N_array,'VariableNames',{'MOL','MOR','SSpL','SSpR'});
names = {'Raw';'Predict';'Match';'NoMatch'};
Ta1 = table(names);
Ta = [Ta1,Ta];
%% load predicted svd
load(fullfile(svd_predict_folder,[mn '_prediction_MO_cutting.mat']));
U_predict = Unew1(1:8:end,1:8:end,:);
U_predict2 = imwarp(U_predict,tform_bwn,'OutputView',imref2d(size(mimg1(1:8:end,1:8:end))));
dV_predict = [zeros(size(V_predict,1),1) diff(V_predict,[],2)];
params.downscale = 8; params.lowpass = 0; params.gsmooth = 0;
rate = 1; freq = [2,8];
[trace2d1,~,tracePhase_predict] = spiralPhaseMap_freq(U_predict2,dV_predict,t1,params,freq,rate);
%% load raw svd
U1d = U1(1:8:end,1:8:end,1:50);
[trace2d_raw,~,tracePhase_raw] = spiralPhaseMap_freq(U1d,dV1,t1,params,freq,rate);
%%
% trace2d2 = imresize(trace2d1(:,:,1001:1008),[560,560]);
% trace2d3 = imwarp(trace2d2,tform_bwn,'OutputView',imref2d(size(mimg1)));
figure;
for i = 1:8
    subplottight(2,8,i);
    C = imfuse(squeeze(trace2d1(:,:,1000+i)),squeeze(trace2d_raw(:,:,1000+i)),...
        'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    imshow(C);
end
%%
Upt = imwarp(Unew1,tform,'OutputView',imref2d(size(projectedAtlas1)));
Upt = Upt(1:8:end,1:8:end,:);
%%
% save([mn '_phase.mat'],'tracePhase_raw','tracePhase_predict');
%%
figure;
hist_bin = 40;
pixSize = 0.01; % mm/pix
pixArea = pixSize^2;
frame_all = numel(t1);
ax1 = subplot(1,1,1);
[unique_spirals,scolor,low_color_bound,high_color_bound] = density_color_plot(spirals_predict,hist_bin);
unique_spirals_unit = unique_spirals(:,3)/(hist_bin*hist_bin*pixArea); % spiral counts/mm^2
unique_spirals_unit = unique_spirals_unit./frame_all*35; % spirals/(mm^2*s)
[ax1,cb1]= plotDesnity(ax1,unique_spirals,unique_spirals_unit,coords,maskPath,BW,st,atlas1);
caxis([0,2]);
axis on;
%%
clear index
sprials_temp = spirals_predict;
index1 = (sprials_temp(:,1)>350 & sprials_temp(:,1)<450 &...
    sprials_temp(:,2)>400 & sprials_temp(:,2)<500 );

index2 = (sprials_temp(:,1)>400 & sprials_temp(:,1)<550 &...
    sprials_temp(:,2)>250 & sprials_temp(:,2)<350);

index3 = (sprials_temp(:,1)>750 & sprials_temp(:,1)<850 &...
    sprials_temp(:,2)>400 & sprials_temp(:,2)<500 );

index4 = (sprials_temp(:,1)>650 & sprials_temp(:,1)<800 &...
    sprials_temp(:,2)>250 & sprials_temp(:,2)<350);
 
index = any([index1,index2,index3,index4],2);

% MO_index = cat(1,area_index{1},area_index{2});
% [index,b] = ismember(spirals_predict(:,1:2),MO_index,'rows');
% spirals_predict1 = spirals_predict(index,:);
%% get craniotomy ROI
[centers1,radius1] = getCraniotomyROI(mimg1);
%%
centers1 = centers1/8;
radius1 = radius1/8;
%%
th2 = 1:5:360; 
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
ax1 = subplot(1,2,1);
frame_select = unique(spirals_predict(:,5));
tracePhase_predict1 = tracePhase_predict(:,:,frame_select);
tracePhase_raw1 = tracePhase_raw(:,:,frame_select);
phase_diff = angdiff(tracePhase_predict1,tracePhase_raw1);
r = circ_r(phase_diff, [], [], 3);
imagesc(1-r);
hold on;
for kk = 1:2
    px1 = centers1(kk,1);
    py1 = centers1(kk,2);
    r = radius1(kk,1);
    cx2 = r*cosd(th2)+px1;
    cy2 = r*sind(th2)+py1;
    hold on;
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
end
axis image; axis off;
% colormap(ax1,'hot');
colormap(ax1,'parula');
caxis([0,0.3]);
colorbar;
title('Angular variance');

ax2 = subplot(1,2,2);
frame_select = unique(spirals_predict(index,5));
tracePhase_predict1 = tracePhase_predict(:,:,frame_select);
tracePhase_raw1 = tracePhase_raw(:,:,frame_select);
phase_diff = angdiff(tracePhase_predict1,tracePhase_raw1);
r = circ_r(phase_diff, [], [], 3);
imagesc(1-r);
hold on;
for kk = 1:2
    px1 = centers1(kk,1);
    py1 = centers1(kk,2);
    r = radius1(kk,1);
    cx2 = r*cosd(th2)+px1;
    cy2 = r*sind(th2)+py1;
    hold on;
    plot([cx2 cx2(1)],[cy2 cy2(1)],'k','LineWidth',1);          % draw the circle at max radius
end
axis image; axis off;
% colormap(ax2,'hot');
colormap(ax2,'parula');
caxis([0,0.3]);
colorbar;
title('Angular variance');
%%
print(h1, [mn '_angular variance'], '-dpdf', '-bestfit', '-painters');
%% registered angular variance
frame_select = unique(spirals_predict(index,5));
tracePhase_predict1 = tracePhase_predict(:,:,frame_select);
tracePhase_raw1 = tracePhase_raw(:,:,frame_select);
phase_diff = angdiff(tracePhase_predict1,tracePhase_raw1);
r1 = circ_r(phase_diff, [], [], 3);
r2 = imresize(r1,[560,560]);
