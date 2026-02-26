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
%%
U_predict = Unew1(1:8:end,1:8:end,:);
dV_predict = [zeros(size(V_predict,1),1) diff(V_predict,[],2)];
params.downscale = 8; params.lowpass = 0; params.gsmooth = 0;
rate = 1; freq = [2,8];
[trace2d1,~,tracePhase_predict] = spiralPhaseMap_freq(U_predict,dV_predict,t1,params,freq,rate);
% load raw svd
merge_folder = 'E:\manipulation_prediction\cutting_prediction\merged_data';
load(fullfile(merge_folder,[mn '_merge_cutting.mat']));
Unew2 = reshape(Unew1,560,560,size(Unew1,2));
U1d = Unew2(1:8:end,1:8:end,1:50); 
V1 = Vnew1(:,size(V0,2)+1:end);
dV1 = [zeros(size(V1,1),1) diff(V1,[],2)];
dV1 = dV1(1:50,:);

% U1d = U_predict;
% U1d = imtranslate(U1d,[0,5]);
% dV1 = dV_predict;
% U1d = U1(1:8:end,1:8:end,1:50);
[trace2d_raw,~,tracePhase_raw] = spiralPhaseMap_freq(U1d,dV1,t1,params,freq,rate);
%%
% U_predict = U1d;
% U_predict = imtranslate(U_predict,[0,5]);
% dV_predict = dV1;
% [trace2d1,~,tracePhase_predict] = spiralPhaseMap_freq(U_predict,dV_predict,t1,params,freq,rate);
%%
% trace2d2 = imresize(trace2d1(:,:,1001:1008),[560,560]);
% trace2d3 = imwarp(trace2d2,tform_bwn,'OutputView',imref2d(size(mimg1)));
figure;
for k = 1:8
    subplottight(2,8,k);
    C = imfuse(squeeze(tracePhase_predict(:,:,2000+k)),squeeze(tracePhase_raw(:,:,2000+k)),...
        'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    imshow(C);
    subplottight(2,8,k+8);
    phase_diff = angdiff(squeeze(tracePhase_predict(:,:,2000+k)),squeeze(tracePhase_raw(:,:,2000+k)));
    imagesc(phase_diff);
    caxis([-pi,pi]);
    axis image; axis off;
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
%%
% mimg1d = imtranslate(mimg1,[0,5*8]);  
mimg1d = mimg1;
[centers1,radius1] = getCraniotomyROI(mimg1d);
centers1 = centers1/8;
radius1 = radius1/8;
%%
figure;
subplot(1,2,1); imagesc(mimg0); axis image;
subplot(1,2,2); imagesc(mimg1); axis image;
%%
th2 = 1:5:360; 
h1 = figure('Renderer', 'painters', 'Position', [100 100 600 300]);
ax1 = subplot(1,3,1);
frame_select = unique(spirals_predict(:,5));
tracePhase_predict1 = tracePhase_predict(:,:,frame_select);
tracePhase_raw1 = tracePhase_raw(:,:,frame_select);
phase_diff = angdiff(tracePhase_predict1,tracePhase_raw1);
ra = circ_r(phase_diff, [], [], 3);
imagesc(1-ra);
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
colormap(ax1,'parula');
% caxis([0,0.6]);
colorbar;
title('Angular variance');

ax1 = subplot(1,3,2);
% frame_select = unique(spirals_predict(index,5));
frame_select = unique(spirals_predict(:,5));
% tracePhase_predict1 = tracePhase_predict(:,:,frame_select);
tracePhase_raw1 = tracePhase_raw(:,:,frame_select);
tracePhase_raw1_shift = imtranslate(tracePhase_raw1,[0,3]);
phase_diff = angdiff(tracePhase_raw1_shift,tracePhase_raw1);
rb = circ_r(phase_diff, [], [], 3);
imagesc(1-rb);
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
colormap(ax1,'parula');
caxis([0,0.4]);
colorbar;
title('Angular variance');
%%
print(h1, [mn '_angular variance5'], '-dpdf', '-bestfit', '-painters');
%%
shift1 = 3;

frame_select = unique(spirals_predict(:,5));
tracePhase_raw1 = tracePhase_raw(:,:,frame_select);
tracePhase_predict1 = tracePhase_predict(:,:,frame_select);
% phase_diff = angdiff(tracePhase_predict1,tracePhase_raw1);
phase_diff = tracePhase_predict1;

[a,b,c] = size(phase_diff);
phase_diff2a = zeros(a,b,c,5);
phase_diff2a(:,:,:,1) = imtranslate(phase_diff,[0,shift1]); 
phase_diff2a(:,:,:,2) = imtranslate(phase_diff,[0,-shift1]); 
phase_diff2a(:,:,:,3) = phase_diff; 
phase_diff2a(:,:,:,4) = imtranslate(phase_diff,[shift1,0]); 
phase_diff2a(:,:,:,5) = imtranslate(phase_diff,[-shift1,0]); 

cir_var1 = circ_r(phase_diff2a,[],[],4);
cir_var2 = squeeze(mean(cir_var1,3));
%%
figure;
imagesc(1-cir_var2);
% caxis([0,0.3]);
colorbar;
%%
shift1 = 3;
frame_select = unique(spirals_predict(:,5));
tracePhase_raw1 = tracePhase_raw(:,:,frame_select);
tracePhase_predict1 = tracePhase_predict(:,:,frame_select);
% phase_diff = angdiff(tracePhase_predict1,tracePhase_raw1);
phase_mean = tracePhase_raw1;

[a,b,c] = size(phase_mean);
phase_mean2a = zeros(a,b,c,5);
phase_mean2a(:,:,:,1) = imtranslate(phase_mean,[0,shift1]); 
phase_mean2a(:,:,:,2) = imtranslate(phase_mean,[0,-shift1]); 
phase_mean2a(:,:,:,3) = phase_mean; 
phase_mean2a(:,:,:,4) = imtranslate(phase_mean,[shift1,0]); 
phase_mean2a(:,:,:,5) = imtranslate(phase_mean,[-shift1,0]); 

cir_mean1 = circ_mean(phase_mean2a,[],4);
cir_mean2 = circ_r(cir_mean1,[],[],3);
%%
figure;
imagesc(cir_mean2);
caxis([0,0.05]);
colorbar;
