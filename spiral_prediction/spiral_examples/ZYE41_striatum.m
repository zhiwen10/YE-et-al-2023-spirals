githubDir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubDir, 'npy-matlab')))
addpath(genpath(fullfile(githubDir, 'Pipelines_2020-12-09')))
addpath(genpath(fullfile(githubDir, 'spikes')))
addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Kilosort2-master')) % path to kilosort folder
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase'))
addpath(genpath('C:\Users\Steinmetz lab\Documents\git\AP_scripts_cortexlab-master'))
addpath(genpath(fullfile(githubDir, 'wheelAnalysis')))
addpath(genpath(fullfile(githubDir, 'PatternDetection')))
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\colorcet'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\plot_arrow'));
addpath(genpath('C:\Users\Steinmetz lab\Documents\MATLAB\Colormaps'));
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%%
dataFolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\AllenProjection\allenData';
[atlas, metaAVGT] = nrrdread(fullfile(dataFolder, 'annotation_50.nrrd'));
atlas = double(atlas);
atlas1 =zeros(264,228);
for k1 = 1:264
    for k2 = 1:228
        if any(atlas(:,k1,k2))
            firstIndex = find(atlas(:,k1,k2)>0,1);
            atlas1(k1,k2) = atlas(firstIndex,k1,k2);
        end
    end
end
dfolder = 'C:\Users\Steinmetz lab\Documents\git\allenCCF\AllenCCFMap';
st = loadStructureTree(fullfile(dfolder,'structure_tree_safe_2017.csv')); % a table of what all the labels mean
root1 = '/997/';
areaName = {'MOp','MOs'};
maskPath{1} = '/997/8/567/688/695/315/500/985/'; % MOp
maskPath{2} = '/997/8/567/688/695/315/500/993/'; % MOs
maskPath{3} = '/997/8/567/688/695/315/31/'; % ACA
maskPath{4} = '/997/8/567/688/695/315/453/378/'; % SS2
maskPath{5} = '/997/8/567/688/695/315/453/322/'; % SSp
maskPath{6} = '/997/8/567/688/695/315/247/'; % AUD
maskPath{7} = '/997/8/567/688/695/315/669/'; % VIS
maskPath{8} = '/997/8/567/688/695/315/254/'; % RSP
maskPath{9} = '/997/8/567/688/695/315/22'; % VISa
maskPath{10} = '/997/8/567/688/695/315/541/'; % TEa
maskPath{11} = '/997/8/567/688/695/315/677/'; % VISC
scale3 = 5;
lineColor = 'k'; lineColor1 = 'k';
hemi = [];
%% get original and predicted dV by kernel regression.
codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\compare_flow1';
T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
kk = 15; % zye41 STR
% kk = 11; % zye60
% kk = 39; % zye67
% kk = 8; % zye58
ops = get_session_info(T,kk);
%%
[U,V,t,mimg] = get_wf_svd(ops.serverRoot);
dV = [zeros(size(V,1),1) diff(V,[],2)];
%%
goodcluster = ops.curation;
[sp] = get_spikes(ops.ksSubFolder,goodcluster);
%% 
[syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT(ops,t);
WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));
dV1 = double(dV(:,~isnan(WF2ephysT1)));
V1 = double(V(1:50,~isnan(WF2ephysT1)));
[MUA_std] = get_MUA_bin(sp,WF2ephysT);
dV1 = double(dV1(1:50,:));
%%
folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\spirals_raw_fftn';
dfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\registration';
fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
load(fullfile(folder,[fname '_spirals_group_fftn.mat']));
tformName = [fname '_tform.mat'];
load(fullfile(dfolder,tformName));
spiral_duration = cellfun(@(x) size(x,1), archiveCell);
indx2 = (spiral_duration>9);
groupedCells = archiveCell(indx2);
spirals_filt1 = cell2mat(groupedCells);
%%
sizeTemplate = [1320,1140];
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
Utransformed = imwarp(U,tform,'OutputView',imref2d(size(projectedAtlas1)));
mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedAtlas1)));
Utransformed1 = Utransformed(1:8:end,1:8:end,:);
mimgtransformed2 = mimgtransformed(1:8:end,1:8:end);
BW = logical(projectedAtlas1);
BW1 = BW(1:8:end,1:8:end);
%% no registration first for flow vector
scale = 4;
Ut = U(1:scale:end,1:scale:end,1:50);
predict_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch\prediction';
load(fullfile(predict_folder,[fname '_dv_predict.mat']));
%% ZYE58
for i = 1:numel(groupedCells)
    % epochs = 7442:7444;
    spiral_temp = groupedCells{i};
    frames = spiral_temp(:,5);
    epochs = frames(1):frames(end);
    %%
    params.downscale = 8;
    params.lowpass = 0;
    params.gsmooth = 0;
    rate = 1;
    frameStart = epochs(1); frameEnd = epochs(end);
    frameTemp = frameStart-35:frameEnd+35; % extra 2*35 frames before filter data 
    dV_raw_epoch = dV(1:50,frameTemp);
    [trace2d1,traceAmp1,tracePhase1] = spiralPhaseMap4(U(:,:,1:50),dV_raw_epoch,t,params,rate);
    trace2d1 = trace2d1(:,:,1+35/rate:end-35/rate)./mimg; 
    tracePhase1 = tracePhase1(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
    dV_predict_epoch = dV_predict(:,frameTemp);
    [trace2d1p,traceAmp1p,tracePhase1p] = spiralPhaseMap4(U(:,:,1:50),dV_predict_epoch,t,params,rate);
    trace2d1p = trace2d1p(:,:,1+35/rate:end-35/rate)./mimg; 
    tracePhase1p = tracePhase1p(:,:,1+35/rate:end-35/rate); % reduce 2*35 frames after filter data 
    
    frame_n = size(tracePhase1,3);
    h1 = figure('Renderer', 'painters', 'Position', [100 100 800 400]);
    for i = 1:frame_n   
        frame_raw = squeeze(tracePhase1(:,:,i));   
        ax1 = subplottight(2,frame_n,i);
        im1 = imagesc(frame_raw);
        colormap(ax1,colorcet('C06'));
        % set(im1, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
        axis image; axis off;
        indx = (spirals_filt1(:,5) == epochs(i));
        spiral_temp = spirals_filt1(indx,:);
        if not(isempty(spiral_temp))
            hold on;
            scatter(spiral_temp(:,1),spiral_temp(:,2),6,'w','filled');
        end
        %%
        frame_predict = squeeze(tracePhase1p(:,:,i));
        ax2 = subplottight(2,frame_n,frame_n+i);
        im2 = imagesc(frame_predict);
        colormap(ax2,colorcet('C06'));
        % set(im2, 'AlphaData', BW2, 'AlphaDataMapping', 'scaled');
        axis image; axis off;
    end
    print(h1,[fname '_' num2str(epochs(1))], '-dpdf', '-bestfit', '-painters');
    close all;
end
%%
dV1 = dV1(:,1:size(dV_predict,2));
% explained_var_all = get_variance_explained2(Ut,dV1(:,1:50000),dV_predict(:,1:50000));
%% angle difference for oscillation period  
epochs = [62950:63050];
frame = 25:35;

% epochs = [110950:111050];
% frame = 40:50;
% epochs = [100100:100200];
% epochs = [110950:111050];
% wf_t = [1131.5, 1135.5];
% wf_t = [1800,1805];
% epoch_real(1) = find(t>=wf_t(1),1,'first');
% epoch_real(2) = find(t>=wf_t(2),1,'first');
% epochs = epoch_real(1):epoch_real(2);
[data_raw] = get_flowfield_structure2(Utransformed1(:,:,1:50),dV1,epochs,mimgtransformed2);
[data_predict] = get_flowfield_structure2(Utransformed1(:,:,1:50),dV_predict,epochs,mimgtransformed2);
[Axy_all] = compare_flowfield_angle(data_raw,data_predict);
%%
fname1 = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_tform'];
regist_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\registration';
load(fullfile(regist_folder, fname1));
point(1,:) = [101,91]*4;
[point_t(:,1),point_t(:,2)] = transformPointsForward(tform,point(:,1),point(:,2));
%%
% color1 = colorcet('C06','N',8);
colorn = size(MUA_std,1);
color1 = gray(colorn);
randi = randperm(colorn);    
colorAll = color1(randi,:);
%%
figure; 
imagesc(mimgtransformed);
colormap(gray);
pp = drawpolygon('LineWidth',2,'Color','cyan');
BW2down = createMask(pp);
BW2 = (BW&BW2down);
save([fname '_mask.mat'],'BW','BW1','BW2');
%%
fname1 = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_tform_4x'];
regist_folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\registration_4x';
load(fullfile(regist_folder, fname1));
explained_var_all(explained_var_all(:)<0) = 0;
mean_var = squeeze(mean(explained_var_all,3));
mean_var_t = imwarp(mean_var,tform,'OutputView',imref2d(size(projectedTemplate1)));
%

% frame = 21:31;
% frame = 40:50;
% frame = 10:19;
% frame = 42:52;
% frame = 15:25;
filename = [ops.mn, '_', num2str(epochs(1)), '_', num2str(epochs(end))];
count1 = plot_single_pixel3(ops,sp,mean_var_t, mimgtransformed, mimg,BW2, point, point_t,...
    data_predict, data_raw,WF2ephysT,epochs,Axy_all,frame,filename,colorAll,maskPath,st,atlas1,hemi,scale3,lineColor);
%%
[data_raw] = get_flowfield_structure2(Utransformed1(:,:,1:50),dV1,epochs,mimgtransformed2);
[data_predict] = get_flowfield_structure2(Utransformed1(:,:,1:50),dV_predict,epochs,mimgtransformed2);
phase_raw = data_raw.tracePhase1;
phase_predict = data_predict.tracePhase1;
phase_raw = permute(phase_raw,[2,3,1]);
phase_predict = permute(phase_predict,[2,3,1]);
%
vx_raw = permute(data_raw.vxRaw,[2,3,1]);
vy_raw = permute(data_raw.vyRaw,[2,3,1]);
vx_predict = permute(data_predict.vxRaw,[2,3,1]);
vy_predict = permute(data_predict.vyRaw,[2,3,1]);
%%
% a = 50; b = 100; % height
% c = 20; d = 70; % width
a = 40; b = 90; % height
c = 20; d = 70; % width
% a = 50; b = 100; % height
% c = 90; d = 140; % width
% a = 20; b = 140; % height
% c = 5; d = 140; % width
rect_roi1 = [a:b]; % height
rect_roi2 = [c:d]; % width
%%
BW2down = BW2(1:8:end,1:8:end);
figure;
ax0 = subplot(1,1,1);
im0 = imagesc(squeeze(phase_raw(:,:,35)));
colormap(ax0,colorcet('C06'));
set(im0, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
axis image;
v1 = [c a;d a; d b;c b];
f1 = [1 2 3 4];
hold on;
patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',2);
%%
mean_var_t2 = mean_var_t(1:8:end,1:8:end);
BW_var = ((mean_var_t2>=0)& BW2down);
indx = find(not(BW_var));
theta = angle(Axy_all);
theta(:,indx) = [];
[S s] = circ_var(theta, [], [], 2);
c1 = 1-S;
%%
clear vx_raw2 vy_raw2 vx_predict2 vy_predict2
skip = 5; zoom_scale = 6;
[vx_raw2,vy_raw2] = flow_vector_scale(vx_raw,vy_raw,skip,zoom_scale);
[vx_predict2,vy_predict2] = flow_vector_scale(vx_predict,vy_predict,skip,zoom_scale);
clear vx_raw3 vy_raw3 vx_predict3 vy_predict3
vx_predict3 = reshape(vx_predict2,[],size(vx_predict2,3));
vx_predict3(not(BW_var),:) = nan; vx_predict3 = reshape(vx_predict3, size(vx_predict2));
vy_predict3 = reshape(vy_predict2,[],size(vy_predict2,3));
vy_predict3(not(BW_var),:) = nan; vy_predict3 = reshape(vy_predict3, size(vy_predict2));
vx_raw3 = reshape(vx_raw2,[],size(vx_raw2,3));
vx_raw3(not(BW_var),:) = nan; vx_raw3 = reshape(vx_raw3, size(vx_raw2));
vy_raw3 = reshape(vy_raw2,[],size(vy_raw2,3));
vy_raw3(not(BW_var),:) = nan; vy_raw3 = reshape(vy_raw3, size(vy_raw2));
%%
scale4 = 5/8;
frame_n = numel(frame);
h1 = figure('Renderer', 'painters', 'Position', [100 100 800 400]);
for i = 1:numel(frame)   
    frame_raw = squeeze(phase_raw(:,:,frame(i)));
    
    ax1 = subplottight(4,frame_n,i);
    im1 = imagesc(frame_raw);
    colormap(ax1,colorcet('C06'));
    set(im1, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);
    
    frame_predict = squeeze(phase_predict(:,:,frame(i)));
    ax2 = subplottight(4,frame_n,frame_n+i);
    im2 = imagesc(frame_predict);
    colormap(ax2,colorcet('C06'));
    set(im2, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
    hold on;
    patch('Faces',f1,'Vertices',v1,'FaceColor','None','EdgeColor','k','lineWidth',1);
    
    ax3(i) = subplottight(4,frame_n,i+2*frame_n);
    ax3(i).Position(1) = ax3(i).Position(1)+0.005;
    ax3(i).Position(2) = ax3(i).Position(2);
    ax3(i).Position(3) = ax3(i).Position(3)-0.01;
    ax3(i).Position(4) = ax3(i).Position(4)-0.025;
    im3 = imagesc(frame_raw);
    colormap(ax3(i),colorcet('C06'));
    set(im3, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
    hold on;
    imH1Raw4 = quiver(squeeze(vx_raw3(:,:,frame(i))),squeeze(vy_raw3(:,:,frame(i))),'k','lineWidth',0.5,'autoScale','off');   
    caxis([-pi,pi]);
    ylim([a,b]);
    xlim([c,d]);
    title(num2str(c1(frame(i))));
    
    ax4(i) = subplottight(4,frame_n,3*frame_n+i);
    ax4(i).Position(1) = ax4(i).Position(1)+0.005;
    ax4(i).Position(2) = ax4(i).Position(2);
    ax4(i).Position(3) = ax4(i).Position(3)-0.01;
    ax4(i).Position(4) = ax4(i).Position(4)-0.025;
    im4 = imagesc(frame_predict);
    colormap(ax4(i),colorcet('C06'));
    set(im4, 'AlphaData', BW2down, 'AlphaDataMapping', 'scaled');
    axis image; axis off;
    hold on;
    plotOutline(maskPath([1:11]),st,atlas1,hemi,scale4,lineColor);
    hold on;
    imH1Raw4 = quiver(squeeze(vx_predict3(:,:,frame(i))),squeeze(vy_predict3(:,:,frame(i))),'k','lineWidth',0.5,'autoScale','off');  
    caxis([-pi,pi]);  
    ylim([a,b]);
    xlim([c,d]);

end
print(h1, [fname '_frame_' num2str(epochs(1))], '-dpdf', '-bestfit', '-painters');