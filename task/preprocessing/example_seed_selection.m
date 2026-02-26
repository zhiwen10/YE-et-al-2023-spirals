%% spiral summary
githubdir = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\git';
addpath(genpath(fullfile(githubdir, 'Pipelines'))); % https://github.com/SteinmetzLab/Pipelines
addpath(genpath(fullfile(githubdir, 'widefield'))); % https://github.com/cortex-lab/widefield
addpath(genpath(fullfile(githubdir, 'npy-matlab'))); % https://github.com/kwikteam/npy-matlab
addpath(genpath('C:\Users\Steinmetz lab\OneDrive - UW\Documents\MATLAB\widefield_DIY\task'));
%% load atlas brain horizontal projection and outline
data_folder1 = 'E:\spiral_data_share\data';   
load(fullfile(data_folder1,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder1,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder1);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
BW = logical(projectedAtlas1);
%%
spiral_folder = 'C:\Users\Steinmetz lab\OneDrive - UW\Documents\MATLAB\widefield_DIY\task\code_20240904\spirals';
freq = [2,8];
freq_folder = [num2str(freq(1)) '_' num2str(freq(2)) 'Hz'];
spiral_folder1 = fullfile(spiral_folder,freq_folder);
fnames = {'ZYE_0085','ZYE_0088','ZYE_0090','ZYE_0091'};
%%
frames_post = 71:85;
framesN = numel(frames_post);
%% sessions and trials with large spirals
id = 1;
mn = fnames{id};
load(fullfile(spiral_folder1,[mn '_spirals_task_sort.mat']));
load(fullfile(spiral_folder,[mn '_task_outcome.mat']));
load(fullfile(spiral_folder,[mn '_task_trial_ID.mat']));
spiral_large_count = zeros(size(spiral_all,1),1);
for i = 1:size(spiral_all,1)
    spiral_temp = cat(1,spiral_all{i,frames_post});
    spiral_large_count(i,1) = sum(spiral_temp(:,3)>=80);
end
index = (spiral_large_count>=5 & T_all.label=="correct" & T_all.left_contrast >= 0.25);
trial_select = trial_all(index,:);
%%
data_folder = 'E:\task';
mn = 'ZYE_0085';
T_session = readtable(fullfile(data_folder,[mn '.xlsx']));
T1 = T_session(T_session.label == "task",:);
T1 = T1((T1.hit_left>0.7 & T1.hit_right>0.7),:);
sessions = size(T1,1);
%%
trial = 6;
trial_select1 = trial_select(trial_select.session == trial,:);
mn = char(T1.MouseID{trial});
tda = T1.date(trial);  
en = T1.folder(trial); 
block_en = T1.block_folder(trial); 
td = datestr(tda,'yyyy-mm-dd');
tdb = datestr(tda,'yyyymmdd');
%% load data and block 
fname = [mn '_' tdb '_' num2str(en)];
session_root = fullfile(data_folder,'task_svd',fname);
[U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
expRoot = expPath(mn, td, block_en);
expDir = dir(fullfile(expRoot,'*_Block.mat'));
load(fullfile(expDir.folder,expDir.name));
% get photodiode time
serverRoot = expPath(mn,td,en);
win = [0,5000];
allPD2 = getPhotodiodeTime(serverRoot,win);
%%
ntrial = numel(block.events.endTrialValues);
norepeatValues = block.events.repeatNumValues(1:ntrial);            % sometimes, last trial don't have a response
response = block.events.responseValues(1:ntrial);
norepeat_indx = (norepeatValues==1)';
allPD2 = allPD2(1:ntrial);
allPD2 = allPD2(norepeat_indx);
%% load atlas registration
downscale = 16;
load(fullfile(data_folder,'rfmap',mn,[fname '.mat']));
sizeTemplate = [1320,1140];
Ut = imwarp(U,tform,'OutputView',imref2d(sizeTemplate));
mimgt = imwarp(mimg,tform,'OutputView',imref2d(sizeTemplate));
mimgt = mimgt(1:downscale:end,1:downscale:end); 
Ut1 = Ut(1:downscale:end,1:downscale:end,1:50); 
dV1 = double(dV(1:50,:));
%% get time index for all events
allPD3 = allPD2(trial_select1.trial);
clear indx
for i = 1:numel(allPD3)
    ta  = t-allPD3(i);
    indx(i,1) = find(ta>0, 1,'first');
    t2(i,1) = t(indx(i,1));
end
%%
% load('mask_ZYE_0091.mat');
lineColor = 'k';
downscale = 16;
downscale1 = 1;
scale1 = 5;
hemi = [];
% BW3 = BW2(1:downscale:end,1:downscale:end);
BW3 = BW(1:downscale:end,1:downscale:end);
downscale = 16;
scale3 = 5/16;
lineColor = 'k'; lineColor1 = 'w';
%%
fig_total = floor(size(indx,1)/3);
for jj = 1:fig_total
    h1 =figure('Renderer', 'painters', 'Position', [100 100 950 950]);
    for kk = 1:3
        kk1 = (jj-1)*3+kk;
        % kk1 = 102+kk;
        indx1 = indx(kk1);
        indx_trial1 = trial_select1.trial(kk1);
        %%
        Ut1a = double(reshape(Ut1,size(Ut1,1)*size(Ut1,2),size(Ut1,3)));
        trace = Ut1a*dV1(:,indx1:indx1+35);
        trace2d = reshape(trace,size(Ut1,1),size(Ut1,2),size(trace,2))./mimgt;
        trace2d = trace2d-trace2d(:,:,1);
        cmax = prctile(trace2d(:),99.9); cmin = prctile(trace2d(:),0.1);
        cmax = 0.03; cmin = -0.03;
        %%
        freq = [2,8];
        Fs = 35;
        trace1 = Ut1a*dV1(:,indx1-70:indx1+70);
        meanTrace = trace1 -mean(trace1,2);
        meanTrace = meanTrace';  
        [f1,f2] = butter(2,freq/(Fs/2), 'bandpass');
        traceFilt = filtfilt(f1,f2,meanTrace);
        traceHilbert =hilbert(traceFilt);
        tracePhase = angle(traceHilbert);
        tracePhase = tracePhase';
        tracePhase = reshape(tracePhase,size(Ut1,1),size(Ut1,2),size(tracePhase,2));
        tracePhase = tracePhase(:,:,71:106);
        
        meanTrace2 = meanTrace';
        meanTrace2 = reshape(meanTrace2,size(Ut1,1),size(Ut1,2),size(meanTrace2,2));
        meanTrace2 = meanTrace2(:,:,71:106);
        cmax = prctile(meanTrace2(:),99.9); cmin = prctile(meanTrace2(:),0.1);
        
        traceFilt2 = traceFilt';
        traceFilt2 = reshape(traceFilt2,size(Ut1,1),size(Ut1,2),size(traceFilt2,2));
        traceFilt2 = traceFilt2(:,:,71:106);
        cmax2 = prctile(traceFilt2(:),99.9); cmin2 = prctile(traceFilt2(:),0.1);
        %%
        for i = 1:18
            ax4 = subplottight(9,18,(kk-1)*18*3+i);
            im1 = imagesc(squeeze(meanTrace2(:,:,i)));
            % im1 = imagesc(squeeze(meanTrace2(:,:,i)));
            hold on;
            plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
            axis image; axis off;
            set(im1, 'AlphaData', BW3, 'AlphaDataMapping', 'scaled');
            colormap(ax4,'parula');
            caxis([cmin,cmax]);
            if i==1
                title(['Trial' num2str(indx_trial1) '_rt_'...
                    'L' num2str(round(trial_select.left_contrast(kk1)*100))],...
                    'interpreter','none');
            end
            
            ax6 = subplottight(9,18,kk*18*3-18*2+i);
            im1 = imagesc(squeeze(traceFilt2(:,:,i)));
            % im1 = imagesc(squeeze(meanTrace2(:,:,i)));
            hold on;
            plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
            axis image; axis off;
            set(im1, 'AlphaData', BW3, 'AlphaDataMapping', 'scaled');
            colormap(ax6,'parula');
            caxis([cmin2,cmax2]);

            ax5 = subplottight(9,18,kk*18*3-18+i);
            im1 = imagesc(squeeze(tracePhase(:,:,i)));
            hold on;
            plotOutline(maskPath([1:11]),st,atlas1,hemi,scale3,lineColor);
            axis image; axis off;
            set(im1, 'AlphaData', BW3, 'AlphaDataMapping', 'scaled');
            colormap(ax5,colorcet('C06'));
            caxis([-pi,pi]);        
        end
    %     cb4 = subplottight(10,18,2*18*trial-18);
    %     im1 = imagesc(squeeze(trace2d(:,:,18)),'visible','off');
    %     colormap(ax4,'parula');
    %     caxis([cmin,cmax]);
    %     axis off;
    %     cb = colorbar;
    %     a =  cb.Position; %gets the positon and size of the color bar
    %     set(cb,'Position',[a(1) a(2)+0.05 a(3) a(4)/8])% To change size
    end
    %%
    print(h1,[mn '_wf_map_mean_trials_spirals_large3-' num2str(jj)],'-dpdf', '-bestfit', '-painters');
    close all;
end