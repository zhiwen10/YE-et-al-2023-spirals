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
%% get original and predicted dV by kernel regression.
codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\compare_flow1';
T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
% kk = 16; % zye41 STR
% kk = 11; % zye60
% kk = 39; % zye67
% kk = 8; % zye58
for kk = 8
    ops = get_session_info(T,kk);
    [U,V,t,mimg] = get_wf_svd(ops.serverRoot);
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    %%
    folder = 'C:\Users\Steinmetz lab\Documents\MATLAB\YE2023\spiral_detection\spirals';
    dfolder = fullfile(folder,ops.mn,ops.td,num2str(ops.en));
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    %% load spike data
    goodcluster = ops.curation;
    [sp] = get_spikes(ops.ksSubFolder,goodcluster);
    %%
    [syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT(ops,t);
    WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));
    dV1 = double(dV(:,~isnan(WF2ephysT1)));
    V1 = double(V(1:50,~isnan(WF2ephysT1)));
    [MUA_std] = get_MUA_bin(sp,WF2ephysT);
    dV1 = double(dV1(1:50,:));
    %% no registration first for flow vector
    scale = 4;
    Ut = U(1:scale:end,1:scale:end,1:50);
    t1 = t(~isnan(WF2ephysT1));
    kernel_t = [-0.5,0.5];
    %% get 10 segment sample index
    sample_length = size(dV1,2);
    epochN = 10;
    epochSize = floor(sample_length/epochN);
    for i = 1:epochN
        epoch_indx(i,1) = 1+epochSize*(i-1);
        epoch_indx(i,2) = epochSize*i;
    end
    interval = epoch_indx(:,2)-epoch_indx(:,1);
    %%
    dV1_odd = [];
    odd_indx = [];
    for i = 1:5
        train_indx1 = [epoch_indx(2*i-1,1):epoch_indx(2*i-1,2)];
        odd_indx = [odd_indx,train_indx1];   
    end
    MUA_std_odd = MUA_std(:,odd_indx);
    dV1_odd = dV1(:,odd_indx);
    even_indx = [];
    MUA_std_even = [];
    dV1_even = [];
    for i = 1:5
        even_indx = [epoch_indx(2*i,1):epoch_indx(2*i,2)];
        MUA_std_even(:,:,i) = MUA_std(:,even_indx);
        dV1_even(:,:,i) = dV1(:,even_indx);
    end
    %%
    [dV_predict_even,k_cv] = mua_prediction_full(dV1_odd,MUA_std_odd,dV1_even,MUA_std_even,kernel_t);
    %%
    dV_raw_even = [];
    even_indx = [];
    for i = 1:5
        train_indx1 = [epoch_indx(2*i,1):epoch_indx(2*i,2)];
        even_indx = [even_indx,train_indx1];   
    end
    MUA_std_even2 = MUA_std(:,even_indx);
    dV1_even2 = dV1(:,even_indx);
    odd_indx = [];
    MUA_std_odd2 = [];
    dV1_odd2 = [];
    for i = 1:5
        odd_indx = [epoch_indx(2*i-1,1):epoch_indx(2*i-1,2)];
        MUA_std_odd2(:,:,i) = MUA_std(:,odd_indx);
        dV1_odd2(:,:,i) = dV1(:,odd_indx);
    end
    [dV_predict_odd,k_cv] = mua_prediction_full(dV1_even2,MUA_std_even2,dV1_odd2,MUA_std_odd2,kernel_t);
    %%
    dV_predict = [];
    for i = 1:10    
        if mod(i,2)==1
            indx = (i-1)/2+1;
            dV_predict = [dV_predict, squeeze(dV_predict_odd(:,:,indx))];
        else
            indx = i/2;
            dV_predict = [dV_predict, squeeze(dV_predict_even(:,:,indx))];
        end
    end
    %%
    epochN = 10;
    dV1 = dV1(:,1:size(dV_predict,2));
    for i = 1:epochN
        dV1_epochs(:,:,i) = dV1(:,epoch_indx(i,1):epoch_indx(i,2));
        dV_predict_epochs(:,:,i) = dV_predict(:,epoch_indx(i,1):epoch_indx(i,2));
    end
    %%
    explained_var_all = get_variance_explained(Ut,dV1_epochs,dV_predict_epochs);
    % cmap1 = flipud(cbrewer2('div', 'RdYlBu', 32,'linear'));
    %%
    % cmap1 = colorcet('L17');
    cmap1 = flipud(inferno);
    h1 = figure; 
    ax1 = subplot(1,1,1);
    imagesc(squeeze(mean(explained_var_all,3)));
    axis image; axis off;
    colormap(ax1,cmap1);
    colorbar;
    % caxis([0,1]);
    print(h1, [fname '_variance'], '-dpdf', '-bestfit', '-painters');
    %%
    save([fname '_dv_predict'],'dV_predict');    
end