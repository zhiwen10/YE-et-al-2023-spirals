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
addpath(genpath(fullfile(githubDir, 'spiralDetection')));
%%
load(fullfile('C:\Users\Steinmetz lab\Documents\git\allenCCF\Browsing Functions', 'ctxOutlines.mat'));
load('C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\atlasMaskArea\segmentationCode\projectedOutlineAtlas.mat');
[projectedAtlas1,projectedTemplate1] = filterProjectedAtlas(projectedAtlas,projectedTemplate);
%% get original and predicted dV by kernel regression.
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
codefolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareSpiral\compareSpiral_batch';
T = readtable(fullfile(codefolder, 'session_list_sorted2.csv'));
T1 = T((contains(T.Area,area(1,[1,2,4])) & T.meanVar>=0.1),:);
%%
rfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\compareFlowField\variance_summary\variance_order';
gfolder = 'C:\Users\Steinmetz lab\Documents\MATLAB\widefield_DIY\phase\reducedRankRegression\registration_4x';
BW = logical(projectedAtlas1);
% for kk = 1:size(T1,1) 
for kk = 28
    clear ops U V t mimg dV sp Ut syncTL syncProbe WF2ephysT1 WF2ephysT dV1 V1 MUA_std
    ops = get_session_info(T1,kk);
    filename = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_explained_var'];
    ffname = fullfile(rfolder,[filename '.mat']);
    if not(isfile(ffname))
        [U,V,t,mimg] = get_wf_svd(ops.serverRoot);
        dV = [zeros(size(V,1),1) diff(V,[],2)];
        % load spike data
        [sp] = get_spikes(ops.ksSubFolder,ops.curation);
        % no registration first
        scale = 4;
        Ut = U(1:scale:end,1:scale:end,1:50);
        % prediction from spiking data
        [syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT(ops,t);
        WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));
        dV1 = double(dV(:,~isnan(WF2ephysT1)));
        V1 = double(V(1:50,~isnan(WF2ephysT1)));
        % depth = T.depth_all(indx(kk));
        [MUA_std] = get_MUA_bin(sp,WF2ephysT);
        spike_n = size(MUA_std,1);
        dV1 = double(dV1(1:50,:));
        %% first poredict using single unit, and order single units by contribution
        explained_var_all = [];
        explained_var = [];
        for n_select = 1:spike_n
            MUA_std1 = MUA_std(n_select,:);
            [dV_raw,dV_predict,epoch_indx] = get_prediction(dV1,MUA_std1);
            explained_var = get_variance_explained(Ut,dV_raw,dV_predict);
            explained_var = mean(explained_var,3,'omitnan');
            explained_var_all(:,:,n_select) = explained_var;
        end        
        %%
        reg_name = [ops.mn '_' ops.tdb '_' num2str(ops.en) '_tform_4x.mat'];
        reg_fullname =fullfile(gfolder,reg_name);
        load(reg_fullname);
        sizeTemplate = size(projectedTemplate1);
        var_reg = imwarp(explained_var_all,tform,'OutputView',imref2d(sizeTemplate));

        var_reg1 = [];
        for i = 1:size(var_reg,3)
            temp = squeeze(var_reg(:,:,i));
            temp(~BW) = nan;
            var_reg1(:,:,i) = temp;
        end
        mean_var = squeeze(mean(var_reg1,[1,2],'omitnan'));

        [B,I] = sort(mean_var,'descend');
        %% increase number of units increamentally by contribution
        explained_var_all1 = [];
        explained_var1 = [];
        count1 = 1;
        for n_select = 1:spike_n
            %% full prediction
            MUA_std1 = MUA_std(I(1:n_select),:);
            [dV_raw,dV_predict,epoch_indx] = get_prediction(dV1,MUA_std1);
            explained_var1 = get_variance_explained(Ut,dV_raw,dV_predict);
            explained_var1 = mean(explained_var1,3,'omitnan');
            explained_var_all1(:,:,n_select) = explained_var1;
        end
        sizeTemplate = size(projectedTemplate1);
        var_reg2 = imwarp(explained_var_all1,tform,'OutputView',imref2d(sizeTemplate));
        var_reg3 = [];
        for i = 1:size(var_reg2,3)
            temp = squeeze(var_reg2(:,:,i));
            temp(~BW) = nan;
            var_reg3(:,:,i) = temp;
        end
        mean_var2 = squeeze(mean(var_reg3,[1,2],'omitnan'));
        save(ffname,'I','B','mean_var','explained_var_all1','mean_var2');
    end
end

