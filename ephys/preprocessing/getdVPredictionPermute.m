function getdVPredictionPermute(T,data_folder,save_folder)
%% get original and predicted dV by kernel regression.
area{1} = 'THAL';
area{2} = 'STR';
area{3} = 'CORTEX';
area{4} = 'MB';
%%
for kk = 1:size(T,1)
    clearvars -except T kk data_folder save_folder
    %% load widefield and spiking data
    ops = get_session_info2(T,kk,data_folder);
    [U,V,t,mimg] = loadUVt1(ops.session_root);                             % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    [sp] = loadKSdir2(ops.session_root);
    [syncTL,syncProbe,WF2ephysT1] = get_wf2ephysT2(ops,t);
    WF2ephysT = WF2ephysT1(~isnan(WF2ephysT1));
    dV1 = double(dV(:,~isnan(WF2ephysT1)));
    V1 = double(V(1:50,~isnan(WF2ephysT1)));
    t1 = t(~isnan(WF2ephysT1));
    [MUA_std] = get_MUA_bin(sp,WF2ephysT);
    dV1 = double(dV1(1:50,:));
    fname = [ops.mn '_' ops.tdb '_' num2str(ops.en)];
    %%
    kernel_t = [-0.5,0.5];    
    scale = 4;
    Ut = U(1:scale:end,1:scale:end,1:50);
    %% get 10 segment sample index
    sample_length = size(dV1,2);
    epochN = 10;
    epochSize = floor(sample_length/epochN);
    for i = 1:epochN
        epoch_indx(i,1) = 1+epochSize*(i-1);
        epoch_indx(i,2) = epochSize*i;
    end
    interval = epoch_indx(:,2)-epoch_indx(:,1);
    %% permute the index of spike clusters
    perm_indx = randperm(size(MUA_std,1));
    %% seperate into 5 odd and 5 even series, and use odd to predict even
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
    MUA_std_even = MUA_std_even(perm_indx,:,:);
    [dV_predict_even,k_cv] = mua_prediction_full(dV1_odd,MUA_std_odd,dV1_even,MUA_std_even,kernel_t);
    %% seperate into 5 odd and 5 even series, and use even to predict odd
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
    MUA_std_odd2 = MUA_std_odd2(perm_indx,:,:);
    [dV_predict_odd,k_cv] = mua_prediction_full(dV1_even2,MUA_std_even2,dV1_odd2,MUA_std_odd2,kernel_t);
    %% fully reconstruct the entire duration from predicted odd and even epochs
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
    %% get variance explained by averaging across the 10 epochs
    epochN = 10;
    dV1 = dV1(:,1:size(dV_predict,2));
    for i = 1:epochN
        dV1_epochs(:,:,i) = dV1(:,epoch_indx(i,1):epoch_indx(i,2));
        dV_predict_epochs(:,:,i) = dV_predict(:,epoch_indx(i,1):epoch_indx(i,2));
    end
    explained_var_all = get_variance_explained(Ut,dV1_epochs,dV_predict_epochs);
    explained_var_all(explained_var_all(:)<0) = 0;
    %% save predicted dV and variance explained
    save(fullfile(save_folder,[fname '_dv_predict.mat']),'dV_predict','explained_var_all'); 
end