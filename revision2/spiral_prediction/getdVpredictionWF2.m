function getdVpredictionWF2(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
spath = string(st.structure_id_path);
%% right SSp index
areaPath(1) = "/997/8/567/688/695/315/453/"; % SS
areaPath(2) = "/997/8/567/688/695/315/247/"; % AUD
areaPath(3) = "/997/8/567/688/695/315/669/"; % VIS
areaPath(4) = "/997/8/567/688/695/315/254/"; % RSP
areaPath(5) = "/997/8/567/688/695/315/22/312782546/"; % VISa
areaPath(6) = "/997/8/567/688/695/315/22/417/"; % VISrl
areaPath(7) = "/997/8/567/688/695/315/541/"; % TEa
areaPath(9) = "/997/8/567/688/695/315/677/"; % VISC
areaPath(10) = "/997/8/567/688/695/1089/"; % HPF
areaPath(11) = "/997/8/567/688/695/315/677/"; % CTXsp
sensoryArea = strcat(areaPath(:));
%%
for kk = 1:15
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    %% load SVD
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];                                % get derivative of V
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));
    %%
    U = U(:,:,1:50);
    V = V(1:50,:);
    %%
    Utransformed = imwarp(U,tform,...
        'OutputView',imref2d(size(projectedTemplate1)));
    mimgtransformed = imwarp(mimg,tform,...
        'OutputView',imref2d(size(projectedTemplate1)));
    mimg1 = mimg/max(mimg(:));
    scale = 8;
    mimgT2 = mimgtransformed(1:scale:end,1:scale:end);
    %% mask and Kernel regression map for SSp and MO
    hemi = 'right';
    sensoryArea = [];
    [indexright,UselectedRight] = select_area(sensoryArea,spath,st,coords,...
        Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
    [Unew_right,Vnew_right] = redoSVD(UselectedRight,V);
    Uright = Unew_right(:,1:50);
    Vright = Vnew_right(1:50,:);
    %% prepare regressor and signal for regression
    Vright1 = zscore(Vright,[],2);
    %%
%     regressor = VSSp1(:,1:30000);
%     Vregressor_GPU = gpuArray(regressor);
%     signal1 = V(:,1:30000);
%     signal_GPU = gpuArray(signal1);
%     kk1 = gather(Vregressor_GPU'\signal_GPU');
%     % prepare test regressor 
%     regressor_test_z = VSSp1(:,30001:58000);
%     V_predict = regressor_test_z'*kk1;
%     V_predict = V_predict';
    %% get 10 segmented sample index
    sample_length = size(V,2);
    epochN = 10;
    epoch_indx = [];
    epochSize = floor(sample_length/epochN);
    for i = 1:epochN
        epoch_indx(i,1) = 1+epochSize*(i-1);
        epoch_indx(i,2) = epochSize*i;
    end
    interval = epoch_indx(:,2)-epoch_indx(:,1);
    %% seperate into 5 odd and 5 even series, and use odd to predict even
    V1_odd = [];
    Vright_odd = [];
    odd_indx = [];
    train_indx1 = [];
    for i = 1:5
        train_indx1 = [epoch_indx(2*i-1,1):epoch_indx(2*i-1,2)];
        odd_indx = [odd_indx,train_indx1];   
    end
    Vright_odd = Vright1(:,odd_indx);
    V1_odd = V(:,odd_indx);
    kk1 = Vright_odd'\V1_odd';
 
    even_indx = [];
    Vright_even_temp = [];
    V1_predict_even = [];
    for i = 1:5
        even_indx = [epoch_indx(2*i,1):epoch_indx(2*i,2)];
        Vright_even_temp = Vright1(:,even_indx);
        V1_predict_even(:,:,i) = Vright_even_temp'*kk1;
    end
    %% seperate into 5 odd and 5 even series, and use even to predict odd
    V1_even2 = [];
    Vright_even2 = [];
    even_indx = [];
    train_indx1 = [];
    for i = 1:5
        train_indx1 = [epoch_indx(2*i,1):epoch_indx(2*i,2)];
        even_indx = [even_indx,train_indx1];   
    end
    Vright_even2 = Vright1(:,even_indx);
    V1_even2 = V(:,even_indx);
    kk1 = Vright_even2'\V1_even2';
    
    odd_indx = [];
    Vright_odd2_temp = [];
    V1_predict_odd = [];
    for i = 1:5
        odd_indx = [epoch_indx(2*i-1,1):epoch_indx(2*i-1,2)];
        Vright_odd2_temp = Vright1(:,odd_indx);
        V1_predict_odd(:,:,i) = Vright_odd2_temp'*kk1;
    end

    %% fully reconstruct the entire duration from predicted odd and even epochs
    V_predict = [];
    for i = 1:10    
        if mod(i,2)==1
            indx = (i-1)/2+1;
            V_predict = [V_predict, squeeze(V1_predict_odd(:,:,indx))'];
        else
            indx = i/2;
            V_predict = [V_predict, squeeze(V1_predict_even(:,:,indx))'];
        end
    end
    %% get variance explained by averaging across the 10 epochs
    clear V_epochs V_predict_epochs
    epochN = 10;
    V2 = V(:,1:size(V_predict,2));
    for i = 1:epochN
        V_epochs(:,:,i) = V2(:,epoch_indx(i,1):epoch_indx(i,2));
        V_predict_epochs(:,:,i) = V_predict(:,epoch_indx(i,1):epoch_indx(i,2));
    end
    %%
    Ut = Utransformed(1:8:end,1:8:end,:);
    explained_var_all = get_variance_explained(Ut,V_epochs,V_predict_epochs);
    explained_var_all(explained_var_all(:)<0) = 0;
    %% save predicted dV and variance explained
    save(fullfile(save_folder,[fname '_v_predict.mat']),'V_predict','explained_var_all'); 
end