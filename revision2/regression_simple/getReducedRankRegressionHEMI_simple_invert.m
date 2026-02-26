function getReducedRankRegressionHEMI_simple_invert(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
spath = string(st.structure_id_path);
%%
point(8,:) = [115 105]; %% VISp
point(7,:) = [97 79]; %% RSP
point(6,:) = [83,93]; %% SSp-tr
point(5,:) = [73,96]; %% SSp-ul
point(4,:) = [65, 105]; %% SSp-ll
point(3,:) = [55,118]; %% SSp-m
point(2,:) = [69,122]; %% SSp-n
point(1,:) = [84,117]; %% SSp-bfd
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll','SSp-ul','SSp-tr','RSP','VISp'};
%% cortex surface outline
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
%%
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
    mimgtransformed = imwarp(mimg,tform,'OutputView',imref2d(size(projectedTemplate1)));
    mimg1 = mimg/max(mimg(:));
    scale = 8;
    mimgT2 = mimgtransformed(1:scale:end,1:scale:end);
    %% mask and Kernel regression map for left and right
    hemi = 'right';
    [indexright,UselectedRight] = select_area(sensoryArea,spath,st,coords,...
        Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
    [Unew_right,Vnew_right] = redoSVD(UselectedRight,V(:,1:60000));
    Uright = Unew_right(:,1:50);
    Vright = Vnew_right(1:50,:);

    hemi = 'left';
    [indexleft,UselectedLeft] = select_area(sensoryArea,spath,st,coords,...
        Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
    [Unew_left,Vnew_left] = redoSVD(UselectedLeft,V(:,1:60000));
    Uleft = Unew_left(:,1:50);
    Vleft = Vnew_left(1:50,:);
    
    BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
    BW_left = BW_empty; BW_left(indexleft) =1;
    BW_right = BW_empty; BW_right(indexright) =1;
    %% prepare regressor and signal for regression  
    regressor1 = zscore(Vright,[],2);
    regressor = regressor1(:,1:20000);
    Vregressor_GPU =gpuArray(regressor);
    Vregressor_GPU = Vregressor_GPU(:,20:end-20);
    signal1 = Uleft*Vleft(:,1:20000);
    shift_t = [-10:2:10];
    signal2 = nan(numel(shift_t),size(signal1,1),size(signal1,2));
    for i = 1:numel(shift_t)
        signal2(i,:,:) = circshift(signal1,shift_t(i),2);
    end
    %%
    signal2 = signal2(:,:,20:end-20);
    signal2 = reshape(signal2,size(signal2,1)*size(signal2,2),size(signal2,3));
    signal2 = gpuArray(signal2);
    kk1 = gather(Vregressor_GPU'\signal2');
    %%
    clear signal2
    %% prepare test regressor 
    data = UselectedRight*V(:,40001:60000);    
    % subtract mean as we did before
    data = data-mean(data,2);
    regressor_test = Unew_right' * data;   
    regressor_test_z = zscore(regressor_test,[],2);
    traceSSp2 = UselectedLeft*V(:,40001:60000);
    traceSSp_predict = regressor_test_z(1:50,:)'*kk1;
    traceSSp_predict = traceSSp_predict';
    traceSSp_predict = reshape(traceSSp_predict,numel(shift_t),size(traceSSp2,1),size(traceSSp2,2));
    %%
    for m = 1:numel(shift_t)
        traceSSp2_shift = circshift(traceSSp2,shift_t(m),2);
        traceSSp2_shift = traceSSp2_shift(:,20:end-20);
        traceSSp_predict_temp = squeeze(traceSSp_predict(m,:,20:end-20));
        explained_var(m,1) = sseExplainedCal(traceSSp2_shift(:)',traceSSp_predict_temp(:)');
    end
    %%
    save(fullfile(save_folder,[fname '-hemi_simple2.mat']),'explained_var','shift_t');
    close all;
end