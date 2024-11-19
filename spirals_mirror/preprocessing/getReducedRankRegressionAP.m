function getReducedRankRegressionAP(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
spath = string(st.structure_id_path);
%%
point(6,:) = [83,93]; %% SSp-tr
point(5,:) = [73,96]; %% SSp-ul
point(4,:) = [65, 105]; %% SSp-ll
point(3,:) = [55,118]; %% SSp-m
point(2,:) = [69,122]; %% SSp-n
point(1,:) = [84,117]; %% SSp-bfd
point(8,:) = [115 105]; %% VISp
point(7,:) = [97 79]; %% RSP
nameList = {'SSp-bfd','SSp-n','SSp-m','SSp-ll','SSp-ul','SSp-tr','RSP','VISp'};
%% cortex surface outline
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
%% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frotalArea = strcat(areaPath(:));
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
    Utransformed = imwarp(U,tform,...
        'OutputView',imref2d(size(projectedTemplate1)));
    mimgtransformed = imwarp(mimg,tform,...
        'OutputView',imref2d(size(projectedTemplate1)));
    mimg1 = mimg/max(mimg(:));
    scale = 8;
    mimgT2 = mimgtransformed(1:scale:end,1:scale:end);
    %% mask and Kernel regression map for SSp and MO
    hemi = 'right';
    [indexSSp,UselectedSSp] = select_area(sensoryArea,spath,st,coords,...
        Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
    [Unew_SSp,Vnew_SSp] = redoSVD(UselectedSSp,V(:,1:58000));
    USSp = Unew_SSp(:,1:50);
    VSSp = Vnew_SSp(1:50,:);

    [indexMO,UselectedMO] = select_area(frotalArea,spath,st,coords,...
        Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
    [Unew_MO,Vnew_MO] = redoSVD(UselectedMO,V(:,1:58000));
    UMO = Unew_MO(:,1:50);
    VMO = Vnew_MO(1:50,:);
    %%
    BW_empty = zeros(size(projectedAtlas1(1:scale:end,1:scale:end)));
    %% prepare regressor and signal for regression
    regressor1 = zscore(VMO,[],2);
    regressor = regressor1(:,1:40000);
    Vregressor_GPU =gpuArray(regressor);

    signal1 = USSp*VSSp(:,1:40000);
    signal_GPU = gpuArray(signal1);
    %% reduced rank regression regression
    [a1, b1, R] = CanonCor2(signal_GPU', Vregressor_GPU');
    a = gather(a1); b = gather(b1); R2 = gather(R);
    kk1 = b(:,1:50) *a(:,1:50)';
    k1_real = UMO*kk1;
    %% prepare test regressor 
    data = UselectedMO*V(:,40001:60000);    
    % subtract mean as we did before
    data = data-mean(data,2);
    regressor_test = Unew_MO' * data;   
    regressor_test_z = zscore(regressor_test,[],2);
    %%
    kernel_full2 = zeros(size(BW_empty,1)*size(BW_empty,2),...
        size(BW_empty,1)*size(BW_empty,2));
    for j = 1:numel(indexSSp)
         kernel_temp2 = zeros(size(BW_empty,1)*size(BW_empty,2),1);
         kernel_temp2(indexMO) = k1_real(:,j);
         kernel_full2(:,indexSSp(j)) = kernel_temp2;
    end
    kernel_full2 = reshape(kernel_full2,[size(mimgT2,1),size(mimgT2,2),...
    size(mimgT2,1),size(mimgT2,2)]);
%%
    traceSSp2 = UselectedSSp*V(:,40001:60000);
    for n = 1:50
        traceSSp_predict = regressor_test_z(1:50,:)' * b(:,1:n) * a(:,1:n)';
        traceSSp_predict = traceSSp_predict';
        explained_var5(n,1) = sseExplainedCal(traceSSp2(:)',traceSSp_predict(:)');
    end
    %%
    save(fullfile(save_folder,[fname '-AP.mat']),...
        'kernel_full2','explained_var5','a','b');
    close all;
end