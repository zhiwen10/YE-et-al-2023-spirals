%% widefield data preprocess
githubdir = 'C:\Users\Steinmetz lab\Documents\git';
addpath(genpath(fullfile(githubdir, 'YE-et-al-2023-spirals'))); 
%%
data_folder = 'E:\spiral_data_share\data';  
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables','horizontal_cortex_template_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
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
% MO index
clear areaPath spath3
areaPath(1) = "/997/8/567/688/695/315/500/"; % MO
areaPath(2) = "/997/8/567/688/695/315/31/"; % ACA
areaPath(3) = "/997/8/567/688/695/315/972/"; % PL
areaPath(4) = "/997/8/567/688/695/315/44/"; % ILA
frontalArea = strcat(areaPath(:));
hemi = [];
%%
data_folder = 'E:\task2';
list_folder = 'C:\Users\Steinmetz lab\Documents\git\YE-et-al-2023-spirals\revision2\cutting';
T1 = readtable(fullfile(list_folder,'cutting_session_list.xlsx'));
%%
T2 = T1(ismember(T1.area,'SSp') & ismember(T1.hemisphere,'bilateral'),:);
% T2 = T1(ismember(T1.area,'SSp') & ismember(T1.hemisphere,'right'),:);
sessions = unique(T2.session_id);
%%
merge_folder = 'E:\manipulation_prediction\cutting_prediction\merged_data';
%%
for i = 1:numel(sessions)
    %%
    Ti = T2(T2.session_id == sessions(i) & strcmp(T2.label,'spontaneous'),:);
    mn = Ti.MouseID{1};
    tda = Ti.date(1);
    en = Ti.folder(1);
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'rfmap',mn,[fname '.mat'])); 
    session_root = fullfile(data_folder,'task_svd',fname);
    [U,V,t,mimg] = loadUVt2(session_root); 
    %%
    load(fullfile(merge_folder,[mn '_merge_cutting.mat']));
    Unew1 = reshape(Unew1,size(mimg,1),size(mimg,2),size(Unew1,2));
    load(fullfile(data_folder,'rfmap',mn,[fname '.mat']));
    %%
    Unew1 = Unew1(:,:,1:50);
    Vnew1 = Vnew1(1:50,:);
    %%
    Utransformed = imwarp(Unew1,tform,...
        'OutputView',imref2d(size(projectedTemplate1)));
    mimgtransformed = imwarp(mimg,tform,...
        'OutputView',imref2d(size(projectedTemplate1)));
    mimg1 = mimg/max(mimg(:));
    scale = 8;
    mimgT2 = mimgtransformed(1:scale:end,1:scale:end);
    %%
    [indexMO,UselectedMO] = select_area(frontalArea,spath,st,...
        coords,Utransformed,projectedAtlas1,projectedTemplate1,hemi,scale);
    %% mask and Kernel regression map before muscimol
    [Unew_MO,Vnew_MO] = redoSVD(UselectedMO,Vnew1);
    UMO = Unew_MO(:,1:50);
    VMO = Vnew_MO(1:50,:);
    VMO1 = zscore(VMO,[],2);
    kk1 = VMO1(:,1:size(V,2))'\Vnew1(:,1:size(V,2))';
    %%
    U1 = Unew1(1:8:end,1:8:end,:);
    Ut = reshape(U1,size(U1,1)*size(U1,2),size(U1,3));
    %%
%     data = UselectedMO*Vnew1(:,size(V,2)+1:end);
%     data = data-mean(data,2);
%     V2 = UMO' * data;
%     V2a = zscore(V2,[],2);
%     V_predict = V2a'*kk1;

    V_predict = VMO1(:,size(V,2)+1:end)'*kk1;  
    V_predict = V_predict';
    %%
    explained_var_all = get_variance_explained(U1,Vnew1(:,size(V,2)+1:end),V_predict);
    % explained_var_all = get_variance_explained(U1,V_predict,Vnew1(:,size(V,2)+1:end));
    % explained_var_all(explained_var_all(:)<0) = 0;
    %%
    save([mn '_prediction_MO_cutting.mat'],'Unew1','V_predict','explained_var_all');
    %%
    figure;
    imagesc(explained_var_all);
    caxis([0,1]);
    colorbar;
    title(fname);
        %%
    % point = [47,52];
    % point = [39,42];
    point = [10,30];
    trace1_merge = squeeze(U1(point(1),point(2),:))'*Vnew1(:,size(V,2)+1:end);
    trace1a = squeeze(U1(point(1),point(2),:))'*V_predict;
    figure;
    plot(1/35:1/35:size(trace1a,2)/35,trace1_merge,'k');
    hold on;
    plot(1/35:1/35:size(trace1a,2)/35,trace1a,'r');
end