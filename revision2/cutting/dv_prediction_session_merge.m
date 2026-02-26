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
%%
data_folder = 'E:\task2';
list_folder = 'C:\Users\Steinmetz lab\Documents\git\YE-et-al-2023-spirals\revision2\cutting';
T1 = readtable(fullfile(list_folder,'cutting_session_list.xlsx'));
%%
T2 = T1(ismember(T1.area,'SSp') & ismember(T1.hemisphere,'bilateral'),:);
sessions = unique(T2.session_id);
%%
% for i = 1:numel(sessions)
for i = 5
    %%
    Ti = T2(T2.session_id == sessions(i) & strcmp(T2.label,'spontaneous'),:);
    mn = Ti.MouseID{1};
    tda = Ti.date(1);
    en = Ti.folder(1);
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'rfmap',mn,[fname '.mat'])); 
    session_root = fullfile(data_folder,'task_svd',fname);
    [U1,V1,t1,mimg1] = loadUVt2(session_root); 
    %%
    Ti = T2(T2.session_id == sessions(i) & strcmp(T2.label,'cutting'),:);
    mn = Ti.MouseID{1};
    tda = Ti.date(1);
    en = Ti.folder(1);
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'rfmap',mn,[fname '.mat'])); 
    session_root = fullfile(data_folder,'task_svd',fname);
    [U2,V2,t2,mimg2] = loadUVt2(session_root); 
    %% outline binary image
    clear cpstruct
    mimg1 = mimg1/max(mimg1(:));
    mimg2 = mimg2/max(mimg2(:));
    h = cpselect(mimg2,mimg1);
    %%
    [movingPoints,fixedPoints] = cpstruct2pairs(cpstruct); 
    tform = fitgeotrans(movingPoints,fixedPoints,'affine');
    %%
    sizeTemplate = size(mimg1);
    mimg2_tf = imwarp(mimg2,tform,'OutputView',imref2d(sizeTemplate));
    %%
    h1 = figure('Renderer', 'painters', 'Position', [100 100 800 600]);
    ax1 = subplot(1,1,1); 
    C = imfuse(mimg1,mimg2_tf,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    imshow(C);
    %%
    h1 = figure('Renderer', 'painters', 'Position', [100 100 600 600]);
    ax1 = subplot(1,1,1); 
    imagesc(mimg2);
    axis image;
    %%
    pointL = [320,124];
    pointR = [320,400];
    figure;
%     trace1L = squeeze(U1(pointL(1),pointL(2),:))'*V1./mimg1(pointL(1),pointL(2));
%     trace1R = squeeze(U1(pointR(1),pointR(2),:))'*V1./mimg1(pointR(1),pointR(2));
%     trace2L = squeeze(U2(pointL(1),pointL(2),:))'*V2./mimg2(pointL(1),pointL(2));
%     trace2R = squeeze(U2(pointR(1),pointR(2),:))'*V2./mimg2(pointR(1),pointR(2));
    trace1L = squeeze(U1(pointL(1),pointL(2),:))'*V1;
    trace1R = squeeze(U1(pointR(1),pointR(2),:))'*V1;
    trace2L = squeeze(U2(pointL(1),pointL(2),:))'*V2;
    trace2R = squeeze(U2(pointR(1),pointR(2),:))'*V2;
    plot(t1,trace1L,'k');
    hold on;
    plot(t1,trace1R,'r');
    hold on;
    plot(t2,trace2L+5000,'k');
    hold on;
    plot(t2,trace2R+5000,'r');
    %%
    U2t = imwarp(U2,tform,'OutputView',imref2d(size(mimg1)));
    U2ta = reshape(U2t,size(U2t,1)*size(U2t,2),size(U2t,3));
    %%
    U1a = reshape(U1,size(U1,1)*size(U1,2),size(U1,3));
    %% merge SVD
    [Unew1,Vnew1,Sv,totalVar] = mergeSVD(U1a,V1,U2ta,V2);
    save([mn '_merge_cutting.mat'],'Unew1','Vnew1','Sv','totalVar');
    %%
    figure;
    imagesc(mimg1);
    axis image;
    %%
    Unew2 = reshape(Unew1,size(mimg1,1),size(mimg1,2),size(Unew1,2));
    % point = [300,450];
    point = [200,200];
    trace1_merge = squeeze(Unew2(point(1),point(2),:))'*Vnew1;
    trace1a = squeeze(U1(point(1),point(2),:))'*V1;
    trace1b = squeeze(U2t(point(1),point(2),:))'*V2;
    trace1 = [trace1a,trace1b];
    figure;
    plot(1/35:1/35:size(trace1,2)/35,trace1_merge,'k');
    hold on;
    plot(1/35:1/35:size(trace1,2)/35,trace1,'r');
    %%
    U_raw = U2t(1:8:end,1:8:end,:);
    Unew2 = reshape(Unew1,560,560,size(Unew1,2));
    Unew2 = Unew2(1:8:end,1:8:end,:);
    U_raw1 = reshape(U_raw,size(U_raw,1)*size(U_raw,2),size(U_raw,3));
    Unew2 = reshape(Unew2,size(Unew2,1)*size(Unew2,2),size(Unew2,3));    
    V_raw = V2;
    V_predict = Vnew1(:,size(V1,2)+1:end);
    traceRaw = U_raw1(:,1:50)*V_raw(1:50,:);
    tracePredict = Unew2(:,1:50)*V_predict(1:50,:);
    explained_var = sseExplainedCal(traceRaw ,tracePredict);
    explained_var = reshape(explained_var,size(U_raw,1),size(U_raw,2));
    %%
    figure;
    imagesc(explained_var);
    caxis([0,1]);
    colorbar
    %%
    U_raw = U1(1:8:end,1:8:end,:);
    Unew2 = reshape(Unew1,560,560,size(Unew1,2));
    Unew2 = Unew2(1:8:end,1:8:end,:);
    U_raw1 = reshape(U_raw,size(U_raw,1)*size(U_raw,2),size(U_raw,3));
    Unew2 = reshape(Unew2,size(Unew2,1)*size(Unew2,2),size(Unew2,3));    
    V_raw = V1;
    V_predict = Vnew1(:,1:size(V1,2));
    traceRaw = U_raw1(:,1:50)*V_raw(1:50,:);
    tracePredict = Unew2(:,1:50)*V_predict(1:50,:);
    explained_var = sseExplainedCal(traceRaw ,tracePredict);
    explained_var = reshape(explained_var,size(U_raw,1),size(U_raw,2));
    %%
    figure;
    imagesc(explained_var)
    caxis([0,1]);
    colorbar
end