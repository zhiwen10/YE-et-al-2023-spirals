function getSpiralSpeed(T,data_folder,save_folder)
%% load atlas brain horizontal projection and outline
load(fullfile(data_folder,'tables','horizontal_cortex_atlas_50um.mat'));
load(fullfile(data_folder,'tables',...
    'isocortex_horizontal_projection_outline.mat'));                       % 10um resolution
BW = logical(projectedAtlas1);
[maskPath,st] = get_cortex_atlas_path(data_folder);                        % get cortical atlas path and tree
root1 = '/997/';
ctx = '/997/8/567/688/';
[row,col] = find(BW);
brain_index = [col,row];
%%
params.lowpass = 0;                                                        % 0 is bandpass filter between 2-8Hz
params.gsmooth = 0;                                                        % no spatial smoothing       
Fs = 35;                                                                   % frame sampling rate
scale1 = 8;
padding = 30;
halfpadding = padding/2;
grid_x = [250,350];
grid_y = [350,550];
for kk = 1:15
    clear filteredSpirals2 filteredSpirals3 tracePhase_raw tracePhase_padded...
        angle_offset_all distance_offset_all indx2 groupedCells
    %%
    mn = T.MouseID{kk};
    tda = T.date(kk);
    en = T.folder(kk);    
    td = datestr(tda,'yyyy-mm-dd');
    tdb = datestr(tda,'yyyymmdd');
    fname = [mn '_' tdb '_' num2str(en)];
    subfolder = [mn '_' tdb '_' num2str(en)];
    session_root = fullfile(data_folder,'spirals','svd',subfolder);
    [U,V,t,mimg] = loadUVt1(session_root);                                 % load U,V, t
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    %% registration
    fname = [mn '_' tdb '_' num2str(en)];
    load(fullfile(data_folder,'spirals','rf_tform',[fname '_tform.mat']));   % load atlas transformation matrix tform;    
    load(fullfile(data_folder,'spirals','spirals_grouping',...
        [fname '_spirals_group_fftn.mat']));
    %%
    indx2 = cellfun(@(x) size(x,1), archiveCell);
    groupedCells = archiveCell(indx2>=2);
    filteredSpirals2 = cell2mat(groupedCells);
    %%
    index1 = (filteredSpirals2(:,2)>grid_x(1) & ...
        filteredSpirals2(:,2)<grid_x(2) & ...
        filteredSpirals2(:,1)>grid_y(1) & ...
        filteredSpirals2(:,1)<grid_y(2));
    filteredSpirals3 = filteredSpirals2(index1,:);
    %%
    U1 = U(1:scale1:end,1:scale1:end,1:50);
    [~,~,tracePhase_raw] = spiralPhaseMap_freq(...
        U1,dV(1:50,:),t,params,freq,rate);
    tracePhase_raw = permute(tracePhase_raw,[2,3,1]);
    [tracePhase_padded] = padZeros(tracePhase_raw,halfpadding);
    %%
    [angle_offset_all,distance_offset_all] = ...
        get_anglular_speed(...
        filteredSpirals3,tracePhase_padded,scale1,halfpadding);
    distance_offset_all1 = cellfun(@(x) mean(x(:,end),'omitnan'),...
        distance_offset_all,'UniformOutput',false);
    distance_offset_all2 = vertcat(distance_offset_all1{:});
    angle_offset_all1 = cellfun(@(x) mean(x(:,end),'omitnan'),...
        angle_offset_all,'UniformOutput',false);
    angle_offset_all2 = vertcat(angle_offset_all1{:});
    mean_angle_offset = cellfun(@(x) mean(x,1),angle_offset_all,...
        'UniformOutput',false);
    mean_angle_all = [mean_angle_offset{:}];
    mean_distance_offset = cellfun(@(x) mean(x,1),distance_offset_all,...
        'UniformOutput',false);
    save(fullfile(save_folder,fname),...
        'angle_offset_all','distance_offset_all');
end